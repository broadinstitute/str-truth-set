import hail as hl
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/gangstr@sha256:0fb6333120be32b715b5edf9b1b0d7d166fe3b3ee6322fd33f97ddccd3784ac3"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

REPEAT_SPECS_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/gangstr/positive_loci.GangSTR.*_of_015.bed"
REPEAT_SPECS_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/gangstr/negative_loci.GangSTR.*_of_015.bed"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/gangstr"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    parser_group.add_argument("--repeat-specs", action="append", help="Path of repeat spec bed file(s) to process")

    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    if args.positive_loci:
        positive_or_negative_loci = "positive_loci"
        repeat_spec_file_stats_list = hl.hadoop_ls(REPEAT_SPECS_POSITIVE_LOCI)
        if len(repeat_spec_file_stats_list) == 0:
            raise ValueError(f"No files found matching {REPEAT_SPECS_POSITIVE_LOCI}")
    elif args.negative_loci:
        positive_or_negative_loci = "negative_loci"
        repeat_spec_file_stats_list = hl.hadoop_ls(REPEAT_SPECS_NEGATIVE_LOCI)
        if len(repeat_spec_file_stats_list) == 0:
            raise ValueError(f"No files found matching {REPEAT_SPECS_NEGATIVE_LOCI}")
    elif args.repeat_specs:
        positive_or_negative_loci = os.path.basename(args.repeat_specs[0]).replace(".bed", "").replace(".gz", "")
        repeat_spec_file_stats_list = [{"path": path} for path in args.repeat_specs]
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: GangSTR: {positive_or_negative_loci}: {bam_path_ending}")
    output_dir = os.path.join(args.output_dir, positive_or_negative_loci)
    if not args.force:
        json_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.json"))
        logging.info(f"Precached {len(json_paths)} json files")

    step1s = []
    step1_output_json_paths = []
    for repeat_spec_i, repeat_spec_file_stats in enumerate(repeat_spec_file_stats_list):
        repeat_spec_path = repeat_spec_file_stats["path"]

        if args.n and repeat_spec_i >= args.n:
            break

        s1 = bp.new_step(f"Run GangSTR #{repeat_spec_i}",
                         arg_suffix=f"gangstr",
                         step_number=1,
                         image=DOCKER_IMAGE,
                         cpu=1,
                         storage="100Gi")
        step1s.append(s1)

        local_fasta = s1.input(args.reference_fasta, localize_by=Localize.COPY)
        if args.reference_fasta_fai:
            s1.input(args.reference_fasta_fai, localize_by=Localize.COPY)

        local_bam = s1.input(args.input_bam, localize_by=Localize.COPY)
        if args.input_bai:
            s1.input(args.input_bai, localize_by=Localize.COPY)

        local_repeat_spec = s1.input(repeat_spec_path)

        output_prefix = re.sub(".bed$", "", local_repeat_spec.filename)
        s1.command(f"echo Genotyping $(cat {local_repeat_spec} | wc -l) loci")
        s1.command("set -ex")
        s1.command(f"""/usr/bin/time --verbose GangSTR \
            --ref {local_fasta} \
            --bam {local_bam} \
            --regions {local_repeat_spec} \
            --out {output_prefix} |& tee {output_prefix}.log""")

        s1.command("ls -lhrt")
        s1.command(f"python3.9 -m str_analysis.convert_gangstr_vcf_to_expansion_hunter_json {output_prefix}.vcf")

        s1.command(f"gzip {output_prefix}.vcf")
        s1.command(f"gzip {output_prefix}.log")
        s1.output(f"{output_prefix}.vcf.gz", output_dir=os.path.join(output_dir, f"vcf"))
        s1.output(f"{output_prefix}.log.gz", output_dir=os.path.join(output_dir, f"log"))
        s1.output(f"{output_prefix}.json", output_dir=os.path.join(output_dir, f"json"))

        step1_output_json_paths.append(os.path.join(output_dir, f"json", f"{output_prefix}.json"))

    # step2: combine json files
    s2 = bp.new_step(name="Combine GangSTR outputs",
                     step_number=2,
                     image=DOCKER_IMAGE,
                     cpu=1,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=output_dir)
    for step1 in step1s:
        s2.depends_on(step1)

    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    for json_path in step1_output_json_paths:
        local_path = s2.input(json_path, localize_by=Localize.COPY)
        s2.command(f"ln -s {local_path}")

    output_prefix = f"combined.{positive_or_negative_loci}"
    s2.command("set -x")
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv --include-extra-gangstr-fields "
               f"--output-prefix {output_prefix}")
    s2.command(f"bgzip {output_prefix}.{len(step1_output_json_paths)}_json_files.bed")
    s2.command(f"tabix {output_prefix}.{len(step1_output_json_paths)}_json_files.bed.gz")
    s2.command("gzip *.tsv")
    s2.command("ls -lhrt")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.variants.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.alleles.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.bed.gz")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.bed.gz.tbi")
    bp.run()


if __name__ == "__main__":
    main()


