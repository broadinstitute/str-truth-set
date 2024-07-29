import hailtop.fs as hfs
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/hipstr@sha256:3c18925e50b35f7bd69362e22f4f1f5fb8c2c8fc79c016120f44459b7eef53c3"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

REGIONS_BED_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/hipstr/positive_loci.*_of_015.bed"
REGIONS_BED_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/hipstr/negative_loci.*_of_015.bed"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/hipstr"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    parser_group.add_argument("--regions-bed", action="append", help="Path of regions bed file(s) to process")

    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    if args.positive_loci:
        positive_or_negative_loci = "positive_loci"
        regions_bed_file_paths = [x.path for x in hfs.ls(REGIONS_BED_POSITIVE_LOCI)]
        if len(regions_bed_file_paths) == 0:
            raise ValueError(f"No files found matching {REGIONS_BED_POSITIVE_LOCI}")
    elif args.negative_loci:
        positive_or_negative_loci = "negative_loci"
        regions_bed_file_paths = [x.path for x in hfs.ls(REGIONS_BED_NEGATIVE_LOCI)]
        if len(regions_bed_file_paths) == 0:
            raise ValueError(f"No files found matching {REGIONS_BED_NEGATIVE_LOCI}")
    elif args.regions_bed:
        positive_or_negative_loci = os.path.basename(args.regions_bed[0]).replace(".bed", "").replace(".gz", "")
        regions_bed_file_paths = [x.path for x in hfs.ls(args.regions_bed)]
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    if args.n:
        regions_bed_file_paths = regions_bed_file_paths[:args.n]

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: HipSTR: {positive_or_negative_loci}: {bam_path_ending}")
    output_dir = os.path.join(args.output_dir, positive_or_negative_loci)
    if not args.force:
        json_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.json"))
        logging.info(f"Precached {len(json_paths)} json files")

    create_hipstr_steps(
        bp,
        reference_fasta=args.reference_fasta,
        input_bam=args.input_bam,
        input_bai=args.input_bai,
        regions_bed_file_paths=regions_bed_file_paths,
        output_dir=output_dir,
        output_prefix=f"combined.{positive_or_negative_loci}",
        reference_fasta_fai=args.reference_fasta_fai,
        male_or_female="female")
    bp.run()

def create_hipstr_steps(bp, *, reference_fasta, input_bam, input_bai, regions_bed_file_paths, output_dir, output_prefix, reference_fasta_fai=None,
                        male_or_female="female"):
    step1s = []
    step1_output_paths = []

    hfs_ls_results = hfs.ls(input_bam)
    if len(hfs_ls_results) == 0:
        raise ValueError(f"No files found matching {input_bam}")
    input_bam_file_stats = hfs_ls_results[0]

    for repeat_spec_i, regions_bed_file_path in enumerate(regions_bed_file_paths):
        s1 = bp.new_step(f"Run HipSTR #{repeat_spec_i}",
                         arg_suffix=f"hipstr",
                         step_number=1,
                         image=DOCKER_IMAGE,
                         cpu=1,
                         localize_by=Localize.COPY,
                         storage=f"{int(input_bam_file_stats.size/10**9) + 25}Gi",
                         output_dir=output_dir)
        step1s.append(s1)

        local_fasta = s1.input(reference_fasta)
        if reference_fasta_fai:
            s1.input(reference_fasta_fai)

        local_bam = s1.input(input_bam)
        if input_bai:
            s1.input(input_bai)

        local_regions_bed = s1.input(regions_bed_file_path)
        input_bam_filename_prefix = re.sub("(.bam|.cram)$", "", os.path.basename(local_bam.filename))
        output_prefix = re.sub(".bed(.gz)?$", "", local_regions_bed.filename)
        s1.command(f"echo Genotyping $(cat {local_regions_bed} | wc -l) loci")
        s1.command("set -ex")
        s1.command(f"""/usr/bin/time --verbose HipSTR \
                --bams {local_bam} \
                --bam-samps {input_bam_filename_prefix} \
                --fasta {local_fasta} \
                --regions {local_regions_bed} \
                --def-stutter-model \
                --min-reads 5 \
                --max-str-len 1000000 \
                --log {output_prefix}.log \
                --str-vcf {output_prefix}.vcf.gz \
                --viz-out {output_prefix}.viz.gz""")

        s1.command("ls -lhrt")

        s1.command(f"python3.9 -m str_analysis.convert_hipstr_vcf_to_expansion_hunter_json {output_prefix}.vcf.gz")
        s1.command(f"gzip {output_prefix}.log")

        s1.output(f"{output_prefix}.vcf.gz", output_dir=os.path.join(output_dir, f"vcf"))
        s1.output(f"{output_prefix}.log.gz", output_dir=os.path.join(output_dir, f"log"))
        s1.output(f"{output_prefix}.viz.gz", output_dir=os.path.join(output_dir, f"viz"))
        s1.output(f"{output_prefix}.json", output_dir=os.path.join(output_dir, f"json"))

        step1_output_paths.append(os.path.join(output_dir, f"json", f"{output_prefix}.json"))

    # step2: combine json files
    s2 = bp.new_step(name=f"Combine HipSTR outputs for {os.path.basename(input_bam)}",
                     step_number=2,
                     image=DOCKER_IMAGE,
                     cpu=2,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=output_dir)
    for step1 in step1s:
        s2.depends_on(step1)

    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    for json_path in step1_output_paths:
        local_path = s2.input(json_path, localize_by=Localize.COPY)
        s2.command(f"ln -s {local_path}")

    s2.command("set -x")
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv --include-extra-hipstr-fields "
               f"--output-prefix {output_prefix}")

    s2.command(f"bgzip {output_prefix}.{len(step1_output_paths)}_json_files.bed")
    s2.command(f"tabix {output_prefix}.{len(step1_output_paths)}_json_files.bed.gz")
    s2.command("ls -lhrt")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.variants.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.alleles.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.bed.gz")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.bed.gz.tbi")


if __name__ == "__main__":
    main()


