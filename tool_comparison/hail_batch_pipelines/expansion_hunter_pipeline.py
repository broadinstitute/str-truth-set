import hail as hl
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/expansion-hunter@sha256:5242b1d8cf477824898f2e634cef49fb8f1430eba93612b478a605fb8215c5d1"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

VARIANT_CATALOG_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/expansion_hunter/positive_loci.EHv5.*_of_308.json"
VARIANT_CATALOG_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/expansion_hunter/negative_loci.EHv5.*_of_305.json"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/expansion_hunter"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    parser.add_argument("--use-illumina-expansion-hunter", action="store_true", help="Go back to using the Illumina "
         "version of ExpansionHunter instead of the optimized version from https://github.com/bw2/ExpansionHunter.git")
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    if args.positive_loci:
        variant_catalog_paths = VARIANT_CATALOG_POSITIVE_LOCI
        positive_or_negative_loci = "positive_loci"
    elif args.negative_loci:
        variant_catalog_paths = VARIANT_CATALOG_NEGATIVE_LOCI
        positive_or_negative_loci = "negative_loci"
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    if args.use_illumina_expansion_hunter:
        tool_exec = "IlluminaExpansionHunter"
        output_dir = os.path.join(args.output_dir, tool_exec, positive_or_negative_loci)
        cache_mates_arg = ""
    else:
        tool_exec = "ExpansionHunter"
        output_dir = os.path.join(args.output_dir, positive_or_negative_loci)
        cache_mates_arg = "--cache-mates "


    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: {tool_exec}: {positive_or_negative_loci}: {bam_path_ending}")
    if not args.force:
        json_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.json"))
        logging.info(f"Precached {len(json_paths)} json files")

    step1s = []
    step1_output_json_paths = []
    for catalog_i, variant_catalog_file_stats in enumerate(hl.hadoop_ls(variant_catalog_paths)):
        variant_catalog_path = variant_catalog_file_stats["path"]

        if args.n and catalog_i >= args.n:
            break

        s1 = bp.new_step(
            f"Run EHv5 #{catalog_i}", arg_suffix=f"eh", step_number=1, image=DOCKER_IMAGE, cpu=2, storage="75Gi")
        step1s.append(s1)

        local_fasta = s1.input(args.reference_fasta, localize_by=Localize.COPY)
        if args.reference_fasta_fai:
            s1.input(args.reference_fasta_fai, localize_by=Localize.COPY)

        local_bam = s1.input(args.input_bam, localize_by=Localize.COPY)
        if args.input_bai:
            s1.input(args.input_bai, localize_by=Localize.COPY)

        local_variant_catalog = s1.input(variant_catalog_path)

        output_prefix = re.sub(".json$", "", local_variant_catalog.filename)
        s1.command(f"echo Genotyping $(cat {local_variant_catalog} | grep LocusId | wc -l) loci")
        s1.command("set -ex")

        s1.command(f"""/usr/bin/time --verbose {tool_exec} {cache_mates_arg} \
            --reference {local_fasta} \
            --reads {local_bam} \
            --variant-catalog {local_variant_catalog} \
            --output-prefix {output_prefix}""")
        s1.command("ls -lhrt")

        s1.output(f"{output_prefix}.json", output_dir=os.path.join(output_dir, f"json"))

        step1_output_json_paths.append(os.path.join(output_dir, f"json", f"{output_prefix}.json"))

    # step2: combine json files
    s2 = bp.new_step(name="Combine EHv5 outputs", step_number=2, image=DOCKER_IMAGE, storage="20Gi", cpu=1,
                     output_dir=output_dir)
    for step1 in step1s:
        s2.depends_on(step1)

    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    for json_path in step1_output_json_paths:
        local_path = s2.input(json_path, localize_by=Localize.COPY)
        s2.command(f"ln -s {local_path}")

    output_prefix = f"combined.{positive_or_negative_loci}"
    s2.command("set -x")
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv --include-extra-expansion-hunter-fields "
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


