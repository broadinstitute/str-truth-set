import hail as hl
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/expansion-hunter:v5"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_BAM_PATH = "gs://str-truth-set/hg38/CHM1_CHM13_2.bam"
CHM1_CHM13_BAI_PATH = "gs://str-truth-set/hg38/CHM1_CHM13_2.bam.bai"

VARIANT_CATALOG_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/expansion_hunter/positive_loci.EHv5.*_of_308.json"
VARIANT_CATALOG_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/expansion_hunter/negative_loci.EHv5.*_of_305.json"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/expansion_hunter"


def main():
    bp = pipeline("STR Truth Set: ExpansionHunter", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    args = bp.parse_known_args()

    if args.positive_loci:
        variant_catalog_paths = VARIANT_CATALOG_POSITIVE_LOCI
        positive_or_negative_loci = "positive_loci"
    elif args.negative_loci:
        variant_catalog_paths = VARIANT_CATALOG_NEGATIVE_LOCI
        positive_or_negative_loci = "negative_loci"
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    output_dir = os.path.join(OUTPUT_BASE_DIR, positive_or_negative_loci)
    if not args.force:
        log_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.log.gz"))
        json_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.json"))
        logging.info(f"Precached {len(json_paths)} json files and {len(log_paths)} log files")

    step1s = []
    step1_output_json_paths = []
    for catalog_i, variant_catalog_file_stats in enumerate(hl.hadoop_ls(variant_catalog_paths)):
        variant_catalog_path = variant_catalog_file_stats["path"]

        s1 = bp.new_step(f"Run EHv5 #{catalog_i}", arg_suffix=f"eh", step_number=1, image=DOCKER_IMAGE, cpu=1)
        step1s.append(s1)

        local_fasta, _ = s1.inputs(REFERENCE_FASTA_PATH, REFERENCE_FASTA_FAI_PATH, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        local_bam, local_bai = s1.inputs(CHM1_CHM13_BAM_PATH, CHM1_CHM13_BAI_PATH, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        local_variant_catalog = s1.input(variant_catalog_path)

        output_prefix = re.sub(".json$", "", local_variant_catalog.filename)
        s1.command("set -ex")

        s1.command(f"""time ExpansionHunter \
                --reference {local_fasta} \
                --reads {local_bam} \
                --variant-catalog {local_variant_catalog} \
                --cache-mates \
                --output-prefix {output_prefix} |& tee {output_prefix}.log""")

        s1.command(f"gzip {output_prefix}.log")
        s1.command("ls -lhrt")

        s1.output(f"{output_prefix}.json", output_dir=os.path.join(output_dir, f"json"))
        s1.output(f"{output_prefix}.log.gz", output_dir=os.path.join(output_dir, f"log"))

        step1_output_json_paths.append(os.path.join(output_dir, f"json", f"{output_prefix}.json"))

    # step2: combine json files
    s2 = bp.new_step(name="Combine EHv5 outputs", step_number=2, image=DOCKER_IMAGE, storage="20Gi", cpu=1)
    for step1 in step1s:
        s2.depends_on(step1)

    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    for json_path in step1_output_json_paths:
        local_path = s2.input(json_path)
        s2.command(f"ln -s {local_path}")

    output_prefix = f"combined.{positive_or_negative_loci}"
    s2.command(f"python3 -m str_analysis.combine_str_json_to_tsv "
               f"--include-extra-expansion-hunter-fields --output-prefix {output_prefix}")
    s2.command("ls -lhrt")

    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.variants.tsv", output_dir=output_dir)
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.alleles.tsv", output_dir=output_dir)

    bp.run()


if __name__ == "__main__":
    main()


