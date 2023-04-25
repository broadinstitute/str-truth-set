import hail as hl
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/str-analysis@sha256:e10c982d1cdc1d8273469745dd35a70ca0fd5e160882063a93491d0ade15f812"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_BAM_PATH = "gs://bw2-delete-after-30-days/CHM1_CHM13_WGS2.downsampled_to_30x.bam"
CHM1_CHM13_BAI_PATH = "gs://bw2-delete-after-30-days/CHM1_CHM13_WGS2.downsampled_to_30x.bam.bai"

VARIANT_CATALOG_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/expansion_hunter/positive_loci.EHv5.*_of_293.json"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results_cost_optimized/expansion_hunter"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--use-illumina-expansion-hunter", action="store_true", help="Go back to using the Illumina "
         "version of ExpansionHunter instead of the optimized version from https://github.com/bw2/ExpansionHunter.git")
    parser.add_argument("--loci-to-exclude", action="append", help="Path of a file containing locus ids (one per line) "
        "to exclude from the variant catalogs before running ExpansionHunter. This can be useful when running with the "
        "--use-illumina-expansion-hunter in order to first filter out the small number of loci that cause "
        "ExpansionHunter to exit with the error message 'Flanks can contain at most 5 characters N but found x Ns.'")
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_BAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_BAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    variant_catalog_paths = VARIANT_CATALOG_POSITIVE_LOCI

    if args.use_illumina_expansion_hunter:
        tool_exec = "IlluminaExpansionHunter"
        cache_mates_arg = ""
        cpu_per_machine = 16
        catalogs_per_cpu = 10
        catalogs_per_machine = 2 * catalogs_per_cpu * cpu_per_machine
    else:
        tool_exec = "ExpansionHunter"
        #cache_mates_arg = "--cache-mates "
        cache_mates_arg = ""
        cpu_per_machine = 16
        catalogs_per_cpu = 10
        catalogs_per_machine = 2 * catalogs_per_cpu * cpu_per_machine

    output_dir = os.path.join(args.output_dir, tool_exec)

    loci_to_exclude = []
    if args.loci_to_exclude:
        for loci_to_exclude_file_path in args.loci_to_exclude:
            with open(loci_to_exclude_file_path, "rt") as f:
                for line in f:
                    loci_to_exclude.append(line.strip())
        print(f"Parsed {len(loci_to_exclude)} locus ids to exclude:", ", ".join(loci_to_exclude[:5]),
              "..." if len(loci_to_exclude) > 5 else "")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    if not args.force:
        existing_json_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.json"))
        logging.info(f"Precached {len(existing_json_paths)} json files")

    variant_catalog_file_stats_list = hl.hadoop_ls(variant_catalog_paths)
    if len(variant_catalog_file_stats_list) == 0:
        raise ValueError(f"No files found matching {variant_catalog_paths}")

    total_catalogs = len(variant_catalog_file_stats_list)
    if args.n:
        total_catalogs = min(total_catalogs, args.n * catalogs_per_machine)

    bp.set_name(f"TR Truth Set: {tool_exec} {cache_mates_arg}: total_cat={total_catalogs:,d}, "
                f"cat/mac={catalogs_per_machine}, cat/cpu={catalogs_per_cpu}, cpu/mac={cpu_per_machine}: "
                f"{bam_path_ending}")

    for batch_i in range(0, len(variant_catalog_file_stats_list), catalogs_per_machine):
        if args.n and batch_i >= args.n:
            break

        current_catalogs = variant_catalog_file_stats_list[batch_i:batch_i + catalogs_per_machine]
        s1 = bp.new_step(
            f"Run EHv5 catalogs #{batch_i}-{batch_i + catalogs_per_machine}", arg_suffix=f"eh", step_number=1,
            image=DOCKER_IMAGE, cpu=cpu_per_machine, storage="75Gi",
            localize_by=Localize.COPY,
            #localize_by=Localize.GSUTIL_COPY,
            delocalize_by=Delocalize.GSUTIL_COPY,
        )
        s1.gcloud_auth_activate_service_account()
        if args.reference_fasta_fai:
            s1.input(args.reference_fasta_fai)

        local_bam = s1.input(args.input_bam)
        if args.input_bai:
            s1.input(args.input_bai)
        local_fasta = s1.input(args.reference_fasta)

        s1.command("set -ex")
        local_variant_catalog_inputs = []
        for catalog_i, variant_catalog_file_stats in enumerate(current_catalogs):
            variant_catalog_path = variant_catalog_file_stats["path"]

            local_variant_catalog = s1.input(variant_catalog_path)
            local_variant_catalog_inputs.append(local_variant_catalog)

        if len(loci_to_exclude) > 0:
            for catalog_i, local_variant_catalog in enumerate(local_variant_catalog_inputs):
                local_variant_catalog_path = f"filtered_catalog.{local_variant_catalog.filename}"

                # add command to filter out excluded loci from variant catalog
                jq_command = f"cat {local_variant_catalog} | "
                jq_command += "jq '.[] | " + " | ".join([f'select(.LocusId != "{i}")' for i in loci_to_exclude]) + "' | "
                jq_command += "jq -s '.' "  # reformat output into proper json
                jq_command += f" > {local_variant_catalog_path}"
                s1.command(jq_command)
            catalog_glob = "filtered_catalog.*.json"
        else:
            for catalog_i, local_variant_catalog in enumerate(local_variant_catalog_inputs):
                local_variant_catalog_path = f"catalog.{local_variant_catalog.filename}"
                s1.command(f"ln -s {local_variant_catalog} {local_variant_catalog_path}")
            catalog_glob = "catalog.*.json"

        s1.command(f"echo Genotyping $(ls {catalog_glob} | wc -l) catalogs: {catalog_glob}")

        s1.command(f"""ls {catalog_glob} | parallel --lb --jobs {catalogs_per_cpu * cpu_per_machine} \
            "{tool_exec} {cache_mates_arg} --reference {local_fasta} --reads {local_bam} --variant-catalog {{}} --output-prefix expansion_hunter.{{/}}"
        """)
        s1.command("ls -lhrt")

        s1.output(f"expansion_hunter*.json", output_dir=os.path.join(output_dir, f"json"))

    bp.run()


if __name__ == "__main__":
    main()


