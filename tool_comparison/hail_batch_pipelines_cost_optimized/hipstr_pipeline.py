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

#REPEAT_SPECS_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/hipstr/positive_loci.*_of_015.bed"
REPEAT_SPECS_POSITIVE_LOCI = "gs://str-truth-set/hg38/ref/other/variant_catalogs_for_at_least_9bp/hipstr/positive_loci.HipSTR.*_of_281.bed"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results_cost_optimized/hipstr"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_BAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_BAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    localize_by=Localize.COPY
    #localize_by=Localize.GSUTIL_COPY
    #localize_by = Localize.HAIL_BATCH_CLOUDFUSE
    memory = "highmem"
    cpu_per_machine = 4
    repeat_specs_per_cpu = 8
    repeat_specs_per_machine = 4 * repeat_specs_per_cpu * cpu_per_machine

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    if not args.force:
        json_paths = bp.precache_file_paths(os.path.join(args.output_dir, f"**/*.json"))
        logging.info(f"Precached {len(json_paths)} json files")

    repeat_spec_file_stats_list = hl.hadoop_ls(REPEAT_SPECS_POSITIVE_LOCI)
    if len(repeat_spec_file_stats_list) == 0:
        raise ValueError(f"No files found matching {REPEAT_SPECS_POSITIVE_LOCI}")
    
    total_repeat_specs = len(repeat_spec_file_stats_list)
    if args.n:
        total_repeat_specs = min(total_repeat_specs, args.n * repeat_specs_per_machine)

    bp.set_name(f"TR Truth Set: HipSTR: total_cat={total_repeat_specs:,d}, "
                f"cat/mac={repeat_specs_per_machine}, cat/cpu={repeat_specs_per_cpu}, cpu/mac={cpu_per_machine}, mem={memory}: "
                f"{bam_path_ending}")

    for batch_i in range(0, len(repeat_spec_file_stats_list), repeat_specs_per_machine):
        if args.n and batch_i >= args.n:
            break

        current_repeat_specs = repeat_spec_file_stats_list[batch_i:batch_i + repeat_specs_per_machine]

        s1 = bp.new_step(
            f"Run HipSTR #{batch_i}-{batch_i + repeat_specs_per_machine}", arg_suffix=f"hipstr", step_number=1,
            image=DOCKER_IMAGE, cpu=cpu_per_machine, memory=memory, storage="75Gi",
            localize_by=localize_by,
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
        for repeat_spec_i, repeat_spec_file_stats in enumerate(current_repeat_specs):
            local_repeat_spec = s1.input(repeat_spec_file_stats["path"])
            s1.command(f"ln -s {local_repeat_spec} repeat_spec.{local_repeat_spec.filename}")
        repeat_spec_glob = "repeat_spec.*.bed"

        s1.command(f"echo Genotyping $(ls {repeat_spec_glob} | wc -l) loci: {repeat_spec_glob}")
        s1.command(f"""ls {repeat_spec_glob} | parallel --halt-on-error 2 --verbose --lb --jobs {repeat_specs_per_cpu * cpu_per_machine} \
            "HipSTR --min-reads 5 --max-str-len 1000000 --def-stutter-model --fasta {local_fasta} --bams {local_bam} --regions {{}} --log {{/}}.log --str-vcf hipstr.{{/}}.vcf.gz && python3 -m str_analysis.convert_hipstr_vcf_to_expansion_hunter_json hipstr.{{/}}.vcf.gz"
        """)
        s1.command("ls -lhrt")

        s1.output(f"hipstr*.json", output_dir=os.path.join(args.output_dir, f"json"))

    bp.run()


if __name__ == "__main__":
    main()


