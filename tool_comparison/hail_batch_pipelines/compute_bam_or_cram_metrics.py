import os
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"


CHM1_CHM13_WGS_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_WES_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WES/CHMI_CHMI3_Nex1.cram"

DOCKER_IMAGE = "docker.io/weisburd/gatk:4.3.0.0"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    parser = bp.get_config_arg_parser()
    parser.add_argument("-R", "--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("-o", "--output-dir", help="Optional output gs:// storage location")
    parser.add_argument("input_bam_or_cram", nargs="*")
    args = bp.parse_known_args()

    if not args.input_bam_or_cram:
        args.input_bam_or_cram = [CHM1_CHM13_WGS_CRAM_PATH, CHM1_CHM13_WES_CRAM_PATH]

    pipeline_name = f"BAM or CRAM metrics: "
    if len(args.input_bam_or_cram) == 1:
        pipeline_name += os.path.basename(args.input_bam_or_cram[0])
    else:
        pipeline_name += f"({len(args.input_bam_or_cram)} files)"
    bp.set_name(pipeline_name)

    for input_bam_or_cram in args.input_bam_or_cram:
        s1 = bp.new_step(pipeline_name, image=DOCKER_IMAGE, cpu=2, output_dir=args.output_dir)  #, storage="250Gi"
        local_fasta = s1.input(args.reference_fasta, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        local_bam_or_cram = s1.input(input_bam_or_cram, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

        output_prefix = re.sub("(.bam|.cram)$", "", local_bam_or_cram.filename)

        s1.command("set -ex")
        s1.command("cd /io/")
        s1.command(f"gatk CollectInsertSizeMetrics "
                   f"STOP_AFTER={5*10**7} "
                   f"I={local_bam_or_cram} "
                   f"O={output_prefix}.insert_size_metrics.txt "
                   f"H={output_prefix}.insert_size_metrics.pdf "
                   f"R={local_fasta}")
        s1.command(f"cat {output_prefix}.insert_size_metrics.txt")

        s1.command(f"gatk CollectWgsMetrics STOP_AFTER={5*10**7} I={local_bam_or_cram} O={output_prefix}.collect_wgs_metrics.txt R={local_fasta}")
        s1.command(f"cat {output_prefix}.collect_wgs_metrics.txt | head -n 8 | tail -n 2")

        if args.output_dir:
            s1.output(f"{output_prefix}.insert_size_metrics.txt")
            s1.output(f"{output_prefix}.insert_size_metrics.pdf")
            s1.output(f"{output_prefix}.collect_wgs_metrics.txt")

    bp.run()


if __name__ == "__main__":
    main()


