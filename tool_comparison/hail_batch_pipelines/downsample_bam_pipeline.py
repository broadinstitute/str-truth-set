import os
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_BAM_PATH = "gs://str-truth-set/hg38/CHM1_CHM13_2.bam"
CHM1_CHM13_BAI_PATH = "gs://str-truth-set/hg38/CHM1_CHM13_2.bam.bai"

CHM1_CHM13_BAM_COVERAGE = 40

DOCKER_IMAGE = "docker.io/weisburd/gatk:4.3.0.0"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai")
    parser.add_argument("--input-bam", default=CHM1_CHM13_BAM_PATH)
    parser.add_argument("--input-bai")
    parser.add_argument("--input-coverage", default=CHM1_CHM13_BAM_COVERAGE, type=float, help="Input bam coverage.")
    parser.add_argument("--target-coverage", default=30, type=float, help="Target coverage. Must be less than the input coverage.")
    parser.add_argument("--output-dir", default=os.path.dirname(CHM1_CHM13_BAM_PATH))
    args = bp.parse_known_args()

    pipeline_name = f"Downsample {os.path.basename(args.input_bam)} to {args.target_coverage}x"
    bp.set_name(pipeline_name)

    if args.input_coverage <= 1:
        parser.error("--input-coverage arg must be > 1")
    if args.target_coverage <= 1:
        parser.error("--target-coverage arg must be > 1")
    if args.target_coverage >= args.input_coverage:
        parser.error(f"--target-coverage arg must be < {args.input_coverage}")

    s1 = bp.new_step(pipeline_name, image=DOCKER_IMAGE, cpu=2, storage="250Gi", output_dir=args.output_dir)
    local_fasta = s1.input(args.reference_fasta, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    if args.reference_fasta_fai:
        s1.input(args.reference_fasta_fai, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    local_bam = s1.input(args.input_bam, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    if args.input_bai:
        s1.input(args.input_bai, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    bam_prefix = re.sub(".bam$", "", local_bam.filename)
    output_bam_filename = f"{bam_prefix}.downsampled_to_{int(args.target_coverage)}x.bam"
    fraction = args.target_coverage / args.input_coverage

    s1.command("set -ex")
    s1.command("cd /io/")
    s1.command(f"time gatk DownsampleSam I={local_bam} O={output_bam_filename} P={fraction:0.3f} CREATE_INDEX=true")
    s1.command(f"mv {output_bam_filename.replace('.bam', '.bai')} {output_bam_filename}.bai")  # rename the .bai file
    s1.command("ls -lh")

    s1.command(f"gatk CollectWgsMetrics STOP_AFTER={5*10**7} I={output_bam_filename} O=metrics.txt R={local_fasta}")
    s1.command(f"cat metrics.txt | head -n 8 | tail -n 2")

    s1.output(output_bam_filename)
    s1.output(f"{output_bam_filename}.bai")

    bp.run()


if __name__ == "__main__":
    main()


