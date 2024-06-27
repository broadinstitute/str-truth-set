import os
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

CHM1_CHM13_CRAM_COVERAGE = 40

DOCKER_IMAGE = "weisburd/gatk@sha256:b9b39f21b51ec9ef17937e6581dc3624a0a232d464616eddffb7678171fba578"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai")
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai")
    #parser.add_argument("--input-coverage", required=True, type=float, help="Input bam coverage.")
    parser.add_argument("--target-coverage", default=30, type=float, help="Target coverage. Must be less than the input coverage.")
    parser.add_argument("--output-dir", required=True, help="Google storage directory for output file")
    args = bp.parse_known_args()

    pipeline_name = f"Downsample {os.path.basename(args.input_bam)} to {args.target_coverage}x"
    bp.set_name(pipeline_name)

    if args.target_coverage <= 1:
        parser.error("--target-coverage arg must be > 1")

    s1 = bp.new_step(pipeline_name, image=DOCKER_IMAGE, cpu=2, memory="highmem", storage="250Gi", output_dir=args.output_dir)
    local_fasta = s1.input(args.reference_fasta, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    if args.reference_fasta_fai:
        s1.input(args.reference_fasta_fai, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    local_bam = s1.input(args.input_bam, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    if args.input_bai:
        s1.input(args.input_bai, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

    bam_or_cram_prefix = re.sub("(.bam|.cram)$", "", local_bam.filename)
    output_bam_filename = f"{bam_or_cram_prefix}.downsampled_to_{int(args.target_coverage)}x.bam"

    s1.command("set -ex")
    s1.command("cd /io/")
    s1.command("wget https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth -O /usr/local/bin/mosdepth")
    s1.command("chmod 777 /usr/local/bin/mosdepth")

    s1.command(f"mosdepth -f {local_fasta} -x coverage {local_bam}")

    s1.command(f"time gatk --java-options '-Xmx11G' DownsampleSam "
               f"--REFERENCE_SEQUENCE {local_fasta} "
               f"-I {local_bam} "
               f"-O {output_bam_filename} "
               f"""-P $(echo "{args.target_coverage} / $(grep total coverage.mosdepth.summary.txt | cut -f 4)" | bc -l | awk '{{printf "%.4f", $0}}') """
               f"--CREATE_INDEX true")

    s1.command(f"samtools calmd -b {output_bam_filename} {local_fasta} > {output_bam_filename}.with_NM_tag.bam")
    s1.command(f"mv {output_bam_filename}.with_NM_tag.bam {output_bam_filename}")
    s1.command(f"samtools index {output_bam_filename}")

    s1.command("ls -lh")

    s1.command(f"mosdepth -f {local_fasta} -x coverage_after_downsampling {output_bam_filename}")
    s1.command(f"cat coverage_after_downsampling.mosdepth.summary.txt | cut -f 4 | tail -n +2 | head -n 23")
    s1.command(f"grep total coverage_after_downsampling.mosdepth.summary.txt > coverage_after_downsampling.total_coverage.txt")
    s1.command(f"cat coverage_after_downsampling.total_coverage.txt")

    #s1.command(f"gatk CollectWgsMetrics STOP_AFTER={5*10**7} I={output_bam_filename} O=metrics.txt R={local_fasta}")
    #s1.command(f"cat metrics.txt | head -n 8 | tail -n 2")

    s1.output(output_bam_filename)
    s1.output(f"{output_bam_filename}.bai")

    bp.run()


if __name__ == "__main__":
    main()


