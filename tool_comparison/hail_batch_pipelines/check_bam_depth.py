import os
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

DOCKER_IMAGE = "weisburd/filter-vcfs@sha256:9ea4e1e648e15bb370efbe59f096d54791c9b6dcb75f1ec3e957a3c897014bac"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai")
    parser.add_argument("--output-dir", required=True, help="Google storage directory for output file")
    parser.add_argument("input_bam", nargs="+")
    args = bp.parse_known_args()


    for input_bam in args.input_bam:
        bam_or_cram_prefix = re.sub("(.bam|.cram)$", "", os.path.basename(input_bam))
        s1 = bp.new_step(f"Coverage: {bam_or_cram_prefix}", image=DOCKER_IMAGE, cpu=1, memory="standard", storage="20Gi", output_dir=args.output_dir)
        local_fasta = s1.input(args.reference_fasta, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        if args.reference_fasta_fai:
            s1.input(args.reference_fasta_fai, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

        local_bam = s1.input(input_bam, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)

        s1.command("set -ex")
        s1.command("wget https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth -O /usr/local/bin/mosdepth")
        s1.command("chmod 777 /usr/local/bin/mosdepth")

        s1.command(f"mosdepth -f {local_fasta} -x {bam_or_cram_prefix}.coverage {local_bam}")
        #s1.output(f"{bam_or_cram_prefix}.coverage.mosdepth.summary.txt")

        s1.command(f"cat {bam_or_cram_prefix}.coverage.mosdepth.summary.txt | cut -f 4 | tail -n +2 | head -n 23")
        s1.command(f"grep total {bam_or_cram_prefix}.coverage.mosdepth.summary.txt > {bam_or_cram_prefix}.total_coverage.txt")

        s1.output(f"{bam_or_cram_prefix}.total_coverage.txt")

    bp.run()


if __name__ == "__main__":
    main()


