import os
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_BAM_PATH = "gs://str-truth-set/hg38/CHM1_CHM13_2.bam"
CHM1_CHM13_BAI_PATH = "gs://str-truth-set/hg38/CHM1_CHM13_2.bam.bai"

DOCKER_IMAGE = "docker.io/weisburd/expansion-hunter-denovo:0.9"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/expansion_hunter_denovo"

MAX_REPEAT_UNIT_LENGTH = 50

def main():
    bp = pipeline("run ExpansionHunterDenovo", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai")
    parser.add_argument("--input-bam", default=CHM1_CHM13_BAM_PATH)
    parser.add_argument("--input-bai")
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    args = bp.parse_known_args()

    s1 = bp.new_step(f"STR Truth Set: ExpansionHunterDenovo", image=DOCKER_IMAGE, cpu=2, output_dir=args.output_dir)
    local_fasta = s1.input(args.reference_fasta, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    if args.reference_fasta_fai:
        s1.input(args.reference_fasta_fai, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    local_bam = s1.input(args.input_bam, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    if args.input_bai:
        s1.input(args.input_bai, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    s1.command("set -ex")

    output_prefix = re.sub(".bam$", "", local_bam.filename)

    s1.command(f"""time ExpansionHunterDenovo profile --reference {local_fasta} --reads {local_bam} \
        --max-unit-len {MAX_REPEAT_UNIT_LENGTH} --output-prefix {output_prefix} |& tee EHdn_command.log""")
    s1.command("ls -lh")

    s1.output(f"{output_prefix}.str_profile.json")
    s1.output(f"{output_prefix}.locus.tsv")
    s1.output(f"{output_prefix}.motif.tsv")
    s1.output("EHdn_command.log")

    bp.run()


if __name__ == "__main__":
    main()


