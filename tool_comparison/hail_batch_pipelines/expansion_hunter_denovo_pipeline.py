import os
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.bai"

EXPANSION_HUNTER_DENOVO_DOCKER_IMAGE = "weisburd/expansion-hunter-denovo@sha256:2b68da96ed3bff60f1d47ec971776083cc5e5a07c673ceff5741a9bc6b8bc0f6"
EXPANSION_HUNTER_DOCKER_IMAGE = "weisburd/expansion-hunter@sha256:a2b3dc962c33733cb0e293e72bbba10449f0eadeaf665adec71debc4b1b7255d"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/expansion_hunter_denovo"

MAX_REPEAT_UNIT_LENGTH = 50

def main():
    bp = pipeline("run ExpansionHunterDenovo", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai")
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai")
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    args = bp.parse_known_args()

    s1 = bp.new_step(image=EXPANSION_HUNTER_DENOVO_DOCKER_IMAGE, step_number=1, cpu=2, output_dir=args.output_dir)
    local_fasta = s1.input(args.reference_fasta, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    if args.reference_fasta_fai:
        s1.input(args.reference_fasta_fai, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    local_bam = s1.input(args.input_bam, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    if args.input_bai:
        s1.input(args.input_bai, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    output_prefix = re.sub("(.bam|.cram)$", "", local_bam.filename)
    s1.name = f"STR Truth Set: ExpansionHunterDenovo: {output_prefix}"

    s1.command("set -ex")
    s1.command(f"""time ExpansionHunterDenovo profile --reference {local_fasta} --reads {local_bam} \
        --max-unit-len {MAX_REPEAT_UNIT_LENGTH} --output-prefix {output_prefix} |& tee EHdn_command.log""")
    s1.command("ls -lh")

    s1.output(f"{output_prefix}.str_profile.json")
    s1.output(f"{output_prefix}.locus.tsv")
    s1.output(f"{output_prefix}.motif.tsv")
    s1.output("EHdn_command.log")

    s2 = bp.new_step(
        "Convert to bed",
        step_number=2,
        image=EXPANSION_HUNTER_DOCKER_IMAGE,
        cpu=2,
        output_dir=args.output_dir,
        depends_on=s1,
    )
    s2.command("set -ex")

    local_input_tsv = s2.input(os.path.join(args.output_dir, f"{output_prefix}.locus.tsv"))
    s2.command(f"python3 -m str_analysis.convert_expansion_hunter_denovo_locus_tsv_to_bed "
               f"{local_input_tsv} "
               f"{output_prefix}.expansion_hunter_denovo.bed")

    s2.command("ls -lh")
    s2.command(f"bgzip {output_prefix}.expansion_hunter_denovo.bed")
    s2.command(f"tabix {output_prefix}.expansion_hunter_denovo.bed.gz")
    s2.output(f"{output_prefix}.expansion_hunter_denovo.bed.gz")
    s2.output(f"{output_prefix}.expansion_hunter_denovo.bed.gz.tbi")

    bp.run()


if __name__ == "__main__":
    main()


