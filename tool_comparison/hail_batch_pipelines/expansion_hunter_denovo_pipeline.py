import hail as hl
import os
import re
from step_pipeline import pipeline, Backend, Localize, Delocalize

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_BAM_PATH = "gs://str-truth-set/hg38/CHM1_CHM13_2.bam"
CHM1_CHM13_BAI_PATH = "gs://str-truth-set/hg38/CHM1_CHM13_2.bam.bai"

CHM1_CHM13_CONFIDENCE_REGIONS = "gs://str-truth-set/hg38/ref/other/full.38.bed.gz"

VARIANT_CATALOG_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/expansion_hunter/positive_loci.EHv5.*_of_289.json"
VARIANT_CATALOG_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/expansion_hunter/negative_loci.EHv5.*_of_289.json"

DOCKER_IMAGE = "docker.io/weisburd/expansion-hunter-denovo:0.9"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/expansion_hunter_denovo"


def main():
    bp = pipeline("run ExpansionHunterDenovo", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    s1 = bp.new_step(f"STR Truth Set: ExpansionHunter", image=DOCKER_IMAGE, cpu=2, output_dir=OUTPUT_BASE_DIR)
    local_fasta, _ = s1.inputs(REFERENCE_FASTA_PATH, REFERENCE_FASTA_FAI_PATH, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    local_bam, local_bai = s1.inputs(CHM1_CHM13_BAM_PATH, CHM1_CHM13_BAI_PATH, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    s1.command("set -ex")

    output_prefix = re.sub(".bam$", "", local_bam.filename)

    s1.command(f"""time ExpansionHunterDenovo profile --reference {local_fasta} --reads {local_bam} \
        --output-prefix {output_prefix} |& tee EHdn_command.log""")
    s1.command("ls -lh")
    s1.output(f"{output_prefix}.str_profile.json")
    s1.output(f"{output_prefix}.locus.tsv")
    s1.output(f"{output_prefix}.motif.tsv")
    s1.output("EHdn_command.log")

    bp.run()

    # download results
    files_to_download_when_done = []
    for remote_path, destination_dir in files_to_download_when_done:
        if not hl.hadoop_exists(remote_path):
            print(f"Output path doesn't exist: {remote_path}. Skipping download..")
            continue
        if not os.path.isdir(destination_dir):
            print(f"Creating local directory: {destination_dir}")
            os.mkdir(destination_dir)
        print(f"Downloading {remote_path} to {destination_dir}/")
        hl.hadoop_copy(remote_path, os.path.join(destination_dir, os.path.basename(remote_path)))


if __name__ == "__main__":
    main()


