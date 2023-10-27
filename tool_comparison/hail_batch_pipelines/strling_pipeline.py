"""
strling version: 0.5.2
strling call

Usage:
  strling call [options] bam bin

Arguments:
  bam              path to bam file
  bin              bin file previously created by `strling extract`

Options:
  -f, --fasta=FASTA          path to fasta file
  -m, --min-support=MIN_SUPPORT
                             minimum number of supporting reads for a locus to be reported (default: 5)
  -c, --min-clip=MIN_CLIP    minimum number of supporting clipped reads for each side of a locus (default: 0)
  -t, --min-clip-total=MIN_CLIP_TOTAL
                             minimum total number of supporting clipped reads for a locus (default: 0)
  -q, --min-mapq=MIN_MAPQ    minimum mapping quality (does not apply to STR reads) (default: 40)
  -l, --loci=LOCI            Annoated bed file specifying additional STR loci to genotype. Format is: chr start stop repeatunit [name]
  -b, --bounds=BOUNDS        STRling -bounds.txt file (usually produced by strling merge) specifying additional STR loci to genotype.
  -o, --output-prefix=OUTPUT_PREFIX
                             prefix for output files (default: strling)
  -v, --verbose
  -h, --help                 Show this help
"""

import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/strling@sha256:e867138fadfa08969c7c0e043fd80999aa60b0592d933f183f3ffb1c5887bfc1"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/strling"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: strling: {bam_path_ending}")
    if not args.force:
        json_paths = bp.precache_file_paths(os.path.join(args.output_dir, f"**/*.json"))
        logging.info(f"Precached {len(json_paths)} json files")

    ### step1: index fasta
    strling_genome_index_path = os.path.join(args.output_dir, "Homo_sapiens_assembly38.fasta.str")
    s1 = bp.new_step(f"index reference", arg_suffix="index-reference",
                     image=DOCKER_IMAGE, cpu=1, memory="lowmem", storage="10Gi",
                     localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
                     delocalize_by=Delocalize.COPY,
                     output_dir=args.output_dir)
    reference_fasta_input, _ = s1.inputs(args.reference_fasta, f"{args.reference_fasta}.fai")
    s1.command(f"""/usr/bin/time -v strling index {reference_fasta_input}""")
    s1.output(os.path.basename(strling_genome_index_path))

    sample_id = re.sub(".bam$", "", os.path.basename(args.input_bam))

    ### step2: extract
    s2 = bp.new_step(f"extract: {sample_id}", arg_suffix="extract",
                     image=DOCKER_IMAGE, cpu=0.5, memory="standard", storage="10Gi",
                     localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
                     delocalize_by=Delocalize.COPY,
                     depends_on=s1,
                     output_dir=f"{args.output_dir}")
    reference_fasta_input, _ = s2.inputs(args.reference_fasta, f"{args.reference_fasta}.fai")
    bam_input, _ = s2.inputs(args.input_bam, f"{args.input_bam}.bai")
    strling_genome_index = s2.input(strling_genome_index_path)

    s2.command("set -euxo pipefail")
    s2.command("/usr/bin/time -v strling extract "
               f"-f {reference_fasta_input} "
               f"-g {strling_genome_index} "
               f"{bam_input} "
               f"{sample_id}.bin ")

    s2.command(f"/usr/bin/time -v strling call "
               f"-f {reference_fasta_input} "
               f"-o {sample_id}.single-sample "
               f"{bam_input} "
               f"{sample_id}.bin")

    s2.command(f"python3 -m str_analysis.convert_strling_calls_to_expansion_hunter_json "
               f"{sample_id}.single-sample-genotype.txt")

    #s2.output(f"{sample_id}.bin")
    #s2.output(f"{sample_id}.single-sample-genotype.txt")
    s2.output(f"{sample_id}.single-sample-genotype.json")

    bp.run()


if __name__ == "__main__":
    main()


