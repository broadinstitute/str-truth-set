"""
BAMLIST=$1
REFERENCE=$2
MARKERS_PER_JOB=$3
"""


import hail as hl
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/popstr@sha256:d05ac1a56baba5d311564781a722bb74497dca0fa27f2fcdc00ba54838b510ba"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

MARKER_INFO_POSITIVE_LOCI = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/positive_loci.popSTR.chr*.markerInfo.gz"
MARKER_INFO_NEGATIVE_LOCI = "gs://str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers/CHM1_CHM13/negative_loci.popSTR.chr*.markerInfo.gz"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/popstr"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    parser_group.add_argument("--marker-info-path-prefix", help="Path of marker info file(s) to process")

    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    args = bp.parse_known_args()

    if args.positive_loci:
        positive_or_negative_loci = "positive_loci"
        output_dir = os.path.join(args.output_dir, positive_or_negative_loci)
        marker_info_file_stats_list = hl.hadoop_ls(MARKER_INFO_POSITIVE_LOCI)
        if len(marker_info_file_stats_list) == 0:
            raise ValueError(f"No files found matching {MARKER_INFO_POSITIVE_LOCI}")
    elif args.negative_loci:
        positive_or_negative_loci = "negative_loci"
        output_dir = os.path.join(args.output_dir, positive_or_negative_loci)
        marker_info_file_stats_list = hl.hadoop_ls(MARKER_INFO_NEGATIVE_LOCI)
        if len(marker_info_file_stats_list) == 0:
            raise ValueError(f"No files found matching {MARKER_INFO_NEGATIVE_LOCI}")
    elif args.marker_info_path_prefix:
        positive_or_negative_loci = os.path.basename(args.marker_info_path_prefix)
        output_dir = args.output_dir
        marker_info_file_stats_list = hl.hadoop_ls(f"{args.marker_info_path_prefix}*.markerInfo.gz")
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: popSTR: {positive_or_negative_loci}: {bam_path_ending}")
    if not args.force:
        json_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.json"))
        logging.info(f"Precached {len(json_paths)} json files")

    step1s = []

    n_cpus = 8
    s1 = bp.new_step(
        f"Run popSTR",
        arg_suffix=f"popSTR",
        step_number=1,
        image=DOCKER_IMAGE,
        cpu=n_cpus,
        memory="highmem",
        storage="100Gi",
        output_dir=output_dir
    )
    step1s.append(s1)

    local_fasta = s1.input(args.reference_fasta, localize_by=Localize.COPY)
    if args.reference_fasta_fai:
        s1.input(args.reference_fasta_fai, localize_by=Localize.COPY)

    local_bam = s1.input(args.input_bam, localize_by=Localize.COPY)
    if args.input_bai:
        s1.input(args.input_bai, localize_by=Localize.COPY)

    local_marker_info_paths = s1.inputs([
        marker_info_file_stats["path"] for marker_info_file_stats in marker_info_file_stats_list])

    sample_id = os.path.basename(args.input_bam).replace(".cram", "").replace(".bam", "")

    s1.command("set -ex")
    s1.command("cd /popSTR-2.0")
    s1.command(f"echo $'{sample_id}\\t{local_bam}' > bamList.txt")
    s1.command(f"cat bamList.txt")
    s1.command(f"ls -lh {local_bam}")
    s1.command(f"ls -lh {local_bam}.bai")
    s1.command("mkdir markerInfo")
    s1.command("cp defaultModel markerInfo/")
    for local_marker_info_path in local_marker_info_paths:
        match = re.search("(chr.{1,2}.markerInfo).gz", str(local_marker_info_path))
        if not match:
            print("WARNING: unexpected marker info path: " + str(local_marker_info_path))
        symlink_filename = match.group(1).replace(".", "")
        #s1.command(f"ln -s {local_marker_info_path} {symlink_filename}.gz")
        s1.command(f"gunzip -c {local_marker_info_path} > markerInfo/{symlink_filename}")
    s1.command("ls -lh .")
    s1.command("echo --------")
    s1.command("ls -lh markerInfo/")
    s1.command(f"./run.sh bamList.txt {local_fasta} {n_cpus}")

    s1.command("ls -lhrt")
    s1.command("ls -lhrt vcfs")
    s1.command(f"tar czf {sample_id}.vcfs.tar.gz vcfs")
    s1.output(f"{sample_id}.vcfs.tar.gz")
    #s1.command(f"python3.9 -m str_analysis.convert_popstr_vcf_to_expansion_hunter_json {output_prefix}.vcf.gz")
    #s1.command(f"gzip {output_prefix}.log")
    #s1.output(f"{output_prefix}.vcf.gz", output_dir=os.path.join(output_dir, f"vcf"))
    #s1.output(f"{output_prefix}.log.gz", output_dir=os.path.join(output_dir, f"log"))
    #s1.output(f"{output_prefix}.viz.gz", output_dir=os.path.join(output_dir, f"viz"))
    #s1.output(f"{output_prefix}.json", output_dir=os.path.join(output_dir, f"json"))

    #step1_output_json_paths.append(os.path.join(output_dir, f"json", f"{output_prefix}.json"))

    bp.run()

    return

    # step2: combine json files
    s2 = bp.new_step(name="Combine popSTR outputs",
                     step_number=2,
                     image=DOCKER_IMAGE,
                     cpu=2,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=output_dir)

    for step1 in step1s:
        s2.depends_on(step1)

    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    for json_path in step1_output_json_paths:
        local_path = s2.input(json_path, localize_by=Localize.COPY)
        s2.command(f"ln -s {local_path}")

    output_prefix = f"combined.{positive_or_negative_loci}"
    s2.command("set -x")
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv --include-extra-popSTR-fields "
               f"--output-prefix {output_prefix}")
    s2.command(f"bgzip {output_prefix}.{len(step1_output_json_paths)}_json_files.bed")
    s2.command(f"tabix {output_prefix}.{len(step1_output_json_paths)}_json_files.bed.gz")
    s2.command("ls -lhrt")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.variants.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.alleles.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.bed.gz")
    s2.output(f"{output_prefix}.{len(step1_output_json_paths)}_json_files.bed.gz.tbi")
    bp.run()


if __name__ == "__main__":
    main()


