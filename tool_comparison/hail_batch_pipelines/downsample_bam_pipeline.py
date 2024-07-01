import os
import re
import hailtop.fs as hfs
from step_pipeline import pipeline, Backend, Localize, Delocalize

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

CHM1_CHM13_CRAM_COVERAGE = 40

DOCKER_IMAGE = "weisburd/gatk@sha256:433406c13c62fb9088bb6cfa842278ff7b6980f540ab8390d224cae20d0f3742"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    parser = bp.get_config_arg_parser()
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--input-index-file", help="Path of BAM or CRAM index file")
    parser.add_argument("-t", "--target-coverage", default=[], action="append", help="Target coverage. Must be less than the input coverage.")
    parser.add_argument("--output-dir", help="Google storage directory for output file. If not specified, it will be based on the input bam")
    parser.add_argument("input_bam_or_cram", nargs="+", default=[CHM1_CHM13_CRAM_PATH], help="Path of input BAM or CRAM file")
    args = bp.parse_known_args()

    if args.input_index_file and len(args.input_bam_or_cram) > 1:
        parser.error("Cannot specify --input-index-file when more than one input BAM or CRAM file is provided. The index files must be in the same directory as the input BAM or CRAM files.")

    bp.set_name(f"Downsample: " + (os.path.basename(args.input_bam_or_cram[0]) if len(args.input_bam_or_cram) == 1 else f"{len(args.input_bam_or_cram)} files"))

    if not args.target_coverage:
        target_coverage = [30]
    else:
        for i, target_coverage in enumerate(args.target_coverage):
            try:
                args.target_coverage[i] = float(target_coverage)
            except Exception as e:
                parser.error(f"Invalid target coverage arg: {target_coverage}")

    if any(t < 2 for t in args.target_coverage):
        parser.error("--target-coverage arg must be >= 2")

    for input_bam_or_cram in args.input_bam_or_cram:
        output_dir = args.output_dir if args.output_dir else os.path.dirname(input_bam_or_cram)

        filename_prefix = re.sub("(.bam|.cram)$", "", os.path.basename(input_bam_or_cram))
        hfs_ls_results = hfs.ls(input_bam_or_cram)
        if len(hfs_ls_results) == 0:
            parser.error(f"Input BAM or CRAM file not found: {input_bam_or_cram}")
        read_data_size = int(hfs_ls_results[0].size/10**9)
        s1 = bp.new_step(f"depth: {os.path.basename(input_bam_or_cram)}", image=DOCKER_IMAGE, arg_suffix="depth", cpu=1, storage=f"{read_data_size + 20}Gi")
        s1.switch_gcloud_auth_to_user_account()
        local_fasta, _ = s1.inputs(args.reference_fasta, args.reference_fasta_fai, localize_by=Localize.COPY)

        if args.input_index_file:
            input_bam_or_cram_index = args.input_index_file
        elif input_bam_or_cram.endswith(".bam"):
            input_bam_or_cram_index = re.sub(".bam$", ".bam.bai", input_bam_or_cram)
        elif input_bam_or_cram.endswith(".cram"):
            input_bam_or_cram_index = re.sub(".cram$", ".cram.crai", input_bam_or_cram)
        else:
            parser.error(f"Input BAM or CRAM file must end with .bam or .cram: {input_bam_or_cram}")

        local_bam, _ = s1.inputs(input_bam_or_cram, input_bam_or_cram_index, localize_by=Localize.GSUTIL_COPY)

        s1.command("set -ex")
        s1.command("curl -L https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth -o /usr/local/bin/mosdepth")
        s1.command("chmod 777 /usr/local/bin/mosdepth")

        s1.command("cd /io/")
        s1.command(f"mosdepth -f {local_fasta} -x {filename_prefix}.coverage {local_bam}")
        #s1.output(f"{filename_prefix}.coverage.mosdepth.summary.txt")

        s1.command(f"cat {filename_prefix}.coverage.mosdepth.summary.txt | cut -f 4 | tail -n +2 | head -n 23")
        s1.command(f"grep total {filename_prefix}.coverage.mosdepth.summary.txt > {filename_prefix}.total_depth.txt")
        s1.command(f"cat {filename_prefix}.total_depth.txt")

        s1.output(f"/io/{filename_prefix}.total_depth.txt", os.path.join(output_dir, f"{filename_prefix}.total_depth.txt"))

        for target_coverage in args.target_coverage:
            s2 = bp.new_step(f"downsample: {os.path.basename(input_bam_or_cram)} to {target_coverage}x",
                             image=DOCKER_IMAGE, cpu=2, memory="highmem", storage=f"{read_data_size+20}Gi",
                             output_dir=output_dir)
            s2.depends_on(s1)
            local_fasta, _ = s2.inputs(args.reference_fasta, args.reference_fasta_fai, localize_by=Localize.COPY)
            local_bam, _ = s2.inputs(input_bam_or_cram, input_bam_or_cram_index, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
            total_depth_file = s2.use_previous_step_outputs_as_inputs(s1, localize_by=Localize.COPY)

            output_bam_filename_prefix = f"{filename_prefix}.downsampled_to_{int(target_coverage)}x"
            s2.command("curl -L https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth -o /usr/local/bin/mosdepth")
            s2.command("chmod 777 /usr/local/bin/mosdepth")

            s2.command("set -ex")
            s2.command("cd /io/")
            s2.command(f"time gatk --java-options '-Xmx11G' DownsampleSam "
                       f"--VALIDATION_STRINGENCY SILENT "
                       f"--REFERENCE_SEQUENCE {local_fasta} "
                       f"-I {local_bam} "
                       f"-O {output_bam_filename_prefix}.bam "
                       f"""-P $(echo "{target_coverage} / $(grep total {total_depth_file} | cut -f 4)" | bc -l | awk '{{printf "%.4f", $0}}') """
                       f"--CREATE_INDEX true")

            s2.command(f"samtools calmd -b {output_bam_filename_prefix}.bam {local_fasta} > {output_bam_filename_prefix}.with_NM_tag.bam")
            s2.command(f"mv {output_bam_filename_prefix}.with_NM_tag.bam {output_bam_filename_prefix}.bam")
            s2.command(f"samtools index {output_bam_filename_prefix}.bam")

            s2.command("ls -lh")

            s2.command(f"mosdepth -f {local_fasta} -x coverage_after_downsampling {output_bam_filename_prefix}.bam")
            s2.command(f"cat coverage_after_downsampling.mosdepth.summary.txt | cut -f 4 | tail -n +2 | head -n 23")
            s2.command(f"grep total coverage_after_downsampling.mosdepth.summary.txt > {output_bam_filename_prefix}.total_depth.txt")
            s2.command(f"cat {output_bam_filename_prefix}.total_depth.txt")

            s2.output(f"{output_bam_filename_prefix}.total_depth.txt")
            s2.output(f"{output_bam_filename_prefix}.bam")
            s2.output(f"{output_bam_filename_prefix}.bam.bai")

    bp.run()


if __name__ == "__main__":
    main()


