import hailtop.fs as hfs
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize

DOCKER_IMAGE = "weisburd/expansion-hunter@sha256:aa315698ca40e3e237a7d01f7c14fc9257b74151dbcad00e25116436f1594a65"

REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

CHM1_CHM13_CRAM_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
CHM1_CHM13_CRAI_PATH = "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

VARIANT_CATALOG_POSITIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/expansion_hunter/positive_loci.EHv5.*_of_293.json"
VARIANT_CATALOG_NEGATIVE_LOCI = "gs://str-truth-set/hg38/variant_catalogs/expansion_hunter/negative_loci.EHv5.*_of_291.json"

OUTPUT_BASE_DIR = "gs://str-truth-set/hg38/tool_results/expansion_hunter"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser_group = parser.add_mutually_exclusive_group(required=True)
    parser_group.add_argument("--positive-loci", action="store_true", help="Genotype truth set loci")
    parser_group.add_argument("--negative-loci", action="store_true", help="Genotype negative (hom-ref) loci")
    parser_group.add_argument("--variant-catalog", help="Path of variant catalog json file(s) to process")

    parser.add_argument("--use-illumina-expansion-hunter", action="store_true", help="Go back to using the Illumina "
         "version of ExpansionHunter instead of the optimized version from https://github.com/bw2/ExpansionHunter.git")
    parser.add_argument("--loci-to-exclude", action="append", help="Path of a file containing locus ids (one per line) "
        "to exclude from the variant catalogs before running ExpansionHunter. This can be useful when running with the "
        "--use-illumina-expansion-hunter in order to first filter out the small number of loci that cause "
        "ExpansionHunter to exit with the error message 'Flanks can contain at most 5 characters N but found x Ns.'")
    parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--use-streaming-mode", action="store_true", help="Run ExpansionHunter with "
        "--analysis-mode streaming. This uses ~100Gb of RAM, but can process all loci at once.")
    parser.add_argument("--input-bam", default=CHM1_CHM13_CRAM_PATH)
    parser.add_argument("--input-bai", default=CHM1_CHM13_CRAI_PATH)
    parser.add_argument("--min-locus-coverage", type=int, help="Sets ExpansionHunter's --min-locus-coverage arg")
    parser.add_argument("--output-dir", default=OUTPUT_BASE_DIR)
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    parser.add_argument("--run-reviewer", action="store_true", help="Generate REViewer read visualizations for all loci")
    args = bp.parse_known_args()

    if args.positive_loci:
        positive_or_negative_loci = "positive_loci"
        variant_catalog_file_paths = [x.path for x in hfs.ls(VARIANT_CATALOG_POSITIVE_LOCI)]
        if len(variant_catalog_file_paths) == 0:
            raise ValueError(f"No files found matching {VARIANT_CATALOG_POSITIVE_LOCI}")
    elif args.negative_loci:
        positive_or_negative_loci = "negative_loci"
        variant_catalog_file_paths = [x.path for x in hfs.ls(VARIANT_CATALOG_NEGATIVE_LOCI)]
        if len(variant_catalog_file_paths) == 0:
            raise ValueError(f"No files found matching {VARIANT_CATALOG_NEGATIVE_LOCI}")
    elif args.variant_catalog:
        positive_or_negative_loci = os.path.basename(args.variant_catalog).replace(".json", "")
        variant_catalog_file_paths = [x.path for x in hfs.ls(args.variant_catalog)]
    else:
        parser.error("Must specify either --positive-loci or --negative-loci")

    if args.use_illumina_expansion_hunter:
        output_dir = os.path.join(args.output_dir, tool_exec, positive_or_negative_loci)
    else:
        output_dir = os.path.join(args.output_dir, positive_or_negative_loci)

    if args.n:
        variant_catalog_file_paths = variant_catalog_file_paths[:args.n]

    bam_path_ending = "/".join(args.input_bam.split("/")[-2:])
    bp.set_name(f"STR Truth Set: {tool_exec}: {positive_or_negative_loci}: {bam_path_ending}")
    if not args.force:
        existing_json_paths = bp.precache_file_paths(os.path.join(output_dir, f"**/*.json"))
        logging.info(f"Precached {len(existing_json_paths)} json files")

    create_expansion_hunter_steps(
        bp,
        reference_fasta=args.reference_fasta,
        input_bam=args.input_bam,
        input_bai=args.input_bai,
        male_or_female="female",
        variant_catalog_file_paths=variant_catalog_file_paths,
        output_dir=output_dir,
        output_prefix=f"combined.{positive_or_negative_loci}",
        reference_fasta_fai=args.reference_fasta_fai,
        use_streaming_mode=args.use_streaming_mode, 
        loci_to_exclude=args.loci_to_exclude,
        min_locus_coverage=args.min_locus_coverage,
        use_illumina_expansion_hunter=args.use_illumina_expansion_hunter,
        run_reviewer=args.run_reviewer)
    bp.run()
        
def create_expansion_hunter_steps(bp, *, reference_fasta, input_bam, input_bai, variant_catalog_file_paths, output_dir, output_prefix, reference_fasta_fai=None, male_or_female="female",
                                  use_streaming_mode=False, loci_to_exclude=None, min_locus_coverage=None, use_illumina_expansion_hunter=False, run_reviewer=False):

    if use_illumina_expansion_hunter:
        tool_exec = "IlluminaExpansionHunter"
        cache_mates_arg = ""
    else:
        tool_exec = "ExpansionHunter"
        cache_mates_arg = "--cache-mates "

    min_locus_coverage_arg = ""
    if min_locus_coverage is not None:
        min_locus_coverage_arg = f"--min-locus-coverage {min_locus_coverage}"

    if loci_to_exclude:
        for loci_to_exclude_file_path in loci_to_exclude:
            with open(loci_to_exclude_file_path, "rt") as f:
                for line in f:
                    loci_to_exclude.append(line.strip())
        print(f"Parsed {len(loci_to_exclude)} locus ids to exclude:", ", ".join(loci_to_exclude[:5]),
              "..." if len(loci_to_exclude) > 5 else "")

    step1s = []
    step1_output_paths = []

    hfs_ls_results = hfs.ls(input_bam)
    if len(hfs_ls_results) == 0:
        raise ValueError(f"No files found matching {input_bam}")
    input_bam_file_stats = hfs_ls_results[0]

    for catalog_i, variant_catalog_path in enumerate(variant_catalog_file_paths):
        if not use_streaming_mode:
            s1 = bp.new_step(
                f"Run EHv5 #{catalog_i} on {os.path.basename(input_bam)} ({os.path.basename(variant_catalog_path)})",
                arg_suffix=f"run-expansion-hunter-step",
                step_number=1,
                image=DOCKER_IMAGE,
                cpu=2,
                localize_by=Localize.COPY,
                storage=f"{int(input_bam_file_stats.size/10**9) + 25}Gi",
                output_dir=output_dir)
        else:
            s1 = bp.new_step(
                f"Run EHv5 #{catalog_i} on {os.path.basename(input_bam)} ({os.path.basename(variant_catalog_path)})",
                arg_suffix=f"run-expansion-hunter-step",
                step_number=1,
                image=DOCKER_IMAGE,
                cpu=16,
                memory="highmem",
                #custom_machine_type="n1-highmem-32",
                #custom_machine_is_preemptible=True,
                storage=f"{int(input_bam_file_stats.size/10**9) + 25}Gi",
                output_dir=output_dir)

        step1s.append(s1)

        local_fasta = s1.input(reference_fasta)
        if reference_fasta_fai:
            s1.input(reference_fasta_fai)

        local_bam = s1.input(input_bam)
        if input_bai:
            s1.input(input_bai)

        s1.command("set -ex")

        local_variant_catalog = s1.input(variant_catalog_path)
        if loci_to_exclude:
            local_variant_catalog_path = f"filtered_{local_variant_catalog.filename}"

            # add command to filter out excluded loci from variant catalog
            jq_command = f"cat {local_variant_catalog} | "
            jq_command += "jq '.[] | " + " | ".join([f'select(.LocusId != "{i}")' for i in loci_to_exclude]) + "' | "
            jq_command += "jq -s '.' "  # reformat output into proper json
            jq_command += f" > {local_variant_catalog_path}"
            s1.command(jq_command)
        else:
            local_variant_catalog_path = str(local_variant_catalog)

        output_prefix = re.sub(".json$", "", local_variant_catalog.filename)
        s1.command(f"echo Genotyping $(cat {local_variant_catalog_path} | grep LocusId | wc -l) loci")

        if not use_streaming_mode:
            s1.command(f"""/usr/bin/time --verbose {tool_exec} {cache_mates_arg} {min_locus_coverage_arg} \
                --reference {local_fasta} \
                --reads {local_bam} \
                --sex {male_or_female} \
                --variant-catalog {local_variant_catalog_path} \
                --output-prefix {output_prefix}""")
        else:
            s1.command(f"""/usr/bin/time --verbose {tool_exec} {min_locus_coverage_arg} \
                --reference {local_fasta} \
                --reads {local_bam} \
                --sex {male_or_female} \
                --variant-catalog {local_variant_catalog_path} \
                --analysis-mode streaming \
                --threads 16 \
                --output-prefix {output_prefix}""")

        s1.command("ls -lhrt")

        s1.output(f"{output_prefix}.json", output_dir=os.path.join(output_dir, f"json"))

        step1_output_paths.append(os.path.join(output_dir, f"json", f"{output_prefix}.json"))

        if run_reviewer:
            reviewer_remote_output_dir = os.path.join(output_dir, f"svg")
            reviewer_output_prefix = re.sub("(.bam|.cram)$", "", local_bam.filename)
            s1.command(f"samtools sort {output_prefix}_realigned.bam -o {output_prefix}_realigned.sorted.bam")
            s1.command(f"samtools index {output_prefix}_realigned.sorted.bam")
            s1.command(f"""/usr/bin/time --verbose REViewer \
                --reference {local_fasta}  \
                --catalog {local_variant_catalog_path} \
                --reads {output_prefix}_realigned.sorted.bam \
                --vcf {output_prefix}.vcf \
                --output-prefix {reviewer_output_prefix}
            """)

            done_file = f"done_generating_reviewer_images_for_{output_prefix}"
            s1.command(f"touch {done_file}")

            s1.output("*.svg", output_dir=reviewer_remote_output_dir, delocalize_by=Delocalize.GSUTIL_COPY)
            s1.output(done_file, output_dir=reviewer_remote_output_dir)

    # step2: combine json files
    s2 = bp.new_step(name=f"Combine EHv5 outputs for {os.path.basename(input_bam)}",
                     step_number=2,
                     arg_suffix=f"combine-expansion-hunter-step",
                     image=DOCKER_IMAGE,
                     cpu=1,
                     memory="highmem",
                     storage="20Gi",
                     output_dir=output_dir)

    for step1 in step1s:
        s2.depends_on(step1)

    s2.command("mkdir /io/run_dir; cd /io/run_dir")
    for json_path in step1_output_paths:
        local_path = s2.input(json_path)
        s2.command(f"ln -s {local_path}")

    s2.command("set -x")
    s2.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv --include-extra-expansion-hunter-fields "
               f"--output-prefix {output_prefix}")
    s2.command(f"bgzip {output_prefix}.{len(step1_output_paths)}_json_files.bed")
    s2.command(f"tabix {output_prefix}.{len(step1_output_paths)}_json_files.bed.gz")
    s2.command("ls -lhrt")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.variants.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.alleles.tsv.gz")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.bed.gz")
    s2.output(f"{output_prefix}.{len(step1_output_paths)}_json_files.bed.gz.tbi")

    return s2

if __name__ == "__main__":
    main()


