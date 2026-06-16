import functools
import hailtop.fs as hfs
import logging
import math
import os
import re
import subprocess

from step_pipeline import pipeline, Backend, Localize, Delocalize


@functools.lru_cache(maxsize=None)
def _count_catalog_loci(catalog_path):
    """Return the number of loci in an ExpansionHunter variant-catalog JSON (one "LocusId" key per locus).

    Used to split a single catalog into evenly sized index ranges for the bw2-fork --start-with/--n-loci shards.
    Reads the whole catalog into memory once -- via `gcloud storage cat` (which retries on transient errors) for
    gs:// paths, or directly off disk for local paths -- instead of hailtop's streaming reader, which times out on
    the ~100MB catalogs. Cached because every EHv5 run for a given sample reuses the same catalog.
    """
    if catalog_path.startswith("gs://"):
        data = subprocess.run(["gcloud", "storage", "cat", catalog_path], capture_output=True, check=True).stdout
    else:
        with open(catalog_path, "rb") as f:
            data = f.read()
    return data.count(b'"LocusId"')

DOCKER_IMAGE = "weisburd/str-analysis-with-expansion-hunter@sha256:3907d127e607a32459b08ba157d322d68476123d6ecb73bfecc5294066884d8b"
DOCKER_IMAGE_DEV = "weisburd/expansion-hunter-dev@sha256:4baa218fdb7bc76af97d6628d13718ff4ad22290d1cc4ffeb2f8e3ea3c5a13b3"

REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

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
    parser.add_argument("--analysis-mode", choices=["seeking", "streaming", "low-mem-streaming", "optimized-streaming"],
                        default="optimized-streaming", help="Run ExpansionHunter with this --analysis-mode.")
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
        tool_exec = "IlluminaExpansionHunter"
        output_dir = os.path.join(args.output_dir, tool_exec, positive_or_negative_loci)
    else:
        tool_exec = "ExpansionHunter"
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
        analysis_mode=args.analysis_mode,
        loci_to_exclude=args.loci_to_exclude,
        min_locus_coverage=args.min_locus_coverage,
        use_illumina_expansion_hunter=args.use_illumina_expansion_hunter,
        run_reviewer=args.run_reviewer)
    bp.run()


def create_expansion_hunter_steps(bp, *, reference_fasta, input_bam, input_bai, variant_catalog_file_paths, output_dir, output_prefix, reference_fasta_fai=None, male_or_female="female",
                                  analysis_mode="seeking", loci_to_exclude=None, min_locus_coverage=None, use_illumina_expansion_hunter=False, run_reviewer=False, num_shards=1,
                                  catalog_prefilter_step=None):

    if use_illumina_expansion_hunter:
        tool_exec = "IlluminaExpansionHunter"
        # the official Illumina ExpansionHunter build only supports seeking/streaming; optimized-streaming and
        # low-mem-streaming are bw2-fork-only modes, so clamp to streaming to avoid a CLI error
        if analysis_mode not in ("seeking", "streaming"):
            analysis_mode = "streaming"
    else:
        tool_exec = "ExpansionHunter"

    min_locus_coverage_arg = ""
    if min_locus_coverage is not None:
        min_locus_coverage_arg = f"--min-locus-coverage {min_locus_coverage}"

    if loci_to_exclude:
        # loci_to_exclude is a list of file paths; read the locus ids out of each into a separate list (mutating the
        # list being iterated would loop forever and then try to open() a locus id as a file)
        locus_ids_to_exclude = []
        for loci_to_exclude_file_path in loci_to_exclude:
            with open(loci_to_exclude_file_path, "rt") as f:
                for line in f:
                    locus_ids_to_exclude.append(line.strip())
        loci_to_exclude = locus_ids_to_exclude
        print(f"Parsed {len(loci_to_exclude)} locus ids to exclude:", ", ".join(loci_to_exclude[:5]),
              "..." if len(loci_to_exclude) > 5 else "")

    step1s = []
    step1_output_paths = []

    hfs_ls_results = hfs.ls(input_bam)
    if len(hfs_ls_results) == 0:
        raise ValueError(f"No files found matching {input_bam}")
    input_bam_file_stats = hfs_ls_results[0]

    # Build the list of genotyping work units: (variant_catalog_path, start_with, n_loci, shard_suffix).
    # EHv5-bw2 streaming modes (low-mem-streaming, optimized-streaming) genotype loci single-threaded, so when
    # num_shards>1 split the single catalog into num_shards index ranges run as parallel 1-cpu jobs via the bw2-fork
    # --start-with/--n-loci flags (start_with=None means "process the whole catalog", the unsharded behavior).
    if (not use_illumina_expansion_hunter) and analysis_mode in ("low-mem-streaming", "optimized-streaming") \
            and num_shards > 1 and len(variant_catalog_file_paths) == 1:
        catalog_path = variant_catalog_file_paths[0]
        total_loci = _count_catalog_loci(catalog_path)
        loci_per_shard = math.ceil(total_loci / num_shards)
        work_units = []
        for i in range(num_shards):
            start = i * loci_per_shard
            if start >= total_loci:
                break
            work_units.append((catalog_path, start, loci_per_shard, f"shard{i:03d}_of_{num_shards:03d}"))
        print(f"EHv5 {analysis_mode}: sharding {total_loci:,} loci into {len(work_units)} parallel 1-cpu jobs of "
              f"{loci_per_shard:,} loci each for {os.path.basename(input_bam)}")
    else:
        work_units = [(p, None, None, None) for p in variant_catalog_file_paths]

    for catalog_i, (variant_catalog_path, start_with, n_loci, shard_suffix) in enumerate(work_units):
        s1 = bp.new_step(
            f"Run EHv5:{analysis_mode} #{catalog_i} on {os.path.basename(input_bam)} "
            f"({os.path.basename(variant_catalog_path)}{(' ' + shard_suffix) if shard_suffix else ''})",
            arg_suffix=f"run-expansion-hunter-step",
            step_number=1,
            image=DOCKER_IMAGE,
            # EHv5-bw2 streaming modes (low-mem-streaming, optimized-streaming) genotype loci single-threaded
            # (measured ~1.1 cores, <=4.3GB RSS at 31x), so cpu=1/highmem (6.5GB > 4.3GB peak) right-sizes them
            # instead of the old 16/lowmem (which paid for 16 cores while using ~1); cpu=1 standard (3.75GB) would
            # OOM at high coverage. --threads only sped up the brief mate-caching. Only the official IlluminaEHv5
            # build (analysis_mode == "streaming") keeps 16/highmem (untested here).
            cpu=16 if analysis_mode == "streaming" else (1 if "streaming" in analysis_mode else 2),
            memory="highmem" if "streaming" in analysis_mode else "standard",
            localize_by=Localize.GSUTIL_COPY,
            storage=f"{int(input_bam_file_stats.size/10**9) + 25}Gi",
            output_dir=output_dir)

        step1s.append(s1)

        # IlluminaEHv5 reads a catalog produced by an upstream prefilter step (drops loci with N-rich flanks), so the
        # genotyping job must wait for that step to write the filtered catalog.
        if catalog_prefilter_step is not None:
            s1.depends_on(catalog_prefilter_step)

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

        # per-catalog-chunk output prefix; kept distinct from the function's output_prefix arg, which names the
        # combined step-2 outputs below
        chunk_prefix = re.sub(".json$", "", local_variant_catalog.filename)
        if shard_suffix:
            chunk_prefix = f"{chunk_prefix}.{shard_suffix}"
        s1.command(f"echo Genotyping $(cat {local_variant_catalog_path} | grep LocusId | wc -l) loci")


        extra_args = ""
        # --cache-mates is a bw2-fork-only flag; the stock Illumina build rejects it
        if analysis_mode == "seeking" and not use_illumina_expansion_hunter: extra_args += "--cache-mates "
        if analysis_mode == "streaming": extra_args += "--threads 16 "       # IlluminaEHv5 build (cpu=16)
        elif "streaming" in analysis_mode: extra_args += "--threads 1 "       # EHv5-bw2 streaming modes (cpu=1)
        # bw2-fork locus slice: this shard genotypes loci [start_with, start_with + n_loci) of the sorted catalog
        if start_with is not None: extra_args += f"--start-with {start_with} --n-loci {n_loci} "

        s1.command(f"""/usr/bin/time --verbose {tool_exec} {extra_args} {min_locus_coverage_arg} \
            --reference {local_fasta} \
            --reads {local_bam} \
            --analysis-mode {analysis_mode} \
            --sex {male_or_female} \
            --variant-catalog {local_variant_catalog_path} \
            --output-prefix {chunk_prefix}""")

        s1.command("ls -lhrt")

        s1.output(f"{chunk_prefix}.json", output_dir=os.path.join(output_dir, f"json"))

        step1_output_paths.append(os.path.join(output_dir, f"json", f"{chunk_prefix}.json"))

        if run_reviewer:
            reviewer_remote_output_dir = os.path.join(output_dir, f"svg")
            reviewer_output_prefix = re.sub("(.bam|.cram)$", "", local_bam.filename)
            s1.command(f"samtools sort {chunk_prefix}_realigned.bam -o {chunk_prefix}_realigned.sorted.bam")
            s1.command(f"samtools index {chunk_prefix}_realigned.sorted.bam")
            s1.command(f"""/usr/bin/time --verbose REViewer \
                --reference {local_fasta}  \
                --catalog {local_variant_catalog_path} \
                --reads {chunk_prefix}_realigned.sorted.bam \
                --vcf {chunk_prefix}.vcf \
                --output-prefix {reviewer_output_prefix}
            """)

            done_file = f"done_generating_reviewer_images_for_{chunk_prefix}"
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
                     localize_by=Localize.GSUTIL_COPY,
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


