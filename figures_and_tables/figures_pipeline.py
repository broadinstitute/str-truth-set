import os

from step_pipeline import pipeline, Backend, Localize, Delocalize
from plot_tool_accuracy_by_allele_size import PURITY_BINS, count_plots_per_purity_bin

DOCKER_IMAGE = "docker.io/weisburd/truth-set-figures@sha256:e65b59b683dcb3b63b2d6ec8b4aa815fbf28282457bf7d103097892703c358be"

# the current pipeline does not stratify accuracy-by-allele-size plots by Q (no tool's Q column is used for plotting)
ACCURACY_BY_ALLELE_SIZE_Q_THRESHOLD = 0


def main():
    bp = pipeline(f"STR Truth Set figures",
                  backend=Backend.HAIL_BATCH_SERVICE,
                  config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("-n", type=int, help="Only run the first n batches for each type of plot. Useful for testing.")
    parser.add_argument("--batch-size", type=int, default=50, help="Number of figures to generate per job")
    parser.add_argument("--output-dir", default="gs://str-truth-set/hg38/figures",
                        help="Google storage output directory for figures.")
    parser.add_argument("--input-table", default="gs://str-truth-set/hg38/combined.results.alleles.tsv.gz")
    args = bp.parse_known_args()

    # generate tool accuracy by allele size. plot_tool_accuracy_by_allele_size.py regenerates the full plot set once
    # per repeat-purity bin, so the shard range covers len(PURITY_BINS) x the per-bin plot count. Compute that count
    # from the plot script (it scales with DEFAULT_TOOLS and the Q setting) rather than a hard-coded literal, so the
    # range never truncates when tools are added; shards past the real total just produce empty jobs (the per-bin
    # count is an upper bound, since the real run skips tools absent from the input table).
    plots_per_purity_bin = count_plots_per_purity_bin(q_threshold=ACCURACY_BY_ALLELE_SIZE_Q_THRESHOLD)
    for i in range(0, plots_per_purity_bin * len(PURITY_BINS), args.batch_size):
        if args.n is not None and i >= args.n:
            break

        s1 = bp.new_step(f"tool accuracy by allele size plots #{i+1}-{i+args.batch_size}",
                         arg_suffix="accuracy-by-allele-size", image=DOCKER_IMAGE, cpu=2, memory="highmem",
                         output_dir=os.path.join(args.output_dir, "accuracy_by_allele_size"))
        local_truth_set_tsv = s1.input(args.input_table, localize_by=Localize.COPY)
        s1.command("set -ex")
        s1.command(f"python3 plot_tool_accuracy_by_allele_size.py --image-type svg --output-dir . --show-title "
                   f"--q-threshold {ACCURACY_BY_ALLELE_SIZE_Q_THRESHOLD} "
                   f"--start-with-plot-i {i} -n {args.batch_size} {local_truth_set_tsv}")
        s1.output(f"*.svg", delocalize_by=Delocalize.GSUTIL_COPY)

    # generate tool accuracy vs Q plots
    for i in range(0, 5*5*4*3*2, args.batch_size):
        if args.n is not None and i >= args.n:
            break

        s2 = bp.new_step(f"tool accuracy vs Q plots #{i+1}-{i+args.batch_size}",
                         arg_suffix="accuracy-vs-Q", image=DOCKER_IMAGE, cpu=2, memory="highmem",
                         output_dir=os.path.join(args.output_dir, "accuracy_vs_Q"))

        local_truth_set_tsv = s2.input(args.input_table, localize_by=Localize.COPY)
        s2.command("set -ex")
        s2.command(f"python3 plot_tool_accuracy_vs_Q.py --output-dir . --show-title --start-with-plot-i {i} -n {args.batch_size} {local_truth_set_tsv}")
        s2.output(f"*.svg", delocalize_by=Delocalize.GSUTIL_COPY)

    bp.run()


if __name__ == "__main__":
    main()


