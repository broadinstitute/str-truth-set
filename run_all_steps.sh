set -ex

# Generate 2 versions of the truth set - one with and one without homopolymers.
# These versions are useful for comparison and for generating stats for the paper.
# However, the final truth set that's used for tool comparisons and for most of
# the paper is the one without homopolymers.

#./run_step_A__create_STR_truth_set.sh --only-high-confidence-regions >& step_A.log  # this version is used as the final truth set
#./run_step_A__create_STR_truth_set.sh --only-high-confidence-regions --include-homopolymers >& step_A_with_homopolymers.log

# Run the rest of the pipeline
#./run_step_B__generate_variant_catalogs_and_upload_to_bucket.sh >& step_B.log
./run_step_C__run_tool_pipelines.sh >& step_C.log
./run_step_D__combine_tool_results.sh >& step_D.log
./run_step_E__generate_plots_and_figures.sh >& step_E.log

# Generate one final version of the truth set including variants outside SynDip high-confidence regions. This is used for computing stats.
./run_step_A__create_STR_truth_set.sh --include-homopolymers >& step_A_raw_with_homopolymers.log

echo Finished running all steps
