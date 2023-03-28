set -ex
#./run_step_A__create_STR_truth_set.sh --only-high-confidence-regions --include-homopolymers >& step_A_with_homopolymers.log
#./run_step_A__create_STR_truth_set.sh --only-high-confidence-regions >& step_A.log

./run_step_B__generate_variant_catalogs_and_upload_to_bucket.sh >& step_B.log
./run_step_C__run_tool_pipelines.sh >& step_C.log
./run_step_D__combine_tool_results.sh >& step_D.log
./run_step_E__generate_plots_and_figures.sh >& step_E.log

echo Finished running all steps
