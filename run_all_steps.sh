set -ex

#./step_A__create_STR_truth_set.sh >& step_A.log
#./step_B__upload_files_to_STR_truth_set_bucket.sh >& step_B.log
./step_C__run_tool_pipelines.sh >& step_C.log
./step_D__combine_tool_results.sh >& step_D.log
./step_E__generate_plots_and_figures.sh >& step_E.log

echo Finished running all steps
