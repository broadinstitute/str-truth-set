set -x
set -e
set -u


cd ./figures_and_tables

python3 -u plot_syndip_indel_size_distribution.py
python3 -u plot_summary_stats.py
python3 -u plot_motif_distribution.py
python3 -u plot_tool_runtime_and_memory.py

#python3 -u plot_tool_comparisons.py --verbose

# upload to gs://str-truth-set/
#gsutil -m cp -n ../tool_comparison/figures/*.svg   gs://str-truth-set/hg38/figures/