set -x
set -e
set -u


cd ./figures_and_tables

python3 -u plot_summary_stats.py
python3 -u plot_tool_comparisons.py --verbose

# upload to gs://str-truth-set/
gsutil -m cp -n ../tool_comparison/figures/*.svg   gs://str-truth-set/hg38/figures/
