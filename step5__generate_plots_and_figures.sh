set -x
set -e
set -u


cd ./figures

python3 -u plot_tool_comparisons.py --verbose

