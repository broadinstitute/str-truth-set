set -ex

python3 plot_syndip_indel_size_distribution.py --width 12 --height 4 --image-type png
python3 plot_tool_accuracy_percent_exactly_right.py --width 11 --height 7 --image-type png
python3 plot_tool_accuracy_percent_exactly_right.py --width 11 --height 7 --exclude-hipstr-no-call-loci --image-type png

python3 plot_tool_runtime_and_memory.py --image-type png

python3 plot_mutation_rates.py --width 8 --height 6 --image-type png
