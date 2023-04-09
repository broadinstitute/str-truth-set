set -ex

# figure 1 panels
python3 plot_syndip_indel_size_distribution.py --width 12 --height 4 --image-type png

# figure 2 panels
python3 plot_summary_stats.py --image-type png --only-plot 2 --width 8 --height 5.5
python3 plot_summary_stats.py --image-type png --only-plot 3 --width 30 --height 6.5
python3 plot_summary_stats.py --image-type png --only-plot 4 --width 16 --height 5.5

# figure 3 panels
python3 plot_motif_distribution.py --width 9 --height 8  --image-type png
python3 plot_motif_distribution.py --width 9 --height 8  --image-type png --only-pure-repeats
python3 plot_summary_stats.py --image-type png --only-plot 6 --width 11 --height 12

# figure 4 panels & supp. figure 2
python3 plot_tool_accuracy_percent_exactly_right.py --width 14 --height 7 --image-type png \
        --min-motif-size 2 --max-motif-size 6 --only-pure-repeats --coverage 30x
python3 plot_tool_accuracy_percent_exactly_right.py --width 17 --height 7 --image-type png \
        --min-motif-size 2 --max-motif-size 6  --only-pure-repeats --by-coverage
python3 plot_tool_accuracy_percent_exactly_right.py --width 14 --height 7 --image-type png \
        --min-motif-size 2 --max-motif-size 6  --only-pure-repeats --exclude-no-call-loci --coverage 30x
python3 plot_tool_accuracy_percent_exactly_right.py --width 17 --height 7 --image-type png \
        --min-motif-size 2 --max-motif-size 6  --only-pure-repeats --exclude-no-call-loci --by-coverage

python3 plot_tool_accuracy_percent_exactly_right.py --width 14 --height 7 --image-type png \
        --min-motif-size 2 --max-motif-size 6 --only-pure-repeats  --excluding-no-call-loci --coverage 30x

# figure 5 panels & supp. figure 3
python3 plot_tool_accuracy_by_allele_size.py --q-threshold 0 --coverage 30x --min-motif-size 2 --max-motif-size 6 \
        --genotype all --only-pure-repeats --show-no-call-loci --tool ExpansionHunter --image-type png --width 20 --height 6.5
python3 plot_tool_accuracy_by_allele_size.py --q-threshold 0 --coverage 30x --min-motif-size 2 --max-motif-size 6 \
        --genotype all --only-pure-repeats --show-no-call-loci --tool GangSTR --image-type png --width 20 --height 6.5
python3 plot_tool_accuracy_by_allele_size.py --q-threshold 0 --coverage 30x --min-motif-size 2 --max-motif-size 6 \
        --genotype all --only-pure-repeats --show-no-call-loci --tool HipSTR --image-type png --width 20 --height 6.5

python3 plot_tool_accuracy_vs_Q.py --image-type png --coverage 30x --min-motif-size 2 --max-motif-size 6 \
       --genotype all --only-pure-repeats --show-no-call-loci --width 9 --height 10
python3 plot_tool_accuracy_vs_Q.py --image-type png --coverage 30x --min-motif-size 2 --max-motif-size 6 \
       --genotype all --only-pure-repeats --hide-no-call-loci --width 9 --height 10

python3 plot_tool_runtime_and_memory.py --image-type png

# figure 6 panels
python3 plot_expansion_hunter_denovo_results.py --width 10 --height 7 --image-type png --only-pure-repeats

# figure 7 panels
python3 plot_summary_stats.py --image-type png --only-plot 5 --only-pure-repeats --width 8 --height 7

python3 plot_gene_constraint_info.py  --image-type png

python3 plot_mutation_rates.py --image-type png --width 8 --height 6


