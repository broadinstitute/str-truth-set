set -ex

# figure 1 panels
python3 plot_syndip_indel_size_distribution.py --width 12 --height 4 --image-type png

# figure 2 panels
python3 plot_summary_stats.py --image-type png --only-plot 2 --width 8 --height 5.5
python3 plot_summary_stats.py --image-type png --only-plot 3 --width 30 --height 6.5
python3 plot_summary_stats.py --image-type png --only-plot 4 --width 16 --height 5.5

# figure 3 panels
python3 plot_summary_stats.py --image-type png --only-plot 6 --width 11 --height 12
python3 plot_motif_distribution.py --width 9 --height 8  --image-type png
python3 plot_motif_distribution.py --width 9 --height 8  --image-type png --only-pure-repeats

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
        --min-motif-size 2 --max-motif-size 6 --only-pure-repeats  --exclude-no-call-loci --coverage 30x

# figure 5 panels & supp. figure 3
 python3 plot_tool_accuracy_by_allele_size.py --q-threshold 0 --coverage 30x --min-motif-size 2 --max-motif-size 6 \
       --genotype all --only-pure-repeats --show-no-call-loci --image-type png \
       --width 20 --height 6.5  # this will generate plots for ExpansionHunter, GangSTR, and HipSTR

python3 plot_tool_accuracy_vs_Q.py --coverage 30x --min-motif-size 2 --max-motif-size 6 --genotype all \
       --image-type png --only-pure-repeats \
       --width 9 --height 10   # this will generate plots with and without no-call loci

python3 plot_tool_runtime_and_memory.py --image-type png

# figure 6 panels
python3 plot_expansion_hunter_denovo_results.py --width 10 --height 7 --image-type png --only-pure-repeats \
  --truth-set-ehdn-input-table ../tool_comparison/results_for_downsampled_30x_bam/expansion_hunter_denovo/CHM1_CHM13_WGS2.downsampled_to_30x.truth_set_EHdn_comparison_table.tsv \
  --ehdn-truth-set-input-table ../tool_comparison/results_for_downsampled_30x_bam/expansion_hunter_denovo/CHM1_CHM13_WGS2.downsampled_to_30x.EHdn_results_table.with_truth_set_concordance.tsv

# figure 7 panels
python3 plot_summary_stats.py --image-type png --only-plot 5 --only-pure-repeats --width 8 --height 7

python3 plot_gene_constraint_info.py  --image-type png

python3 plot_mutation_rates.py --image-type png --width 8.5 --height 8.5

# supp. figure 1 A, B
python3 ../figures_and_tables/plot_summary_stats.py --only-plot 7 --image-type png

# supp. figure 2
cd ..
python3 paper_generate_table3_and_supp_fig2_repeat_catalogs.py

