set -x
set -e
set -u


python3 -u ./figures_and_tables/numbers1_for_intro_and_multiallelic_variants.py
python3 -u ./figures_and_tables/numbers2_for_validation_diagram.py
python3 -u ./figures_and_tables/numbers3_genomic_regions.py
python3 -u ./figures_and_tables/numbers4_variant_catalogs.py
python3 -u ./figures_and_tables/numbers5_tool_comparisons.py


cd ./figures_and_tables


python3 -u plot_syndip_indel_size_distribution.py
python3 -u plot_summary_stats.py
python3 -u plot_motif_distribution.py
python3 -u plot_tool_runtime_and_memory.py
python3 -u plot_expansion_hunter_denovo_results.py
python3 -u plot_gene_constraint_info.py

python3 -u plot_tool_comparisons_percent_exactly_right.py
python3 -u plot_tool_comparisons_by_motif_size.py
python3 -u plot_tool_comparisons_by_num_repeats.py --verbose

# upload to gs://str-truth-set/
gsutil -m cp *.png *.svg   gs://str-truth-set/hg38/figures/
gsutil -m cp ../tool_comparison/figures/*.svg   gs://str-truth-set/hg38/figures/
