set -x
set -e
set -u


#python3 -u ./figures_and_tables/numbers1_for_intro_and_multiallelic_variants.py
#python3 -u ./figures_and_tables/numbers2_for_validation_diagram.py
#python3 -u ./figures_and_tables/numbers3_genomic_regions.py
#python3 -u ./figures_and_tables/numbers4_variant_catalogs.py
#python3 -u ./figures_and_tables/numbers5_tool_comparisons.py


gsutil -m cp -n  ./figures_and_tables/*.png   gs://str-truth-set/hg38/figures/

output_dir=svg_figures
mkdir -p ./figures_and_tables/${output_dir}
cd ./figures_and_tables/


python3 -u plot_syndip_indel_size_distribution.py --output-dir ${output_dir}
python3 -u plot_summary_stats.py --output-dir ${output_dir}
python3 -u plot_motif_distribution.py --output-dir ${output_dir}
python3 -u plot_tool_runtime_and_memory.py --output-dir ${output_dir}
python3 -u plot_expansion_hunter_denovo_results.py --output-dir ${output_dir}
python3 -u plot_gene_constraint_info.py --output-dir ${output_dir}

python3 -u plot_tool_accuracy_percent_exactly_right.py --show-title --output-dir ${output_dir}
python3 -u plot_tool_accuracy_percent_exactly_right.py --show-title --exclude-hipstr-no-call-loci --output-dir ${output_dir}
python3 -u plot_tool_accuracy_by_motif_size.py --output-dir ${output_dir}
python3 -u plot_mutation_rates.py --output-dir ${output_dir}


# upload figures to gs://str-truth-set/
gsutil -m cp ${output_dir}/*.svg   gs://str-truth-set/hg38/figures/${output_dir}/

# generate tool accuracy vs Q and tool accuracy by num repeats plots
python3 figures_pipeline.py --force --batch-size 25

gsutil -m cp -r gs://str-truth-set/hg38/figures/accuracy_vs_Q .
gsutil -m cp -r gs://str-truth-set/hg38/figures/accuracy_by_allele_size  .

#python3 -u plot_tool_accuracy_vs_Q.py --verbose
#python3 -u plot_tool_accuracy_by_allele_size.py --verbose

./generate_figure_panels_for_paper.sh

echo Done with step E
