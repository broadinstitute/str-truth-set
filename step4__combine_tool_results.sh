set -x
set -e
set -u

# EHdn
python3 ./tool_comparison/scripts/intersect_expansion_hunter_denovo_results_with_truth_set.py ./tool_comparison/results*/expansion_hunter_denovo/CHM1_CHM13_WGS2.*locus.tsv

# ExpansionHunter and GangSTR
for results_folder in results_for_exome  results  results_for_downsampled_30x_bam  results_for_downsampled_20x_bam  results_for_downsampled_10x_bam  results_for_downsampled_5x_bam;
do
  echo Processing $results_folder ...
  python3 tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py \
    --output-dir ./tool_comparison/${results_folder}/  \
    STR_truthset.v1.variants.tsv.gz \
    ./tool_comparison/variant_catalogs/negative_loci.tsv.gz

  # ExpansionHunter
  python3 tool_comparison/scripts/add_tool_results_columns.py \
      --tool ExpansionHunter \
     ./tool_comparison/${results_folder}/expansion_hunter/positive_loci/combined.positive_loci.*_json_files.variants.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truthset.v1.variants.for_comparison.tsv.gz

  python3 tool_comparison/scripts/add_tool_results_columns.py \
      --tool ExpansionHunter \
      ./tool_comparison/${results_folder}/expansion_hunter/negative_loci/combined.negative_loci.*_json_files.variants.tsv.gz \
      ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz

  mv ./tool_comparison/${results_folder}/STR_truthset.v1.variants.for_comparison.with_ExpansionHunter_results.tsv.gz ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.tsv.gz
  mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_ExpansionHunter_results.tsv.gz ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz

  # GangSTR
  python3 tool_comparison/scripts/add_tool_results_columns.py \
    --tool GangSTR \
    ./tool_comparison/${results_folder}/gangstr/positive_loci/combined.positive_loci.*_json_files.variants.tsv.gz \
    ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.tsv.gz

  python3 tool_comparison/scripts/add_tool_results_columns.py \
    --tool GangSTR \
    ./tool_comparison/${results_folder}/gangstr/negative_loci/combined.negative_loci.*_json_files.variants.tsv.gz \
    ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz

  mv ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.with_GangSTR_results.tsv.gz ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.tsv.gz
  mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_GangSTR_results.tsv.gz ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz

  # concordance
  python3 tool_comparison/scripts/add_concordance_columns.py \
    ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.tsv.gz

  python3 tool_comparison/scripts/add_concordance_columns.py \
    ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz

  mv ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.with_concordance.tsv.gz ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.variants.tsv.gz
  mv ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.with_concordance.alleles.tsv.gz ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.alleles.tsv.gz
  mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_concordance.tsv.gz ./tool_comparison/${results_folder}/negative_loci.for_comparison.variants.tsv.gz
  mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_concordance.alleles.tsv.gz ./tool_comparison/${results_folder}/negative_loci.for_comparison.alleles.tsv.gz

  gunzip -f ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.variants.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/negative_loci.for_comparison.variants.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/STR_truthset.v1.for_comparison.alleles.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/negative_loci.for_comparison.alleles.tsv.gz

done

python3 ./tool_comparison/scripts/combine_all_results_tables.py
