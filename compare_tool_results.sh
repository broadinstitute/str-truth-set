set -x
set -e
set -u


python3 tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py \
	--output-dir ./tool_comparison/results/  \
	STR_truthset.v1.variants.tsv.gz \
	./tool_comparison/variant_catalogs/negative_loci.tsv.gz

# ExpansionHunter
python3 tool_comparison/scripts/add_tool_results_columns.py \
    --tool ExpansionHunter \
   ./tool_comparison/results/expansion_hunter_positive_loci/positive_loci.expansion_hunter.289_json_files.variants.tsv.gz \
   ./tool_comparison/results/STR_truthset.v1.variants.for_comparison.tsv.gz

python3 tool_comparison/scripts/add_tool_results_columns.py \
    --tool ExpansionHunter \
    ./tool_comparison/results/expansion_hunter_negative_loci/negative_loci.expansion_hunter.289_json_files.variants.tsv.gz \
    ./tool_comparison/results/negative_loci.for_comparison.tsv.gz

# GangSTR
python3 tool_comparison/scripts/add_tool_results_columns.py \
  --tool GangSTR \
  ./tool_comparison/results/gangstr_positive_loci/positive_loci.gangstr.15_json_files.15_json_files.variants.tsv.gz \
  ./tool_comparison/results/STR_truthset.v1.variants.for_comparison.with_ExpansionHunter_results.tsv.gz

python3 tool_comparison/scripts/add_tool_results_columns.py \
  --tool GangSTR \
  ./tool_comparison/results/gangstr_negative_loci/negative_loci.gangstr.15_json_files.15_json_files.variants.tsv.gz \
  ./tool_comparison/results/negative_loci.for_comparison.with_ExpansionHunter_results.tsv.gz

# concordance
python3 tool_comparison/scripts/add_concordance_columns.py \
  ./tool_comparison/results/STR_truthset.v1.variants.for_comparison.with_ExpansionHunter_results.with_GangSTR_results.tsv.gz

python3 tool_comparison/scripts/add_concordance_columns.py \
  ./tool_comparison/results/negative_loci.for_comparison.with_ExpansionHunter_results.with_GangSTR_results.tsv.gz
