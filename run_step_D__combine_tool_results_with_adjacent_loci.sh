set -x
set -e
set -u

# ExpansionHunter, GangSTR, HipSTR
for results_folder in \
    results_with_adjacent_loci__max_dist_50bp \
    results_with_adjacent_loci__max_dist_24bp \
    results_with_adjacent_loci__max_dist_10bp
do
  # compute STR_truth_set.v1.variants.for_comparison.tsv.gz
  python3 -u tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py \
      --output-dir ./tool_comparison/${results_folder}  \
      STR_truth_set.v1.variants.tsv.gz

  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.variants.for_comparison.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz
  
  #tool=$(echo $results_folder | sed 's/results_folder/ExpansionHunter/g' | sed 's/results_with/ExpansionHunter_with/g')
  tool=ExpansionHunter

  set +x
  echo ============================================
  echo Processing $results_folder $tool results
  set -x

  python3 -u tool_comparison/scripts/add_tool_results_columns.py \
    --tool $tool \
    ./tool_comparison/${results_folder}/combined.positive_loci.*_json_files.variants.tsv.gz \
    ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz

  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.with_${tool}_results.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz


  # add tool vs. truth set concordance columns
  python3 -u tool_comparison/scripts/add_concordance_columns.py \
    --concordance-pair-value1 ${tool} \
    --concordance-pair-value2 Truth \
    ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz

  # rename files
  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.with_concordance.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.variants.tsv.gz
  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.with_concordance.alleles.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.alleles.tsv.gz

  gunzip -f ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.variants.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.alleles.tsv.gz

  # clean up intermediate files
  rm ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz
done

# combine all
python3 -u ./tool_comparison/scripts/combine_all_results_tables_with_adjacent_loci.py

#gsutil -m cp ./tool_comparison/combined.results.*.tsv.gz  gs://str-truth-set/hg38/

echo Done with step D
