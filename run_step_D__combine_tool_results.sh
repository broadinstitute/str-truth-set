set -x
set -e
set -u

# ExpansionHunter, GangSTR, HipSTR
for results_folder in \
  results_for_exome \
  results \
  results_for_downsampled_30x_bam \
  results_for_downsampled_20x_bam \
  results_for_downsampled_10x_bam \
  results_for_downsampled_5x_bam
do
  # compute STR_truth_set.v1.variants.for_comparison.tsv.gz
  python3 -u tool_comparison/scripts/compute_truth_set_tsv_for_comparisons.py \
      --output-dir ./tool_comparison/${results_folder}/  \
      STR_truth_set.v1.variants.tsv.gz \
      ./tool_comparison/variant_catalogs/negative_loci.tsv.gz

  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.variants.for_comparison.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz
  
  for i in 1 2 3;
  do
      if [ $i == 1 ]; then
        tool=ExpansionHunter
        subfolder=expansion_hunter
      elif [ $i == 2 ]; then
        tool=GangSTR
        subfolder=gangstr
      elif [ $i == 3 ]; then
        tool=HipSTR
        subfolder=hipstr
      else
        echo ERROR: unexpected i == $i
        exit 1;
      fi


      set +x
      echo ============================================
      echo Processing $results_folder $tool results
      set -x
      
      python3 -u tool_comparison/scripts/add_tool_results_columns.py \
	      --tool $tool \
	      ./tool_comparison/${results_folder}/${subfolder}/positive_loci/combined.positive_loci.*_json_files.variants.tsv.gz \
	      ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz

      python3 -u tool_comparison/scripts/add_tool_results_columns.py \
	      --tool $tool \
	      ./tool_comparison/${results_folder}/${subfolder}/negative_loci/combined.negative_loci.*_json_files.variants.tsv.gz \
	      ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz
      
      mv ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.with_${tool}_results.tsv.gz \
         ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz
      mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_${tool}_results.tsv.gz \
         ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz
  done

  # add tool vs. truth set concordance columns
  python3 -u tool_comparison/scripts/add_concordance_columns.py \
    ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz

  python3 -u tool_comparison/scripts/add_concordance_columns.py \
    ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz

  # rename files
  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.with_concordance.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.variants.tsv.gz
  mv ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.with_concordance.alleles.tsv.gz \
     ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.alleles.tsv.gz
  mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_concordance.tsv.gz \
     ./tool_comparison/${results_folder}/negative_loci.for_comparison.variants.tsv.gz
  mv ./tool_comparison/${results_folder}/negative_loci.for_comparison.with_concordance.alleles.tsv.gz \
     ./tool_comparison/${results_folder}/negative_loci.for_comparison.alleles.tsv.gz

  gunzip -f ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.variants.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/negative_loci.for_comparison.variants.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.alleles.tsv.gz
  gunzip -f ./tool_comparison/${results_folder}/negative_loci.for_comparison.alleles.tsv.gz

  # clean up intermediate files
  rm ./tool_comparison/${results_folder}/STR_truth_set.v1.for_comparison.tsv.gz
  rm ./tool_comparison/${results_folder}/negative_loci.for_comparison.tsv.gz
done

# EHdn
python3 -u ./tool_comparison/scripts/intersect_expansion_hunter_denovo_results_with_truth_set.py \
  --window-size 600 \
  ./tool_comparison/results*/expansion_hunter_denovo/CHM1_CHM13_WGS2.*locus.tsv

# combine all
python3 -u ./tool_comparison/scripts/combine_all_results_tables.py
