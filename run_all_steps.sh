set -ex

# Generate 4 different versions of the truth set by toggling whether to include homopolymers and whether
# to always extend locus coordinates to include interruptions (rather than only using pure repeats when possible)
# these versions are useful for comparison and for generating different stats for the paper.
# However, the final truth set that's used for tool comparisons and for most of the paper is the 4th one - which
# doesn't includes homopolymers and doesn't extend locus coordinates to include interruptions.

./run_step_A__create_STR_truth_set.sh --only-high-confidence-regions >& step_A.log  # this version is used as the final truth set
./run_step_A__create_STR_truth_set.sh --only-high-confidence-regions --include-homopolymers >& step_A_with_homopolymers.log
./run_step_A__create_STR_truth_set.sh --only-high-confidence-regions --always-extend-locus-coordinates-to-include-interruptions >& step_A_extended_with_interruptions.log
./run_step_A__create_STR_truth_set.sh --only-high-confidence-regions --include-homopolymers --always-extend-locus-coordinates-to-include-interruptions >& step_A_extended_with_interruptions_step_A_with_homopolymers.log


for _including_homopolymer in ""  "s_including_homopolymer"; do
  for table_type in "variants" "alleles"; do
    python3 scripts/join_truth_set_tables_generated_using_different_options.py \
      --table1 STR${_including_homopolymer}_truth_set.v1.${table_type}.tsv.gz \
      --table2 STR_loci_extended_with_interruption${_including_homopolymer}_truth_set.v1.${table_type}.tsv.gz \
      --table2-suffix ExtendedWithInterruptions \
      --output-table STR_truth_set.v1.${table_type}.joined.tsv.gz \
      -k Chrom -k VcfPos -k VcfRef -k VcfAlt   -c Locus -c SummaryString -c NumRepeatsInReference -c Motif \
      | python3 -u scripts/add_prefix_to_stdout.py "i${_including_homopolymer}s:${table_type}:       "
  done
done


# Run the rest of the pipeline
./run_step_B__generate_variant_catalogs_and_upload_to_bucket.sh >& step_B.log
./run_step_C__run_tool_pipelines.sh >& step_C.log
./run_step_D__combine_tool_results.sh >& step_D.log
./run_step_E__generate_plots_and_figures.sh >& step_E.log

echo Finished running all steps
