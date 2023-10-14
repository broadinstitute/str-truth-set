set -x
set -e
set -u


cd ./figures_and_tables/
python3 -u plot_tool_accuracy_percent_exactly_right.py --show-title --by-adjacent-loci ../tool_comparison/combined.results.with_adjacent_loci.alleles.tsv.gz

#python3 -u plot_tool_accuracy_by_allele_size.py --compare-loci-with-adjacent-repeats --only-pure-repeats --show-no-call-loci ../tool_comparison/combined.results.with_adjacent_loci.alleles.tsv.gz

echo Done with step E
