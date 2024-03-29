# This script does a fast comparison using bedtools. This is a sanity check on the more nuanced
# comparison done by compute_overlap_with_other_catalogs.py

truthset_bed=${1:-./STR_truth_set.v1.variants.bed.gz}
prefix=$2

# total # of records in the truth set
total=$(gzcat $truthset_bed | grep -v ':0=' | wc -l)
# NOTE: grep -v ':0=' excludes truth set STR variants that are novel (ie. don't have matching STR in their reference
# genome context) since it doesn't make sense to check whether they overlap loci in another catalog.

function subtract {
    label=$1
    catalog_bed=$2


    # t = total records in catalog
    t=$(gzcat $catalog_bed | wc -l)
    # c = records in truth set that don't overlap anything in the catalog.
    c=$(bedtools subtract -a <(gzcat $truthset_bed | grep -v ':0=') -b $catalog_bed -A | wc -l)
    echo "${label} This catalog misses $(printf "%'d" ${c}) out of $(printf "%'d" ${total}) ($(echo scale=1\; 100*"${c}"/"${total}" | bc)%) truth set records. The catalog contains $(printf "%'d" ${t}) records."
}

# process each catalog
subtract "${prefix}             Illumina:" ./ref/other/illumina_variant_catalog.sorted.bed.gz
subtract "${prefix}              GangSTR:" ./ref/other/hg38_ver17.adjusted.bed.gz
subtract "${prefix}               HipSTR:" ./ref/other/hg38.hipstr_reference.adjusted.bed.gz
subtract "${prefix}               PopSTR:" ./ref/other/popstr_catalog_v2.bed.gz
subtract "${prefix}                 TRGT:" ./ref/other/trgt_repeat_catalog.hg38.reformatted_to_motif_only.bed.gz
subtract "${prefix}               adotto:" ./ref/other/adotto_tr_catalog_v1.2.bed.gz

for i in $(seq 6 3 30); do
    subtract "${prefix}   PureRepeats >= ${i}bp, including homopolymers:"  ./ref/other/repeat_specs_GRCh38_without_mismatches.including_homopolymers.sorted.at_least_${i}bp.bed.gz
done
for i in $(seq 6 3 30); do
    subtract "${prefix}   PureRepeats >= ${i}bp:"  ./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_${i}bp.bed.gz
done

