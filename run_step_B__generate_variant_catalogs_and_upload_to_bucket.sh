set -ex

# generate variant catalogs of all STR truth set loci
python3 -u tool_comparison/scripts/convert_truth_set_to_variant_catalogs.py \
  --expansion-hunter-loci-per-run 500 \
  --gangstr-loci-per-run 10000 \
  --output-dir ./tool_comparison/variant_catalogs \
  --syndip-high-confidence-regions-bed ./ref/full.38.bed.gz \
  --syndip-indels-vcf ./ref/full.38.INDELs.vcf.gz \
  --all-hg38-repeats-bed ./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz \
  STR_truth_set.v1.variants.tsv.gz

# upload catalogs to bucket
gsutil -q -m cp tool_comparison/variant_catalogs/positive_loci.bed.gz* tool_comparison/variant_catalogs/negative_loci.bed.gz* gs://str-truth-set/hg38/

gsutil -q -m cp -r tool_comparison/variant_catalogs/expansion_hunter  gs://str-truth-set/hg38/variant_catalogs/
gsutil -q -m cp -r tool_comparison/variant_catalogs/gangstr           gs://str-truth-set/hg38/variant_catalogs/
gsutil -q -m cp -r tool_comparison/variant_catalogs/hipstr            gs://str-truth-set/hg38/variant_catalogs/


set +x

gsutil -q -m cp STR_truth_set.v1.*.bed.gz* STR_truth_set.v1.*.tsv* STR_truth_set.v1.vcf*   gs://str-truth-set/hg38/

for i in ref/full.38.* \
       ref/chm13v2.0.fa* \
       ref/hg38.fa* \
       ref/*.chain \
       ref/other/GRCh38GenomicSuperDup.without_decoys.sorted.bed.gz \
       ref/other/GRCh38GenomicSuperDup.without_decoys.sorted.bed.gz.tbi \
       ref/other/MANE.GRCh38.v1.0.ensembl_genomic.gtf.gz \
       ref/other/MANE.v1.0.ensembl_genomic.sorted.gtf.gz \
       ref/other/MANE.v1.0.ensembl_genomic.sorted.gtf.gz.tbi \
       ref/other/gencode.v42.annotation.gtf.gz \
       ref/other/gencode.v42.annotation.sorted.gtf.gz \
       ref/other/gencode.v42.annotation.sorted.gtf.gz.tbi \
       ref/other/hg38_ver17.adjusted.bed.gz \
       ref/other/hg38_ver17.adjusted.bed.gz.tbi \
       ref/other/illumina_variant_catalog.sorted.bed.gz \
       ref/other/illumina_variant_catalog.sorted.bed.gz.tbi \
       ref/other/known_disease_associated_STR_loci.GRCh38.bed.gz \
       ref/other/known_disease_associated_STR_loci.GRCh38.bed.gz.tbi \
       ref/other/repeat_specs_GRCh38_allowing_mismatches.sorted.trimmed.at_least_6bp.bed.gz \
       ref/other/repeat_specs_GRCh38_allowing_mismatches.sorted.trimmed.at_least_6bp.bed.gz.tbi \
       ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_6bp.bed.gz \
       ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_6bp.bed.gz.tbi \
       ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz \
       ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz.tbi \
       ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_12bp.bed.gz \
       ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_12bp.bed.gz.tbi \
       ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_15bp.bed.gz \
       ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_15bp.bed.gz.tbi;
do
    echo ${i}
    gsutil -q -m cp -n ${i} gs://str-truth-set/hg38/${i}
done
