#!/usr/bin/env bash

# Required command-line tools:
#
#   bedtools
#   gatk
#   python3
#   bgzip
#   tabix
#   curl

# Required data (will be downloaded if it's not available)
#    reference_data_bundle.tar.gz from github.com/
#    hg38.fa
#    chm13v2.0.fa  (T2T reference fasta)

set -u
set -o pipefail

version=v1

syndip_truth_vcf=./ref/full.38.vcf.gz
syndip_confidence_regions_bed=./ref/full.38.bed.gz

hg38_fasta_path=./ref/hg38.fa
t2t_fasta_path=./ref/chm13v2.0.fa

hg38_t2t_chain_path=./ref/hg38-chm13v2.chain
t2t_hg38_chain_path=./ref/chm13v2-hg38.chain

# Install python dependencies
set -x
python3 -m pip install -r requirements.txt -qq
set +x

# Download the reference data bundle if necessary
if [ ! -d ./ref ] || [ ! -d ./ref/other ] || [ ! -f ${hg38_t2t_chain_path} ] || [ ! -f ${t2t_hg38_chain_path} ]
then
  echo "Downloading the STR truth set reference data bundle..."
  set -x
  mkdir -p ./ref
  curl -L https://github.com/broadinstitute/str-truth-set/releases/download/v1/STR_truth_set_reference_bundle.v1.tar.gz -o STR_truth_set_reference_bundle.v1.tar.gz
  tar xzf STR_truth_set_reference_bundle.v1.tar.gz
  rm STR_truth_set_reference_bundle.v1.tar.gz
  set +x
fi

# If it doesn't exist yet, the hg38 truth VCF and confidence regions from the Synthetic Diploid Benchmark (SynDip) [Li 2018]
#     Additional information about SynDip is @ https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6341484 and
#     https://github.com/lh3/CHM-eval
if [ ! -f ${syndip_truth_vcf} ] || [ ! -f ${syndip_confidence_regions_bed} ]
then
  echo "Downloading the SynDip truth vcf and confidence regions..."
  set -x
  curl -L https://github.com/lh3/CHM-eval/releases/download/v0.4/CHM-evalkit-20180221.tar \
    | tar xf - CHM-eval.kit/full.38.vcf.gz CHM-eval.kit/full.38.bed.gz
  mv CHM-eval.kit/full.38.vcf.gz ${syndip_truth_vcf}
  mv CHM-eval.kit/full.38.bed.gz ${syndip_confidence_regions_bed}

  rmdir CHM-eval.kit/

  tabix ${syndip_truth_vcf}
  tabix ${syndip_confidence_regions_bed}

  set +x
fi

# Download the hg38 and T2T reference fastas if necessary
if [ ! -f "${hg38_fasta_path}" ]
then
    echo "$Downloading hg38 reference genome: {hg38_fasta_path}"
    set -x
    curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz \
      | gunzip -c - > "${hg38_fasta_path}"
    gunzip -c
    samtools faidx "${hg38_fasta_path}"
    gatk CreateSequenceDictionary R=${hg38_fasta_path} O=${hg38_fasta_path}.dict
    set +x
fi

if [ ! -f "${t2t_fasta_path}" ]
then
    echo "Downloading T2T reference genome: ${t2t_fasta_path}"
    set -x
    curl -L https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz \
      | gunzip -c - > "${t2t_fasta_path}"
    samtools faidx "${t2t_fasta_path}"
    gatk CreateSequenceDictionary R=${t2t_fasta_path} O=${t2t_fasta_path}.dict
    set +x
fi


function print_input_stats {
  local input_vcf="$1"
  local step_description="$2"
  echo "$step_description"
  echo "        " INPUT: $input_vcf  "      " $(gunzip -c $input_vcf | grep -v ^# | wc -l) variants
  #python3 scripts/vcf_stats.py --label "the input vcf" $input_vcf
}


function print_output_stats {
  local input_vcf="$1"
  local output_vcf="$2"
  echo "        " OUTPUT: $output_vcf "       " $(gunzip -c $output_vcf | grep -v ^# | wc -l) variants
  python3 scripts/vcf_stats.py --min-percent 0 --label "the input vcf" $input_vcf
  python3 scripts/vcf_stats.py --min-percent 0 --label "the output vcf" $output_vcf
}

function print_liftover_output_stats {
  print_output_stats "$1" "$2"
  local output_failed_liftover_vcf="$3"
  python3 scripts/vcf_stats.py --min-percent 0 "$output_failed_liftover_vcf"

  echo Reasons for Liftover failure:

  # print more stats (disable pipefail in case vcf is empty)
  set +eo pipefail
  gunzip -c "$output_failed_liftover_vcf" | grep -v ^# | cut -f 7  | sort | uniq -c
  set -eo pipefail
}


set -euo pipefail

echo ===============
input_vcf=${syndip_truth_vcf}
output_vcf=step1.high_confidence_regions.vcf.gz

print_input_stats $input_vcf "STEP #1: Filter SynDip truth set to SynDip high confidence regions"
set -x

bedtools intersect -header -f 1 -wa -u \
    -a $input_vcf  \
    -b ${syndip_confidence_regions_bed} \
    | bgzip > ${output_vcf}

set +x
print_output_stats ${input_vcf} ${output_vcf}

echo ===============
input_vcf="${output_vcf}"
output_prefix=step2.STRs
output_vcf="${output_prefix}.vcf.gz"

# NOTE: these parameter values were chosen so that the subset of known disease-associated STR loci that have
# non-reference genotypes in the SynDip truth set would be added to the STR truth set and would have the expected
# reference start and end coordinates despite some of the repeat interruptions in the hg38 reference sequence at some
# of these loci

min_str_length=9
min_str_repeats=3

print_input_stats $input_vcf "STEP #2: Filter variants to the subset that are actually STR expansion or contractions"
set -x

python3 -u -m str_analysis.filter_vcf_to_STR_variants \
  -R "${hg38_fasta_path}" \
  --write-bed-file \
  --allow-interruptions \
  --min-str-length "${min_str_length}" \
  --min-str-repeats "${min_str_repeats}" \
  --output-prefix "${output_prefix}" \
  "${input_vcf}"

#   -n 10000 \


# generate overlap statistics before liftover
python3 -u scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs ${output_prefix}.variants.tsv.gz  ${output_prefix}.variants.with_overlap_columns.tsv.gz -c IlluminaSTRCatalog -c GangSTRCatalog17 -c KnownDiseaseAssociatedSTRs

set +x
print_output_stats $input_vcf $output_vcf


echo ===============
input_vcf=$output_vcf
output_vcf=step3.STRs.lifted_to_chm13v2.vcf.gz
output_failed_liftover1_vcf=step3.lifted_to_chm13v2_rejected.vcf.gz

print_input_stats $input_vcf "STEP #3: Liftover variants from hg38 to the T2T reference (chm13v2.0)"
set -x

gatk --java-options "-Xmx7g" LiftoverVcf \
    -R "${t2t_fasta_path}" \
    -C "${hg38_t2t_chain_path}" \
    -I       "${input_vcf}" \
    -O       "${output_vcf}" \
    --REJECT "${output_failed_liftover1_vcf}" \
    --VALIDATION_STRINGENCY LENIENT \
    --RECOVER_SWAPPED_REF_ALT \
    --WARN_ON_MISSING_CONTIG \
    --LIFTOVER_MIN_MATCH 0.00001 \
    --WRITE_ORIGINAL_POSITION \
    --WRITE_ORIGINAL_ALLELES \
    --ALLOW_MISSING_FIELDS_IN_HEADER

set +x
print_liftover_output_stats $input_vcf $output_vcf $output_failed_liftover1_vcf

echo ===============
input_vcf=$output_vcf
output_vcf=step4.STRs.found_in_chm13v2.vcf.gz
print_input_stats ${input_vcf} "STEP #4: Filter out variants that are true homozygous alt or multi-allelic after liftover since this means neither allele matches the T2T reference (chm13v2.0)"
set -x

python3 -u scripts/filter_out_discordant_variants_after_liftover.py --reference-fasta ${t2t_fasta_path} ${input_vcf} ${output_vcf}

set +x
print_output_stats ${input_vcf} ${output_vcf}

echo ===============
input_vcf=${output_vcf}
output_vcf=step5.STRs.lifted_back_to_38.vcf.gz
output_failed_liftover2_vcf=step5.lifted_back_to_38_rejected.vcf.gz

print_input_stats $input_vcf "STEP #5: Liftover the truth variants that passed all checks back to hg38"
set -x

gatk --java-options "-Xmx7g" LiftoverVcf \
     -R "${hg38_fasta_path}" \
     -C "${t2t_hg38_chain_path}" \
     -I       ${input_vcf} \
     -O       ${output_vcf} \
     --REJECT $output_failed_liftover2_vcf \
     --VALIDATION_STRINGENCY LENIENT \
     --RECOVER_SWAPPED_REF_ALT \
     --WARN_ON_MISSING_CONTIG \
     --LIFTOVER_MIN_MATCH 0.00001 \
     --WRITE_ORIGINAL_POSITION \
     --WRITE_ORIGINAL_ALLELES \
     --ALLOW_MISSING_FIELDS_IN_HEADER

set +x
print_liftover_output_stats $input_vcf $output_vcf $output_failed_liftover2_vcf

echo ===============
input_vcf=$output_vcf
output_vcf=step6.STRs.restored_dels_that_failed_liftover.vcf.gz

print_input_stats $input_vcf "STEP #6: Restore deletions that were dropped in the 1st liftover (hg38 => chm13v2.0) due to IndelStraddlesMultipleIntevals flag in picard LiftoverVcf. See https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/vcf/LiftoverVcf.java#L424-L431 for more details on IndelStraddlesMultipleIntevals"
set -x

gunzip -c $input_vcf > temp.vcf
gunzip -c $output_failed_liftover1_vcf | grep -v '^#' | grep IndelStraddlesMultipleIntevals | awk '{ if ( length($4) > length($5) || index($5, ",") > 0) { print $0 } }' >> temp.vcf

gatk FixVcfHeader -I temp.vcf -O temp.fixed.vcf
gatk SortVcf -I temp.fixed.vcf -O $output_vcf

rm temp.vcf temp.fixed.vcf*

set +x
print_output_stats $input_vcf $output_vcf

echo ===============
input_vcf=$output_vcf
output_vcf=step7.STRs.passed_liftover_checks.vcf.gz

print_input_stats $input_vcf "STEP #7: Check VCF positions before vs. after liftover to confirm concordance."
set -x

python3 -u scripts/check_vcf_concordance_before_vs_after_liftover.py \
  -o $output_vcf \
  step2.STRs.vcf.gz \
  step6.STRs.restored_dels_that_failed_liftover.vcf.gz

set +x
print_output_stats $input_vcf $output_vcf

echo ===============
input_vcf=$output_vcf
output_prefix=step7.filtered_STRs

print_input_stats $input_vcf "step8: Print stats and compute overlap with other catalogs."
set -x

# even though all the variants in $input_vcf are already STRs, run the str_analysis.filter_vcf_to_STR_variants script on
# it just to generate the .variants.tsv and .alleles.tsv tables and print summary stats.
python3 -u -m str_analysis.filter_vcf_to_STR_variants \
  -R "${hg38_fasta_path}" \
  --write-bed-file \
  --allow-interruptions \
  --min-str-length "${min_str_length}" \
  --min-str-repeats "${min_str_repeats}" \
  --output-prefix "${output_prefix}" \
  "${input_vcf}"


# compute overlap with various reference annotations
suffix=with_overlap_columns
python3 -u scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs ${output_prefix}.variants.tsv.gz  ${output_prefix}.variants.${suffix}.tsv.gz &
python3 -u scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs ${output_prefix}.alleles.tsv.gz  ${output_prefix}.alleles.${suffix}.tsv.gz &
wait
mv ${output_prefix}.variants.${suffix}.tsv.gz ${output_prefix}.variants.tsv.gz
mv ${output_prefix}.alleles.${suffix}.tsv.gz  ${output_prefix}.alleles.tsv.gz

suffix=with_gencode_v42_columns
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/gencode.v42.annotation.gtf.gz  ${output_prefix}.variants.tsv.gz  ${output_prefix}.variants.${suffix}.tsv.gz &
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/gencode.v42.annotation.gtf.gz  ${output_prefix}.alleles.tsv.gz   ${output_prefix}.alleles.${suffix}.tsv.gz &
wait
mv ${output_prefix}.variants.${suffix}.tsv.gz ${output_prefix}.variants.tsv.gz
mv ${output_prefix}.alleles.${suffix}.tsv.gz ${output_prefix}.alleles.tsv.gz

suffix=with_MANE_columns
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/MANE.v1.0.ensembl_genomic.gtf.gz  ${output_prefix}.variants.tsv.gz ${output_prefix}.variants.${suffix}.tsv.gz &
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/MANE.v1.0.ensembl_genomic.gtf.gz  ${output_prefix}.alleles.tsv.gz  ${output_prefix}.alleles.${suffix}.tsv.gz &
wait
mv ${output_prefix}.variants.${suffix}.tsv.gz ${output_prefix}.variants.tsv.gz
mv ${output_prefix}.alleles.${suffix}.tsv.gz ${output_prefix}.alleles.tsv.gz


# move files to final output filenames
final_output_prefix=STR_truthset.${version}

mv ${output_prefix}.variants.tsv.gz ${final_output_prefix}.variants.tsv.gz
mv ${output_prefix}.alleles.tsv.gz  ${final_output_prefix}.alleles.tsv.gz

mv ${output_prefix}.variants.bed.gz     ${final_output_prefix}.variants.bed.gz
mv ${output_prefix}.variants.bed.gz.tbi ${final_output_prefix}.variants.bed.gz.tbi

mv ${output_prefix}.vcf.gz      ${final_output_prefix}.vcf.gz
mv ${output_prefix}.vcf.gz.tbi  ${final_output_prefix}.vcf.gz.tbi

python3 tool_comparison/scripts/convert_truth_set_to_variant_catalogs.py \
	--expansion-hunter-loci-per-run 500 \
	--gangstr-loci-per-run 10000 \
	--output-dir ./tool_comparison/variant_catalogs \
	--high-confidence-regions-bed ./ref/full.38.bed.gz \
	--all-repeats-bed ./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz \
	STR_truthset.v1.variants.tsv.gz


# compute overlap with various reference annotations for negative loci
negative_loci_bed_path=tool_comparison/variant_catalogs/negative_loci.bed.gz
output_prefix=tool_comparison/variant_catalogs/negative_loci

suffix=with_overlap_columns
python3 -u scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs ${negative_loci_bed_path}  ${output_prefix}.${suffix}.tsv.gz
mv ${output_prefix}.${suffix}.tsv.gz  ${output_prefix}.tsv.gz

suffix=with_gencode_v42_columns
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/gencode.v42.annotation.gtf.gz  ${output_prefix}.tsv.gz  ${output_prefix}.${suffix}.tsv.gz 
mv ${output_prefix}.${suffix}.tsv.gz  ${output_prefix}.tsv.gz

suffix=with_MANE_columns
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/MANE.v1.0.ensembl_genomic.gtf.gz  ${output_prefix}.tsv.gz ${output_prefix}.${suffix}.tsv.gz
mv ${output_prefix}.${suffix}.tsv.gz  ${output_prefix}.tsv.gz

set +x

echo ===============

