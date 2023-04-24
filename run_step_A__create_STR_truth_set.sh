#!/usr/bin/env bash

# Optional command-line args:
#   --only-high-confidence-variants: only include variants that are within SynDip high confidence regions
#   --include-homopolymers: include homopolymer variants in the truth set
#   --only-pure-repeats: only include STRs that don't contain interruptions
#   --always-extend-locus-coordinates-to-include-interruptions This option is mutually-exclusive with
#       --only-pure-repeats and tells the STR filter to extend locus coordinates to include the interruptions even
#       when the variant passes the filter based on only pure repeats. By default, the STR filter will only looks for
#       pure repeats in the reference context in these situations.


allow_interruptions="only-if-pure-repeats-not-found"
always_look_for_interrupted_repeats="false"
only_high_confidence_regions="false"
include_homopolymers="false"

if [[ " $@ " =~ " --only-pure-repeats " ]]; then
    allow_interruptions="no"
elif [[ " $@ " =~ " --always-extend-locus-coordinates-to-include-interruptions " ]]; then
    allow_interruptions="always"
fi

if [[ " $@ " =~ " --only-high-confidence-regions " ]]; then
    only_high_confidence_regions="true"
fi
if [[ " $@ " =~ " --include-homopolymers " ]]; then
    include_homopolymers="true"
fi


# Check that required command-line tools are installed
for command in "python3" "curl" "bgzip" "tabix" "bedtools" "gatk"
do 
    if ! command -v "${command}" &> /dev/null; then
	echo "${command} command not found. Please install it and add it to your PATH.";
	exit 1
    fi
done


# Required data (will be downloaded if it's not available)
#    reference_data_bundle.tar.gz from github.com/
#    hg38.fa
#    chm13v2.0.fa  (T2T reference fasta)

set -e
set -u
set -o pipefail

version=v1

# Reference data paths
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
  curl --silent -L https://storage.googleapis.com/str-truth-set/STR_truth_set_reference_bundle.v1.tar.gz -o STR_truth_set_reference_bundle.v1.tar.gz
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
  curl --silent -L https://github.com/lh3/CHM-eval/releases/download/v0.5/CHM-evalkit-20180222.tar \
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
    echo "Downloading hg38 reference genome: ${hg38_fasta_path}"
    set -x
    curl --silent -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz \
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
    curl --silent -L https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz \
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
  local prefix="$3"
  echo "        " OUTPUT: $output_vcf "       " $(gunzip -c $output_vcf | grep -v ^# | wc -l) variants
  python3 scripts/vcf_stats.py --prefix "${prefix}:input:" --min-percent 0 --label "the input vcf" $input_vcf
  python3 scripts/vcf_stats.py --prefix "${prefix}:output:" --min-percent 0 --label "the output vcf" $output_vcf
}

function print_liftover_output_stats {
  print_output_stats "$1" "$2" "$4"
  local output_failed_liftover_vcf="$3"
  local prefix="$4"
  python3 scripts/vcf_stats.py --count-by-filter --prefix "${prefix}:failed-liftover:" --min-percent 0 "$output_failed_liftover_vcf"

  echo "${prefix}:failed-liftover:          Reasons for Liftover failure:"

  # print more stats (disable pipefail in case vcf is empty)
  set +eo pipefail

  OLDIFS="$IFS"
  IFS=$'\n'
  for l in $(gunzip -c "$output_failed_liftover_vcf" | grep -v ^# | cut -f 7  | sort | uniq -c | sort -n -r); do
      echo "${prefix}:failed-liftover:          ${l}"
  done
  IFS="$OLDIFS"

  set -eo pipefail
}


set -euo pipefail


# NOTE: these parameter values were chosen so that the subset of known disease-associated STR loci that have
# non-reference genotypes in the SynDip truth set would be added to the STR truth set and would have the expected
# reference start and end coordinates despite some of the repeat interruptions in the hg38 reference sequence at some
# of these loci

min_str_length=9
min_str_repeats=3
min_repeat_unit_length=2
max_repeat_unit_length=50


if [ $allow_interruptions == "no" ]; then
  STR_type="pure_STR"
elif [ $allow_interruptions == "only-if-pure-repeats-not-found" ]; then
  STR_type="STR"
elif [ $allow_interruptions == "always" ]; then
  STR_type="STR_loci_extended_with_interruption"
fi

allow_interruptions_arg="--allow-interruptions ${allow_interruptions}"

if [ $include_homopolymers == "true" ]; then
  STR_type="${STR_type}s_including_homopolymer"
  min_repeat_unit_length=1
fi

if [ $only_high_confidence_regions == "false" ]; then
  STR_type="raw_${STR_type}"
fi


echo ===============
input_vcf=${syndip_truth_vcf}

# generate stats and summary files for the raw SynDip VCF before filtering to high-confidence regions. These files aren't used in downstream analysis.
python3 scripts/get_indels_from_vcf.py ${input_vcf}  | python3 -u scripts/add_prefix_to_stdout.py "step0:${STR_type}:  "

if [ $only_high_confidence_regions == "true" ]
then
  output_vcf=step1.high_confidence_regions.vcf.gz

  echo ===============
  print_input_stats $input_vcf "STEP #1: Filter SynDip truth set to SynDip high confidence regions"
  set -x

  tabix -f ${syndip_confidence_regions_bed}

  bedtools intersect -header -f 1 -wa -u \
      -a $input_vcf  \
      -b ${syndip_confidence_regions_bed} \
      | bgzip > ${output_vcf}

  tabix -f ${output_vcf}

  set +x

  print_output_stats ${input_vcf} ${output_vcf} "step1"

  python3 scripts/get_indels_from_vcf.py ${output_vcf}  | python3 -u scripts/add_prefix_to_stdout.py "step1:${STR_type}:  "

  input_vcf=${output_vcf}
fi

output_prefix=step2.${STR_type}s
output_vcf="${output_prefix}.vcf.gz"

echo ===============
echo Starting to process ${STR_type}s: ${allow_interruptions_arg} \
  --min-repeat-unit-length ${min_repeat_unit_length} --max-repeat-unit-length ${max_repeat_unit_length}


print_input_stats $input_vcf "STEP #2: Filter variants to the subset that are actually STR expansion or contractions"
set -x

python3 -u -m str_analysis.filter_vcf_to_STR_variants ${allow_interruptions_arg} \
  -R "${hg38_fasta_path}" \
  --write-bed-file \
  --min-str-length "${min_str_length}" \
  --min-str-repeats "${min_str_repeats}" \
  --min-repeat-unit-length ${min_repeat_unit_length} \
  --max-repeat-unit-length ${max_repeat_unit_length} \
  --output-prefix "${output_prefix}" \
  --write-vcf-with-filtered-out-variants \
  --verbose \
  "${input_vcf}" \
  | python3 -u scripts/add_prefix_to_stdout.py "step2:${STR_type}:  "

#   -n 10000 \

# generate overlap statistics before liftover
#python3 -u scripts/compute_overlap_with_other_catalogs.py --all-repeats \
#  --all-motifs ${output_prefix}.variants.tsv.gz \
#  ${output_prefix}.variants.with_overlap_columns.tsv.gz \
#  -c IlluminaSTRCatalog -c GangSTRCatalog17 -c HipSTRCatalog -c KnownDiseaseAssociatedSTRs \
#  | python3 -u scripts/add_prefix_to_stdout.py "step2: "

set +x
print_output_stats $input_vcf $output_vcf "step2:${STR_type}"


echo ===============
input_vcf=$output_vcf
converted_vcf=step3.${STR_type}s.with_dels_converted_for_liftover.vcf.gz
output_vcf=step3.${STR_type}s.lifted_to_chm13v2.vcf.gz
output_failed_liftover1_vcf=step3.${STR_type}s.lifted_to_chm13v2_rejected.vcf.gz

print_input_stats $input_vcf "STEP #3: Liftover variants from hg38 to the T2T reference (chm13v2.0)"
set -x

python3 scripts/convert_monoallelic_deletions_to_snvs_for_liftover.py \
  ${input_vcf} \
  ${converted_vcf} \
  | python3 -u scripts/add_prefix_to_stdout.py "step3:${STR_type}:  "

gatk --java-options "-Xmx7g" LiftoverVcf \
    -R "${t2t_fasta_path}" \
    -C "${hg38_t2t_chain_path}" \
    -I       "${converted_vcf}" \
    -O       "${output_vcf}" \
    --REJECT "${output_failed_liftover1_vcf}" \
    --VALIDATION_STRINGENCY LENIENT \
    --RECOVER_SWAPPED_REF_ALT \
    --WARN_ON_MISSING_CONTIG \
    --LIFTOVER_MIN_MATCH 0.00001 \
    --WRITE_ORIGINAL_POSITION \
    --WRITE_ORIGINAL_ALLELES \
    --ALLOW_MISSING_FIELDS_IN_HEADER \
    | python3 -u scripts/add_prefix_to_stdout.py "step3:${STR_type}:  "

set +x
print_liftover_output_stats $input_vcf $output_vcf $output_failed_liftover1_vcf  "step3:${STR_type}"

echo ===============
input_vcf=$output_vcf
output_vcf=step4.${STR_type}s.matched_t2t.vcf.gz
didnt_match_t2t_vcf=step4.${STR_type}s.didnt_match_t2t.vcf.gz
print_input_stats ${input_vcf} "STEP #4: Filter out variants where neither allele matches the T2T reference (chm13v2.0)"
set -x

python3 -u scripts/filter_out_discordant_variants_after_liftover.py --reference-fasta ${t2t_fasta_path} ${input_vcf} ${output_vcf} ${didnt_match_t2t_vcf} \
    | python3 -u scripts/add_prefix_to_stdout.py "step4:${STR_type}: "

python3 -u scripts/convert_vcf_to_tsv.py ${didnt_match_t2t_vcf} \
    | python3 -u scripts/add_prefix_to_stdout.py "step4:${STR_type}: "

set +x
print_output_stats ${input_vcf} ${output_vcf} "step4:${STR_type}"

echo ===============
input_vcf=${output_vcf}
output_vcf=step5.${STR_type}s.lifted_back_to_38.vcf.gz
output_failed_liftover2_vcf=step5.${STR_type}s.lifted_back_to_38_rejected.vcf.gz

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
     --ALLOW_MISSING_FIELDS_IN_HEADER \
     | python3 -u scripts/add_prefix_to_stdout.py "step5:${STR_type}: "

set +x
print_liftover_output_stats $input_vcf $output_vcf $output_failed_liftover2_vcf  "step5:${STR_type}"

echo ===============
input_vcf=$output_vcf
output_vcf=step6.${STR_type}s.restored_dels_that_failed_liftover.vcf.gz

print_input_stats $input_vcf "STEP #6: Restore variants that were dropped in the 1st liftover (hg38 => chm13v2.0) due to IndelStraddlesMultipleIntevals flag in picard LiftoverVcf. See https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/vcf/LiftoverVcf.java#L424-L431 for more details on IndelStraddlesMultipleIntevals"
set -x

gunzip -c $input_vcf > temp.vcf
gunzip -c $output_failed_liftover1_vcf | grep -v '^#' | grep IndelStraddlesMultipleIntevals >> temp.vcf

gatk FixVcfHeader -I temp.vcf -O temp.fixed.vcf
gatk SortVcf -I temp.fixed.vcf -O $output_vcf

rm temp.vcf temp.fixed.vcf*

set +x
print_output_stats $input_vcf $output_vcf "step6:${STR_type}"

echo ===============
input_vcf=$output_vcf
output_vcf=step7.${STR_type}s.passed_liftover_checks.vcf.gz

print_input_stats $input_vcf "STEP #7: Check VCF positions before vs. after liftover to confirm concordance."
set -x

python3 -u scripts/check_vcf_concordance_before_vs_after_liftover.py \
  --log-prefix "step7:${STR_type}" \
  -o $output_vcf \
  step2.${STR_type}s.vcf.gz \
  step6.${STR_type}s.restored_dels_that_failed_liftover.vcf.gz

set +x
print_output_stats $input_vcf $output_vcf "step7:${STR_type}"

echo ===============
input_vcf=$output_vcf
output_prefix=step7.filtered_${STR_type}s

print_input_stats $input_vcf "step8: Print stats and compute overlap with other catalogs."
set -x

# even though all the variants in $input_vcf are already STRs, run the str_analysis.filter_vcf_to_STR_variants script on
# it just to generate the .variants.tsv and .alleles.tsv tables and print summary stats.
python3 -u -m str_analysis.filter_vcf_to_STR_variants ${allow_interruptions_arg} \
  -R "${hg38_fasta_path}" \
  --write-bed-file \
  --min-str-length "${min_str_length}" \
  --min-str-repeats "${min_str_repeats}" \
  --min-repeat-unit-length ${min_repeat_unit_length} \
  --max-repeat-unit-length ${max_repeat_unit_length} \
  --copy-info-field-keys-to-tsv SkippedValidation \
  --copy-info-field-keys-to-tsv NumRepeatsInT2T \
  --output-prefix "${output_prefix}" \
  --verbose \
  "${input_vcf}" | python3 -u scripts/add_prefix_to_stdout.py "step8:${STR_type}:  "


# compute overlap with various reference annotations
./scripts/compute_overlap_with_other_catalogs_using_bedtools.sh ${output_prefix}.variants.bed.gz "step8:overlap:${STR_type}"

suffix=with_overlap_columns
python3 -u scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs ${output_prefix}.variants.tsv.gz  ${output_prefix}.variants.${suffix}.tsv.gz | python3 -u scripts/add_prefix_to_stdout.py "step8:${STR_type}:  "
python3 -u scripts/compute_overlap_with_other_catalogs.py --all-repeats --all-motifs ${output_prefix}.alleles.tsv.gz  ${output_prefix}.alleles.${suffix}.tsv.gz | python3 -u scripts/add_prefix_to_stdout.py "step8:${STR_type}: "
mv ${output_prefix}.variants.${suffix}.tsv.gz ${output_prefix}.variants.tsv.gz
mv ${output_prefix}.alleles.${suffix}.tsv.gz  ${output_prefix}.alleles.tsv.gz

suffix=with_gencode_v42_columns
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/gencode.v42.annotation.sorted.gtf.gz  ${output_prefix}.variants.tsv.gz  ${output_prefix}.variants.${suffix}.tsv.gz | python3 -u scripts/add_prefix_to_stdout.py "step8:${STR_type}:  "
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/gencode.v42.annotation.sorted.gtf.gz  ${output_prefix}.alleles.tsv.gz   ${output_prefix}.alleles.${suffix}.tsv.gz | python3 -u scripts/add_prefix_to_stdout.py "step8:${STR_type}:  "
mv ${output_prefix}.variants.${suffix}.tsv.gz ${output_prefix}.variants.tsv.gz
mv ${output_prefix}.alleles.${suffix}.tsv.gz ${output_prefix}.alleles.tsv.gz

suffix=with_MANE_columns
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/MANE.v1.0.ensembl_genomic.sorted.gtf.gz  ${output_prefix}.variants.tsv.gz ${output_prefix}.variants.${suffix}.tsv.gz | python3 -u scripts/add_prefix_to_stdout.py "step8:${STR_type}:  "
python3 -u scripts/compute_overlap_with_gene_models.py ./ref/other/MANE.v1.0.ensembl_genomic.sorted.gtf.gz  ${output_prefix}.alleles.tsv.gz  ${output_prefix}.alleles.${suffix}.tsv.gz | python3 -u scripts/add_prefix_to_stdout.py "step8:${STR_type}:  "
mv ${output_prefix}.variants.${suffix}.tsv.gz ${output_prefix}.variants.tsv.gz
mv ${output_prefix}.alleles.${suffix}.tsv.gz ${output_prefix}.alleles.tsv.gz

# move files to final output filenames
final_output_prefix=${STR_type}_truth_set.${version}

mv ${output_prefix}.variants.tsv.gz ${final_output_prefix}.variants.tsv.gz
mv ${output_prefix}.alleles.tsv.gz  ${final_output_prefix}.alleles.tsv.gz

mv ${output_prefix}.variants.bed.gz     ${final_output_prefix}.variants.bed.gz
mv ${output_prefix}.variants.bed.gz.tbi ${final_output_prefix}.variants.bed.gz.tbi

mv ${output_prefix}.vcf.gz      ${final_output_prefix}.vcf.gz
mv ${output_prefix}.vcf.gz.tbi  ${final_output_prefix}.vcf.gz.tbi

set +x

echo ===============
echo Done with step A
