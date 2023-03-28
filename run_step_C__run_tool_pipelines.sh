# Run tool pipelines using Hail Batch


force=""
#force="--force"

debug=""
#debug="echo"
if [ -z "$debug" ]; then
  set -ex
fi

function run_pipelines {
    input_bam=$1
    input_bai=$2
    output_dir=$3
    local_dir=$4
    run_illumina_expansion_hunter=$5
    min_locus_coverage_arg=$6
    
    # ExpansionHunterDenovo
    $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_denovo_pipeline.py --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter_denovo &

    # ExpansionHunter
    $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --positive-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter ${min_locus_coverage_arg} ${force} &
    $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --negative-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter ${min_locus_coverage_arg} ${force} &

    if [ "$run_illumina_expansion_hunter" == "yes" ]; then
	      $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --use-illumina-expansion-hunter --loci-to-exclude ./tool_comparison/hail_batch_pipelines/truth_set_loci_that_cause_illumina_expansion_hunter_error.txt \
            --positive-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter ${min_locus_coverage_arg}  ${force} &
	      $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --use-illumina-expansion-hunter --loci-to-exclude ./tool_comparison/hail_batch_pipelines/negative_loci_that_cause_illumina_expansion_hunter_error.txt \
	          --negative-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/expansion_hunter ${min_locus_coverage_arg}  ${force} &
    fi

    # GangSTR
    $debug python3 ./tool_comparison/hail_batch_pipelines/gangstr_pipeline.py  --positive-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/gangstr ${force} &
    $debug python3 ./tool_comparison/hail_batch_pipelines/gangstr_pipeline.py  --negative-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/gangstr ${force} &

    # HipSTR
    $debug python3 ./tool_comparison/hail_batch_pipelines/hipstr_pipeline.py  --positive-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/hipstr ${force} &
    $debug python3 ./tool_comparison/hail_batch_pipelines/hipstr_pipeline.py  --negative-loci --input-bam ${input_bam} --input-bai ${input_bai} --output-dir ${output_dir}/hipstr ${force} &

    wait

    # download results
    $debug mkdir -p ${local_dir}/expansion_hunter_denovo
    $debug gsutil -m cp ${output_dir}/expansion_hunter_denovo/CHM*.locus.tsv ${local_dir}/expansion_hunter_denovo/

    $debug mkdir -p ${local_dir}/expansion_hunter/positive_loci/      ${local_dir}/expansion_hunter/negative_loci/
    $debug gsutil -m cp ${output_dir}/expansion_hunter/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/expansion_hunter/positive_loci/
    $debug gsutil -m cp ${output_dir}/expansion_hunter/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/expansion_hunter/negative_loci/

    $debug mkdir -p ${local_dir}/gangstr/positive_loci/      ${local_dir}/gangstr/negative_loci/
    $debug gsutil -m cp ${output_dir}/gangstr/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/gangstr/positive_loci/
    $debug gsutil -m cp ${output_dir}/gangstr/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/gangstr/negative_loci/

    $debug mkdir -p ${local_dir}/hipstr/positive_loci/      ${local_dir}/hipstr/negative_loci/
    $debug gsutil -m cp ${output_dir}/hipstr/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/hipstr/positive_loci/
    $debug gsutil -m cp ${output_dir}/hipstr/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/hipstr/negative_loci/
}

# Downsample coverage
for coverage in 30 20 10 5; do 
    $debug python3 ./tool_comparison/hail_batch_pipelines/downsample_bam_pipeline.py --verbose --target-coverage ${coverage} \
	--output-dir gs://bw2-delete-after-30-days/ --input-bam gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram &
done


# Run pipelines on original bam
run_pipelines \
  "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram" \
  "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai" \
  "gs://str-truth-set/hg38/tool_results" \
  "./tool_comparison/results" \
  "yes" \
  ""


# Run pipelines on exome data
run_pipelines \
  "gs://broad-public-datasets/CHM1_CHM13_WES/CHMI_CHMI3_Nex1.cram" \
  "gs://broad-public-datasets/CHM1_CHM13_WES/CHMI_CHMI3_Nex1.cram.crai" \
  "gs://str-truth-set/hg38/tool_results_for_exome" \
  "./tool_comparison/results_for_exome" \
  "no" \
  ""

wait   # wait for downsampling to finish

# Process other coverage levels
for coverage in 30 20 10 5; do    
    if [ "$coverage" == "10" ] || [ "$coverage" == "5" ]; then
	min_locus_coverage="--min-locus-coverage 3"
    else
	min_locus_coverage=""
    fi
    
    # Rerun pipelines on downsampled bam
    run_pipelines \
      "gs://bw2-delete-after-30-days/CHM1_CHM13_WGS2.downsampled_to_${coverage}x.bam" \
      "gs://bw2-delete-after-30-days/CHM1_CHM13_WGS2.downsampled_to_${coverage}x.bam.bai" \
      "gs://str-truth-set/hg38/tool_results_for_downsampled_${coverage}x_bam" \
      "./tool_comparison/results_for_downsampled_${coverage}x_bam" \
      "no" \
      "${min_locus_coverage}"
done

# Run REViewer
python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --positive-loci --run-reviewer \
	--input-bam gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram --input-bai gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai \
	--output-dir gs://str-truth-set/hg38/tool_results/expansion_hunter &
python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --negative-loci --run-reviewer \
	--input-bam gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram --input-bai gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai \
	--output-dir gs://str-truth-set/hg38/tool_results/expansion_hunter &
wait

python3 ./tool_comparison/scripts/add_reviewer_image_url_to_bed.py -i ./STR_truth_set.v1.variants.bed.gz
python3 ./tool_comparison/scripts/add_reviewer_image_url_to_bed.py -i ./tool_comparison/variant_catalogs/negative_loci.bed.gz

gsutil -m cp ./STR_truth_set.v1.variants.with_reviewer_image_urls.bed.gz* gs://str-truth-set/hg38/
gsutil -m cp ./tool_comparison/variant_catalogs/negative_loci.with_reviewer_image_urls.bed.gz* gs://str-truth-set/hg38/tool_comparison/variant_catalogs/

echo Done with step C
