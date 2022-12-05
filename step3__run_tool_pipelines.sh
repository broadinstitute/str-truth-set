# Run tool pipelines using Hail Batch

set -ex

function run_pipelines {
    input_bam=$1
    output_dir=$2
    local_dir=$3

    # ExpansionHunterDenovo
    python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_denovo_pipeline.py --input-bam ${input_bam} --output-dir ${output_dir}/expansion_hunter_denovo &
    
    # ExpansionHunter
    python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --positive-loci --input-bam ${input_bam} --output-dir ${output_dir}/expansion_hunter &
    python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --negative-loci --input-bam ${input_bam} --output-dir ${output_dir}/expansion_hunter &

    # GangSTR
    python3 ./tool_comparison/hail_batch_pipelines/gangstr_pipeline.py  --positive-loci --input-bam ${input_bam} --output-dir ${output_dir}/gangstr &
    python3 ./tool_comparison/hail_batch_pipelines/gangstr_pipeline.py  --negative-loci --input-bam ${input_bam} --output-dir ${output_dir}/gangstr &

    wait

    # download results
    mkdir -p ${local_dir}/expansion_hunter_denovo
    gsutil -m cp ${output_dir}/expansion_hunter_denovo/CHM1_CHM13_2*.locus.tsv ${local_dir}/expansion_hunter_denovo/

    mkdir -p ${local_dir}/expansion_hunter/positive_loci/      ${local_dir}/expansion_hunter/negative_loci/
    gsutil -m cp ${output_dir}/expansion_hunter/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/expansion_hunter/positive_loci/
    gsutil -m cp ${output_dir}/expansion_hunter/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/expansion_hunter/negative_loci/

    mkdir -p ${local_dir}/gangstr/positive_loci/      ${local_dir}/gangstr/negative_loci/
    gsutil -m cp ${output_dir}/gangstr/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/gangstr/positive_loci/
    gsutil -m cp ${output_dir}/gangstr/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/gangstr/negative_loci/
}


# Run pipelines on original bam
run_pipelines \
  "gs://str-truth-set/hg38/CHM1_CHM13_2.bam" \
  "gs://str-truth-set/hg38/tool_results" \
  "./tool_comparison/results"


# Downsample to 30x
python3 ./tool_comparison/hail_batch_pipelines/downsample_bam_pipeline.py --verbose --target-coverage 30

# Rerun pipelines on downsampled bam
run_pipelines \
  "gs://str-truth-set/hg38/CHM1_CHM13_2.downsampled_to_30x.bam" \
  "gs://str-truth-set/hg38/tool_results_for_downsampled_30x_bam" \
  "./tool_comparison/results_for_downsampled_30x_bam"


# Downsample to 20x
python3 ./tool_comparison/hail_batch_pipelines/downsample_bam_pipeline.py --verbose --target-coverage 20

# Rerun pipelines on downsampled bam
run_pipelines \
  "gs://str-truth-set/hg38/CHM1_CHM13_2.downsampled_to_20x.bam" \
  "gs://str-truth-set/hg38/tool_results_for_downsampled_20x_bam" \
  "./tool_comparison/results_for_downsampled_20x_bam"


# Run pipelines on exome data
run_pipelines \
    "gs://bw-proj/CHMI_CHMI3_Nex1.cram" \
    "gs://str-truth-set/hg38/tool_results_for_exome" \
    "./tool_comparison/results_for_exome"

