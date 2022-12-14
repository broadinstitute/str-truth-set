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

    # HipSTR
    python3 ./tool_comparison/hail_batch_pipelines/hipstr_pipeline.py  --positive-loci --input-bam ${input_bam} --output-dir ${output_dir}/hipstr &
    python3 ./tool_comparison/hail_batch_pipelines/hipstr_pipeline.py  --negative-loci --input-bam ${input_bam} --output-dir ${output_dir}/hipstr &    

    wait

    # download results
    mkdir -p ${local_dir}/expansion_hunter_denovo
    gsutil -m cp ${output_dir}/expansion_hunter_denovo/CHM*.locus.tsv ${local_dir}/expansion_hunter_denovo/

    mkdir -p ${local_dir}/expansion_hunter/positive_loci/      ${local_dir}/expansion_hunter/negative_loci/
    gsutil -m cp ${output_dir}/expansion_hunter/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/expansion_hunter/positive_loci/
    gsutil -m cp ${output_dir}/expansion_hunter/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/expansion_hunter/negative_loci/

    mkdir -p ${local_dir}/gangstr/positive_loci/      ${local_dir}/gangstr/negative_loci/
    gsutil -m cp ${output_dir}/gangstr/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/gangstr/positive_loci/
    gsutil -m cp ${output_dir}/gangstr/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/gangstr/negative_loci/

    mkdir -p ${local_dir}/hipstr/positive_loci/      ${local_dir}/hipstr/negative_loci/
    gsutil -m cp ${output_dir}/hipstr/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz ${local_dir}/hipstr/positive_loci/
    gsutil -m cp ${output_dir}/hipstr/negative_loci/combined.negative_loci.*_json_files.*.tsv.gz ${local_dir}/hipstr/negative_loci/    
}

    
# Run pipelines on original bam
run_pipelines \
  "gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram" \
  "gs://str-truth-set/hg38/tool_results" \
  "./tool_comparison/results"


# Run pipelines on exome data
run_pipelines \
    "gs://broad-public-datasets/CHM1_CHM13_WES/CHMI_CHMI3_Nex1.cram" \
    "gs://str-truth-set/hg38/tool_results_for_exome" \
    "./tool_comparison/results_for_exome"



# Downsample
for coverage in 30 20 10 5; do 
    python3 ./tool_comparison/hail_batch_pipelines/downsample_bam_pipeline.py --verbose --target-coverage ${coverage} \
	--output-dir gs://bw2-delete-after-5-days/ --input-bam gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram

    # Rerun pipelines on downsampled bam
    run_pipelines \
	"gs://bw2-delete-after-5-days/CHM1_CHM13_WGS2.downsampled_to_${coverage}x.bam" \
	"gs://str-truth-set/hg38/tool_results_for_downsampled_${coverage}x_bam" \
	"./tool_comparison/results_for_downsampled_${coverage}x_bam"
done
