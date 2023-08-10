# Run tool pipelines using Hail Batch


force=""
#force="--force"

debug=""
#debug="echo"
if [ -z "$debug" ]; then
  set -ex
fi

input_bam="gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram"
input_bai="gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai"

for max_dist in 50 24 10; do
    output_dir="gs://str-truth-set/hg38/tool_results/expansion_hunter_with_adjacent_loci__max_dist_${max_dist}bp"
    variant_catalog_paths="gs://str-truth-set/hg38/variant_catalogs/expansion_hunter_with_adjacent_loci__max_dist_${max_dist}bp/positive*.with_adjacent_loci.json"
    run_reviewer=""
    if [ "${max_dist}" == "50" ]; then
	echo Running reviewer for max_dist = ${max_dist}
	run_reviewer="--run-reviewer"
    fi
    
    $debug python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py \
       --input-bam ${input_bam} \
       --input-bai ${input_bai} \
       --variant-catalog-paths ${variant_catalog_paths} \
       --output-dir ${output_dir} \
       ${run_reviewer} \
       --verbose \
       ${min_locus_coverage_arg} \
       ${force} 

    local_dir=./tool_comparison/results_with_adjacent_loci__max_dist_${max_dist}bp
    
    $debug mkdir -p ${local_dir}
    $debug gsutil -m cp ${output_dir}/positive_loci/combined.positive_loci.*_json_files.*.tsv.gz  ${local_dir}/
done

# Run REViewer
python3 ./tool_comparison/hail_batch_pipelines/expansion_hunter_pipeline.py  --positive-loci --run-reviewer \
	--input-bam gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram --input-bai gs://broad-public-datasets/CHM1_CHM13_WGS2/CHM1_CHM13_WGS2.cram.crai \
        --output-dir gs://str-truth-set/hg38/tool_results/expansion_hunter &

echo Done
