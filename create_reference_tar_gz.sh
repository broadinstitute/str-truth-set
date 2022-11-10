set -x
output_path=STR_truth_set_reference_bundle.v1.tar.gz
tar czf ${output_path} ref/*.chain ref/other/*.*
echo Generated $(pwd)/${output_path}
ls -lh ${output_path}
