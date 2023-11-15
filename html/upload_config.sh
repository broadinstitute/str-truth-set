set -ex

python3 check_tgg_viewer_config_paths.py tgg-viewer-config.json
gsutil -m cp tgg-viewer-config.json gs://str-truth-set/hg38/
echo Uploaded gs://str-truth-set/hg38/tgg-viewer-config.json

python3 check_tgg_viewer_config_paths.py tgg-viewer-config-catalogs-only.json
gsutil -m cp tgg-viewer-config-catalogs-only.json gs://str-truth-set/hg38/
echo Uploaded gs://str-truth-set/hg38/tgg-viewer-config-catalogs-only.json 
