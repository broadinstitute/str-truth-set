set -ex

python3 check_tgg_viewer_config_paths.py tgg-viewer-config.json

gsutil -m cp tgg-viewer-config.json gs://str-truth-set/hg38/
