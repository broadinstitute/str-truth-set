"""This script checks that all file paths references in tgg-viewer-config.json are valid."""

import argparse
import hail as hl
import json
import os
import sys


def is_url(path):
    """Takes a string and returns True if it starts with gs://, http://, or https://."""

    return path.startswith("gs://") or path.startswith("http://") or path.startswith("https://")


def get_all_urls(json_data):
    """Recursively traverse the given JSON data structure and yield all gs:// and http paths found at the leaves."""
    if isinstance(json_data, str) and is_url(json_data):
        path = json_data
        yield path

    elif isinstance(json_data, dict):
        for value in json_data.values():
            for path in get_all_urls(value):
                yield path

    elif isinstance(json_data, (list, tuple)):
        for value in json_data:
            for path in get_all_urls(value):
                yield path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("json_config_path", nargs="?", default="tgg-viewer-config.json",
                        help="Path to tgg-viewer-config.json")
    args = parser.parse_args()

    if not os.path.exists(args.json_config_path):
        parser.error(f"File not found: {args.json_config_path}")

    with open(args.json_config_path, "rt") as f:
        config = json.load(f)

    urls = list(get_all_urls(config))
    if len(urls) == 0:
        parser.error(f"No URLs found in {args.json_config_path}")

    hl.init(log="/dev/null")

    missing_file_counter = 0
    total_file_counter = 0
    for url in urls:
        total_file_counter += 1
        if not hl.hadoop_is_file(url):
            print(f"File not found: {url}")
            missing_file_counter += 1

        # check that index file also exists
        for suffix, idx_suffix in (
            (".bam", ".bai"),
            (".cram", ".crai"),
            (".bed.gz", ".bed.gz.tbi"), (".bed.bgz", ".bed.gz.tbi"),
            (".vcf.gz", ".vcf.gz.tbi"), (".vcf.bgz", ".vcf.bgz.tbi"),
            (".gtf.gz", ".gtf.gz.tbi"), (".gtf.bgz", ".gtf.bgz.tbi"),
            (".json", None),
            (".bw", None),  # BigWig files don't have an index file
        ):
            if url.endswith(suffix):
                if idx_suffix is not None:
                    total_file_counter += 1
                    if not hl.hadoop_is_file(f"{url}{idx_suffix}") and \
                       not hl.hadoop_is_file(f"{url.replace(suffix, idx_suffix)}"):
                        print(f"File not found: {url}{idx_suffix}")
                        missing_file_counter += 1
                break
        else:
            print("WARNING: unexpected file extension:", url)

    if missing_file_counter:
        print(f"ERROR: {missing_file_counter:,d} out of {total_file_counter:,d} files are missing.")
        sys.exit(1)

    print(f"Done checking urls. All {total_file_counter:,d} file paths are valid.")


if __name__ == "__main__":
    main()