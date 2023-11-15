"""This script checks that all file paths references in tgg-viewer-config.json are valid."""

import argparse
import hailtop.fs as hfs
import json
import os
import requests
import sys
import tqdm

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

def check_url(url):
    if url.startswith("http"):
        r = requests.get(url, stream=True, verify=False)
        
        if r.status_code == 200:
            content = next(r.iter_content(10))
            return True
        
        return False
    elif url.startswith("gs://"):
        return hfs.is_file(url)
    else:
        print(f"ERROR: Unexpected url prefix: {url}")
        return False
        
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

    missing_file_counter = 0
    total_file_counter = 0
    for url in tqdm.tqdm(urls, unit=" urls"):
        total_file_counter += 1
        if not check_url(url):
            print(f"File not found: {url}")
            missing_file_counter += 1
            continue
        
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
                    if not check_url(f"{url}{idx_suffix}") and not check_url(f"{url.replace(suffix, idx_suffix)}"):
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
