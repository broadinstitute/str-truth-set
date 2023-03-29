"""This script checks that all file paths references in tgg-viewer-config.json are valid."""

import argparse
import hail as hl
import json
import os
import sys


def is_url(path):
    return path.startswith("gs://") or path.startswith("http://") or path.startswith("https://")


def get_all_urls(json_data):
    """Recursively traverse the given JSON data structure and yield all URLs found at the leaves."""
    if isinstance(json_data, str) and is_url(json_data):
        yield json_data
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
    for url in urls:
        if not hl.hadoop_is_file(url):
            print(f"File not found: {url}")
            missing_file_counter += 1

    if missing_file_counter:
        print(f"ERROR: {missing_file_counter} out of {len(url)} files are missing.")
        sys.exit(1)

    print(f"Done checking urls. All {len(url)} file paths are valid.")


if __name__ == "__main__":
    main()