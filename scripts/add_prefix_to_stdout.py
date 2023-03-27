"""This script reads lines from stdin and adds a prefix to each line.
Useful for adding a prefix to the stdout of a command.
"""

import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("prefix")
    args = p.parse_args()

    for line in sys.stdin:
        print(args.prefix + line, end="")


if __name__ == "__main__":
    main()
