#!/usr/bin/env python3
"""Annotate every Zuck repeat/spacer and search candidate viral protospacers offline."""
import argparse
from pathlib import Path
import sys
import time
import signal

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from zuck_annotations import annotate


def main():
    def cancel(signum, frame):
        raise KeyboardInterrupt("Annotation cancelled")
    signal.signal(signal.SIGTERM, cancel)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--catalogue", type=Path, help="CRISPRloci snapshot with exact-search.sqlite3")
    source.add_argument("--array-database", type=Path, help="Limited native array TSV projection")
    parser.add_argument("--viral-reference", type=Path)
    parser.add_argument("--viral-sha256")
    parser.add_argument("--crisprspacer", default="crisprspacer")
    parser.add_argument("--blastn", default="blastn")
    parser.add_argument("--mismatches", type=int, choices=(0, 1, 2), default=2)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=300)
    args = vars(parser.parse_args())
    if not 1 <= args["threads"] <= 8 or not 1 <= args["timeout"] <= 86400:
        parser.error("threads must be 1–8 and timeout 1–86400 seconds")
    args["destination"] = args.pop("output")
    deadline = time.monotonic() + args["timeout"]
    def check():
        if time.monotonic() > deadline:
            raise TimeoutError("Annotation deadline exceeded; no completed annotation exists")
    annotate(**args, checkpoint=check)


if __name__ == "__main__":
    main()
