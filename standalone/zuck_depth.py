#!/usr/bin/env python3
"""Estimate metagenomic depth from FASTQ and optionally route a Zuck 6.1 run."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
ROUTED_FLAGS = ("read-first-shallow-rescue", "read-first-family-admission")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--input-files", type=Path, nargs="+", required=True,
                        metavar="FASTQ", help="One FASTQ or two mate files, plain or gzip.")
    common.add_argument("--scanner", type=Path, default=HERE / "zuck-kmer-spectrum",
                        help="Path to the k-mer spectrum executable.")
    commands.add_parser("estimate", parents=[common], help="Write the coverage estimate as JSON.")
    run = argparse.ArgumentParser(add_help=False)
    run.add_argument("--zuck", type=Path, default=HERE / "Zuck", help="Zuck 6.1 executable.")
    run.add_argument("--settings", type=Path, default=HERE / "v6.1.conf",
                     help="Native preset. The route controls the two read-first switches below.")
    run.add_argument("--output-folder", type=Path, required=True,
                     help="A new or empty result directory.")
    for name in ROUTED_FLAGS:
        run.add_argument("--" + name, choices=("true", "false"),
                         help="Override this switch after coverage routing.")
    commands.add_parser("run", parents=[common, run], help="Use the original 6.1 shallow-only routing policy.")
    commands.add_parser("adaptive", parents=[common, run], help="Run native sparse-junction recovery at low predicted depth.")
    return parser.parse_args()


def executable(path):
    path = path.resolve(strict=True)
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(f"Cannot execute {path}")
    return str(path)


def estimate(scanner, reads):
    from zuck_depth_model import predict

    result = subprocess.run([scanner, *map(str, reads)], stdout=subprocess.PIPE, check=True)
    spectrum = json.loads(result.stdout)
    prediction = predict(spectrum)
    if prediction["regime"] not in {"low", "intermediate", "high", "mixed", "uncertain"}:
        raise ValueError(f"Unrecognized coverage regime: {prediction['regime']}")
    return {**prediction, "spectrum": spectrum}


def json_text(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def prepare_output(folder):
    folder.mkdir(parents=True, exist_ok=True)
    if any(folder.iterdir()):
        raise ValueError(f"Choose a new or empty result directory: {folder}")
    return folder / "coverage_routing.json"


def exit_status(returncode):
    return 128 - returncode if returncode < 0 else returncode


def run(args, reads, prediction):
    binary = executable(args.zuck)
    settings = args.settings.resolve(strict=True)
    settings_bytes = settings.read_bytes()
    # Native presets can redirect result files outside --output-folder.
    # A routed run places all new outputs inside the requested directory.
    for line in settings_bytes.decode().splitlines():
        key = line.partition("=")[0].strip().replace("_", "-")
        if key in {"output-file", "graph-folder", "cycles-folder"}:
            raise ValueError(f"Remove {key} from the preset to keep results inside --output-folder.")
    folder = args.output_folder.resolve()
    low = prediction["regime"] == "low"
    switches = dict(zip(ROUTED_FLAGS, (True, not low)))
    mode = "shallow_only" if low else "full_6.1"
    policy = "shallow_only_if_low_v1"
    if args.command == "adaptive":
        switches["read-first-family-admission"] = True
        switches["read-first-sparse-junctions"] = low
        mode = "sparse_junctions" if low else "full_6.1"
        policy = "native_sparse_junctions_if_low_v1"
    overrides = {}
    for flag in ROUTED_FLAGS:
        value = getattr(args, flag.replace("-", "_"))
        if value is not None:
            switches[flag] = overrides[flag] = value == "true"
    command = [binary, "--settings", str(settings), "--input-files", *map(str, reads),
               "--output-folder", str(folder)]
    for flag, value in switches.items():
        command.extend(("--" + flag, str(value).lower()))

    routing = {"policy": policy, "selected_mode": mode,
               "effective_switches": switches, "explicit_overrides": overrides, "argv": command,
               "settings_sha256": hashlib.sha256(settings_bytes).hexdigest(), "native_exit_status": None}
    report = {**prediction, "routing": routing}
    receipt = prepare_output(folder)
    # Exclusive creation also reserves an initially empty directory between launchers.
    with receipt.open("x", encoding="utf-8") as stream:
        stream.write(json_text(report))
    try:
        result = subprocess.run(command, stdout=sys.stderr)
        status = exit_status(result.returncode)
    except OSError as error:
        routing["launch_error"] = str(error)
        status = 126
    routing["native_exit_status"] = status
    rendered = json_text(report)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=folder,
                                     prefix=".coverage_routing.", delete=False) as stream:
        stream.write(rendered)
        staging = stream.name
    os.replace(staging, receipt)
    sys.stdout.write(rendered)
    return status


def main():
    args = arguments()
    if len(args.input_files) not in (1, 2):
        raise ValueError("--input-files accepts one FASTQ or a pair of mates")
    # Preserve user-supplied aliases so a blinded scanner sees opaque filenames.
    reads = [path.absolute() for path in args.input_files]
    if not all(path.is_file() for path in reads):
        raise ValueError("Each FASTQ must be a regular file")
    if len(reads) == 2 and reads[0].samefile(reads[1]):
        raise ValueError("Two mate inputs must identify different files")
    prediction = estimate(executable(args.scanner), reads)
    if args.command in {"run", "adaptive"}:
        return run(args, reads, prediction)
    sys.stdout.write(json_text(prediction))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except subprocess.CalledProcessError as error:
        print(f"Coverage scanner failed with status {error.returncode}", file=sys.stderr)
        sys.exit(exit_status(error.returncode))
    except (OSError, ValueError, KeyError) as error:
        print(f"Coverage routing failed: {error}", file=sys.stderr)
        sys.exit(2)
