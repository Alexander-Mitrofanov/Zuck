#!/usr/bin/env bash
set -euo pipefail

# Replay a retained graph with different workers and fixed logical partitions.
usage() {
  cat <<EOF
Usage: $0 --binary PATH --input GRAPH_DIR --threads LIST [--runs N] -o OUTCSV

  --binary PATH     Zuck executable
  --input PATH      Retained graph directory (the --graph input to Zuck)
  --threads LIST    Comma-separated positive worker counts (e.g. 1,4,8)
  --runs N          Positive runs per worker count (default: 3)
  -o OUTCSV         CSV output (overwritten)
EOF
}
fail() { echo "Error: $*" >&2; exit 1; }

binary="" graph="" thread_list="" runs=3 output=""
while (($#)); do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --binary|--input|--threads|--runs|-o)
      (($# >= 2)) || fail "$1 requires a value"
      case "$1" in
        --binary) binary=$2 ;;
        --input) graph=$2 ;;
        --threads) thread_list=$2 ;;
        --runs) runs=$2 ;;
        -o) output=$2 ;;
      esac
      shift 2 ;;
    *) fail "unknown argument: $1" ;;
  esac
done
[[ -n "$binary" && -n "$graph" && -n "$thread_list" && -n "$output" ]] || { usage >&2; exit 1; }
[[ -x "$binary" ]] || fail "binary is not executable: $binary"
[[ -d "$graph" ]] || fail "graph directory does not exist: $graph"
[[ "$runs" =~ ^[1-9][0-9]*$ ]] || fail "--runs must be a positive integer"
[[ ${#runs} -le 9 ]] || fail "--runs is too large (maximum 999999999)"
[[ "$thread_list" =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]] || fail "--threads must contain positive integers"
[[ -x /usr/bin/time ]] || fail "/usr/bin/time is required"
command -v python3 >/dev/null || fail "python3 is required to read run manifests"

# Each invocation owns its scratch directory, including after failed runs.
scratch=$(mktemp -d "${TMPDIR:-/tmp}/zuck-benchmark.XXXXXX")
trap 'rm -rf -- "$scratch"' EXIT
mkdir -p -- "$(dirname -- "$output")"
printf 'threads,run,elapsed_s,max_rss_kb,cycles\n' > "$output"
IFS=',' read -r -a workers <<< "$thread_list"
for threads in "${workers[@]}"; do
  for ((run=1; run<=runs; run++)); do
    result="$scratch/result-$threads-$run"
    echo "Running threads=$threads run=$run"
    /usr/bin/time -f '%e %M' -o "$scratch/time.txt" \
      "$binary" --graph "$graph" --threads "$threads" --output-folder "$result" \
      > "$scratch/stdout.txt" 2> "$scratch/stderr.txt" || {
        cat "$scratch/stderr.txt" >&2
        fail "Zuck failed for threads=$threads run=$run"
      }
    read -r elapsed rss_kb < "$scratch/time.txt"
    cycles=$(python3 - "$result/run_manifest.json" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as manifest:
    print(json.load(manifest)["counts"]["loaded_cycle_records"])
PY
)
    printf '%s,%s,%s,%s,%s\n' "$threads" "$run" "$elapsed" "$rss_kb" "$cycles" >> "$output"
  done
done
echo "Done. Results at: $output"
