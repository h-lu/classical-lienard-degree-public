#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
unset PYTHONOPTIMIZE
export LC_ALL=C PYTHONDONTWRITEBYTECODE=1
replay_dir=$(mktemp -d)
trap 'rm -rf -- "$replay_dir"' EXIT
./build.sh replay > "$replay_dir/manuscript.txt"
cmp REPLAY.txt "$replay_dir/manuscript.txt"
bash verification/fixed-degree12/replay.sh > "$replay_dir/fixed-six.txt"
cmp verification/fixed-degree12-replay.txt "$replay_dir/fixed-six.txt"
python3 -B verification/check_fixed_inputs.py > "$replay_dir/fixed-inputs.txt"
cmp verification/fixed-inputs.txt "$replay_dir/fixed-inputs.txt"
echo 'All 13 manuscript programs, 4 fixed-field programs and the coefficient/escape check reproduced their three transcripts byte for byte.'
