#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

printf '%s  %s\n' \
  '124145efc96913b16fbb4360a767d850909dcfd36610257b8301e1f05c1fd5cd' \
  'reference/checks/verify_seed.py' | sha256sum --check

python3 -B -u outer.py
python3 -B -u compact.py
python3 -B -u shooting.py
python3 -B -u fixed_outer.py
