#!/usr/bin/env bash
# Rebuild both papers and run the interval and algebraic certificates in a temporary directory.
set -euo pipefail
cd -- "$(dirname -- "$0")"
export LC_ALL=C TZ=UTC PYTHONDONTWRITEBYTECODE=1
unset PYTHONOPTIMIZE
sha256sum --quiet -c SHA256SUMS
echo 'Scientific file checksums verified.'
scratch_dir=$(mktemp -d)
trap 'rm -rf -- "$scratch_dir"' EXIT
cp -R -- source "$scratch_dir/source"
(cd "$scratch_dir/source" && ./build.sh pdfs && python3 -B make-replay-archive.py)
cmp main.pdf "$scratch_dir/source/main.pdf"
cmp low-degree-supplement.pdf "$scratch_dir/source/low-degree-supplement.pdf"
cmp source/certificates.tar.gz "$scratch_dir/source/certificates.tar.gz"
echo 'Both PDFs and the regenerated certificate archive match byte for byte.'
tar -xzf "$scratch_dir/source/certificates.tar.gz" -C "$scratch_dir"
(
    cd "$scratch_dir/certificates"
    sha256sum --quiet -c CERTIFICATE-SHA256SUMS
    ./replay-all.sh
    python3 -B verification/check_action_coordinates.py > "$scratch_dir/coordinates.txt"
    cmp verification/action-coordinates.txt "$scratch_dir/coordinates.txt"
    cat "$scratch_dir/coordinates.txt"
    sha256sum --quiet -c CERTIFICATE-SHA256SUMS
)
echo 'The coordinate test reproduced its fourth transcript.'
sha256sum source/REPLAY.txt source/verification/fixed-degree12-replay.txt \
    source/verification/fixed-inputs.txt source/verification/action-coordinates.txt source/certificates.tar.gz
sha256sum --quiet -c SHA256SUMS
echo 'Scientific file checksums verified after replay.'
