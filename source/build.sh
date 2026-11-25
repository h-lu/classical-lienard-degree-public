#!/usr/bin/env bash
set -euo pipefail

cd -- "$(dirname -- "$0")"
export LC_ALL=C TZ=UTC FORCE_SOURCE_DATE=1 SOURCE_DATE_EPOCH=1795392000
export PYTHONDONTWRITEBYTECODE=1
unset PYTHONOPTIMIZE

mode=${1:-all}
case "$mode" in
    main|supplement|pdfs|replay|all) ;;
    *) echo "Usage: ./build.sh [main|supplement|pdfs|replay|all]" >&2; exit 2 ;;
esac

stage=$(mktemp -d "./.build.XXXXXXXX")
trap 'rm -rf -- "$stage"' EXIT

build_pdf() {
    local stem=$1
    local pass
    for pass in 1 2 3 4 5 6; do
        if ! pdflatex -no-shell-escape -halt-on-error -file-line-error \
            -interaction=nonstopmode -output-directory="$stage" "$stem.tex" \
            > "$stage/$stem.stdout"; then
            tail -n 40 "$stage/$stem.stdout" >&2
            return 1
        fi
        if [[ $pass == 1 ]]; then
            if ! (cd "$stage" && BIBINPUTS=..: bibtex "$stem") > "$stage/$stem.bibstdout"; then
                cat "$stage/$stem.bibstdout" >&2
                return 1
            fi
            if grep -n 'Warning--' "$stage/$stem.blg"; then
                echo "Unresolved bibliography diagnostic in $stem" >&2
                return 1
            fi
        fi
        if [[ $pass -ge 3 ]] && ! grep -qE 'Rerun to get|Rerun to get outlines|Label\(s\) may have changed|undefined references' "$stage/$stem.log"; then
            break
        fi
    done
    if grep -nE 'Warning|Overfull|Underfull|Undefined|multiply defined|Missing character:' "$stage/$stem.log"; then
        echo "Unresolved typesetting diagnostic in $stem" >&2
        return 1
    fi
    mv -f -- "$stage/$stem.pdf" "$stem.pdf"
    sha256sum "$stem.pdf"
}

if [[ $mode == main || $mode == pdfs || $mode == all ]]; then
    build_pdf main
fi
if [[ $mode == supplement || $mode == pdfs || $mode == all ]]; then
    build_pdf low-degree-supplement
fi
if [[ $mode == replay || $mode == all ]]; then
    for program in checks/verify_*.py; do
        echo "$program"
        python3 -B "$program"
    done
fi
