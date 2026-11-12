#!/bin/bash
# LD-clump the prepared GWAS variants against the LD reference with plink2.
# Usage: run_clump.sh <PLINK2> <ref pfile prefix> <clump input> <out prefix> <P1> <P2> <R2> <KB> <threads>
set -euo pipefail

PLINK2=${1:?}; REF=${2:?}; INPUT=${3:?}; OUT=${4:?}
P1=${5:?}; P2=${6:?}; R2=${7:?}; KB=${8:?}; THREADS=${9:?}

mkdir -p "$(dirname "$OUT")"
# Record the flags this plink2 build supports (the HPC build is 2023-11-23).
"$PLINK2" --help clump > "$OUT.plink2_help_clump.txt" 2>&1 || true

# ID and P name the columns written by prep_clump_input.py. The default r² is phased;
# the 1000G hg38 panel is phased.
"$PLINK2" --pfile "$REF" \
    --clump "$INPUT" \
    --clump-id-field ID \
    --clump-p-field P \
    --clump-p1 "$P1" \
    --clump-p2 "$P2" \
    --clump-r2 "$R2" \
    --clump-kb "$KB" \
    --threads "$THREADS" \
    --out "$OUT"

if [[ ! -s $OUT.clumps ]]; then
    echo "[clump] no $OUT.clumps written (no index variants below P1=$P1?); see $OUT.log" >&2
    exit 1
fi
echo "[clump] $(($(wc -l < "$OUT.clumps") - 1)) clumps -> $OUT.clumps"
