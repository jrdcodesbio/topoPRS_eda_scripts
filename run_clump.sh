#!/bin/bash
# LD-clump the prepared GWAS variants against the LD reference with plink2.
# Usage: run_clump.sh <PLINK2> <ref pfile prefix> <clump input> <out prefix> <P1> <P2> <R2> <KB> <threads>
set -euo pipefail

PLINK2=${1:?}; REF=${2:?}; INPUT=${3:?}; OUT=${4:?}
P1=${5:?}; P2=${6:?}; R2=${7:?}; KB=${8:?}; THREADS=${9:?}
UNPHASED=${UNPHASED:-1}   # 1: --clump-unphased (default); 0: plink2's phased r²

mkdir -p "$(dirname "$OUT")"
# Record the flags this plink2 build supports (the HPC build is 2023-11-23).
"$PLINK2" --help clump > "$OUT.plink2_help_clump.txt" 2>&1 || true

# Unphased r² is the default because the phase status of the reference is not established
# (plink2's own default is phased r², which is only meaningful for phased data).
EXTRA=()
if [[ $UNPHASED == 1 ]]; then
    grep -q -- '--clump-unphased' "$OUT.plink2_help_clump.txt" \
        || { echo "[clump] this plink2 build lacks --clump-unphased; see $OUT.plink2_help_clump.txt" >&2; exit 1; }
    EXTRA+=(--clump-unphased)
fi

# ID and P name the columns written by prep_clump_input.py.
"$PLINK2" --pfile "$REF" \
    --clump "$INPUT" \
    --clump-id-field ID \
    --clump-p-field P \
    --clump-p1 "$P1" \
    --clump-p2 "$P2" \
    --clump-r2 "$R2" \
    --clump-kb "$KB" \
    ${EXTRA[@]+"${EXTRA[@]}"} \
    --threads "$THREADS" \
    --out "$OUT"

if [[ ! -s $OUT.clumps ]]; then
    echo "[clump] no $OUT.clumps written (no index variants below P1=$P1?); see $OUT.log" >&2
    exit 1
fi
echo "[clump] $(($(wc -l < "$OUT.clumps") - 1)) clumps -> $OUT.clumps"
