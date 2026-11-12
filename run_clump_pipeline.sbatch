#!/bin/bash

#SBATCH --job-name=topoEDA_clump
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=06:00:00
#SBATCH --output=/fslustre/labs/ext_rosen_joseph_mayo_edu/topoEDA/slurmLogs/clump_%j.out
#SBATCH --error=/fslustre/labs/ext_rosen_joseph_mayo_edu/topoEDA/slurmLogs/clump_%j.err
#SBATCH --partition=lg-n64-256g

# GWAS QC -> match to the LD reference -> plink2 clump -> locus tables.
# LD reference: MBB Extracted_PRS_MBB_gwide (hg38, European, relatives already removed per the user).
#
# Submit from the topoEDA root. slurmLogs/ must exist first, because Slurm opens the log
# files before this script runs:
#   mkdir -p slurmLogs && sbatch scripts/run_clump_pipeline.sbatch
# Override thresholds without editing:
#   sbatch --export=ALL,R2=0.01,KB=1000 scripts/run_clump_pipeline.sbatch
# Other overrides: REF=<pfile prefix>, UNPHASED=0 (phased r² instead of the default unphased).
# Each threshold set writes to its own results/clump/<tag>/ directory, so runs don't overwrite each other.
set -euo pipefail

TOPO=${TOPO:-/fslustre/labs/ext_rosen_joseph_mayo_edu/topoEDA}
PLINK2=${PLINK2:-/usr/local/biotools/plink2/v2.00a6LM.2023.1123/plink2}
PYTHON=${PYTHON:-python3}          # needs pandas + numpy (see requirements.txt)
THREADS=${SLURM_NTASKS:-16}

P1=${P1:-5e-8}                     # index variant P threshold
P2=${P2:-5e-8}                     # clump member P threshold
R2=${R2:-0.1}
KB=${KB:-250}
UNPHASED=${UNPHASED:-1}            # 1: unphased r² (phase status of the reference is unverified)
export UNPHASED

module load python/3.10.11

GWAS=$TOPO/sourceData/GWAS/allen3wayGWAS.h.tsv.gz
SCRIPTS=$TOPO/scripts
REF=${REF:-/fslustre/labs/ext_rosen_joseph_mayo_edu/MBB_CTD_PRS_GWAS/Extracted_PRS_MBB_gwide}
TAG=p1_${P1}_p2_${P2}_r2_${R2}_kb_${KB}
RES=$TOPO/results/clump/$TAG
QC=$TOPO/QC_OUT
mkdir -p "$RES" "$QC/gwas_schema" "$QC/clump/$TAG" "$TOPO/slurmLogs"

step() { echo "[$(date '+%F %T')] === $* ==="; }

step "config"
echo "TOPO=$TOPO PLINK2=$PLINK2 PYTHON=$($PYTHON --version 2>&1) THREADS=$THREADS"
echo "P1=$P1 P2=$P2 R2=$R2 KB=$KB UNPHASED=$UNPHASED TAG=$TAG"
echo "REF=$REF"
"$PLINK2" --version
$PYTHON -c "import pandas, numpy; print('pandas', pandas.__version__, 'numpy', numpy.__version__)"
[[ -s $GWAS ]] || { echo "GWAS file not found: $GWAS" >&2; exit 1; }
for ext in pgen pvar psam; do
    [[ -s $REF.$ext ]] || { echo "LD reference file not found: $REF.$ext" >&2; exit 1; }
done

# Record exactly which reference this run used.
PROV=$QC/clump/$TAG/ref_provenance
mkdir -p "$PROV"
{
    echo -e "ref	$REF"
    echo -e "samples	$(($(wc -l < "$REF.psam") - 1))"
    echo -e "variants	$(grep -vc '^#' "$REF.pvar")"
} | tee "$PROV/ref_counts.tsv"
head -5 "$REF.psam" > "$PROV/psam_head.txt"
grep '^##' "$REF.pvar" | head -50 > "$PROV/pvar_meta_head.txt" || true
# head closes the pipe early; under pipefail that SIGPIPE must not abort the job.
(grep -v '^#' "$REF.pvar" | head -3 | tee "$PROV/pvar_first_records.txt") || true
[[ -f $REF.log ]] && cp "$REF.log" "$PROV/"

step "1/4 GWAS schema QC"
$PYTHON "$SCRIPTS/inspect_gwas.py" --gwas "$GWAS" --out-dir "$QC/gwas_schema"

step "2/4 match significant GWAS variants to reference"
P_MAX=$($PYTHON -c "print(max(float('$P1'), float('$P2')))")
$PYTHON "$SCRIPTS/prep_clump_input.py" \
    --gwas "$GWAS" --ref-pvar "$REF.pvar" --p-max "$P_MAX" \
    --out-clump "$RES/clump_input.tsv" --qc-dir "$QC/clump/$TAG"

step "3/4 plink2 clump"
bash "$SCRIPTS/run_clump.sh" "$PLINK2" "$REF" "$RES/clump_input.tsv" "$RES/allen_ipf" \
    "$P1" "$P2" "$R2" "$KB" "$THREADS"
cp "$RES/allen_ipf.log" "$QC/clump/$TAG/plink2_clump.log"

step "4/4 locus tables"
$PYTHON "$SCRIPTS/summarize_clumps.py" \
    --clumps "$RES/allen_ipf.clumps" --clump-input "$RES/clump_input.tsv" \
    --unmatched "$QC/clump/$TAG/unmatched_sig_snps.tsv" --kb "$KB" \
    --out-dir "$RES" --qc-dir "$QC/clump/$TAG"

step "done: results in $RES, QC in $QC"
