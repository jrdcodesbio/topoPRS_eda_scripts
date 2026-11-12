#!/bin/bash

#SBATCH --job-name=topoEDA_clump
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=06:00:00
#SBATCH --output=/fslustre/labs/ext_rosen_joseph_mayo_edu/topoEDA/slurmLogs/clump_%j.out
#SBATCH --error=/fslustre/labs/ext_rosen_joseph_mayo_edu/topoEDA/slurmLogs/clump_%j.err
#SBATCH --partition=lg-n64-256g

# GWAS QC -> LD reference -> match to reference -> plink2 clump -> locus tables.
#
# Submit from the topoEDA root. slurmLogs/ must exist first, because Slurm opens the log
# files before this script runs:
#   mkdir -p slurmLogs && sbatch scripts/run_clump_pipeline.sbatch
# Override thresholds without editing:
#   sbatch --export=ALL,R2=0.01,KB=1000 scripts/run_clump_pipeline.sbatch
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
REBUILD_LDREF=${REBUILD_LDREF:-0}

module load python/3.10.11

GWAS=$TOPO/sourceData/GWAS/allen3wayGWAS.h.tsv.gz
SCRIPTS=$TOPO/scripts
REF=$TOPO/refData/1000G_hg38_EUR/EUR_hg38
TAG=p1_${P1}_p2_${P2}_r2_${R2}_kb_${KB}
RES=$TOPO/results/clump/$TAG
QC=$TOPO/QC_OUT
mkdir -p "$RES" "$QC/gwas_schema" "$QC/clump/$TAG" "$TOPO/slurmLogs"

step() { echo "[$(date '+%F %T')] === $* ==="; }

step "config"
echo "TOPO=$TOPO PLINK2=$PLINK2 PYTHON=$($PYTHON --version 2>&1) THREADS=$THREADS"
echo "P1=$P1 P2=$P2 R2=$R2 KB=$KB TAG=$TAG"
"$PLINK2" --version
$PYTHON -c "import pandas, numpy; print('pandas', pandas.__version__, 'numpy', numpy.__version__)"
[[ -s $GWAS ]] || { echo "GWAS file not found: $GWAS" >&2; exit 1; }

step "1/5 GWAS schema QC"
$PYTHON "$SCRIPTS/inspect_gwas.py" --gwas "$GWAS" --out-dir "$QC/gwas_schema"

step "2/5 LD reference"
if [[ -s $REF.pgen && $REBUILD_LDREF != 1 ]]; then
    echo "using existing $REF.pgen (set REBUILD_LDREF=1 to rebuild)"
else
    # Probe the host the downloads actually come from. Any HTTP status (even 403/404) means
    # it is reachable; only "000" (no response) means no connectivity.
    code=$(curl -s --max-time 20 -o /dev/null -w '%{http_code}' https://www.dropbox.com/ || true)
    echo "connectivity probe https://www.dropbox.com/ -> HTTP $code"
    if [[ $code == 000 ]]; then
        echo "No internet from $(hostname). Build the reference on a login/transfer node first:" >&2
        echo "  bash $SCRIPTS/build_ldref.sh $TOPO $PLINK2 4" >&2
        exit 1
    fi
    bash "$SCRIPTS/build_ldref.sh" "$TOPO" "$PLINK2" "$THREADS"
fi

step "3/5 match significant GWAS variants to reference"
P_MAX=$($PYTHON -c "print(max(float('$P1'), float('$P2')))")
$PYTHON "$SCRIPTS/prep_clump_input.py" \
    --gwas "$GWAS" --ref-pvar "$REF.pvar" --p-max "$P_MAX" \
    --out-clump "$RES/clump_input.tsv" --qc-dir "$QC/clump/$TAG"

step "4/5 plink2 clump"
bash "$SCRIPTS/run_clump.sh" "$PLINK2" "$REF" "$RES/clump_input.tsv" "$RES/allen_ipf" \
    "$P1" "$P2" "$R2" "$KB" "$THREADS"
cp "$RES/allen_ipf.log" "$QC/clump/$TAG/plink2_clump.log"

step "5/5 locus tables"
$PYTHON "$SCRIPTS/summarize_clumps.py" \
    --clumps "$RES/allen_ipf.clumps" --clump-input "$RES/clump_input.tsv" \
    --unmatched "$QC/clump/$TAG/unmatched_sig_snps.tsv" --kb "$KB" \
    --out-dir "$RES" --qc-dir "$QC/clump/$TAG"

step "done: results in $RES, QC in $QC"
