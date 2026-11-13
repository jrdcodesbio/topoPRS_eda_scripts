#!/bin/bash

#SBATCH --job-name=topoEDA_annotate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --output=/fslustre/labs/ext_rosen_joseph_mayo_edu/topoEDA/slurmLogs/annotate_%j.out
#SBATCH --error=/fslustre/labs/ext_rosen_joseph_mayo_edu/topoEDA/slurmLogs/annotate_%j.err
#SBATCH --partition=lg-n64-256g

# Annotate the clumped loci against loops, contact domains and compartment scores (lung + T cell).
# Run after run_clump_pipeline.sbatch. Submit from the topoEDA root:
#   sbatch scripts/run_annotate.sbatch
# Use the loci of another clump parameter set, or change the 'near' window:
#   sbatch --export=ALL,TAG=p1_5e-8_p2_5e-8_r2_0.01_kb_1000,NEAR_KB=25 scripts/run_annotate.sbatch
# Needs pyBigWig in the python module env: pip install --user pyBigWig
set -euo pipefail

TOPO=${TOPO:-/fslustre/labs/ext_rosen_joseph_mayo_edu/topoEDA}
PYTHON=${PYTHON:-python3}
TAG=${TAG:-p1_5e-8_p2_5e-8_r2_0.1_kb_250}   # the run_clump_pipeline.sbatch default tag
NEAR_KB=${NEAR_KB:-10}

module load python/3.10.11

SCRIPTS=$TOPO/scripts
CLUMP=$TOPO/results/clump/$TAG
GEN=$TOPO/sourceData/genFiles
OUT=$TOPO/results/annotate/$TAG
QC=$TOPO/QC_OUT/annotate/$TAG
mkdir -p "$OUT" "$QC"

echo "[$(date '+%F %T')] TAG=$TAG NEAR_KB=$NEAR_KB PYTHON=$($PYTHON --version 2>&1)"

if ! $PYTHON -c "import pyBigWig" >/dev/null 2>&1; then
    echo "[$(date '+%F %T')] pyBigWig not found; installing..."
    $PYTHON -m pip install --user pyBigWig || {
        echo "Failed to install pyBigWig" >&2
        exit 1
    }
fi

$PYTHON -c "import pandas, numpy, pyBigWig; print('pandas', pandas.__version__, 'numpy', numpy.__version__, 'pyBigWig', getattr(pyBigWig, '__version__', 'unknown'))"
for f in "$CLUMP/loci.tsv" "$CLUMP/locus_members.tsv"; do
    [[ -s $f ]] || { echo "missing $f; run run_clump_pipeline.sbatch first (or set TAG)" >&2; exit 1; }
done

$PYTHON "$SCRIPTS/annotate_loci.py" \
    --loci "$CLUMP/loci.tsv" --members "$CLUMP/locus_members.tsv" \
    --reference "lung=$GEN/lung" --reference "Tcell=$GEN/Tcell" \
    --near-kb "$NEAR_KB" --out-dir "$OUT" --qc-dir "$QC"

echo "[$(date '+%F %T')] done: results in $OUT, QC in $QC"
