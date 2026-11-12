#!/bin/bash
# Build the EUR LD reference from the 1000 Genomes hg38 high-coverage panel
# (Byrska-Bishop et al. 2022, 3202 samples), using the plink2-formatted files from
#   https://www.cog-genomics.org/plink/2.0/resources  (release "2022-08-04", option b38b)
#
# Usage: build_ldref.sh <TOPO root> <PLINK2 binary> <threads>
# Output: $TOPO/refData/1000G_hg38_EUR/EUR_hg38.{pgen,pvar,psam}; QC logs in $TOPO/QC_OUT/ldref/
#
# Needs internet access for the downloads (~3.7 GiB). If compute nodes are offline,
# run this once on a login/transfer node: bash scripts/build_ldref.sh <TOPO> <PLINK2> 4
set -euo pipefail

TOPO=${1:?TOPO root required}
PLINK2=${2:?PLINK2 path required}
THREADS=${3:-4}

REF=$TOPO/refData/1000G_hg38_EUR
RAW=$REF/raw
QC=$TOPO/QC_OUT/ldref
mkdir -p "$RAW" "$QC"

# URLs copied from the plink2 resources page (checked 2026-11-12).
# The rsID/no-annotation pvar is used. IDs are reset below anyway; the rsIDs remain in raw/ for reference.
PGEN_URL='https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1'
PVAR_URL='https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1'
PSAM_URL='https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1'
# 1st+2nd-degree relatives for this release (--king-cutoff 0.0884, 629 samples).
DEG2_URL='https://www.dropbox.com/s/4zhmxpk5oclfplp/deg2_hg38.king.cutoff.out.id?dl=1'

fetch() {  # fetch <url> <dest>: skip if already present; download to .part so a failed run leaves no partial file
    local url=$1 dest=$2
    if [[ -s $dest ]]; then echo "[ldref] exists: $dest"; return; fi
    echo "[ldref] downloading $(basename "$dest")"
    curl -L --fail --retry 3 -o "$dest.part" "$url"
    mv "$dest.part" "$dest"
}

fetch "$PGEN_URL" "$RAW/all_hg38.pgen.zst"
fetch "$PVAR_URL" "$RAW/all_hg38.pvar.zst"   # renamed from all_hg38_rs_noannot.pvar.zst, as the page instructs
fetch "$PSAM_URL" "$RAW/all_hg38.psam"       # renamed from hg38_corrected.psam
fetch "$DEG2_URL" "$RAW/deg2_hg38.king.cutoff.out.id"

if [[ ! -s $RAW/all_hg38.pgen ]]; then
    "$PLINK2" --zst-decompress "$RAW/all_hg38.pgen.zst" "$RAW/all_hg38.pgen"
fi

# Select EUR by the SuperPop column; its position is read from the header, not hard-coded.
echo "[ldref] psam header: $(head -1 "$RAW/all_hg38.psam")"
awk -F'\t' '
    NR==1 { for (i=1;i<=NF;i++) { if ($i=="#IID") iid=i; if ($i=="SuperPop") sp=i }
            if (!iid || !sp) { print "psam lacks #IID or SuperPop column" > "/dev/stderr"; exit 1 }
            print "#IID"; next }
    $sp=="EUR" { print $iid }' "$RAW/all_hg38.psam" > "$REF/EUR.keep"
echo "[ldref] EUR samples in psam: $(($(wc -l < "$REF/EUR.keep") - 1))"

# Pass 1: unrelated EUR, autosomal biallelic ACGT SNPs, MAF >= 1% in EUR, positional IDs.
# The MAF filter is applied after sample filtering, so it uses EUR frequencies.
"$PLINK2" --pfile "$RAW/all_hg38" vzs \
    --keep "$REF/EUR.keep" \
    --remove "$RAW/deg2_hg38.king.cutoff.out.id" \
    --chr 1-22 \
    --snps-only just-acgt \
    --max-alleles 2 \
    --maf 0.01 \
    --set-all-var-ids '@:#:$r:$a' \
    --threads "$THREADS" \
    --make-pgen --out "$REF/tmp_EUR_hg38"

# Pass 2: drop every copy of a duplicated chr:pos:ref:alt. This is a separate pass so
# --rm-dup sees the new IDs whatever plink2's internal order of operations is.
"$PLINK2" --pfile "$REF/tmp_EUR_hg38" \
    --rm-dup exclude-all list \
    --threads "$THREADS" \
    --make-pgen --out "$REF/EUR_hg38"

rm -f "$REF"/tmp_EUR_hg38.{pgen,pvar,psam}
cp "$REF"/EUR_hg38.log "$REF"/tmp_EUR_hg38.log "$QC"/ 2>/dev/null || true
[[ -f $REF/EUR_hg38.rmdup.list ]] && cp "$REF/EUR_hg38.rmdup.list" "$QC/"

{
    echo -e "metric\tvalue"
    echo -e "samples\t$(($(wc -l < "$REF/EUR_hg38.psam") - 1))"
    echo -e "variants\t$(grep -vc '^#' "$REF/EUR_hg38.pvar")"
    echo -e "pvar_first_chrom\t$(grep -v '^#' "$REF/EUR_hg38.pvar" | head -1 | cut -f1)"
} | tee "$QC/ldref_counts.tsv"
echo "[ldref] done: $REF/EUR_hg38"
