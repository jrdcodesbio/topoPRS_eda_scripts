# topoEDA: IPF GWAS x chromatin architecture

## Purpose and current scope

Exploratory integration of the **Allen 3-way IPF GWAS (hg38/GRCh38)** with two ENCODE Hi-C references:

- **Lung:** `ENCSR305GRJ`.
- **Primary CD8 T cells:** `ENCSR321BHC`.\

Explore where associated variants lie, which regions/promoters they contact, and how their structural context differs between references. This is **EDA and technique development**, not TopoPRS training, a new GWAS, or patient-level prediction. These reference contexts are not IPF cases versus controls. Large differences between maps are expected; document them without blocking exploration or repeatedly focusing on loop counts.

## Execution environment and paths

**Write code for deployment and execution on the HPC system where the data reside.**  The dev environment cannot access these paths. I will pusn to a git repo and pull the modified code onto the HPC from there. 

Project root:

```text
/fslustre/labs/ext_rosen_joseph_mayo_edu/topoEDA/
├── scripts/
│   └── qc_parser.py
└── sourceData/
    ├── GWAS/
    │   └── allen3wayGWAS.h.tsv.gz
    ├── genFiles/
    │   ├── lung/
    │   │   ├── compartments.bigWig
    │   │   ├── contactDom.bedpe.gz
    │   │   ├── contactMat.hic
    │   │   └── loops.bedpe.gz
    │   └── Tcell/
    │       ├── compartments.bigWig
    │       ├── contactDom.bedpe.gz
    │       ├── ENCFF009ONH.hic
    │       ├── ENCFF009ONH.hic.1
    │       └── loops.bedpe.gz
    └── phenoFiles/                 # Contents/use not established; not needed for current EDA.
```

Created by the clumping pipeline (not present before the first run):

```text
├── results/clump/<tag>/            # clump input, plink2 output, loci tables
├── QC_OUT/                         # gwas_schema/, clump/<tag>/ (incl. ref_provenance/)
└── slurmLogs/                      # must exist before sbatch (Slurm opens logs first)
```

Paths are case-sensitive.

Keep scripts under `scripts/`; create `results/` for new outputs. Keep `sourceData/` unchanged. Existing QC outputs were written outside the project to `/fslustre/labs/ext_rosen_joseph_mayo_edu/hic_qc/`, that needs ot be changed. from now on all scripts should publish qc to a specific dir "QC_OUT" within topoEDA.

## Inputs and parsing

| Input | Meaning / handling |
|---|---|
| GWAS TSV.gz | Summary association statistics, not individual genotypes. hg38 per user; inspect actual columns, position convention, alleles, missing values, and P-value representation before coding joins. |
| `loops.bedpe.gz` | Called focal contacts: two anchor intervals plus contact/background/significance statistics. Use for variant-to-anchor and partner-region mapping. |
| `contactDom.bedpe.gz` | Arrowhead contact domains represented as diagonal interval pairs; annotate domain membership, not additional loop edges. |
| `compartments.bigWig` | Numerical compartment track. Establish score orientation before assigning A/B labels. |
| `.hic` | Multiresolution contact matrices underlying the structural calls. Query selected regions/resolutions; never load the entire genome as a dense matrix. |

BEDPE headers start with `#chr1`; subsequent `# juicer_tools ...` lines are comments. Read the header explicitly, preserve metadata, and make duplicate names unique (domain files have two `score` columns). Do not discard the header and infer names from the first data row.

First six columns: `chr1 x1 x2 chr2 y1 y2`. BEDPE coordinates are zero-based, half-open. Verify the GWAS coordinate convention before conversion and standardize chromosome names. Preserve `highRes_*`/`local*` fields without silently substituting their coordinates. Confirm Hi-C assembly from metadata, not the filenames.

## Completed work

| Reference | Loops | Domains |
|---|---:|---:|
| Lung | 729 | 945 |
| T cell | 74,287 | 1,562 |

Loop tables have 33 columns; domain tables have 16. No exact duplicate pairs; all listed contacts are cis, and all domains are diagonal. No additional significance filter was applied. Neither matrices nor bigWigs have been inspected. Zero trans calls in BEDPE does not establish absence of trans contacts in `.hic`.

## GWAS clumping pipeline (implemented; synthetic tests only)

LD-clumps genome-wide-significant variants into loci for the Hi-C annotation step. The 6 unit tests in `tests/test_clump_pipeline.py` pass on synthetic data (`python -m unittest tests/test_clump_pipeline.py`). **HPC status:** step 1 (`inspect_gwas.py`) ran on the real GWAS (10,767,487 rows; 3,553 with p<5e-8, 3,978 with p<1e-6, 4,819 with p<1e-5; no p=0). The clump steps have not yet been run with the MBB reference.

| Step | Script | What it does |
|---|---|---|
| 1 | `inspect_gwas.py` | Chunked schema/missing-value/P-value QC of the GWAS; changes nothing. |
| 2 | `prep_clump_input.py` | Matches significant GWAS variants to the reference `.pvar` on chromosome, position and unordered allele pair; writes the clump input using reference IDs. Unmatched and incomplete rows are listed in QC, not dropped silently. Stops if a matched reference ID is `.`/empty or occurs more than once in the pvar (`duplicate_ref_ids.tsv`), and warns if under 50% of significant variants match (check genome build/coverage). |
| 3 | `run_clump.sh` | plink2 `--clump` (reference IDs, `P` column). Uses **unphased** r² (`--clump-unphased`) by default because the phase status of the reference is not established; `UNPHASED=0` restores plink2's phased default. |
| 4 | `summarize_clumps.py` | Writes `loci.tsv` (one row per clump) and `locus_members.tsv` (one row per locus/variant, original GWAS columns joined); flags unmatched significant SNPs lying within `KB` of a locus. |

**LD reference:** `MBB_CTD_PRS_GWAS/Extracted_PRS_MBB_gwide.{pgen,pvar,psam}` under `/fslustre/labs/ext_rosen_joseph_mayo_edu/` (per the user: European, hg38, relatives already removed with a plink KING matrix; not re-verified here). Nothing is built or downloaded; override the location with `REF=<pfile prefix>`. Each run copies the reference's sample/variant counts, `.psam` head, pvar header/first records and `.log` to `QC_OUT/clump/<tag>/ref_provenance/`.

`gwas_common.py` holds shared readers (all fields read as strings; explicit missing tokens). Dependencies: `requirements.txt` (pandas, numpy) plus plink2.

Run (from the topoEDA root):

```bash
mkdir -p slurmLogs && sbatch scripts/run_clump_pipeline.sbatch
sbatch --export=ALL,R2=0.01,KB=1000 scripts/run_clump_pipeline.sbatch   # override thresholds
sbatch --export=ALL,UNPHASED=0 scripts/run_clump_pipeline.sbatch        # phased r2 instead
```

Defaults: `P1=P2=5e-8`, `R2=0.1`, `KB=250`, `UNPHASED=1`. Each parameter set writes to `results/clump/p1_<P1>_p2_<P2>_r2_<R2>_kb_<KB>/`, with QC in `QC_OUT/`.

Caveats: loci are LD-based on the MBB samples, so they depend on how well that cohort's LD matches the GWAS population (both European per the user; not verified here). The reference may be a variant subset ('Extracted'); significant variants absent from it cannot be clumped, so check `unmatched_sig_snps.tsv`, `unmatched_near_loci.tsv` and `frac_sig_matched` in `prep_summary.json`. Whether `--clump-unphased` is accepted by the HPC plink2 build (2023-11-23) is checked at run time against its saved `--help clump`.

## Next implementation priorities

1. ~~Inspect and harmonize the GWAS schema~~ Implemented (`inspect_gwas.py`, configurable thresholds, original statistics retained); pending a run on HPC data. Loci tables are now available as input to step 2.
2. Join variants to loop anchors on either endpoint, domains, and compartment scores. Distinguish anchor overlap from being between anchors or merely nearby. Add promoter/gene mapping only after obtaining a compatible annotation.
3. Produce interpretable tables/plots of locus connectivity and reference-specific/shared annotations. Distinguish direct edges, shared partners, and connected components; significant SNPs in LD are not independent loci.
4. Inspect selected GWAS regions in both matrices at common coordinates/resolutions. Explore trans contacts separately from the cis loop list. Missing loop calls are not proof of absent contact.

Use modular Python CLI scripts with explicit paths/output directories, documented dependencies, logging, and small test fixtures. Prefer chunked GWAS reading and indexed genomic queries. Report assumptions and distinguish synthetic tests from HPC-data execution. Do not invent input columns, gene annotations, file provenance, or completed analyses.

## Documentation

- [ENCODE Hi-C outputs](https://www.encodeproject.org/hic/)
- [HiCCUPS loop fields](https://github.com/aidenlab/juicer/wiki/HiCCUPS)
- [Arrowhead domains](https://github.com/aidenlab/juicer/wiki/Arrowhead)

Project status and paths above are based on the user's HPC listings and completed QC run.


## sbatch format: 
# Working example: `run_clump_pipeline.sbatch`. The template below is the original reference.
# example below- adjust paths to project. slurmLogs dir will be created within TopoEDA

#!/bin/bash

#SBATCH --job-name=DOSE_XTRACT_PHEWAS       # job name (shows up in squeue)
#SBATCH --nodes=1                          # all plink2 calls run on 1 node
#SBATCH --ntasks=16                         # one task (serial / multithreaded)
#SBATCH --time=20:00:00                    # walltime (HH:MM:SS), adjust as needed
#SBATCH --output=/fslustre/labs/ext_rosen_joseph_mayo_edu/phewasHYE/slurmLogs/slurm_MBBPHE_dose_%j.out
#SBATCH --error=/fslustre/labs/ext_rosen_joseph_mayo_edu/phewasHYE/slurmLogs/slurm_MBBPHE_dose_%j.err
#SBATCH --partition=lg-n64-256g

# plink2 version- what we will use for clumping and general PGEN/sumstats manipulation
PLINK2=/usr/local/biotools/plink2/v2.00a6LM.2023.1123/plink2

