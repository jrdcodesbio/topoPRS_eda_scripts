# topoEDA: IPF GWAS x chromatin architecture

## Purpose and current scope

Exploratory integration of the **Allen 3-way IPF GWAS (hg38/GRCh38)** with two ENCODE Hi-C references:

- **Lung:** `ENCSR305GRJ`.
- **Primary CD8 T cells:** `ENCSR321BHC`. Do not label this as the previously considered activated-CD8 experiment.

Explore where associated variants lie, which regions/promoters they contact, and how their structural context differs between references. This is **EDA and technique development**, not TopoPRS training, a new GWAS, or patient-level prediction. These reference contexts are not IPF cases versus controls. Large differences between maps are expected; document them without blocking exploration or repeatedly focusing on loop counts.

## Execution environment and paths

**Write code for deployment and execution on the HPC system where the data reside.**  The dev environment cannot access these paths. I will pusn to a git repo and pull the modified code on the HPC from there. 

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

Paths are case-sensitive. GWAS is approximately 394 MB compressed; lung matrix approximately 46 GB. The T-cell `.hic` and `.hic.1` are approximately 1.6 GB and 34 GB, respectively. **Determine which download is complete using expected file size/checksum; do not select by filename/size alone or delete either automatically.**

Keep scripts under `scripts/`; create `results/` for new outputs. Keep `sourceData/` unchanged. Existing QC outputs were written outside the project to `/fslustre/labs/ext_rosen_joseph_mayo_edu/hic_qc/`.

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

`scripts/qc_parser.py` successfully parsed all four BEDPEs:

| Reference | Loops | Domains |
|---|---:|---:|
| Lung | 729 | 945 |
| T cell | 74,287 | 1,562 |

Loop tables have 33 columns; domain tables have 16. No exact duplicate pairs; all listed contacts are cis, and all domains are diagonal. No additional significance filter was applied. Neither matrices nor bigWigs have been inspected. Zero trans calls in BEDPE does not establish absence of trans contacts in `.hic`.

## Next implementation priorities

1. Inspect and harmonize the GWAS schema; report usable variants and association distributions. Make significance thresholds configurable and retain original statistics.
2. Join variants to loop anchors on either endpoint, domains, and compartment scores. Distinguish anchor overlap from being between anchors or merely nearby. Add promoter/gene mapping only after obtaining a compatible annotation.
3. Produce interpretable tables/plots of locus connectivity and reference-specific/shared annotations. Distinguish direct edges, shared partners, and connected components; significant SNPs in LD are not independent loci.
4. Inspect selected GWAS regions in both matrices at common coordinates/resolutions. Explore trans contacts separately from the cis loop list. Missing loop calls are not proof of absent contact.

Use modular Python CLI scripts with explicit paths/output directories, documented dependencies, logging, and small test fixtures. Prefer chunked GWAS reading and indexed genomic queries. Report assumptions and distinguish synthetic tests from HPC-data execution. Do not invent input columns, gene annotations, file provenance, or completed analyses.

## Documentation

- [ENCODE Hi-C outputs](https://www.encodeproject.org/hic/)
- [HiCCUPS loop fields](https://github.com/aidenlab/juicer/wiki/HiCCUPS)
- [Arrowhead domains](https://github.com/aidenlab/juicer/wiki/Arrowhead)

Project status and paths above are based on the user's HPC listings and completed QC run.
