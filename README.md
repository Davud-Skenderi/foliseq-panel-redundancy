# Spatial panel redundancy and Foli-seq

Which genes should future Foli-seq panels target? This repository compares the gene panels of commercial
spatial transcriptomics platforms (10x Xenium, NanoString CosMx and GeoMx, Vizgen MERSCOPE, plus Illumina
AmpliSeq and IDT xGen for reference) with the current Foli-seq panels (mouse mPS3, human hPS2; Huang et al.,
*Nature Biotechnology* 2025). It looks for **consistency, redundancy and gaps** at the pathway level, then uses
public expression data and Foli-seq stool results to rank genes to keep, review and add.

Start with `results/Panel_redundancy_results.xlsx`, sheet **foli_candidate_genes**. It has one row per gene,
with all evidence columns and a Keep / Review / Add call.

## Main findings

- Vendors choose mostly different genes (Jaccard mostly 0.01–0.35) to probe largely the same pathways
  (pathway-profile similarity 0.5–0.9). A pathway can be tracked with about 7–8 well-chosen genes
  (held-out r ≈ 0.87). About 300 covariance-chosen genes reconstruct gut/IBD pathway scores as well as
  1,000–10,000-gene panels.
- Foli-seq mPS3/hPS2 cover cytokine signalling about as deeply as the largest spatial panels (83–91% of
  KEGG cytokine–receptor genes) and carry little internal redundancy (15–16% of genes predictable from the rest).
- 279 of 535 mPS3 genes changed in Foli-seq stool data in at least one colitis model (paper Tables S3–S7).
  Whether vendors include a gene does not predict this (55% response for genes on no spatial panel, 60% for 8–12).
- The shared gap is gut epithelium: enterocyte, goblet, enteroendocrine and stem programs are thinly covered by
  every panel, Foli-seq included. Candidate additions are grouped by cell type in the workbook.

Figures are in `figures/`.

## Layout

```
inputs/vendor_panels/      vendor gene lists as downloaded (MERSCOPE human IO files are empty)
inputs/foli_supplement/    Huang et al. 2025 Supplementary Tables S1–S9
data/                      harmonised panels (gene_by_panel.csv, panels_long.csv), Foli-seq panels, stool DEGs
results/                   all result tables + Panel_redundancy_results.xlsx
figures/                   fig1–fig9
scripts/                   pipeline (run in numeric order)
```

## Running it

Requires Python 3.11+ with pandas, numpy, scipy, scikit-learn, statsmodels, scanpy, anndata, pyarrow,
openpyxl, matplotlib and rdata.

```bash
pip install pandas numpy scipy scikit-learn statsmodels scanpy anndata pyarrow openpyxl matplotlib rdata
bash scripts/fetch_public_data.sh          # ~400 MB of public data into external/
for s in 00_foli_inputs 01_harmonize_panels 02_genesets 03_pathway_coverage 04_expression_metacells \
         05_covariance_redundancy 05b_pathway_core_heldout 06_optimize_panel; do python scripts/$s.py; done
FOLI_PANEL=FOLI_mPS3 python scripts/06b_foli_focus.py
FOLI_PANEL=FOLI_hPS2 SKIP_DENOVO=1 python scripts/06b_foli_focus.py
for s in 06c_foli_proxies 09_candidates 07_figures 08_workbook; do python scripts/$s.py; done
```

Steps 5 and 6 take roughly 30–60 minutes on a laptop. Step 5 reuses earlier vendor-panel results in `results/`
if present, so delete `results/panel_redundancy_by_dataset.csv` for a full recompute.

| Step | Script | What it does |
| --- | --- | --- |
| 0 | `00_foli_inputs.py` | Foli-seq panels and stool-responsive genes from the paper's supplement |
| 1 | `01_harmonize_panels.py` | parses every gene list; maps mouse/dog to human (HGNC + MGI/AGR orthologs) |
| 2 | `02_genesets.py` | MSigDB 7.5.1 gene sets plus a gut/IBD lens (221 sets) |
| 3 | `03_pathway_coverage.py` | pathway coverage and enrichment per panel; gene and pathway overlap between panels |
| 4 | `04_expression_metacells.py` | public scRNA-seq into metacells; HPA + GTEx tissue profiles |
| 5 | `05_covariance_redundancy.py` | within-panel redundancy, pathway fidelity, consensus gene–gene correlation |
| 5b | `05b_pathway_core_heldout.py` | split-half validated core genes per pathway |
| 6 | `06_optimize_panel.py` | covariance-optimised panel (generic objective) and vendor benchmark |
| 6b | `06b_foli_focus.py` | gut/IBD objective: unique contribution and ranked additions for a Foli-seq panel |
| 6c | `06c_foli_proxies.py` | Foli-seq genes on spatial panels, directly or via co-expressed proxies |
| 9 | `09_candidates.py` | per-gene evidence table and Keep / Review / Add calls |
| 7, 8 | `07_figures.py`, `08_workbook.py` | figures and results workbook |

Paths resolve relative to the repository; override with `PANEL_EXT`, `PANEL_RAW`, `PANEL_SUP` or `PANEL_PROJ`.

## Data sources

- Huang Y, et al. Fecal exfoliome sequencing captures immune dynamics of the healthy and inflamed gut.
  *Nature Biotechnology* (2025). doi:10.1038/s41587-025-02894-4 (Supplementary Tables S1–S9)
- Vendor panel gene lists from 10x Genomics, NanoString, Vizgen, Illumina and IDT (some are form-gated
  downloads; keep this repository private unless redistribution is cleared)
- MSigDB 7.5.1 via [msigdbr](https://github.com/cran/msigdbr); orthologs via [babelgene](https://github.com/cran/babelgene);
  symbol updates via [HGNChelper](https://github.com/cran/HGNChelper)
- HPA and GTEx tissue expression via [TissueEnrich](https://github.com/bioc/TissueEnrich)
- Public single-cell data from [YosefLab/scVI-data](https://github.com/YosefLab/scVI-data): 10x PBMC 10k,
  mouse spleen and lymph node, and a Human Heart Cell Atlas subsample

## Caveats

- Co-expression comes from blood, spleen/lymph node, heart and tissue-atlas data, not gut. Epithelial
  redundancy may differ; Foli-seq count data or a colon single-cell atlas would be the better source.
- Stool evidence is mouse-only, so human-only hPS2 genes have none.
- MSigDB intestinal cell-type signatures are approximate; some include immune genes.
