# Spatial panel comparison for Foli-seq panel design

This repo compares the gene panels of commercial spatial transcriptomics platforms (10x Xenium, NanoString CosMx and
GeoMx, Vizgen MERSCOPE) with the current Foli-seq panels (mouse mPS3, human hPS2; Huang et al., *Nature
Biotechnology* 2025) to inform which genes future Foli-seq panels should target.

Each step answers one plain question with an established method, and every output is a gene list, a count or an
accuracy. The evaluation of panels follows the published framework of **Spapros** (Kuemmerle et al. 2024, *Nature
Methods*, "Probe set selection for targeted spatial transcriptomics") and uses its code.

| Question | Method | Reference | Main output |
| --- | --- | --- | --- |
| Q1 Consistency: which genes do vendors agree on? | Count vendors per gene; UpSet plot of set overlaps | Lex et al. 2014 (UpSet) | `results/v2/q1_consensus_genes.csv` |
| Q2 Gaps: which cell types and pathways does each panel cover? | Count signature/marker genes on each panel; Hallmark over-representation (hypergeometric, BH FDR) | Busslinger et al. 2021; Gao et al. 2018; Liberzon et al. 2015; Wolf et al. 2018 | `results/v2/q2*_*.csv` |
| Q3 Redundancy: which panel genes carry the same information? | Pairwise Pearson correlation; pairs with \|r\| ≥ 0.8 | Spapros `gene_corr` | `results/v2/q3_correlated_pairs.csv` |
| Q4 Function: can the panel identify cell types? | XGBoost classifier on panel genes, 5-fold cross-validation, per-cell-type accuracy | Spapros `forest_clfs`; Chen & Guestrin 2016 | `results/v2/q4_accuracy_by_celltype.csv` |
| Q4b Do the consensus genes help? | Same classifier on mPS3 + the vendor-consensus genes it lacks, against mPS3 + the same number of random genes | Spapros `forest_clfs` | `results/v2/q4b_addition_test.csv` |
| Q5 Stool relevance | Foli-seq stool DEG tables | Huang et al. 2025, Tables S3–S7 | `data/foli_stool_degs.csv` |
| Q6 Recommendations | Fixed rules over Q1–Q5 evidence (no weighted score) | — | `results/v2/q6_recommendations.csv` |

Figures are in `figures/v2/`. A full methods guide (each step's rationale, citation, result and limits) is in the Claude Doc "Foli-seq panel design: methods guide (v2)". The first version of this analysis (cosine similarity of enrichment profiles, PCA–ridge
reconstruction, weighted priority score) is kept in `archive/v1/` for the record and is superseded.

## Reference data and its limits

- **Q3/Q4 currently use immune references**: mouse spleen + lymph node CITE-seq with 26 annotated cell types
  (Gayoso et al. 2021, *Nature Methods*) and human PBMC 10k (10x Genomics; redundancy only, unannotated).
- **Gut references still needed.** The same scripts accept any annotated `.h5ad`. Planned: the mouse small-intestinal
  epithelium atlas of Haber et al. 2017 (*Nature*; GEO GSE92332, `GSE92332_atlas_UMIcounts.txt.gz`) and Foli-seq
  count data. Until then, epithelial conclusions rest on Q2 (gene lists) only.
- Stool evidence (Q5) is mouse-only.

## Running

```bash
pip install pandas numpy scipy statsmodels scanpy anndata pyarrow openpyxl matplotlib rdata upsetplot xgboost
pip install --no-deps spapros && pip install rich venndata
bash scripts/fetch_public_data.sh
python scripts/00_foli_inputs.py && python scripts/01_harmonize_panels.py && python scripts/02_genesets.py
python scripts/10_consistency.py && python scripts/12_spapros_eval.py && python scripts/11_coverage.py
python scripts/15_addition_test.py && python scripts/13_recommendations.py && python scripts/14_figures_v2.py
```

`12_spapros_eval.py` takes about an hour on a laptop (cached per panel, so it can be resumed).

## Layout

```
inputs/vendor_panels/      vendor gene lists as received (MERSCOPE human IO files are empty)
inputs/foli_supplement/    Huang et al. 2025 Supplementary Tables S1–S9 (Foli-seq panels, stool DEGs)
data/                      harmonised gene x panel table, Foli-seq panels, stool DEGs
scripts/                   00–02 inputs; 10–15 analysis (Q1–Q6) and figures
results/v2/, figures/v2/   current results
archive/v1/                superseded first version
```

## Key references

- Huang Y, et al. Fecal exfoliome sequencing captures immune dynamics of the healthy and inflamed gut. *Nat Biotechnol* (2025).
- Kuemmerle LB, et al. Probe set selection for targeted spatial transcriptomics. *Nat Methods* (2024).
- Lex A, et al. UpSet: visualization of intersecting sets. *IEEE TVCG* (2014).
- Liberzon A, et al. The Molecular Signatures Database hallmark gene set collection. *Cell Syst* (2015).
- Busslinger GA, et al. Human gastrointestinal epithelia of the esophagus, stomach, and duodenum resolved at single-cell resolution. *Cell Rep* (2021).
- Gao S, et al. Tracing the temporal-spatial transcriptome landscapes of the human fetal digestive tract using single-cell RNA-sequencing. *Nat Cell Biol* (2018).
- Gayoso A, et al. Joint probabilistic modeling of single-cell multi-omic data with totalVI. *Nat Methods* (2021).
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* (2018).
- Soneson C, Robinson MD. Bias, robustness and scalability in single-cell differential expression analysis. *Nat Methods* (2018).
- Chen T, Guestrin C. XGBoost: a scalable tree boosting system. *KDD* (2016).
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B* (1995).
- Haber AL, et al. A single-cell survey of the small intestinal epithelium. *Nature* (2017).
