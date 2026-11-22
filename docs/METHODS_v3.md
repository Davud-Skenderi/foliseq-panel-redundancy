# What this analysis is, in plain terms (v3)

This file explains the *why* behind every step, for someone picking the project up cold.
It documents the v3 additions: gut references, the accuracy histogram, the corrected redundancy
figure, and the robustness test. For the v2 methods see the README table.

---

## 1. The problem the project exists to solve

Imaging-based spatial transcriptomics (10x Xenium, Vizgen MERSCOPE, NanoString/Bruker CosMx) measures
mRNA *in situ* by binding probes to chosen transcripts and reading them out optically. The number of
transcripts that can be distinguished is bounded by the optical encoding, so the genes must be chosen
in advance. **Genes left off a panel are never measured and cannot be recovered after the fact.**

Foli-seq (Huang et al., *Nat Biotechnol* 2025) sequences the fecal exfoliome — host cells shed into
stool — as a non-invasive readout of gut state. Deciding which genes a future Foli-seq panel should
target is the question this repo serves. Commercial panels are used as evidence: where vendors agree,
that is weak evidence of informativeness; where none of them look, that is a candidate gap.

## 2. The three ingredients, and why they are kept separate

Every number in this repo is produced by combining exactly three things:

1. **Panels** — the gene lists being compared (`data/gene_by_panel.csv`, a gene x panel 0/1 matrix).
2. **Metrics** — the definitions of "redundant" and "identifiable" (Spapros, Kuemmerle et al. 2024).
3. **Reference** — annotated single-cell data the metrics are computed *on*.

Holding (1) and (2) fixed while varying (3) is what turns a result into a **robustness test**. This
matters because a metric computed on one reference is not a property of the panel alone — it is a
property of the panel *as seen through that reference*.

## 3. Why a different reference shifts the whole distribution

A reference with more cell types, or with finer distinctions between them, is a harder classification
problem. Chance accuracy for `k` balanced classes is `1/k`, so moving from 26 immune cell types to 72
intestinal ones lowers every panel's score. **The entire distribution shifts left.** That shift is
expected and is not evidence that any panel got worse.

The question that matters is whether the **ordering** of panels survives the shift:

- **Ordering preserved** -> the comparison reflects the panels. The result is portable.
- **Ordering changes** -> the metric was partly measuring the reference. Any claim about which panel
  is "better" is only valid for the reference it was computed on.

A panel that moves *against* the general shift is the interesting case: it indicates content specific
to one tissue rather than general-purpose cell-typing power. For a gut assay that is exactly what you
want to detect.

This is quantified in `results/v2/q7_reference_rank_correlation.csv` (Spearman rank correlation of
per-panel mean accuracy between each pair of references) and drawn in
`figures/v3/v3_fig7_reference_robustness.png`.

## 4. The references

| Tag | Source | Tissue / state | Cell types |
| --- | --- | --- | --- |
| `MM_SLN` | Gayoso et al. 2021 CITE-seq | mouse spleen + lymph node | 26 immune |
| `HS_PBMC` | 10x Genomics PBMC 10k | human blood (unannotated) | redundancy only |
| `GUT_HEALTHY` | Gut Cell Atlas, Elmentaite et al. 2021 | human intestine, healthy | 72 |
| `GUT_IBD` | metaplasia in inflammatory gut disease | human intestine, inflamed | 54 |
| `GUT_CRC` | colorectal cancer immune hubs | human intestine, malignant | 45 |

The three gut references are built by `scripts/16_build_gut_refs.py` from the CELLxGENE Census
(release 2025-11-08), restricted to `tissue_general in {large intestine, small intestine}` and to
primary data (`is_primary_data == True`) so the same cells are not counted twice across datasets.
Cells are subsampled to 250 per annotated cell type; types with fewer than 40 cells are dropped
(matching Spapros `n_cells_min`). Raw counts are kept — normalisation happens downstream.

**The full measured transcriptome is retained, not just panel genes.** This matters for the random
baselines: a random gene set drawn from the union of commercial panels is not a null, because the
large panels *are* most of that union. Drawing from the transcriptome makes the baseline meaningful.

Census `var.feature_name` already carries HGNC symbols, so no ortholog mapping is applied to these
references and `common.py` is not imported by `16_`/`17_`. Mouse panels continue to be compared in
human-ortholog space via `data/gene_by_panel.csv`, which was harmonised upstream by `01_`.

## 5. The metrics, and what each one actually asks

### Q4 — cell-type identification (`spapros.evaluation.metrics.xgboost_forest_classification`)

*Question:* using only this panel's genes, can a classifier tell the reference's cell types apart?

Gradient-boosted trees, 5-fold cross-validation, confusion matrix normalised per true cell type so
that common and rare types count equally. **Per-cell-type accuracy is the diagonal of that matrix.**
A panel yields one accuracy per cell type — 72 numbers on `GUT_HEALTHY`.

*Why a histogram and not a bar of the mean.* v2's fig5 reduced those numbers to their mean. Two panels
with the same mean can behave completely differently: uniformly mediocre, versus excellent on most
types and failing on a few. The second is actionable — you can name the cell types that fail and add
markers for them. The first is not. `figures/v3/v3_fig5h_accuracy_distribution_*.png` plots the
distribution, with the mean drawn as a dashed line so the v2 quantity is still visible.

On `MM_SLN` this immediately shows **Foli-seq hPS2 is bimodal** — one group of cell types near 55%,
another near 85% — which its 70% mean concealed entirely.

### Q3 — redundancy (`spapros.evaluation.metrics.correlation_matrix`)

*Question:* how much of the panel is spent on genes that carry information another panel gene already
carries? Pairwise Pearson correlation on log-normalised counts; a pair is redundant at |r| >= 0.8.

*The bug this version fixes.* v2's fig4 plotted `n_pairs_r_ge_0_8`, the **raw count** of redundant
pairs. The number of *possible* pairs is `n(n-1)/2`, so a 4,686-gene panel has roughly 340x more
opportunities to produce a redundant pair than a 330-gene panel. The raw count therefore ranks panels
mostly by size. The code comment above that figure already said it was plotting "% genes without a
|r| >= 0.8 partner" — the Spapros definition, which `q3_redundancy_summary.csv` stores as
`pct_genes_nonredundant`. The figure did not plot it.

**Correcting this inverts the conclusion.** On mouse spleen/lymph node, CosMx Human 6K has the *most*
redundant pairs by raw count (31) and nearly the *fewest* per gene (0.47% of its genes have a
redundant partner). On human PBMC the CosMx IO 100-plex add-on has 18 pairs — unremarkable by count —
but 18.2% of its genes have a redundant partner, far the highest of any panel. That is the sensible
result: a small focused immuno-oncology panel is dense in co-regulated genes, while a discovery panel
spreads across the transcriptome.

`figures/v3/v3_fig4c_redundancy_normalised.png` plots both size-normalised forms: % of genes with a
redundant partner, and redundant pairs per 1,000 possible pairs.

## 6. Baselines

Random gene sets of matched size are the control: they say what a panel of that size achieves with no
curation at all. v2 used 535 genes x 3 draws. v3 adds 330 x 3 (the measured size of Foli-seq mPS3), so
each Foli-seq panel has a size-matched null rather than one borrowed from a different size.


## 9. A confound in the current reference set, and how to remove it

`MM_SLN` is **mouse** tissue. All three gut references are **human**. So the contrast currently
labelled "immune vs gut" is also "mouse vs human", and the two cannot be separated with these four
references alone.

The data shows why this matters. Across the 9 shared panels the mean shift from immune to gut is
small (-4.5 to +4.2 points), but the two Foli-seq panels move in **opposite directions**:

| Panel | MM_SLN (mouse immune) | mean of 3 human gut refs | shift |
| --- | --- | --- | --- |
| Foli-seq mPS3 (mouse panel) | 76.9% | 72.8% | **-4.1** |
| Foli-seq hPS2 (human panel) | 69.5% | 73.7% | **+4.2** |

Each panel does better on the reference matching its own species. That is exactly the signature of a
species effect rather than a tissue effect, and it is the single most likely alternative explanation
for any tissue-specific claim made from these four references.

**Adding Haber et al. 2017 (GSE92332, mouse small intestinal epithelium) resolves it.** With a mouse
gut reference the design becomes 2x2 — mouse/human x immune/gut — and the tissue effect can be
estimated with species held fixed. The repo README already named Haber as planned work; this analysis
turns that from a nice-to-have into the next necessary step.

## 10. Curation beats size-matched random, cleanly

With random baselines drawn from the **transcriptome** (not from the union of commercial panels),
a 330-gene random set is far behind the curated panels of comparable size:

| Reference | random 330 (mean of 3) | Foli-seq mPS3 | gap |
| --- | --- | --- | --- |
| Gut Cell Atlas, healthy | 43.2% | 69.7% | +26.5 |
| Inflamed gut | 45.2% | 70.9% | +25.7 |
| Colorectal cancer | 46.7% | 77.9% | +31.1 |

This is worth stating carefully. A null drawn from the union of commercial panel genes would be a far
weaker test, because the large panels make up most of that union — a random draw from it is mostly
panel genes. Drawing from the transcriptome is what makes the comparison meaningful.

## 7. Known limitations

- **Subsampling.** 250 cells per type (Spapros then uses up to 200). Adequate for per-type accuracy,
  but rare types rest on fewer cells and their accuracies are noisier.
- **One seed.** v2 set `n_seeds=1` instead of the Spapros default of 5, for compute. Per-cell-type
  accuracies therefore carry run-to-run variance that is not quantified here.
- **Cross-species comparison.** Mouse panels are mapped to human orthologs. Orthology is imperfect and
  one-to-many mappings are resolved by priority (MGI, then AGR, then babelgene support count).
- **Cell-type labels are the references' own.** "72 cell types" reflects the Gut Cell Atlas authors'
  annotation granularity, not a biological constant. A reference annotated more coarsely would raise
  every score.
- **The |r| >= 0.8 threshold is a convention**, taken from Spapros. `q3_correlated_pairs_*.csv` stores
  only pairs above it, so the distribution of r below 0.8 cannot be recovered from those files without
  recomputing the correlation matrices.
- **Species is confounded with tissue** in the current reference set; see section 9.
- **Cache keys.** `12_spapros_eval.py` caches confusion matrices as `q4_confusion_{panel}.csv` with no
  reference in the filename. Running a second reference through that script would silently reuse the
  first reference's results. `17_multiref_eval.py` keys its cache `q4_confusion_{REF}__{panel}.csv`.

## 8. What this does not yet answer

The gut references make epithelial conclusions possible for the first time, but nothing here yet
tests the host-microbiome interface specifically — mucins, antimicrobial peptides, pattern-recognition
receptors, tight junctions, IgA transport. That gene set is the natural next analysis: quantify how
little of it the commercial panels carry, which is the strongest available argument for Foli-seq
targeting content no vendor panel provides.
