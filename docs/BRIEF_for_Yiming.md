# Panel comparison — update for Yiming

Four things, in the order you asked for them. Figures are in `figures/v3/`.

---

## 1. The accuracy figure, as a distribution

**Figure:** `v3_fig5h_accuracy_distribution_MM_SLN.png`

The old fig5 plotted each panel's *mean* accuracy across cell types. A mean can't distinguish a panel
that is uniformly mediocre from one that is excellent on most cell types and fails on a few — and
those mean very different things for design.

Plotted as a distribution, **Foli-seq hPS2 is bimodal on mouse spleen/lymph node**: one group of cell
types near 55%, another near 85%. Its "70%" was an average of two behaviours. The second case is
actionable — you can name the failing cell types and add markers.

---

## 2. Gut single-cell references, and whether the comparison holds up

**Figure:** `v3_fig7_reference_robustness.png`

The panel comparison was running on immune tissue only — mouse spleen/lymph node and human PBMC — so
nothing epithelial had ever been tested against real gut data.

Added, all from the CELLxGENE Census (release 2025-11-08), primary data only:

| Reference | Tissue | Cell types |
| --- | --- | --- |
| Gut Cell Atlas (Elmentaite 2021) | healthy human intestine | 72 |
| Metaplasia / inflammatory gut disease | inflamed | 54 |
| Colorectal cancer immune hubs | malignant | 45 |
| Haber 2017 (GSE92332) | mouse small intestine | 15 |

**The panel ordering survives the change of reference**: all pairwise Spearman correlations of
per-panel mean accuracy fall between **0.92 and 0.97**. Absolute scores shift — a reference with more
cell types is a harder problem — but which panel beats which does not.

Panels are scored as **lift over a size-matched random gene set measured on the same reference**,
which is how panel size is normalised. The baselines confirm the references differ in difficulty:
random 535 genes score 57.9% on mouse immune, 52.2% on mouse gut epithelium, 48.2% on the 72-type
human gut.

Curation clearly works: a 330-gene random set scores 43–47% on gut where Foli-seq mPS3 scores 70–78%.

---

## 3. What the panels carry at the host-microbiome interface

**Data:** `results/v2/q11_interface_coverage.csv`, `q11_interface_orphans.csv`

This is the part that speaks directly to what Foli-seq should target. It uses no classifier and no
normalisation — it is counting genes against a 154-gene interface set (antimicrobial peptides,
mucins, pattern-recognition receptors, tight junctions, IgA transport, trefoil/goblet, Paneth
markers), denominated by what the CosMx Whole Transcriptome panel can see.

| Panel | size | interface genes | **% of its own gene budget** |
| --- | --- | --- | --- |
| Xenium Mouse 5K | 4,949 | 57 | 1.15 |
| CosMx Human 6K | 6,338 | 57 | 0.90 |
| Xenium Human 5K | 5,000 | 52 | 1.04 |
| **Foli-seq hPS2** | **543** | **42** | **7.73** |
| GeoMx CTA | 1,890 | 34 | 1.80 |
| **Foli-seq mPS3** | **535** | **29** | **5.42** |
| MERSCOPE PanNeuro | 812 | 3 | 0.37 |

**The Foli-seq panels spend 5–8% of their budget on interface biology; every commercial panel spends
0.4–2.9%** — roughly sevenfold higher density. hPS2 carries 42 interface genes in 543 slots where
Xenium Human 5K carries 52 in 5,000.

**55 of the 154 interface genes are on no targeted panel at all**, concentrated in:
- 30 of 45 antimicrobial peptides — nearly the whole beta-defensin family
- 13 of 21 mucins — including **MUC3A, MUC12, MUC13, MUC17**, the membrane-tethered mucins of the
  intestinal brush border

**One gap in our own panels:** hPS2 carries **1 of 34** tight junction / barrier genes; mPS3 carries
4; Xenium Mouse 5K carries 20. For an assay reading cells shed from an epithelial barrier, that looks
like the most actionable item for the next revision.

---

## 4. A correction to fig4 that should be made before it goes anywhere

**Figure:** `v3_fig4c_redundancy_normalised.png`

The old fig4 plotted the **raw count** of gene pairs with |r| ≥ 0.8. Possible pairs grow as
n(n−1)/2, so a 4,686-gene panel has roughly 340× more chances to produce one than a 330-gene panel —
the count ranks panels mostly by size. The code comment above the figure already described the
size-normalised Spapros definition, and `q3_redundancy_summary.csv` already stores it as
`pct_genes_nonredundant`; it just wasn't the quantity plotted.

**Normalised, the ordering inverts.** CosMx Human 6K has the most redundant pairs by raw count (31)
and nearly the fewest per gene (0.47%). The CosMx IO 100-plex add-on looks unremarkable at 18 pairs
but has **18.2%** of its genes carrying a redundant partner — the highest of any panel. Which makes
sense: a small focused immuno-oncology panel is dense in co-regulated genes; a discovery panel
spreads thin.

---

## Questions

1. **Which pathway database** for the Q2-style analysis — GO, Reactome, KEGG, MSigDB Hallmark?
2. **Is the 154-gene interface set the right one?** It is my construction and is defined explicitly
   in `scripts/25_interface_coverage.py`. Changing it changes the 55-gene gap list.
3. **What resolution does Foli-seq need?** Fifteen-way epithelial subtyping may be finer than the
   assay requires, which would mean the identification metric is calibrated too strictly.
4. **Mouse or human as primary?** It matters — see the reserve note below.

---

## Held in reserve (only if asked)

- **"How do you know that isn't just the reference you chose?"** The design is now a 2×2 —
  mouse/human × immune/gut — so each effect is estimated twice with the other factor fixed. The
  tissue effect replicates (r = +0.97, sign agreeing 7/9). A general species effect does not
  (r = −0.29).
- **"Is that difference real?"** Classifiers now run at 5 seeds on three of the four grid cells;
  standard errors are 1.4–2.6 points. Treat anything under ~3 points as noise.
- **Species transfer.** mPS3 beats hPS2 by +4.2 ± 2.2 points on mouse gut epithelium and loses by
  −7.2 ± 2.6 on human. An 11.4-point swing, ~3 combined SE. **mPS3's mouse validation should not be
  assumed to transfer to human samples.**
- **Two panels score below random on gut epithelium.** CosMx IO −14.1 ± 1.7 (mouse) and −19.7 ± 2.4
  (human); MERSCOPE Brain similar. Not broken panels — specialised ones outside their domain. CosMx
  IO scores 71.5% on human immune tissue with 94 genes.
- **A targeted panel gets most of the way to measuring everything.** Against the CosMx Whole
  Transcriptome ceiling, hPS2 reaches **90.3%** of the achievable range on human gut epithelium with
  543 genes — 85.2% against the whole transcriptome's 87.5%.
- **Cross-species caution.** Ortholog mapping collapses 19 mouse alpha-defensins into a single human
  DEFA1, and MHC class I 15-to-1. Mucins, PRRs, tight junctions and IgA transport show no collapses,
  so those comparisons are safe and antimicrobial-peptide comparisons are not.
- **Two bugs in existing code.** `12_spapros_eval.py` caches confusion matrices without the reference
  in the key, so a second reference run through it would silently reload the first one's results; and
  it discards the classifier's standard-deviation output, which is why nothing had error bars.
