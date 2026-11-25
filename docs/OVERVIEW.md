# Foli-seq panel design: what this project is, and everything in it

This document assumes no prior context. It starts from why the project exists and builds up to what
has been found and what is still open. `docs/METHODS_v3.md` holds the technical detail; this is the
explanation.

---

## 1. The biology, from the start

Every cell in the body carries the same DNA but uses different parts of it. A gene being "expressed"
means the cell is actively transcribing it into mRNA. Measuring which mRNAs are present, and how
many, tells you what kind of cell you are looking at and what it is doing. A neuron and a gut lining
cell have identical genomes and completely different mRNA profiles.

**Single-cell sequencing** dissociates a tissue into individual cells and measures the mRNA in each
one. You learn which cell types are present and what each is doing. The cost is that dissociation
destroys the tissue architecture: you know a cell was in the sample, but not where it sat or what it
was touching.

**Spatial transcriptomics** keeps the position. Two families of method exist, and the difference
between them is why this project exists.

*Capture-based* methods lay tissue on a surface covered in spatially barcoded capture spots. mRNA
binds, gets sequenced, and the barcode records where it came from. You get the whole transcriptome
with no genes chosen in advance. The cost is resolution (a spot may cover several cells) and capture
efficiency. Visium, Stereo-seq and Curio Seeker work this way.

*Imaging-based* methods are the opposite trade. Short DNA probes bind specific mRNA molecules, carry
fluorescent tags, and the tissue is photographed over many rounds. Each gene gets a distinct optical
barcode — a pattern of which rounds and colours it lights up in. You see individual mRNA molecules at
their real position inside individual cells. Xenium (10x Genomics), MERSCOPE (Vizgen) and CosMx
(Bruker) work this way.

**The constraint that drives everything below:** the number of distinguishable optical barcodes is
bounded by how many imaging rounds and colour channels the instrument can run. So the genes must be
chosen before the experiment. That chosen list is a **panel**, and genes left off it are never
measured and cannot be recovered afterwards. Panel design is one of the few decisions in genomics
that is genuinely irreversible.

## 2. What Foli-seq is

The intestine continuously replaces its lining. The discarded cells — the **exfoliome** — leave the
body in stool. Foli-seq sequences those cells, giving a readout of gut state **without a biopsy**:
no colonoscopy, no anaesthesia, just a stool sample.

If Foli-seq is to be a targeted assay, someone must choose its gene list. That is the question this
repository serves.

The approach is to learn from what commercial vendors already chose. Vendors have spent real money
deciding which genes are worth a probe slot. Where they agree, that is weak evidence those genes are
informative. Where none of them look, that is a candidate gap.

## 3. The three ingredients

Every number here comes from combining exactly three things:

1. **Panels** — the gene lists being compared (`data/gene_by_panel.csv`, a gene x panel 0/1 matrix).
2. **Metrics** — the definitions of "redundant" and "identifiable", taken from Spapros
   (Kuemmerle et al. 2024), a published framework for evaluating targeted panels.
3. **Reference** — annotated single-cell data the metrics are computed *on*.

The crucial point: **a metric computed on one reference is not a property of the panel.** It is a
property of the panel *as seen through that reference*. "This panel achieves 77% accuracy" really
means "77% at distinguishing the 26 cell types annotated in this particular dataset."

Holding panels and metrics fixed while varying the reference turns a result into a robustness test.
That is what most of the work below is doing.

## 4. What the metrics actually ask

**Q4, cell-type identification.** Hide every gene except the panel's, then train a classifier
(gradient-boosted trees, 5-fold cross-validation) to predict each cell's annotated type. The
confusion matrix is normalised per true cell type so common and rare types count equally, and the
diagonal gives one accuracy per cell type.

This is **not** proposing that Foli-seq do cell typing. It is a measuring instrument for information
content. The question you want answered — do these genes carry the biological signal that
distinguishes these cells? — is abstract. The classifier makes it concrete: if cell identity can be
recovered from the panel alone, the information is present; if not, the panel is blind to that
variation and no downstream analysis will recover what was never measured.

In a reference dataset the annotations are known — that is what makes it a reference, and they serve
as an answer key for scoring. In a real experiment you have no answer key, and the classifier is
measuring whether inference would be possible.

**Q3, redundancy.** Pairwise correlation between genes within a panel. If two genes always rise and
fall together, the second is buying little. A pair counts as redundant at |r| >= 0.8.

## 5. Why raw accuracy cannot be compared across references, and what replaces it

Each reference poses a different problem: 26 immune cell types, 37 immune, 15 epithelial, 12
epithelial, 72 whole-tissue. Guessing among 26 options is harder than among 12, so a panel scoring
70% on one and 65% on another tells you nothing — the second may simply be harder.

The fix is to measure every panel against a **random gene set of the same size, evaluated on the same
reference**:

    lift = panel accuracy - random-535 accuracy, on the same reference

Lift asks "how much did curating these genes buy, on this problem", and subtracts the problem's
difficulty from both sides. The baselines confirm the difficulties genuinely differ: random 535 genes
score 57.9% on mouse immune, 52.2% on mouse gut epithelium, 63.2% on human gut epithelium, 48.2% on
the 72-type human gut.

**Lift is imperfect and we measured how.** Across references, lift correlates **+0.78** with the
number of cell types — more classes leaves more headroom for curation to fill. A headroom-normalised
alternative still correlates +0.63. Neither fully removes the dependence, which is why §8 runs an
explicit control.

**Where a ceiling exists, there is a better scale still.** The CosMx Whole Transcriptome panels
measure essentially everything the platform can, so their accuracy bounds what any panel on that
instrument could reach. With a floor and a ceiling:

    ceiling fraction = (panel - random) / (whole transcriptome - random)

0 means no better than random; 100 means as good as measuring everything. Unlike lift, this has a
meaningful maximum.

## 6. The references

| Tag | Source | Species | Tissue | Cell types | Seeds |
| --- | --- | --- | --- | --- | --- |
| `MM_SLN` | Gayoso et al. CITE-seq | Mouse | Immune (spleen/lymph node) | 26 | 5 |
| `HUMAN_IMMUNE` | Human Immune Health Atlas | Human | Immune (blood) | 37 | 1 |
| `HABER_GUT` | Haber et al. 2017, GSE92332 | Mouse | Gut epithelium | 15 | 5 |
| `GUT_HEALTHY_EPI` | Gut Cell Atlas, epithelial subset | Human | Gut epithelium | 12 | 5 |
| `GUT_HEALTHY` | Gut Cell Atlas | Human | Gut, whole tissue | 72 | 1 |
| `GUT_IBD` | metaplasia in inflammatory gut disease | Human | Gut, inflamed | 54 | 1 |
| `GUT_CRC` | colorectal cancer immune hubs | Human | Gut, malignant | 45 | 1 |

All but Haber come from the CELLxGENE Census (release 2025-11-08), restricted to primary data so no
cell is counted twice. Each is subsampled to 250 cells per annotated cell type, types under 40 cells
dropped, raw counts kept, and **the full measured transcriptome retained** — the last point matters,
because a random baseline drawn from the union of commercial panel genes is not a null when the large
panels make up most of that union.

Mouse references are mapped to human symbols with the repo's own `common.to_human`, the same mapping
applied throughout, so references differ by biology rather than by mapping method.

## 7. Why the grid shape matters

The first version of this work compared a mouse immune reference against human gut references. That
cannot work, because every immune reference was mouse and every gut reference human: "immune vs gut"
and "mouse vs human" moved together, and no amount of analysis separates two things that never varied
independently.

Haber (mouse, gut) and the Human Immune Health Atlas (human, immune) complete the grid:

|  | Immune | Gut epithelium |
| --- | --- | --- |
| **Mouse** | MM_SLN (26) | HABER_GUT (15) |
| **Human** | HUMAN_IMMUNE (37) | GUT_HEALTHY_EPI (12) |

Now each effect can be estimated **twice**, holding the other factor fixed:

- tissue effect in mouse = HABER_GUT − MM_SLN
- tissue effect in human = GUT_HEALTHY_EPI − HUMAN_IMMUNE
- species effect on immune = HUMAN_IMMUNE − MM_SLN
- species effect on gut = GUT_HEALTHY_EPI − HABER_GUT

**And replication is a sharper test than either estimate.** The class-count artifact pushes in a
direction set by whether class count rises or falls. For the two tissue estimates it falls both times
(26→15, 37→12), so the artifact could drive both. For the two species estimates it **rises** once
(26→37) and falls once (15→12) — opposite directions. So a species effect with the same sign in both
rows could not be manufactured by the artifact.

## 8. What was found

### Panel ordering is robust to the choice of reference
All pairwise Spearman correlations of per-panel mean accuracy across references fall between
**0.92 and 0.97**. Changing the reference moves everyone's absolute score and barely touches the
ranking.

### Curated panels beat size-matched random, clearly
A 330-gene random set scores 43–47% on gut references where Foli-seq mPS3 scores 70–78%. Gaps of
+26 to +31 points.

### The tissue effect is real, replicates, and is absent for the Foli-seq panels
Moving from immune tissue to gut epithelium, panels lose lift. The two independent estimates — one
computed entirely within mouse, one entirely within human — correlate **r = +0.97**, with the sign
agreeing for 7 of 9 panels. This is not an artifact of one reference.

Re-running MM_SLN at 5 seeds shrank several magnitudes, so the single-run numbers had been somewhat
inflated: CosMx IO from -20.4 to -18.2, CosMx Human UCC from -10.1 to -7.2, CosMx Mouse Neuro from
-3.1 to -0.2. **Both Foli-seq panels land at or above zero** (mPS3 -2.5 in mouse and -0.2 in human;
hPS2 +0.2 and +6.5), meaning they pay no epithelium penalty at all and hPS2 is actually better on
gut epithelium than on immune tissue.

But magnitude matters, and a control was needed. Because lift correlates with class count, and the
tissue comparison reduces class count, the effect's direction was also the artifact's direction.
`scripts/22_celltype_count_control.py` re-ran MM_SLN restricted to a random 15 of its 26 cell types
over 3 draws, matching Haber exactly. The direction survived for all nine panels, but about a third
of the magnitude was artifact:

| Panel | uncontrolled | **controlled** | sd |
| --- | --- | --- | --- |
| CosMx IO (100) | −21.0 | **−21.2** | 2.98 |
| MERSCOPE Brain (128) | −12.6 | **−12.3** | 0.58 |
| MERSCOPE Pan-Mouse (815) | −12.8 | **−10.1** | 1.80 |
| Foli-seq mPS3 | −5.3 | **−2.5** | 1.54 |
| Foli-seq hPS2 | −2.0 | **−0.3** | 2.30 |

**The defensible statement is therefore not "all panels are immune-weighted" but "panels curated for
other tissues lose power on gut epithelium, and the Foli-seq panels do not."** hPS2's interval
comfortably includes zero.

### Two panels perform worse than random on gut epithelium
CosMx IO 100-plex: lift **−14.1 ± 1.7** on mouse and **−19.7 ± 2.4** on human gut epithelium.
MERSCOPE Brain: −12.2 ± 1.6 and −13.2 ± 2.2. Negative lift means a random gene set of the same size
does better. This is not a broken panel but a specialised one outside its domain — CosMx IO scores
**71.5%** on human immune tissue with 94 genes, matching the Foli-seq panels on their own ground.

### The Foli-seq panels each favour their own species — on gut tissue
With 5 seeds and error bars:

| | mouse gut epithelium | human gut epithelium |
| --- | --- | --- |
| mPS3 − hPS2 | **+4.2 ± 2.2** | **−7.2 ± 2.6** |

An 11.4-point swing, roughly 3 combined standard errors. Each panel does better on the species it was
built for.

**But this does not generalise to panels at large.** The two species estimates correlate **r = −0.28**
— they do not replicate. So there is no general species effect; what exists is specific to the
Foli-seq pair on gut tissue, where it is well measured. On immune tissue the same comparison gives
+7.4 (mouse) and −0.6 (human): the direction is consistent but weaker.

The practical consequence stands: **mPS3's mouse validation should not be assumed to transfer to
human samples.**

### A targeted panel gets most of the way to measuring everything
With the whole-transcriptome ceiling: **Foli-seq hPS2 reaches 90.3% of the ceiling on human gut
epithelium with 543 genes** — 85.2% accuracy against the whole transcriptome's 87.5%, a 2.3-point gap
for roughly a 26-fold reduction in genes measured.

### Redundancy ranks by panel size unless normalised
The v2 figure plotted the raw count of gene pairs with |r| >= 0.8. Possible pairs grow as n(n−1)/2, so
a 4,686-gene panel has ~340x more chances to produce one than a 330-gene panel. Normalising **inverts
the ordering**: CosMx 6K is near the least redundant per gene (0.47%) and the CosMx IO 100-plex the
most (18.2% on human PBMC). Small focused panels are dense in co-regulated genes; discovery panels
spread thin.

### Where the panels are thin at the host-microbiome interface
Against a 154-gene interface set (antimicrobial peptides, mucins, pattern-recognition receptors,
tight junctions, IgA transport, trefoil/goblet, Paneth markers), denominated by what CosMx WTx can
see:

- **The Foli-seq panels spend 5–8% of their gene budget there; every commercial panel spends
  0.4–2.9%.** hPS2 carries 42 interface genes in 543 slots; Xenium Human 5K carries 52 in 5,000.
- **55 of the 154 genes are on no targeted panel at all** — most beta-defensins, and the
  membrane-tethered intestinal mucins **MUC3A, MUC12, MUC13, MUC17**.
- **hPS2 carries 1 of 34 tight junction / barrier genes** (mPS3 carries 4; Xenium Mouse 5K carries
  20). For an assay reading cells shed from a barrier, this is the most actionable gap.

This analysis uses no classifier, no normalisation and no metric choice. It is counting genes, which
is why it is the hardest result here to argue with.

### Where ortholog mapping distorts, precisely
194 human symbols absorb two or more mouse genes on Haber, 383 extra genes in all, 3.1% of those
mapped. Modest in aggregate, heavily concentrated: **19 mouse alpha-defensins collapse into a single
human DEFA1**, and MHC class I collapses 15-to-1. Mucins, pattern-recognition receptors, tight
junctions and IgA transport show **zero** collapses.

So cross-species claims about those four categories are safe; claims about antimicrobial peptides
are not.

### Haber is not batch-confounded
Cell type and mouse of origin are nearly independent: Cramér's V = 0.104, normalised MI = 0.020. No
cell type draws more than 28% of its cells from a single mouse across 10 mice (median 23%, even split
would be 10%). The concern that the classifier might be reading donor rather than cell identity does
not hold.

## 9. Mistakes made and corrected

Recorded because they shaped the conclusions.

- **"5K/6K panels are no better than random" — retracted.** The random sets had been drawn from the
  union of commercial panel genes, and the 6K panel is 64% of that pool, so it competed against
  itself. Baselines are now drawn from the transcriptome.
- **A figure with mislabelled rows.** Two sub-plots sorted differently shared a y-axis, so one
  column's labels overwrote the other's. Caught only by looking at the rendered image.
- **A silently failed patch.** A "5-seed" run was configured but the call site still read
  `n_seeds=1`, so results were bit-identical 1-seed numbers written under a 5-seed filename. Caught by
  noticing the timing had not changed. Those caches were deleted and the run redone.
- **A silent out-of-memory kill.** The redundancy step builds an n x n correlation matrix; at
  whole-transcriptome size it exhausted memory and killed the process **after** the classification had
  been cached, with no error in the log. Now skipped above `CORR_MAX_GENES = 8000`.
- **Answering the wrong question first.** Cell-type analysis was built when pathway analysis had been
  asked for — and the repo had already done pathway analysis in Q2.

Two bugs were found in pre-existing code: `12_spapros_eval.py` caches confusion matrices as
`q4_confusion_{panel}.csv` with no reference in the key, so a second reference run through it would
silently reload the first one's results; and it discards `out[1]` from the Spapros classifier, which
is the standard-deviation matrix — the reason no number in this repo had an error bar until now.

## 10. What is still open

**Unquantified.** `HUMAN_IMMUNE` and the three whole-tissue gut references are still at `n_seeds=1`;
`MM_SLN`, `HABER_GUT` and `GUT_HEALTHY_EPI` are at 5. So the mouse row of the grid is internally
consistent and the human row is not yet. The whole-transcriptome
ceiling is a single run, and the ceiling fraction divides by it.

**Incomplete runs.** `GEOMX_H_CTA`, `GEOMX_DOG_CCA`, `XEN_M5K`, `XEN_H5K` and `COSMX_H6K` have not run
on the human gut references; the WTx ceiling exists only on the two epithelial references. Everything
caches per (reference, panel), so runs resume.

**Design questions only the lab can settle.** Whether the 154-gene interface set is the right one;
what resolution Foli-seq actually needs (15-way epithelial subtyping may be far finer than the assay
requires); whether mouse or human is primary, which matters because the species-match effect is real;
and whether |r| >= 0.8 is the intended redundancy threshold.

**Known residual confounds.** `GUT_HEALTHY_EPI` is a subset of `GUT_HEALTHY`, so the compartment
contrast is within-dataset. Haber is small intestine while most human references are colon-weighted.
Class counts are not matched across the grid, so each single effect estimate carries the +0.78
artifact — which is why replication, not magnitude, carries the weight.
