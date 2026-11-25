# Panel comparison — full brief

Written to be read before presenting. Sections 1–3 are the concepts; section 4 walks every figure
with what it shows and how to say it out loud; section 5 has the questions to ask.

---

## 0. The thirty-second version

> We'd only ever tested the panels against immune tissue, so nothing epithelial had been checked
> against real gut data. I added four gut references, converted the accuracy figure from a bar of
> means to a distribution, and normalised the panel comparison for panel size. The ranking of panels
> holds up across every reference. And counting genes directly, our panels spend 5–8% of their gene
> budget on host-microbiome interface biology where commercial panels spend 0.4–2.9% — with 55
> interface genes, including the brush-border mucins, on no targeted panel at all.

---

## 1. What a "reference" is, and what each one is for

A **reference** is an annotated single-cell dataset: someone dissociated a tissue, sequenced each
cell, and labelled every cell with its type. Those labels are the answer key.

We use a reference to **test a panel**. The procedure: throw away every gene except the panel's,
hand what's left to a classifier, and see whether it can recover each cell's known label. If it can,
the panel carries the information that distinguishes those cells. If it can't, the panel is blind to
that variation — and nothing downstream will recover what was never measured.

Two consequences worth holding onto:

- **This is not proposing Foli-seq do cell typing.** The classifier is a measuring instrument for
  information content. We already know the labels; we're asking whether they're *recoverable*.
- **A score is never a property of the panel alone.** "77% accuracy" means "77% at separating the 26
  cell types annotated in this one dataset." Change the dataset and the number changes.

That second point is why there are several references rather than one. Each has a distinct job:

| Reference | What it is | Cell types | Why it's in the study |
| --- | --- | --- | --- |
| `MM_SLN` | mouse spleen / lymph node, CITE-seq | 26 | The original reference. Now serves as the **immune contrast** — the comparison that lets us say panels are immune-weighted. Not a design target; no gene should be chosen because it works on spleen. |
| `HUMAN_IMMUNE` | Human Immune Health Atlas, blood | 37 | Fills the human-immune slot so the species comparison can be checked on immune tissue as well as gut. |
| `HABER_GUT` | Haber 2017, mouse small intestine (GSE92332) | 15 | The only well-annotated **mouse gut** data. Lets us compare immune vs gut *within mouse*, so species can't explain the difference. |
| `GUT_HEALTHY_EPI` | Gut Cell Atlas, epithelial cells only | 12 | Haber is epithelium-only, so a human epithelium-only reference is needed to compare the two **like for like**. |
| `GUT_HEALTHY` | Gut Cell Atlas, whole tissue | 72 | Healthy human gut, all compartments: epithelial, immune, stromal. |
| `GUT_IBD` | metaplasia in inflammatory gut disease | 54 | Inflamed human gut. Tests whether findings survive disease. |
| `GUT_CRC` | colorectal cancer immune hubs | 45 | Malignant human gut. Same purpose. |

**Why so many gut ones?** A result from one dataset might be about that dataset. A result that holds
across healthy, inflamed and malignant gut is about panel design.

**Why the epithelium-only ones matter.** The exfoliome — what Foli-seq reads — is epithelium. The gut
wall also contains immune and stromal cells, which is why whole-tissue references have 45–72 cell
types and the epithelial ones have 12–15. Comparing an epithelium-only reference against a
whole-tissue one mixes up "different species" with "different compartment," so both versions exist.

All references except Haber come from the CELLxGENE Census, release 2025-11-08, primary data only so
no cell is counted twice. Each is subsampled to 250 cells per cell type with types under 40 cells
dropped, so common types don't swamp rare ones.

---

## 2. The two metrics

**Cell-type identification (Q4).** The classifier procedure above. Gradient-boosted trees, 5-fold
cross-validation. The confusion matrix is normalised per true cell type so a rare type counts as much
as a common one, and its diagonal gives **one accuracy per cell type** — 15 numbers on Haber, 72 on
the whole-tissue Gut Cell Atlas.

**Redundancy (Q3).** Correlate every pair of genes within a panel. If two genes always rise and fall
together, the second buys little. A pair counts as redundant at |r| ≥ 0.8.


---

## 2b. The classifier, in detail — and what "removing the other genes" means

### Which genes get removed

The universe is **every gene measured in that reference dataset**, after a light filter (a gene must
appear in at least 10 cells). Not a fixed list, and not "the transcriptome" in the abstract — whatever
that particular experiment detected.

| Reference | genes measured | mPS3 genes present | hidden from the classifier |
| --- | --- | --- | --- |
| `HABER_GUT` | 12,090 | 304 | **11,786** |
| `GUT_HEALTHY` / `GUT_HEALTHY_EPI` | 20,620 | 492 | **20,128** |
| `HUMAN_IMMUNE` | 15,932 | 333 | **15,599** |

On Haber the dataset measures 12,090 genes, mPS3 holds 304 of them, and the other 11,786 columns are
deleted from the matrix. The classifier cannot see them, cannot use them, and does not know they
exist.

Two details that explain numbers elsewhere in the repo:

- Haber's raw file has 13,487 gene columns but only 12,090 survive. That is the mouse-to-human symbol
  mapping: about a thousand mouse genes have no human counterpart and are dropped, and several
  hundred collapse together (19 mouse alpha-defensins become one human DEFA1).
- mPS3 is a 535-gene panel but only 304 appear in Haber; the rest were never detected in that
  experiment. This is why "panel size" and "genes actually used" differ, and why random baselines are
  matched on **measured** size rather than nominal size.

### What the classifier is

A classifier learns to predict a label from measurements. Here the measurements are expression levels
of the panel's genes in one cell; the label is that cell's annotated type.

**XGBoost, gradient-boosted decision trees.** A decision tree asks a chain of yes/no questions about
the numbers — "is LGR5 above 2.1? if so, is MUC2 below 0.4?" — until it reaches a guess. One tree is
weak. Boosting builds many in sequence, each trained on the mistakes the previous ones made, and the
final prediction is their weighted vote. Trees suit single-cell data because expression is sparse,
zero-inflated and non-linear, and trees need no scaling or distributional assumptions.

**Five-fold cross-validation.** The cells are split into five groups; the model trains on four and
predicts the held-out one, rotating until every cell has been predicted by a model that never saw it.
Without this the score would measure memorisation.

**The confusion matrix.** Rows are the true cell type, columns the predicted one, each row normalised
to sum to 1 — so a row reads "of the cells that really were goblet cells, what fraction did the model
call goblet, what fraction enterocyte, and so on." Normalising per row is deliberate: a type with 60
cells counts as much as one with 3,000.

**The diagonal is what gets reported** — the fraction of each cell type correctly identified. One
number per cell type: 15 on Haber, 72 on the whole-tissue Gut Cell Atlas. The old fig5 averaged them
into one number per panel; the histogram shows all of them.

**Five seeds.** The whole procedure is repeated with five random initialisations; the average is the
result and the spread is the error bar. That spread is 1.4–2.6 points, which is why differences under
about 3 points should not be called differences.

**What it is not.** It is not a proposal for how Foli-seq data would be analysed. It is a probe for
whether information is present. The labels are already known — they are the dataset's annotations.
The question is whether a panel's genes are sufficient to *recover* them.

### Important: this is not new work

`scripts/12_spapros_eval.py` predates this round, and **`v2_fig5_identification_summary.png` and
`v2_fig6_celltype_identification.png` are both its outputs**. The classifier was already in what was
presented before — fig5 as the bar chart of mean accuracy, fig6 as the per-cell-type heatmap. Present
it as the thing those figures were already doing, now explained and run against gut tissue rather
than immune tissue only.

---

## 2c. What changed since the last update

**Already in place:** Q1–Q6, six v2 figures, the classifier, the redundancy analysis, the MSigDB
Hallmark pathway enrichment, the UpSet overlap plots, and both Foli-seq candidate panels.

**New in this round:**

- **Six references.** Three whole-tissue human gut (healthy, inflamed, malignant), Haber mouse gut, a
  human epithelium-only subset, and a human immune atlas. Previously everything ran on immune tissue.
- **Two panels** added to the matrix: the CosMx whole-transcriptome human and mouse lists.
- **fig5 replotted as a distribution** (`v3_fig5h`) rather than a bar of means.
- **fig4 corrected** (`v3_fig4c`): it had plotted raw redundant-pair counts, which rank by panel size.
- **Five new analyses:** reference robustness, species versus tissue, the cell-type-count control, the
  whole-transcriptome ceiling, and host-microbiome interface coverage.
- **Error bars**, for the first time, from 5-seed runs.
- **Two bugs found in pre-existing code** (see section 6).

---

## 3. Why accuracies can't be compared across references, and the two fixes

Each reference is a different difficulty of problem. Guessing among 72 cell types is harder than
among 12 — chance alone is 1/72 versus 1/12. So a panel scoring 70% on one reference and 85% on
another tells you nothing on its own.

**Fix one: lift.** Score every panel against a **random gene set of the same size, run on the same
reference**:

> **lift = panel accuracy − random-535 accuracy, on the same reference**

That subtracts the problem's difficulty from both sides, and it is also how panel size is normalised
— a 5,000-gene panel is compared against 5,000 random genes, not against a 300-gene panel.

The baselines show the difficulties really do differ: random 535 genes score **57.9%** on mouse
immune, **52.2%** on mouse gut epithelium, **63.2%** on human gut epithelium, **48.2%** on the
72-type human gut. Note these don't simply track cell-type count, which is why an empirical baseline
is needed instead of a formula.

Lift is imperfect and we measured how: across references it correlates **+0.78** with the number of
cell types, because a harder problem leaves more headroom for curation to fill. Figure 9 (section 4)
is the control for that.

**Fix two: ceiling fraction**, where a ceiling exists. The CosMx Whole Transcriptome panel measures
essentially everything the instrument can, so its accuracy bounds what *any* panel on that platform
could reach. With a floor and a ceiling:

> **ceiling fraction = (panel − random) / (whole transcriptome − random)**

0 means no better than random genes; 100 means as good as measuring everything. Unlike lift, this has
a meaningful maximum.

---

## 4. Every figure

### `v3_fig5h_accuracy_distribution_MM_SLN.png` — the histogram you asked for

**What you're looking at.** One small panel per gene panel. X-axis is per-cell-type accuracy, 0–100%.
Y-axis is how many cell types fall in each bin. Grey = random size-matched panels for comparison;
blue = commercial panel; red = Foli-seq. The dashed vertical line is the mean — the single number the
old fig5 plotted.

**What it shows.** A mean hides shape. **Foli-seq hPS2 is bimodal here**: one group of cell types near
55%, another near 85%. Its "70%" was an average of two different behaviours.

**How to say it.** *"The bar chart gave one number per panel. Plotting the distribution shows hPS2
isn't uniformly 70% — it's good on about half the cell types and poor on the other half. That's
actionable in a way a mean isn't, because you can name the failing types and add markers for them."*

**Companion figures** `..._GUT_HEALTHY.png`, `..._GUT_IBD.png`, `..._GUT_CRC.png` are the same plot
on the three whole-tissue gut references. Worth mentioning that hPS2's bimodality is much weaker
there — 11–24% of cell types below 60%, versus 38% on mouse immune — so the bimodality was partly
about testing a human panel on mouse immune tissue.

---

### `v3_fig7_reference_robustness.png` — does the comparison survive changing the reference?

**Three panels.**

- **(a)** One violin per reference, showing the spread of per-cell-type accuracies pooled over the
  shared panels. Reading down, the distributions **shift left** as references get harder.
- **(b)** One line per reference. Panels along the y-axis ordered by mean accuracy; each line shows
  where that reference puts them. Lines running roughly parallel = ordering preserved.
- **(c)** Spearman rank correlation of per-panel mean accuracy between each pair of references:

|  | MM_SLN | GUT_HEALTHY | GUT_IBD | GUT_CRC |
| --- | --- | --- | --- | --- |
| MM_SLN | — | 0.93 | 0.93 | 0.95 |
| GUT_HEALTHY | 0.93 | — | 0.95 | 0.92 |
| GUT_IBD | 0.93 | 0.95 | — | 0.97 |
| GUT_CRC | 0.95 | 0.92 | 0.97 | — |

**What it shows.** Absolute scores move a lot between references; the **ranking barely moves**, 0.92
to 0.97 everywhere.

**How to say it.** *"Changing the reference shifts every panel's score, because a reference with more
cell types is a harder problem. But which panel beats which stays the same across immune, healthy
gut, inflamed gut and cancer. So the comparison is about the panels, not about which dataset I
happened to pick."*

**Likely question — "then why bother with four references?"** Because that was the thing to check,
and it could have come out the other way. If the ordering had scrambled, every panel conclusion would
only have been valid for one dataset.

---

### `v3_fig4c_redundancy_normalised.png` — the correction to fig4

**What you're looking at.** Two columns, one per reference (mouse spleen/lymph node, human PBMC). Two
rows: the top plots **% of a panel's genes that have a redundant partner**, the bottom plots
**redundant pairs per 1,000 possible pairs**. Each bar is annotated with the raw pair count and the
panel's gene count, so the old and new quantities are both visible.

**What it shows.** The old fig4 plotted the **raw count** of redundant pairs. Possible pairs grow as
n(n−1)/2, so a 4,686-gene panel has ~340× more chances to produce one than a 330-gene panel — the
count ranks panels mostly by size. Normalised, **the ordering inverts**: CosMx Human 6K has the most
redundant pairs by raw count (31) and nearly the fewest per gene (0.47%), while the CosMx IO
100-plex has 18 pairs but **18.2%** of its genes carrying a redundant partner, the highest of any
panel.

**How to say it.** *"The comment above that figure describes the Spapros definition — percent of genes
with a redundant partner — and `q3_redundancy_summary.csv` already stores it. The figure was plotting
raw pair counts instead, which scale with panel size. Normalised, the conclusion flips: the big
discovery panels are the least redundant per gene and the small focused panels the most. Which makes
biological sense — a 100-gene immuno-oncology panel is packed with co-regulated immune genes."*

**Tone note.** Frame it as something you found while checking, not as a criticism. The right metric
was already computed; only the plotted column was wrong.

---

### `v3_fig10_ceiling_fraction.png` — how close a panel gets to measuring everything

**What you're looking at.** One horizontal bar pair per panel. Blue = mouse gut epithelium, green =
human gut epithelium. X-axis is **% of the range from a random 535-gene set to the whole
transcriptome**. The dashed line at 100 is the whole transcriptome. Foli-seq labels are in red.
Negative bars mean the panel does worse than random genes.

**What it shows.** On human gut epithelium the floor is 63.2% and the ceiling 87.5%, a 24.3-point
usable range. **Foli-seq hPS2 reaches 90.3% of it with 543 genes** — 85.2% accuracy against the whole
transcriptome's 87.5%. On mouse gut epithelium (floor 52.2%, ceiling 72.5%) the 5K/6K panels reach
83–84% and mPS3 reaches 67.4%.

**How to say it.** *"Adding the whole-transcriptome list gives an upper bound — nothing on that
instrument beats measuring everything. On human gut epithelium hPS2 gets 90% of the way to that
ceiling with 543 genes. Measuring the other 13,000 buys 2.3 accuracy points."*

---

### `v3_fig9_celltype_count_control.png` — the control for cell-type count

**What you're looking at.** One bar pair per panel. Grey = tissue effect as originally measured (26
immune cell types vs Haber's 15). Blue = the same effect with cell-type count **matched** at 15 on
both sides. Error bars are the spread across 3 independent draws. Negative means the panel loses
power moving from immune tissue to gut epithelium.

**Why it exists.** Lift correlates +0.78 with cell-type count. The tissue comparison *reduces* cell
types from 26 to 15, so the effect pointed in the same direction the artifact would. To separate
them, `MM_SLN` was re-run restricted to a random 15 of its 26 types, matching Haber exactly.

**What it shows.** The direction survives for all nine panels, but about a third of the magnitude was
artifact. Panels built for other tissues still lose a lot — CosMx IO −21.2, MERSCOPE Brain −12.3 —
while **the Foli-seq panels fall to near zero**: mPS3 −2.5 and hPS2 −0.3, whose interval includes
zero.

**How to say it.** *"The honest version isn't 'all panels are immune-weighted.' It's that panels
curated for other tissues lose power on gut epithelium and ours don't."*

---

### `v3_fig8_species_vs_tissue.png` — separating species from tissue

**Two panels.** **(a)** One line per panel, walking left to right through references so that **each
step changes exactly one thing** — tissue, then species, then compartment, marked by the dotted
vertical lines. Y-axis is lift. Foli-seq mPS3 is the thick red circle line, hPS2 the thick orange
square line. **(b)** Per panel, the tissue effect and the species effect as separate bars.

**Why it exists.** Originally every immune reference was mouse and every gut reference human, so
"immune vs gut" and "mouse vs human" moved together and neither could be measured. Haber and the
epithelium-matched human reference break that tie.

**What it shows.** The crossing between red and orange at the species step: **each Foli-seq panel does
better on its own species.** mPS3 beats hPS2 by +4.2 ± 2.2 on mouse gut epithelium and loses by
−7.2 ± 2.6 on human — an 11.4-point swing.

**How to say it.** *"Our two candidate panels swap places when you change species, by about 11 points.
So mPS3's mouse validation shouldn't be assumed to carry over to human samples."*

---

### `v3_fig11_complete_grid.png` — the completed design and the replication check

**Two panels.** **(a)** The four cells of the grid along the x-axis — mouse immune (26 types), human
immune (37), mouse gut epithelium (15), human gut epithelium (12) — with lift on the y-axis and error
bars where 5 seeds were run. Zero means no better than random. **(b)** A scatter: each effect
estimated **twice**, once per row or column, plotted against itself. Points near the dotted diagonal
replicate.

**What it shows.** The **tissue effect replicates** — r = +0.97 between the mouse-only and human-only
estimates, sign agreeing for 7 of 9 panels. A **general species effect does not** — r = −0.29.

**How to say it.** *"With all four cells filled, each effect can be estimated twice with the other
factor held fixed. The tissue effect replicates at 0.97. A species effect across panels in general
does not — so the species result is specific to our two panels on gut tissue, not a general
property."*

**Only raise this if asked** how we know a result isn't an artifact of the chosen reference.

---

### Interface coverage — no figure, a table

`results/v2/q11_interface_coverage.csv`. **This is the most on-target result and the hardest to argue
with, because it involves no classifier, no normalisation and no metric choice. It is counting
genes.**

A 154-gene host-microbiome interface set — antimicrobial peptides, mucins, pattern-recognition
receptors, tight junctions, IgA transport, trefoil/goblet, Paneth markers — denominated by what the
CosMx Whole Transcriptome panel can see.

| Panel | size | interface genes | **% of its own gene budget** |
| --- | --- | --- | --- |
| Xenium Mouse 5K | 4,949 | 57 | 1.15 |
| CosMx Human 6K | 6,338 | 57 | 0.90 |
| Xenium Human 5K | 5,000 | 52 | 1.04 |
| **Foli-seq hPS2** | **543** | **42** | **7.73** |
| GeoMx CTA | 1,890 | 34 | 1.80 |
| **Foli-seq mPS3** | **535** | **29** | **5.42** |
| MERSCOPE PanNeuro | 812 | 3 | 0.37 |

**55 of the 154 genes are on no targeted panel at all** — 30 of 45 antimicrobial peptides (nearly the
whole beta-defensin family) and 13 of 21 mucins, including **MUC3A, MUC12, MUC13, MUC17**, the
membrane-tethered mucins of the intestinal brush border.

**And a gap in our own panels:** hPS2 carries **1 of 34** tight junction / barrier genes. mPS3 carries
4. Xenium Mouse 5K carries 20.

**How to say it.** *"hPS2 carries 42 interface genes in 543 slots; Xenium Human 5K carries 52 in
5,000. Our panels are about sevenfold denser in this biology. But 55 interface genes are on no
targeted panel anywhere — including the brush-border mucins — and hPS2 has only one of 34 tight
junction genes, which for an assay reading shed barrier cells is probably the thing to fix next."*

---

## 5. Questions to ask

1. **Which pathway database** for the Q2-style analysis — GO, Reactome, KEGG, MSigDB Hallmark?
2. **Is the 154-gene interface set right?** It's my construction, defined explicitly in
   `scripts/25_interface_coverage.py`. Changing it changes the 55-gene gap list.
3. **What resolution does Foli-seq need?** Fifteen-way epithelial subtyping may be finer than the
   assay requires, which would mean the identification metric is calibrated too strictly.
4. **Mouse or human as primary?** The species result means this choice has consequences.
5. **Is |r| ≥ 0.8 the intended redundancy threshold?** It's the Spapros convention, and only pairs
   above it are stored, so the distribution below can't be recovered without recomputing.

---

## 6. If he asks why there's more here than he asked for

The honest answer is also the good one: **checking the original claim kept revealing confounds, and
each had to be resolved before the claim could be stated.** The immune references were all mouse while
the gut references were all human, so "immune vs gut" couldn't be separated from "mouse vs human."
Fixing that revealed that Haber is epithelium-only while the human references are whole tissue.
Fixing *that* revealed that the lift metric doesn't fully normalise for cell-type count. Scope grew
because verification grew.

Two bugs in existing code also turned up and are worth mentioning: `12_spapros_eval.py` caches
confusion matrices without the reference in the key, so running a second reference through it would
silently reload the first one's results; and it discards the classifier's standard-deviation output,
which is why no number in the repo had an error bar until now.

---

## 7. Numbers worth having memorised

- Panel ordering across references: **Spearman 0.92–0.97**
- Random 330-gene baseline on gut: **43–47%**; Foli-seq mPS3: **70–78%**
- Standard errors on 5-seed references: **1.4–2.6 points** — treat anything under ~3 points as noise
- hPS2 reaches **90.3%** of the whole-transcriptome ceiling with **543** genes
- Interface budget share: Foli-seq **5.4–7.7%**, commercial **0.4–2.9%**
- **55 of 154** interface genes on no targeted panel; hPS2 has **1 of 34** tight junction genes
- Foli-seq species swing: **11.4 points**, roughly 3 combined standard errors
