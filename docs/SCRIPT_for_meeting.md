# Meeting script

Lines in **bold** are meant to be said more or less as written. Square brackets are stage directions.
Target is 10–12 minutes of talking, leaving room for him to interrupt — and let him interrupt, the
interruptions are the useful part.

**Have open, in this order:** `v3_fig5h_accuracy_distribution_MM_SLN.png`,
`v3_fig7_reference_robustness.png`, `results/v2/q11_interface_coverage.csv`,
`v3_fig4c_redundancy_normalised.png`. Everything else stays closed unless asked for.

---

## Opening — about 45 seconds

**"So the main thing I did was get the panel comparison onto gut tissue. Everything in Q3 and Q4 had
been running against immune references — mouse spleen and lymph node, and human PBMC — because those
were the well-annotated datasets we had. Which means every epithelial conclusion was resting on gene
lists alone and had never been checked against real gut data. I added gut references and re-ran it.**

**I also converted fig5 into a distribution like you asked, and while I was in there I found a
problem with fig4 that I want to show you.**

**Four things total. Can I go through them?"**

[Stop. Let him steer. If he wants a different order, take it — the pieces are independent.]

---

## 1. The histogram — about 2 minutes

[Put up `v3_fig5h_accuracy_distribution_MM_SLN.png`]

**"This is fig5 replotted. Each small panel is one gene panel. The x-axis is per-cell-type accuracy,
the y-axis is how many cell types land in each bin. Grey behind is random gene sets of the same size,
for scale. The dashed line is the mean — that's the single number the old fig5 was plotting."**

**"The reason to plot it this way: a mean can't tell you the shape. A panel that's uniformly 70%
across every cell type and a panel that's 85% on most and 50% on the rest give you the same mean, and
they mean completely different things for design."**

**"And that second case is what we actually have. Look at hPS2 — it's bimodal. One cluster of cell
types around 55%, another around 85%. The 70% was an average of two different behaviours."**

**"Which is more useful, because now it's actionable. You can go and name the cell types in that
left cluster and add markers for them. You can't do anything with a mean."**

[If he asks whether it's bimodal everywhere:]

**"No, and that's interesting. On the human gut references the left tail is much smaller — about 11
to 24 percent of cell types below 60, against 38 percent on mouse immune. So some of that bimodality
was a human panel being tested on mouse immune tissue, not a property of the panel."**

---

## 2. The gut references — about 3 minutes

**"So the references. I added four, all from the CELLxGENE Census except one."**

**"Three human gut: the Gut Cell Atlas for healthy, a metaplasia and inflammatory disease dataset,
and a colorectal cancer one. 72, 54 and 45 annotated cell types. Three rather than one because a
result from a single dataset might just be about that dataset — if it holds across healthy, inflamed
and malignant, it's about panel design."**

**"And Haber 2017, mouse small intestine, 15 cell types, pulled from GEO. That one I added for a
specific reason I'll come to if it's useful."**

[Put up `v3_fig7_reference_robustness.png`]

**"The question this answers is whether the panel comparison survives changing the reference at all.
Panel (a) is the spread of per-cell-type accuracies on each reference — you can see everything shifts
left as references get harder. Panel (b) is each panel's mean accuracy, one line per reference.
Panel (c) is the rank correlation between every pair of references."**

**"The absolute numbers move a lot. The ranking barely moves — Spearman 0.92 to 0.97 across all
four. So which panel beats which is a property of the panels, not of which dataset I picked."**

**"On normalising for panel size: every panel is scored against a random gene set of the same size,
run on the same reference. So a 5,000-gene panel competes against 5,000 random genes, not against a
300-gene panel. And the baselines confirm the references genuinely differ in difficulty — random 535
genes get 58% on mouse immune, 52% on mouse gut, 48% on the 72-type human gut."**

**"Curation is clearly doing real work, by the way. A 330-gene random set gets 43 to 47 percent on
the gut references. mPS3 gets 70 to 78."**

---

## 3. Interface coverage — about 3 minutes, this is the important one

[Put up `results/v2/q11_interface_coverage.csv`]

**"This one I think is the most relevant to what gene set Foli-seq should target, and it's also the
one that doesn't depend on any of the methodology — there's no classifier here, no normalisation, no
metric choice. It's counting genes."**

**"I built a 154-gene host-microbiome interface set — antimicrobial peptides, mucins,
pattern-recognition receptors, tight junction genes, IgA transport, the trefoil and goblet
programme, Paneth markers. And I used the CosMx whole-transcriptome list as the denominator, so
'missing' means missing relative to what the platform could have included, not relative to the
genome."**

**"Our panels spend 5 to 8 percent of their gene budget on that biology. Every commercial panel
spends 0.4 to 2.9 percent. hPS2 carries 42 interface genes in 543 slots; Xenium Human 5K carries 52
in five thousand. So we're about sevenfold denser."**

[Pause here. This is the headline — let it land.]

**"But the part I think matters more: 55 of those 154 genes are on no targeted panel at all. Thirty
of 45 antimicrobial peptides, basically the whole beta-defensin family. And 13 of 21 mucins —
including MUC3A, MUC12, MUC13 and MUC17, which are the membrane-tethered mucins of the brush border.
That's core barrier biology and nobody has it."**

**"And one thing in our own panels I wanted to flag. hPS2 has one of 34 tight junction and barrier
genes. mPS3 has four. Xenium Mouse 5K has twenty. For an assay that's reading cells shed off an
epithelial barrier, being that thin on barrier genes seems like the thing to fix in the next
revision."**

---

## 4. The fig4 problem — about 2 minutes

[Put up `v3_fig4c_redundancy_normalised.png`]

[Tone matters here. You found it while checking. It's not a criticism.]

**"Last thing. While I was going through the code I noticed fig4 is plotting the raw count of gene
pairs above the 0.8 correlation threshold. The number of possible pairs goes as n squared over two,
so a 4,686-gene panel has about 340 times more chances to produce a redundant pair than a 330-gene
one. So that figure is mostly ranking panels by size."**

**"And the comment above it in the code actually describes the right metric — percent of genes with
a redundant partner, the Spapros definition — and q3_redundancy_summary already stores it as
pct_genes_nonredundant. It just wasn't the column being plotted."**

**"Normalised, the ordering flips. CosMx 6K has the most redundant pairs by raw count, 31, and almost
the fewest per gene — 0.47 percent. And the CosMx IO 100-plex looks unremarkable at 18 pairs but has
18.2 percent of its genes carrying a redundant partner, which is the highest of anything."**

**"Which makes sense biologically — a small focused immuno-oncology panel is packed with co-regulated
immune genes, and a big discovery panel spreads thin."**

---

## 5. What I'm not sure about — about 2 minutes

[Say this section. It's the part that makes the rest credible, and these are real open questions, not
false modesty.]

**"A few things I want to flag that I'm genuinely unsure about."**

**"First — scope. You asked me to look at different references and convert a figure to a histogram.
I ended up doing a fair bit more than that, because every time I checked something I found a confound
that had to be dealt with first. The immune references were all mouse and the gut ones all human, so
'immune versus gut' couldn't be separated from 'mouse versus human'. Fixing that turned up that Haber
is epithelium-only while the human references are whole tissue. And so on. I think each step was
necessary but I'd rather you tell me if I went too far down that path."**

**"Second — the interface gene set is mine. I built it from gene families I thought were relevant.
It's defined explicitly in the script so it's auditable, but if you'd define it differently the
55-gene gap list changes. I'd rather use your version."**

**"Third, and this one I keep going back and forth on — the whole Q4 metric assumes what we care
about is telling cell types apart. That's what Spapros optimises for. But if what Foli-seq actually
needs to report is disease state, or inflammation, then cell-type discriminability might be a proxy
for the wrong thing. A panel could be great at separating goblet cells from tuft cells and useless
at detecting inflammation. I don't know which we want."**

**"And related — the gut references annotate 12 to 15 epithelial cell types. I don't know if
Foli-seq needs that resolution. If it really only needs to distinguish a few coarse states, then
I've been calibrating against a much harder problem than the real one."**

---

## 6. Questions to put to him

1. **"Which pathway database do you want for the Q2 side — GO, Reactome, KEGG, Hallmark?"**
2. **"Is the interface gene set right, or would you define it differently?"**
3. **"What does Foli-seq actually need to output — cell types, or state?"**
4. **"Mouse or human as the primary target?"** [This one matters — see reserve below.]
5. **"Is 0.8 the threshold you want for redundancy?"**

---

## Reserve — only if he asks

**If he asks "how do you know that isn't just an artifact of the reference you chose":**

**"That's why I added Haber. With a mouse gut reference the design becomes a 2x2 — mouse and human,
immune and gut — so each effect can be estimated twice with the other factor held fixed. The tissue
effect replicates, correlation 0.97 between the two estimates. A general species effect doesn't,
minus 0.29. So the species thing is specific to our two panels on gut tissue, not a general
property."**

**If he asks "is that difference real" about any number:**

**"The classifiers run at five seeds on three of the four references now, and standard errors are
1.4 to 2.6 points. So I wouldn't call anything under about three points a difference."**

**If he asks about mouse versus human transfer:**

**"mPS3 beats hPS2 by 4.2 points on mouse gut epithelium and loses by 7.2 on human — an 11-point
swing, about three combined standard errors. Each panel does better on the species it was built for.
So I don't think mPS3's mouse validation should be assumed to carry over to human samples."**

**If he asks whether the big panels are better:**

**"Per gene, no. Against the whole-transcriptome ceiling, hPS2 gets 90 percent of the achievable
range on human gut epithelium with 543 genes — 85.2 against the whole transcriptome's 87.5. And two
panels actually score below random on gut epithelium: CosMx IO and MERSCOPE Brain. Not because
they're bad, they're just outside their domain — CosMx IO gets 71.5 percent on human immune tissue
with 94 genes."**

**If he asks about the classifier itself:**

**"It's the same one fig5 and fig6 were already using — boosted decision trees, five-fold
cross-validation, confusion matrix normalised per true cell type. On Haber it sees 304 genes and the
other 11,786 measured genes are deleted from the matrix. The diagonal of the confusion matrix is the
per-cell-type accuracy, which is what the histogram is showing."**

**If he asks about bugs:**

**"Two in the existing code. 12_spapros_eval caches confusion matrices with the panel name but not
the reference, so if you ran a second reference through it you'd silently get the first one's
results back. And it takes element zero from the Spapros classifier and drops element one, which is
the standard deviation — that's why nothing in the repo had error bars."**

---

## Closing

**"That's what I have. The whole thing's on a branch with a PR open, and there's a document in
docs/ that explains every figure. What do you want me to do next?"**

[Then stop talking.]
