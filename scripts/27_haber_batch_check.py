"""Is the Haber mouse-gut classifier reading cell type, or reading which mouse the cell came from?

WHY ASK. Haber's cells come from several mice. The cell-type labels are what we want the classifier
to learn. If a given cell type happens to come mostly from one mouse, then anything that distinguishes
mice -- strain, batch, capture day, sequencing run -- is also a cue for that cell type. The classifier
will happily use it, accuracy goes up, and the number no longer means what we think it means.

This does not require re-running classifiers. It is a property of the experimental design, measurable
from the labels alone:

  1. Cross-tabulate cell type against mouse of origin.
  2. Measure how far that table is from independence (Cramer's V, and normalised mutual information).
  3. Report, per cell type, the share of its cells coming from its single most common mouse.

If every cell type is drawn evenly from every mouse, batch cannot be a cue and the accuracies stand.
If some cell types sit almost entirely in one mouse, those types' accuracies are suspect and should
be read with that in mind.
"""
import os, sys
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ, EXT
from scipy.stats import chi2_contingency
from sklearn.metrics import normalized_mutual_info_score

def main():
    a = ad.read_h5ad(EXT + "expr/haber_gut.h5ad", backed="r")
    obs = a.obs[["celltype", "batch"]].astype(str)
    ct = pd.crosstab(obs.celltype, obs.batch)
    print(f"Haber: {len(obs)} cells, {ct.shape[0]} cell types, {ct.shape[1]} mice\n")
    print("cells per mouse:", dict(obs.batch.value_counts()), "\n")

    chi2, p, dof, _ = chi2_contingency(ct)
    n = ct.values.sum()
    v = np.sqrt(chi2 / (n * (min(ct.shape) - 1)))
    nmi = normalized_mutual_info_score(obs.celltype, obs.batch)
    print(f"cell type vs mouse:  Cramer's V = {v:.3f}   normalised MI = {nmi:.3f}   p = {p:.2e}")
    print("  (V = 0 means cell type and mouse are independent; 1 means one determines the other)\n")

    frac = ct.div(ct.sum(axis=1), axis=0)
    top = frac.max(axis=1).sort_values(ascending=False)
    even = 1.0 / ct.shape[1]
    print(f"share of each cell type's cells coming from its single most common mouse")
    print(f"  (even split across {ct.shape[1]} mice would be {even:.2f})\n")
    for t, f in top.items():
        flag = "  <-- concentrated" if f > 0.6 else ""
        print(f"   {t:34s} {f:.2f}   n={int(ct.loc[t].sum()):5d}{flag}")

    out = pd.DataFrame({"celltype": top.index, "max_mouse_share": top.values,
                        "n_cells": ct.sum(axis=1).reindex(top.index).values,
                        "n_mice_present": (ct > 0).sum(axis=1).reindex(top.index).values})
    out.to_csv(PROJ + "results/v2/q13_haber_batch_structure.csv", index=False)
    conc = out[out.max_mouse_share > 0.6]
    print(f"\n{len(conc)} of {len(out)} cell types draw >60% of their cells from one mouse")
    if len(conc):
        print("   affected:", ", ".join(conc.celltype))
    print(f"\nmedian share from top mouse: {out.max_mouse_share.median():.2f}")

if __name__ == "__main__":
    main()
