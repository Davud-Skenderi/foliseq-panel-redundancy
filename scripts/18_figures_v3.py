"""v3 figures: a corrected redundancy figure, accuracy as a DISTRIBUTION, and a reference-robustness panel.

Why these three:
  F4c  v2's fig4 plotted the raw NUMBER of gene pairs with |r| >= 0.8 per panel. The number of possible
       pairs grows as n(n-1)/2, so a 4,686-gene panel has ~340x more chances to produce a pair than a
       330-gene panel. The raw count therefore ranks panels largely by size. The figure's own code comment
       says it shows "% genes without a |r| >= 0.8 partner", which is the Spapros definition and is already
       computed in q3_redundancy_summary.csv as pct_genes_nonredundant. This plots that, and the pair rate
       per 1,000 possible pairs beside it.
  F5h  v2's fig5 plotted each panel's MEAN accuracy over cell types. A mean hides shape: a panel that is
       uniformly mediocre and one that is excellent on most types and fails on a few can share a mean.
       This plots the distribution of per-cell-type accuracies.
  F7   Robustness. Panels and metrics are held fixed while the reference is varied. Absolute accuracy is
       expected to shift (a reference with more and finer cell types is harder). What should survive is the
       ORDERING of panels. Spearman rank correlation between references tests exactly that.
"""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

O, F = PROJ + "results/v2/", PROJ + "figures/v3/"
os.makedirs(F, exist_ok=True)
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
BLUE, RED, GREY = "#2a78d6", "#e34948", "#9a9995"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 200,
                     "savefig.bbox": "tight", "figure.facecolor": "white"})
SHORT = {"FOLI_mPS3": "Foli-seq mPS3", "FOLI_hPS2": "Foli-seq hPS2", "COSMX_HIO": "CosMx IO add-on (100)",
         "MER_M_BRAINTRAIN": "MERSCOPE Brain (128)", "MER_M_PAN": "MERSCOPE Pan-Mouse (815)",
         "MER_M_PANNEURO": "MERSCOPE PanNeuro (815)", "COSMX_MNEURO": "CosMx Mouse Neuro (1K)",
         "COSMX_MUCC": "CosMx Mouse UCC (1K)", "COSMX_HUCC": "CosMx Human UCC (1K)",
         "GEOMX_H_CTA": "GeoMx CTA (1.8K)", "GEOMX_DOG_CCA": "GeoMx Canine (1.9K)",
         "XEN_M5K": "Xenium Mouse 5K", "XEN_H5K": "Xenium Human 5K", "COSMX_H6K": "CosMx Human 6K"}
REF_TITLE = {"MM_SLN": "Mouse spleen / lymph node (26 immune types)",
             "HS_PBMC": "Human PBMC", "GUT_HEALTHY": "Gut Cell Atlas — healthy",
             "GUT_IBD": "Inflamed gut (metaplasia/IBD)", "GUT_CRC": "Colorectal cancer"}
short = lambda p: SHORT.get(p, p)
save = lambda fig, n: (fig.savefig(F + n), plt.close(fig))
is_rand = lambda p: p.startswith("RANDOM")
colour = lambda p: RED if p.startswith("FOLI") else (GREY if is_rand(p) else BLUE)



def accuracy_table(ref):
    """Per-cell-type accuracy for a reference.

    Prefers the aggregated CSV written by 17_multiref_eval.py. If that run is still in progress,
    rebuilds the table from the per-panel cached confusion matrices so partial results are usable.
    """
    agg = O + (f"q4_accuracy_by_celltype_{ref}.csv" if ref != "MM_SLN" else "q4_accuracy_by_celltype.csv")
    aggA = pd.read_csv(agg, index_col=0) * 100 if os.path.exists(agg) else None
    cols = {}
    for f in sorted(glob.glob(O + f"q4_confusion_{ref}__*.csv")):
        panel = os.path.basename(f).split("__", 1)[1][:-4]
        cm = pd.read_csv(f, index_col=0)
        cols[panel] = pd.Series(np.diag(cm.values), index=cm.index)
    cacheA = pd.DataFrame(cols) * 100 if cols else None
    # A subset run (--panels) rewrites the aggregate with only that subset, so prefer whichever
    # source carries more panels; the per-panel caches are never overwritten.
    if aggA is not None and (cacheA is None or aggA.shape[1] >= cacheA.shape[1]):
        return aggA, "complete"
    if cacheA is None:
        return None, "missing"
    return cacheA, "partial"


# ------------------------------------------------------------------ F4 corrected
def fig4_corrected():
    R = pd.read_csv(O + "q3_redundancy_summary.csv")
    R["pct_genes_redundant"] = 100 - R.pct_genes_nonredundant
    R["pair_rate_per_1k"] = 1000 * R.n_pairs_r_ge_0_8 / (R.n_genes_measured * (R.n_genes_measured - 1) / 2)
    datasets = [d for d in ("MM_SLN", "HS_PBMC") if d in set(R.dataset)]
    # NOTE: no sharey — each column is sorted by its own reference, so a shared y-axis would
    # overwrite one column's row labels with the other's ordering.
    fig, axs = plt.subplots(2, len(datasets), figsize=(6.2 * len(datasets), 8.6), sharey=False)
    axs = np.atleast_2d(axs)
    for col, ds in enumerate(datasets):
        r0 = R[R.dataset == ds]
        for row, (metric, lab) in enumerate([
                ("pct_genes_redundant", "% of panel genes with a |r| ≥ 0.8 partner"),
                ("pair_rate_per_1k", "Redundant pairs per 1,000 possible pairs")]):
            r = r0.sort_values(metric).set_index("panel_id")
            ax = axs[row, col]
            ax.barh(range(len(r)), r[metric], color=[colour(i) for i in r.index], height=.68)
            for i, (v, g, c) in enumerate(zip(r[metric], r.n_genes_measured, r.n_pairs_r_ge_0_8)):
                ax.text(v + r[metric].max() * .015, i, f"{v:.2f}   ({int(c)} pairs, {int(g)} genes)",
                        va="center", fontsize=6.5, color=INK2)
            ax.set_yticks(range(len(r))); ax.set_yticklabels([short(i) for i in r.index], fontsize=7)
            ax.set_xlim(0, r[metric].max() * 1.45)
            ax.set_xlabel(lab, fontsize=8)
            if row == 0:
                ax.set_title(REF_TITLE.get(ds, ds), fontsize=9.5, loc="left", color=INK)
            ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
    fig.suptitle("Redundancy, normalised for panel size  (v2 fig4 plotted raw pair counts, which rank by size)",
                 fontsize=10, color=INK, x=.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, .975])
    save(fig, "v3_fig4c_redundancy_normalised.png")
    return R


# ------------------------------------------------------------------ F5 histogram
def fig5_histogram(ref):
    A, status = accuracy_table(ref)
    if A is None:
        return None
    rand_cols = [c for c in A.columns if is_rand(c)]
    panels = [c for c in A.columns if not is_rand(c)]
    order = A[panels].mean().sort_values(ascending=False).index.tolist()
    rand_all = A[rand_cols].values.ravel() if rand_cols else None
    bins = np.linspace(0, 100, 21)
    ncol = 3
    nrow = int(np.ceil(len(order) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.5 * ncol, 1.85 * nrow), sharex=True, sharey=True)
    for ax, p in zip(axs.ravel(), order):
        v = A[p].dropna().values
        if rand_all is not None:
            ax.hist(rand_all, bins=bins, weights=np.full(len(rand_all), 1 / len(rand_cols)),
                    color=GREY, alpha=.45, label="random panels")
        ax.hist(v, bins=bins, color=colour(p), alpha=.85, label=short(p))
        ax.axvline(v.mean(), color=INK, lw=1.1, ls="--")
        ax.set_title(f"{short(p)}   mean {v.mean():.0f}%, median {np.median(v):.0f}%",
                     fontsize=7.8, loc="left", color=INK)
        ax.grid(axis="y", color=GRID); ax.set_axisbelow(True)
    for ax in axs.ravel()[len(order):]:
        ax.set_visible(False)
    for ax in axs[-1]:
        ax.set_xlabel("per-cell-type accuracy (%)", fontsize=8)
    for ax in axs[:, 0]:
        ax.set_ylabel("cell types", fontsize=8)
    h = [plt.Rectangle((0, 0), 1, 1, color=GREY, alpha=.45), plt.Rectangle((0, 0), 1, 1, color=BLUE, alpha=.85),
         plt.Rectangle((0, 0), 1, 1, color=RED, alpha=.85)]
    fig.legend(h, ["random size-matched panels", "commercial panel", "Foli-seq panel"],
               loc="lower right", frameon=False, fontsize=8, ncol=3, bbox_to_anchor=(.99, -.012))
    fig.suptitle(f"Distribution of per-cell-type accuracy — {REF_TITLE.get(ref, ref)}   "
                 f"({A.shape[0]} cell types{'' if status == 'complete' else '; PARTIAL RUN'})   "
                 f"dashed line = mean (what v2 fig5 plotted)",
                 fontsize=9.5, color=INK, x=.01, ha="left")
    fig.tight_layout(rect=[0, .015, 1, .975])
    save(fig, f"v3_fig5h_accuracy_distribution_{ref}.png")
    return A


# ------------------------------------------------------------------ F7 robustness
def fig_robustness():
    frames, status = {}, {}
    for r in ["MM_SLN", "GUT_HEALTHY", "GUT_IBD", "GUT_CRC"]:
        A, st = accuracy_table(r)
        if A is not None and A.shape[1] >= 4:
            frames[r], status[r] = A, st
    if len(frames) < 2:
        return None
    refs = list(frames)
    means = pd.DataFrame({r: frames[r].mean() for r in refs})
    shared = [p for p in means.index if not is_rand(p) and means.loc[p].notna().all()]

    fig = plt.figure(figsize=(14.0, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 1.0, 1.0], wspace=.32)

    # (a) distributions shift, per reference.
    # Pooled over the SHARED panel set only: pooling different panel sets per reference would
    # confound the reference shift with a change in which panels contribute.
    ax = fig.add_subplot(gs[0, 0])
    for i, r in enumerate(refs):
        vals = frames[r][shared].values.ravel()
        vals = vals[np.isfinite(vals)]
        parts = ax.violinplot([vals], positions=[i], vert=False, showmeans=True, widths=.8)
        for b in parts["bodies"]:
            b.set_facecolor(BLUE); b.set_alpha(.45)
        for k in ("cbars", "cmins", "cmaxes", "cmeans"):
            parts[k].set_color(INK2); parts[k].set_linewidth(1)
    ax.set_yticks(range(len(refs)))
    ax.set_yticklabels([f"{REF_TITLE.get(r, r)}\n({frames[r].shape[0]} cell types, "
                        f"{len(shared)} shared panels{'' if status[r]=='complete' else ', partial run'})"
                        for r in refs], fontsize=7)
    ax.set_xlabel(f"per-cell-type accuracy (%), pooled over the {len(shared)} shared panels", fontsize=8)
    ax.set_title("a  Changing the reference shifts the whole distribution", fontsize=9, loc="left", color=INK)
    ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)

    # (b) per-panel mean across references
    ax = fig.add_subplot(gs[0, 1])
    ordr = means.loc[shared].mean(axis=1).sort_values().index
    for r in refs:
        ax.plot(means.loc[ordr, r], range(len(ordr)), "o-", ms=3.5, lw=1, label=REF_TITLE.get(r, r), alpha=.85)
    ax.set_yticks(range(len(ordr))); ax.set_yticklabels([short(p) for p in ordr], fontsize=7)
    ax.set_xlabel("mean accuracy (%)", fontsize=8)
    ax.set_title("b  Ordering is largely preserved", fontsize=9, loc="left", color=INK)
    ax.legend(fontsize=6.2, frameon=False, loc="lower right")
    ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)

    # (c) rank correlation between references
    ax = fig.add_subplot(gs[0, 2])
    n = len(refs)
    S = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            S[i, j] = spearmanr(means.loc[shared, refs[i]], means.loc[shared, refs[j]]).statistic
    im = ax.imshow(S, cmap="BuPu", vmin=0, vmax=1)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{S[i, j]:.2f}", ha="center", va="center", fontsize=7.5,
                    color="white" if S[i, j] > .6 else INK2)
    ax.set_xticks(range(n)); ax.set_xticklabels(refs, rotation=40, ha="right", fontsize=7)
    ax.set_yticks(range(n)); ax.set_yticklabels(refs, fontsize=7)
    ax.set_title("c  Spearman rank correlation of panel ordering", fontsize=9, loc="left", color=INK)
    ax.spines[:].set_visible(False)
    fig.colorbar(im, ax=ax, fraction=.04, pad=.02)
    fig.suptitle("Robustness: panels and metrics fixed, reference varied", fontsize=10.5, color=INK, x=.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, .95])
    save(fig, "v3_fig7_reference_robustness.png")
    means.to_csv(O + "q7_mean_accuracy_by_reference.csv")
    pd.DataFrame(S, index=refs, columns=refs).to_csv(O + "q7_reference_rank_correlation.csv")
    return means


if __name__ == "__main__":
    R = fig4_corrected(); print("fig4 corrected written")
    for ref in ["MM_SLN", "GUT_HEALTHY", "GUT_IBD", "GUT_CRC"]:
        if fig5_histogram(ref) is not None:
            print("histogram written:", ref)
    m = fig_robustness()
    print("robustness figure written" if m is not None else "robustness: needs >=2 references")
