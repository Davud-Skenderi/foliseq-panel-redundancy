"""Figures for the v2 analysis (results/v2 -> figures/v2). Every figure shows counts or accuracies, nothing derived."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ  # noqa: E402
import os, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

O = PROJ + "results/v2/"; F = PROJ + "figures/v2/"; os.makedirs(F, exist_ok=True)
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
BLUE = "#2a78d6"; RED = "#e34948"
SEQ = LinearSegmentedColormap.from_list("blue", ["#f4f3f0", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 200, "savefig.bbox": "tight", "figure.facecolor": "white"})
SHORT = {"FOLI_mPS3": "Foli-seq mPS3", "FOLI_hPS2": "Foli-seq hPS2", "COSMX_HIO": "CosMx IO add-on (100)",
         "MER_M_BRAINTRAIN": "MERSCOPE Brain (128)", "MER_M_PAN": "MERSCOPE Pan-Mouse (815)",
         "MER_M_PANNEURO": "MERSCOPE PanNeuro (815)", "COSMX_MNEURO": "CosMx Mouse Neuro (1K)",
         "COSMX_MUCC": "CosMx Mouse UCC (1K)", "COSMX_HUCC": "CosMx Human UCC (1K)", "GEOMX_H_CTA": "GeoMx CTA (1.8K)",
         "GEOMX_DOG_CCA": "GeoMx Canine (1.9K)", "XEN_M5K": "Xenium Mouse 5K", "XEN_H5K": "Xenium Human 5K",
         "COSMX_H6K": "CosMx Human 6K", "RANDOM_535": "Random 535 genes", "ALL_GENES": "All genes"}


def save(fig, n):
    fig.savefig(F + n); plt.close(fig)


def heat(ax, data, rowlab, collab, vmax=None, fmt="{:.0f}", cmap=SEQ):
    im = ax.imshow(data, cmap=cmap, aspect="auto", vmin=0, vmax=vmax if vmax is not None else np.nanmax(data))
    hi = (vmax if vmax is not None else np.nanmax(data)) * 0.6
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            v = data[i, j]
            if np.isfinite(v):
                ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=6.5, color="white" if v > hi else INK2)
    ax.set_yticks(range(len(rowlab))); ax.set_yticklabels(rowlab, fontsize=7.5)
    ax.set_xticks(range(len(collab))); ax.set_xticklabels(collab, rotation=50, ha="right", fontsize=7.5)
    ax.spines[:].set_visible(False)
    return im


# F1 UpSet of vendors + Foli-seq (drawn directly: bars = intersection sizes, dots = which sets)
V = pd.read_csv(O + "q1_gene_vendor_matrix.csv", index_col=0)
cols = ["10x Xenium", "NanoString CosMx", "Vizgen MERSCOPE", "NanoString GeoMx", "on_FoliSeq"]
names = ["10x Xenium", "NanoString CosMx", "Vizgen MERSCOPE", "NanoString GeoMx", "Foli-seq (mPS3 + hPS2)"]
X = V[cols].astype(int)
inter = X.value_counts().reset_index(name="n").sort_values("n", ascending=False)
inter = inter[inter.n >= 60].reset_index(drop=True)
fig = plt.figure(figsize=(10.5, 5.2))
gsp = fig.add_gridspec(2, 2, height_ratios=[2.2, 1.3], width_ratios=[1.0, 5], hspace=0.05, wspace=0.42)
axb = fig.add_subplot(gsp[0, 1]); axd = fig.add_subplot(gsp[1, 1], sharex=axb); axs_ = fig.add_subplot(gsp[1, 0], sharey=axd)
xs = np.arange(len(inter))
has_foli = inter.on_FoliSeq.values == 1
axb.bar(xs, inter.n, color=[RED if f else BLUE for f in has_foli], width=.66)
for x_, n_ in zip(xs, inter.n):
    axb.text(x_, n_ + 40, f"{n_:,}", ha="center", fontsize=6.8, color=INK2)
axb.set_ylabel("genes in exactly this combination"); axb.tick_params(axis="x", bottom=False, labelbottom=False)
axb.grid(axis="y", color=GRID); axb.set_axisbelow(True)
for k, c in enumerate(cols):
    yy = len(cols) - 1 - k
    on_ = inter[c].values == 1
    axd.scatter(xs[~on_], [yy] * (~on_).sum(), color=GRID, s=40, zorder=2)
    axd.scatter(xs[on_], [yy] * on_.sum(), color=INK, s=40, zorder=3)
for x_ in xs:
    ys = [len(cols) - 1 - k for k, c in enumerate(cols) if inter.loc[x_, c] == 1]
    if len(ys) > 1:
        axd.plot([x_, x_], [min(ys), max(ys)], color=INK, lw=1.5, zorder=2)
axd.set_yticks(range(len(cols))); axd.set_yticklabels([names[len(cols) - 1 - y] for y in range(len(cols))], fontsize=7.5)
axd.tick_params(axis="y", length=0); axd.set_xticks([]); axd.spines[:].set_visible(False)
tot = X.sum()
axs_.barh([len(cols) - 1 - k for k in range(len(cols))], tot.values, color=MUTED, height=.5)
axs_.invert_xaxis(); axs_.tick_params(axis="y", left=False, labelleft=False)
axs_.set_xlabel("genes per set", fontsize=7.5); axs_.spines["left"].set_visible(False)
save(fig, "v2_fig1_upset_vendors_foliseq.png")

# F2 epithelial signature coverage (count of signature genes on panel)
EP = pd.read_csv(O + "q2a_epithelial_signature_coverage.csv")
pan = ["FOLI_mPS3", "FOLI_hPS2", "COSMX_HIO", "MER_M_PAN", "COSMX_MUCC", "COSMX_HUCC", "GEOMX_H_CTA", "XEN_M5K", "XEN_H5K", "COSMX_H6K"]
lab = [f"{c} — {s.split()[0]} (n={n})" for c, s, n in zip(EP.cell_type, EP.source, EP.n_signature_genes)]
data = EP[[f"{p}_n" for p in pan]].values.astype(float)
fig, ax = plt.subplots(figsize=(8.5, 7.5))
im = heat(ax, data, lab, [SHORT[p] for p in pan])
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01).set_label("signature genes on the panel")
save(fig, "v2_fig2_epithelial_signature_coverage.png")

# F3 immune marker coverage (of top-20 markers)
IM = pd.read_csv(O + "q2b_immune_marker_coverage.csv")
data = IM[[f"{p}_n" for p in pan]].values.astype(float)
fig, ax = plt.subplots(figsize=(8.5, 7.5))
im = heat(ax, data, [f"{c} (n={n})" for c, n in zip(IM.cell_type, IM.n_markers)], [SHORT[p] for p in pan], vmax=20)
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01).set_label("of the cell type's top-20 markers, number on the panel")
save(fig, "v2_fig3_immune_marker_coverage.png")

# F5 summary + F6 per-cell-type heatmap: cell-type identification accuracy (Spapros XGBoost, 5-fold CV)
A = pd.read_csv(O + "q4_accuracy_by_celltype.csv", index_col=0)
S = pd.read_csv(O + "q4_celltype_identification.csv")
S["base"] = S.panel_id.str.replace(r"_\d$", "", regex=True)
Sm = S.groupby("base").agg(acc=("mean_celltype_accuracy", "mean"), n=("n_genes_measured", "mean"), pct=("pct_celltypes_acc_gt_0_8", "mean"))
A["RANDOM_535"] = A[[c for c in A.columns if c.startswith("RANDOM")]].mean(axis=1)
order_cols = list(Sm.sort_values("acc").index)
A = A[order_cols]
A = A.loc[A.mean(axis=1).sort_values().index]
fig, ax = plt.subplots(figsize=(10, 8))
im = heat(ax, A.values * 100, list(A.index), [f"{SHORT.get(c, c)}\n{int(Sm.loc[c, 'n'])} genes" for c in A.columns], vmax=100)
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01).set_label("cross-validated accuracy for this cell type (%)")
save(fig, "v2_fig6_celltype_identification.png")

fig, ax = plt.subplots(figsize=(6.8, 4.6))
Ss = Sm.sort_values("acc")
cols_ = [RED if i.startswith("FOLI") else (MUTED if i in ("RANDOM_535", "ALL_GENES") else BLUE) for i in Ss.index]
ax.barh(range(len(Ss)), Ss.acc * 100, color=cols_, height=.68)
for i, (a_, p_) in enumerate(zip(Ss.acc, Ss.pct)):
    ax.text(a_ * 100 + 0.8, i, f"{a_*100:.0f}%  ({p_:.0f}% of types > 80%)", va="center", fontsize=7, color=INK2)
ax.set_yticks(range(len(Ss))); ax.set_yticklabels([f"{SHORT.get(i, i)} — {int(n)} genes" for i, n in zip(Ss.index, Ss.n)], fontsize=7.5)
ax.set_xlim(0, 118); ax.set_xlabel("Mean accuracy identifying 26 mouse immune cell types (%)")
ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
save(fig, "v2_fig5_identification_summary.png")

# F4 redundancy (% genes without a |r| >= 0.8 partner; Spapros definition)
R = pd.read_csv(O + "q3_redundancy_summary.csv")
fig, axs = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True)
order = R[R.dataset == "MM_SLN"].sort_values("n_pairs_r_ge_0_8").panel_id.tolist()
for ax, ds, t in zip(axs, ["MM_SLN", "HS_PBMC"], ["Mouse spleen/lymph node", "Human PBMC"]):
    r = R[R.dataset == ds].set_index("panel_id").reindex(order)
    ax.barh(range(len(r)), r.n_pairs_r_ge_0_8, color=[RED if i.startswith("FOLI") else BLUE for i in r.index], height=.68)
    for i, (n_, g_) in enumerate(zip(r.n_pairs_r_ge_0_8, r.n_genes_measured)):
        ax.text(n_ + 0.3, i, f"{int(n_)}  (of {int(g_)} genes)", va="center", fontsize=7, color=INK2)
    ax.set_title(t, fontsize=9, loc="left", color=INK); ax.set_xlabel("Gene pairs with |Pearson r| ≥ 0.8")
    ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
axs[0].set_yticks(range(len(order))); axs[0].set_yticklabels([SHORT[p] for p in order], fontsize=7.5)
fig.tight_layout(); save(fig, "v2_fig4_redundant_pairs.png")
print(sorted(os.listdir(F)))
