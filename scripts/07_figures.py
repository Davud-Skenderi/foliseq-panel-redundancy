"""Step 7: figures (PNG, light theme) + results workbook."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import numpy as np, pandas as pd, matplotlib, json, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list

D = PROJ + "data/"; O = PROJ + "results/"; F = PROJ + "figures/"
os.makedirs(F, exist_ok=True)
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#fcfcfb"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQ = LinearSegmentedColormap.from_list("blue", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#f2a3a3", "#c23434"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.dpi": 200, "savefig.bbox": "tight"})
meta = pd.read_csv(D + "panel_meta.csv")
PLAT_COL = {"10x Xenium": CAT[0], "NanoString CosMx": CAT[1], "Vizgen MERSCOPE": CAT[2], "NanoString GeoMx": CAT[6],
            "Illumina AmpliSeq": CAT[3], "IDT xGen": MUTED, "Foli-seq (Wang Lab)": CAT[7]}
SHORT = {"XEN_H5K": "Xenium 5K (Hs)", "XEN_M5K": "Xenium 5K (Mm)", "COSMX_H6K": "CosMx 6K (Hs)", "COSMX_HUCC": "CosMx UCC 1K (Hs)",
         "COSMX_MUCC": "CosMx UCC 1K (Mm)", "COSMX_HIO": "CosMx IO add-on (Hs)", "COSMX_MNEURO": "CosMx Neuro 1K (Mm)",
         "COSMX_MNEURO_PROT": "CosMx Neuro protein (Mm)", "MER_M_PAN": "MERSCOPE Pan-Mouse 815", "MER_M_PANNEURO": "MERSCOPE PanNeuro 815",
         "MER_M_BRAINTRAIN": "MERSCOPE Brain 128", "GEOMX_H_WTA": "GeoMx WTA (Hs)", "GEOMX_M_WTA": "GeoMx WTA (Mm)",
         "GEOMX_H_CTA": "GeoMx CTA (Hs)", "GEOMX_DOG_CCA": "GeoMx Canine CA", "AMPLI_H_TX": "AmpliSeq Transcriptome",
         "AMPLI_COMP_V3": "AmpliSeq Comp. v3 (DNA)", "IDT_XGEN_PANCANCER": "xGen Pan-Cancer (DNA)", "FOLI_mPS3": "Foli-seq mPS3 (Mm)", "FOLI_hPS2": "Foli-seq hPS2 (Hs)"}
plat = dict(zip(meta.panel_id, meta.platform))
TARGETED_RNA = [p for p in SHORT if p not in ("GEOMX_H_WTA", "GEOMX_M_WTA", "AMPLI_H_TX", "AMPLI_COMP_V3", "IDT_XGEN_PANCANCER", "COSMX_MNEURO_PROT")]


def save(fig, name):
    fig.savefig(F + name); plt.close(fig)


# ---------- Fig 1: panel size vs pathway breadth
s = pd.read_csv(O + "panel_scope_summary.csv")
fig, ax = plt.subplots(figsize=(7.2, 4.2))
for _, r in s.iterrows():
    ax.scatter(r.n_genes_in_universe, r.gut_lens_scorable_frac * 100, s=46, color=PLAT_COL[r.platform],
               edgecolor="white", linewidth=1.5, zorder=3, marker="D" if r.modality == "DNA" else "o")
lab = {"FOLI_mPS3": (-8, -9), "FOLI_hPS2": (-8, 3), "COSMX_HIO": (6, 2), "MER_M_PAN": (-8, -6), "COSMX_HUCC": (-8, 6), "GEOMX_H_CTA": (6, -6),
       "XEN_H5K": (-8, -3), "COSMX_H6K": (-8, 2), "GEOMX_H_WTA": (-8, -12), "MER_M_PANNEURO": (6, -4), "COSMX_MNEURO": (6, 2),
       "AMPLI_COMP_V3": (6, 2), "IDT_XGEN_PANCANCER": (6, -8)}
for p, (dx, dy) in lab.items():
    r = s[s.panel_id == p].iloc[0]
    ax.annotate(SHORT[p], (r.n_genes_in_universe, r.gut_lens_scorable_frac * 100), xytext=(dx, dy),
                textcoords="offset points", fontsize=7.5, color=INK2, ha="right" if dx < 0 else "left")
ax.set_xscale("log"); ax.set_xlabel("Genes on panel (in annotated universe, log scale)")
ax.set_ylabel("Gut/IBD-lens pathways scorable (%)\n(≥5 member genes on panel)")
ax.grid(color=GRID, linewidth=0.6); ax.set_axisbelow(True); ax.set_ylim(-3, 103)
handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=k, markersize=6) for k, c in PLAT_COL.items()]
ax.legend(handles=handles, frameon=False, fontsize=7.5, loc="upper left")
save(fig, "fig1_size_vs_breadth.png")

# ---------- Fig 2: gene Jaccard + pathway-scope similarity (targeted panels)
J = pd.read_csv(O + "panel_jaccard.csv", index_col=0).loc[TARGETED_RNA, TARGETED_RNA]
S = pd.read_csv(O + "panel_pathway_similarity.csv", index_col=0).loc[TARGETED_RNA, TARGETED_RNA]
order = leaves_list(linkage(1 - S.values[np.triu_indices(len(S), 1)], "average"))
lbl = [SHORT[TARGETED_RNA[i]] for i in order]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
for ax, Mx, title, cmap, vmin in [(axs[0], J, "Shared genes (Jaccard)", SEQ, 0), (axs[1], S, "Shared biology (cosine of pathway-enrichment profiles)", SEQ, 0)]:
    V = Mx.values[np.ix_(order, order)]
    im = ax.imshow(V, cmap=cmap, vmin=vmin, vmax=1)
    ax.set_xticks(range(len(lbl))); ax.set_xticklabels(lbl, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(lbl))); ax.set_yticklabels(lbl, fontsize=7)
    for i in range(len(lbl)):
        for j in range(len(lbl)):
            if i != j:
                ax.text(j, i, f"{V[i, j]:.2f}".lstrip("0").replace("-0", "-"), ha="center", va="center", fontsize=5.2,
                        color="white" if abs(V[i, j]) > 0.6 else INK2)
    ax.set_title(title, fontsize=9, color=INK, loc="left"); ax.spines[:].set_visible(False)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
fig.tight_layout(); save(fig, "fig2_panel_redundancy_matrices.png")

# ---------- Fig 3: gut/IBD pathway coverage heatmap
cov = pd.read_csv(O + "coverage_long.csv.gz")
pick = ["HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_IL2_STAT5_SIGNALING",
        "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT", "HALLMARK_APOPTOSIS", "HALLMARK_HYPOXIA",
        "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_NOTCH_SIGNALING",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_XENOBIOTIC_METABOLISM",
        "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM", "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
        "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION", "KEGG_TIGHT_JUNCTION", "KEGG_JAK_STAT_SIGNALING_PATHWAY",
        "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY", "REACTOME_INTERLEUKIN_10_SIGNALING", "REACTOME_INTERLEUKIN_17_SIGNALING",
        "REACTOME_O_LINKED_GLYCOSYLATION_OF_MUCINS", "GAO_LARGE_INTESTINE_24W_C10_ENTEROCYTE", "GAO_LARGE_INTESTINE_24W_C8_GOBLET_CELL",
        "GAO_LARGE_INTESTINE_ADULT_CE_OLFM4HIGH_STEM_CELL", "GAO_LARGE_INTESTINE_ADULT_CA_ENTEROENDOCRINE_CELLS",
        "BUSSLINGER_DUODENAL_PANETH_CELLS", "BUSSLINGER_DUODENAL_TUFT_CELLS", "GAO_LARGE_INTESTINE_ADULT_CJ_IMMUNE_CELLS",
        "GAO_LARGE_INTESTINE_ADULT_CI_MESENCHYMAL_CELLS"]
pick = [p for p in pick if p in set(cov.gene_set)]
cols = ["FOLI_mPS3", "FOLI_hPS2", "COSMX_HIO", "MER_M_BRAINTRAIN", "MER_M_PANNEURO", "MER_M_PAN", "COSMX_MNEURO", "COSMX_HUCC", "COSMX_MUCC",
        "GEOMX_H_CTA", "GEOMX_DOG_CCA", "XEN_H5K", "XEN_M5K", "COSMX_H6K", "GEOMX_H_WTA"]
H = cov[cov.gene_set.isin(pick) & cov.panel_id.isin(cols)].pivot(index="gene_set", columns="panel_id", values="coverage").loc[pick, cols]
size = cov[cov.panel_id == cols[0]].set_index("gene_set").set_size.reindex(pick)
fig, ax = plt.subplots(figsize=(9.5, 9))
im = ax.imshow(H.values * 100, cmap=SEQ, vmin=0, vmax=100, aspect="auto")
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        v = H.values[i, j] * 100
        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.3, color="white" if v > 55 else INK2)
pretty = [p.replace("HALLMARK_", "H: ").replace("KEGG_", "KEGG: ").replace("REACTOME_", "R: ").replace("GAO_LARGE_INTESTINE_", "Colon: ")
          .replace("BUSSLINGER_DUODENAL_", "Duodenum: ").replace("_", " ").lower().capitalize() + f"  (n={size[p]})" for p in pick]
ax.set_yticks(range(len(pick))); ax.set_yticklabels(pretty, fontsize=7)
ax.set_xticks(range(len(cols))); ax.set_xticklabels([SHORT[c] for c in cols], rotation=55, ha="right", fontsize=7.5)
ax.spines[:].set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label("% of pathway genes on panel")
save(fig, "fig3_gut_pathway_coverage.png")

# ---------- Fig 4: how many targeted spatial panels carry each Foli-seq gene
P = pd.read_csv(O + "foli_gene_platform_presence.csv")
bins = [-1, 0, 3, 7, 12]; labels_ = ["0 panels", "1-3", "4-7", "8-12"]
fig, ax = plt.subplots(figsize=(6.4, 3.2))
w = 0.38
for k, (fp, col) in enumerate([("FOLI_mPS3", CAT[7]), ("FOLI_hPS2", CAT[0])]):
    sub = P[P.foli_panel == fp]
    cnt = pd.cut(sub.n_spatial_panels, bins, labels=labels_).value_counts().reindex(labels_)
    x_ = np.arange(len(labels_)) + (k - .5) * w
    ax.bar(x_, cnt.values, width=w - .04, color=col, label=f"{SHORT[fp]} (n={len(sub)})")
    for xx, v in zip(x_, cnt.values):
        ax.text(xx, v + 4, int(v), ha="center", fontsize=7, color=INK2)
ax.set_xticks(range(len(labels_))); ax.set_xticklabels(labels_)
ax.set_xlabel("Targeted spatial panels (of 12) that include the gene"); ax.set_ylabel("Foli-seq genes")
ax.grid(axis="y", color=GRID); ax.set_axisbelow(True); ax.legend(frameon=False, fontsize=7.5)
save(fig, "fig4_foli_genes_on_platforms.png")

# ---------- Fig 5: within-panel redundancy vs random
if os.path.exists(O + "panel_redundancy_by_dataset.csv"):
    R = pd.read_csv(O + "panel_redundancy_by_dataset.csv")
    R = R[R.panel_id.isin(TARGETED_RNA + ["GEOMX_H_WTA"]) & (R.n_measured >= 10)]
    agg = R.groupby("panel_id")[["frac_self_predictable", "random_frac_self_predictable", "eff_dim_frac", "random_eff_dim_frac",
                                 "tx_recon_R2", "random_tx_recon_R2", "n_measured"]].mean()
    agg = agg.sort_values("frac_self_predictable")
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
    y = np.arange(len(agg))
    for ax, a, b, xl in [(axs[0], "frac_self_predictable", "random_frac_self_predictable", "Genes recoverable from rest of panel (held-out R² ≥ 0.7), %"),
                         (axs[1], "tx_recon_R2", "random_tx_recon_R2", "Off-panel transcriptome reconstructed (mean held-out R²)")]:
        sc_ = 100 if a.startswith("frac") else 1
        for i, (_, r) in enumerate(agg.iterrows()):
            ax.plot([r[b] * sc_, r[a] * sc_], [i, i], color=GRID, lw=2, zorder=1)
        ax.scatter(agg[b] * sc_, y, color=MUTED, s=30, label="size-matched random genes", zorder=2, edgecolor="white")
        ax.scatter(agg[a] * sc_, y, color=CAT[0], s=40, label="vendor panel", zorder=3, edgecolor="white")
        ax.set_xlabel(xl, fontsize=8); ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
    axs[0].set_yticks(y); axs[0].set_yticklabels([SHORT[p] for p in agg.index], fontsize=7.5)
    axs[0].legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.tight_layout(); save(fig, "fig5_within_panel_redundancy.png")
    agg.to_csv(O + "panel_redundancy_summary.csv")

# ---------- Fig 6: gut/IBD-focused optimisation curve vs vendor panels
if os.path.exists(O + "gut_focused_curve.csv"):
    C = pd.read_csv(O + "gut_focused_curve.csv"); V = pd.read_csv(O + "vendor_gut_lens_R2.csv")
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.plot(C["size"], C.optimized_R2, color=CAT[0], lw=2, marker="o", ms=4, label="Covariance-optimised (gut/IBD lens genes)")
    ax.plot(C["size"], C.random_lens_R2, color=MUTED, lw=2, ls="--", label="Random gut/IBD lens genes")
    V = V[V.panel_id.isin(TARGETED_RNA + ["GEOMX_H_WTA"])]
    show = {"FOLI_mPS3": ("Foli-seq mPS3", (-5, -12)), "FOLI_hPS2": ("Foli-seq hPS2", (-5, 6)), "COSMX_HIO": ("CosMx IO add-on", (5, -10)),
            "MER_M_BRAINTRAIN": ("MERSCOPE Brain 128", (5, -3)), "MER_M_PANNEURO": ("MERSCOPE PanNeuro", (5, -3)),
            "MER_M_PAN": ("MERSCOPE Pan-Mouse", (5, -10)), "COSMX_HUCC": ("CosMx UCC 1K (Hs, Mm), Neuro 1K", (5, -14))}
    for _, r in V.iterrows():
        ax.scatter(r.n_in_pool, r.gut_lens_R2, color=PLAT_COL[plat[r.panel_id]], s=40, zorder=3, edgecolor="white")
        if r.panel_id in show:
            t, (dx, dy) = show[r.panel_id]
            ax.annotate(t, (r.n_in_pool, r.gut_lens_R2), xytext=(dx, dy), textcoords="offset points",
                        fontsize=7, color=INK2, ha="right" if dx < 0 else "left")
    big = V[V.n_in_pool >= 1000]
    ax.annotate("GeoMx CTA / Canine, Xenium 5K,\nCosMx 6K, GeoMx WTA", (big.n_in_pool.min(), big.gut_lens_R2.max()),
                xytext=(-6, 9), textcoords="offset points", fontsize=7, color=INK2)
    ax.set_ylim(top=0.85)
    ax.axhline(big.gut_lens_R2.median(), color=GRID, lw=1, ls=":")
    ax.set_xscale("log"); ax.set_xlabel("Genes on panel that are measured in the public data (log scale)")
    ax.set_ylabel("Held-out R² of gut/IBD pathway &\ncell-type scores")
    ax.grid(color=GRID, lw=.6); ax.set_axisbelow(True); ax.legend(frameon=False, fontsize=8, loc="lower right")
    save(fig, "fig6_optimisation_curve.png")

# ---------- Fig 7: pathway redundancy (split-half validated core size)
if os.path.exists(O + "pathway_core_summary.csv"):
    W = pd.read_csv(O + "pathway_core_summary.csv")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    cmap_ = {"Hallmark": CAT[0], "KEGG": CAT[1], "Reactome": CAT[2], "C8 cell types": CAT[6]}
    for c, col in cmap_.items():
        w = W[W.collection == c]
        ax.scatter(w.n_members, w.core_size, s=22, color=col, alpha=.85, edgecolor="white", linewidth=.6, label=f"{c} (n={len(w)})")
    xs = np.geomspace(10, 330, 80)
    ax.plot(xs, xs, color=MUTED, lw=1, ls=":"); ax.text(20.5, 18.6, "y = all members", fontsize=7, color=MUTED)
    ax.plot(xs, xs * .1, color=MUTED, lw=1, ls="--"); ax.text(205, 18.6, "y = 10% of members", fontsize=7, color=MUTED, ha="right")
    ax.set_xscale("log"); ax.set_ylim(0, 20); ax.set_xlim(9, 360)
    ax.set_xlabel("Pathway / cell-type genes measured (median over datasets, log scale)")
    ax.set_ylabel("Core genes that reproduce the full score\n(held-out r median 0.87)")
    ax.grid(color=GRID, lw=.6); ax.set_axisbelow(True); ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1))
    save(fig, "fig7_pathway_redundancy.png")

# ---------- Fig 8: Foli-seq augmentation curves + recommendation counts
if os.path.exists(O + "foli_focus_augmentation_curve_FOLI_mPS3.csv"):
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.2), gridspec_kw=dict(width_ratios=[1, 1.1]))
    ax = axs[0]
    for fp, col in [("FOLI_mPS3", CAT[7]), ("FOLI_hPS2", CAT[0])]:
        f = O + f"foli_focus_augmentation_curve_{fp}.csv"
        if not os.path.exists(f):
            continue
        A = pd.read_csv(f); a = A[A.baseline.isna()]; b = A[A.baseline.notna()]
        ax.plot(a.n_added, a.R2, color=col, lw=2, marker="o", ms=4, label=f"{SHORT[fp]}: covariance-ranked")
        ax.plot(b.n_added, b.R2, color=col, lw=1.5, ls="--", alpha=.7, label=f"{SHORT[fp]}: random lens genes")
    ax.set_xlabel("Genes added to the current panel"); ax.set_ylabel("Held-out R², gut/IBD pathway scores")
    ax.grid(color=GRID, lw=.6); ax.set_axisbelow(True); ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax = axs[1]
    if os.path.exists(O + "foli_candidate_genes.csv"):
        T = pd.read_csv(O + "foli_candidate_genes.csv")
        cats = ["Keep: responds in Foli-seq stool data", "Keep: unique or core signal", "Keep: supporting",
                "Keep: no evidence either way", "Review: redundant with a panel partner",
                "Add: gut epithelial marker", "Add: covariance-ranked", "Add: spatial consensus"]
        cc = [CAT[2], CAT[0], "#86b6ef", "#cde2fb", CAT[7], CAT[6], "#9085e9", "#c3bdf3"]
        n = T.recommendation.value_counts().reindex(cats).fillna(0)
        ax.barh(range(len(cats)), n.values, color=cc, height=.66)
        for i_, v in enumerate(n.values):
            ax.text(v + 3, i_, int(v), va="center", fontsize=7.5, color=INK2)
        ax.set_yticks(range(len(cats))); ax.set_yticklabels(cats, fontsize=7.5); ax.invert_yaxis()
        ax.set_xlabel("Genes"); ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
    fig.tight_layout(); save(fig, "fig8_foli_redundancy_and_additions.png")

# ---------- Fig 9: can spatial panels read out Foli-seq genes?
if os.path.exists(O + "foli_spatial_proxy_summary.csv"):
    SS = pd.read_csv(O + "foli_spatial_proxy_summary.csv")
    SS["proxy_mid"] = SS.proxy_r04 - SS.proxy_r06
    parts = [("direct", "On the panel", CAT[0]), ("proxy_r06", "Proxy |r| ≥ 0.6", "#86b6ef"), ("proxy_mid", "Proxy 0.4–0.6", "#cde2fb"), ("no_coverage", "No usable proxy", "#e6e5e1")]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), sharey=False)
    for ax, fp in zip(axs, ["FOLI_mPS3", "FOLI_hPS2"]):
        S = SS[SS.foli_panel == fp].set_index("panel_id").sort_values(["direct", "proxy_r06"])
        left = np.zeros(len(S))
        for col, lab_, c in parts:
            ax.barh(range(len(S)), S[col], left=left, color=c, height=.66, label=lab_, edgecolor="white", linewidth=1.2)
            for i_, (l, v) in enumerate(zip(left, S[col])):
                if v >= 30:
                    ax.text(l + v / 2, i_, int(v), ha="center", va="center", fontsize=6.8, color="white" if c == CAT[0] else INK2)
            left += S[col].values
        ax.set_yticks(range(len(S))); ax.set_yticklabels([SHORT[p] for p in S.index], fontsize=7.5)
        ax.set_xlabel(f"{SHORT[fp]} genes (of {int(S.n_foli_genes.iloc[0])})"); ax.spines["left"].set_visible(False)
    axs[0].legend(frameon=False, fontsize=7.5, ncol=4, loc="upper left", bbox_to_anchor=(0.0, 1.16))
    fig.tight_layout(); save(fig, "fig9_foli_spatial_readout.png")
print("figures written:", sorted(os.listdir(F)))
