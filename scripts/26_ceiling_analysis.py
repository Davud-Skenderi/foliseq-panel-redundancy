"""Express each targeted panel as a fraction of what measuring everything achieves.

Lift says how far a panel beats a random gene set of the same size. It does not say how far it is
from the best possible. The CosMx Whole Transcriptome panels supply that upper bound: they measure
essentially everything the platform can, so their accuracy is the ceiling for any panel on it.

With a floor (size-matched random) and a ceiling (whole transcriptome), a panel's score becomes
interpretable on an absolute scale:

    ceiling fraction = (panel - random) / (whole transcriptome - random)

0 means the panel does no better than picking genes at random. 1 means it does as well as measuring
the entire transcriptome. This is the number to quote when asking "how much is a targeted panel
giving up", because unlike lift it has a meaningful maximum.

Only the two compartment-matched references are used, since those are where the ceiling was run.
"""
import os, sys, importlib.util
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ

spec = importlib.util.spec_from_file_location(
    "f3", os.path.join(os.path.dirname(os.path.abspath(__file__)), "18_figures_v3.py"))
f3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f3)
O, F = PROJ + "results/v2/", PROJ + "figures/v3/"
REFS = {"HABER_GUT": "Mouse gut epithelium (Haber, 15 types)",
        "GUT_HEALTHY_EPI": "Human gut epithelium (Gut Cell Atlas, 12 types)"}
WTX = ["COSMX_H_WTX", "COSMX_M_WTX"]

rows = []
for ref, title in REFS.items():
    A, st = f3.accuracy_table(ref)
    if A is None:
        print(f"{ref}: missing"); continue
    rand = [c for c in A.columns if c.startswith("RANDOM_535")]
    wtx = [c for c in WTX if c in A.columns]
    if not rand or not wtx:
        print(f"{ref}: need both baseline and ceiling; have rand={len(rand)} wtx={len(wtx)}"); continue
    floor = A[rand].values.mean()
    ceil = A[wtx].mean().mean()
    print(f"\n{ref}: floor (random 535) {floor:.1f}%   ceiling (whole transcriptome) {ceil:.1f}%"
          f"   usable range {ceil-floor:.1f} pp")
    for p in A.columns:
        if p.startswith("RANDOM") or p in WTX:
            continue
        acc = A[p].mean()
        rows.append(dict(reference=ref, panel_id=p, panel_acc=acc, floor=floor, ceiling=ceil,
                         lift=acc - floor,
                         ceiling_fraction=(acc - floor) / (ceil - floor) if ceil > floor else np.nan))
R = pd.DataFrame(rows)
R["ceiling_pct"] = (R.ceiling_fraction * 100).round(1)
R.to_csv(O + "q12_ceiling_fraction.csv", index=False)

for ref in REFS:
    sub = R[R.reference == ref].sort_values("ceiling_fraction", ascending=False)
    if not len(sub): continue
    print(f"\n=== {REFS[ref]} — % of the random-to-whole-transcriptome range captured ===")
    print(sub[["panel_id", "panel_acc", "lift", "ceiling_pct"]].round(1).to_string(index=False))

# ---- figure ----
piv = R.pivot_table(index="panel_id", columns="reference", values="ceiling_pct").dropna()
if len(piv):
    INK, MUTED, GRID, RED, BLUE, GREEN = "#0b0b0b", "#52514e", "#e6e5e1", "#e34948", "#2a78d6", "#3f9e6a"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 200, "savefig.bbox": "tight",
                         "figure.facecolor": "white"})
    piv = piv.sort_values("GUT_HEALTHY_EPI")
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    y = np.arange(len(piv))
    ax.barh(y - .2, piv.HABER_GUT, height=.38, color=BLUE, label=REFS["HABER_GUT"])
    ax.barh(y + .2, piv.GUT_HEALTHY_EPI, height=.38, color=GREEN, label=REFS["GUT_HEALTHY_EPI"])
    ax.axvline(0, color=MUTED, lw=1)
    ax.axvline(100, color=MUTED, lw=1, ls="--")
    ax.text(100, len(piv) - .3, " whole transcriptome", fontsize=7.5, color=MUTED, va="center")
    ax.set_yticks(y)
    labs = ax.set_yticklabels([f3.short(p) for p in piv.index], fontsize=8)
    for lab, p in zip(labs, piv.index):
        if p.startswith("FOLI"):
            lab.set_color(RED); lab.set_fontweight("bold")
    ax.set_xlabel("% of the range from random-535 to whole transcriptome", fontsize=8.5)
    ax.set_title("How close does a targeted panel get to measuring everything?",
                 fontsize=10, loc="left", color=INK)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(F + "v3_fig10_ceiling_fraction.png"); plt.close(fig)
    print("\nwrote figures/v3/v3_fig10_ceiling_fraction.png")
