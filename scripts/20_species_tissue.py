"""Separate the species effect from the tissue effect, using the 2x2 the Haber reference completes.

Without a mouse gut reference, "immune vs gut" and "mouse vs human" move together and neither can be
estimated. HABER_GUT is mouse AND gut, which breaks the tie:

    tissue effect, species held fixed (mouse)   =  HABER_GUT  -  MM_SLN
    species effect, tissue held fixed (gut)     =  human gut  -  HABER_GUT

RAW ACCURACY IS NOT COMPARABLE ACROSS REFERENCES. Each reference poses a different problem: 26 immune
cell types, 15 epithelial ones, 72/54/45 intestinal ones. Chance alone moves with the number of
classes. So every panel is scored as LIFT OVER THE SIZE-MATCHED RANDOM BASELINE measured on the same
reference (RANDOM_535, the baseline common to every run):

    lift(panel, ref) = mean_accuracy(panel, ref) - mean_accuracy(RANDOM_535, ref)

Lift asks "how much does curating these genes buy, on this problem" and is comparable across
references in a way raw accuracy is not.
"""
import os, sys, importlib.util
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ

spec = importlib.util.spec_from_file_location("f3", os.path.join(os.path.dirname(os.path.abspath(__file__)), "18_figures_v3.py"))
f3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f3)
O, F = PROJ + "results/v2/", PROJ + "figures/v3/"
os.makedirs(F, exist_ok=True)

CELL = {"MM_SLN": ("Mouse", "Immune"), "HABER_GUT": ("Mouse", "Gut epithelium"),
        "GUT_HEALTHY_EPI": ("Human", "Gut epithelium"),
        "GUT_HEALTHY": ("Human", "Gut whole"), "GUT_IBD": ("Human", "Gut whole"),
        "GUT_CRC": ("Human", "Gut whole")}
HUMAN_GUT = ["GUT_HEALTHY", "GUT_IBD", "GUT_CRC"]
# compartment-matched human counterpart to HABER_GUT; see 21_build_epi_matched_ref.py
HUMAN_EPI = "GUT_HEALTHY_EPI"

rows = []
tables = {}
for ref in CELL:
    A, st = f3.accuracy_table(ref)
    if A is None:
        print(f"{ref}: missing, skipped"); continue
    base_cols = [c for c in A.columns if c.startswith("RANDOM_535")]
    if not base_cols:
        print(f"{ref}: no RANDOM_535 baseline, skipped"); continue
    base = A[base_cols].values.mean()
    tables[ref] = A
    sp_, ti = CELL[ref]
    for p in A.columns:
        if p.startswith("RANDOM"): continue
        rows.append(dict(reference=ref, species=sp_, tissue=ti, panel_id=p,
                         mean_acc=A[p].mean(), random535=base, lift=A[p].mean() - base,
                         n_celltypes=A.shape[0], status=st))
L = pd.DataFrame(rows)
L.to_csv(O + "q8_lift_by_reference.csv", index=False)
print(f"references used: {sorted(tables)}")
print("\nrandom-535 baseline per reference (accuracy %):")
for r in sorted(tables):
    print(f"  {r:13s} {L[L.reference==r].random535.iloc[0]:5.1f}%   ({L[L.reference==r].n_celltypes.iloc[0]} cell types)")

piv = L.pivot_table(index="panel_id", columns="reference", values="lift")
shared = piv.dropna().index.tolist()
piv = piv.loc[shared]
piv["TISSUE_mouse"] = piv["HABER_GUT"] - piv["MM_SLN"]
piv["SPECIES_gut"] = piv[[c for c in HUMAN_GUT if c in piv.columns]].mean(axis=1) - piv["HABER_GUT"]
if HUMAN_EPI in piv.columns:
    # cleaner species estimate: both sides are gut EPITHELIUM, so compartment is held fixed too
    piv["SPECIES_epi_matched"] = piv[HUMAN_EPI] - piv["HABER_GUT"]
    piv["COMPARTMENT_human"] = piv[[c for c in HUMAN_GUT if c in piv.columns]].mean(axis=1) - piv[HUMAN_EPI]
piv.round(1).to_csv(O + "q8_species_tissue_effects.csv")

print("\nLIFT OVER RANDOM-535 (percentage points), and the two decomposed effects\n")
cols = (["MM_SLN", "HABER_GUT"] + ([HUMAN_EPI] if HUMAN_EPI in piv.columns else [])
        + [c for c in HUMAN_GUT if c in piv.columns] + ["TISSUE_mouse", "SPECIES_gut"]
        + (["SPECIES_epi_matched", "COMPARTMENT_human"] if HUMAN_EPI in piv.columns else []))
print(piv[cols].round(1).sort_values("SPECIES_gut").to_string())

print("\nFoli-seq panels — which reference each one favours (lift, pp):")
for p in ["FOLI_mPS3", "FOLI_hPS2"]:
    if p in piv.index:
        r = piv.loc[p]
        hg = np.mean([r[c] for c in HUMAN_GUT if c in piv.columns])
        print(f"  {p:10s} mouse immune {r.MM_SLN:5.1f} | mouse gut {r.HABER_GUT:5.1f} | human gut {hg:5.1f}")
if {"FOLI_mPS3", "FOLI_hPS2"} <= set(piv.index):
    d_mouse = piv.loc["FOLI_mPS3", "HABER_GUT"] - piv.loc["FOLI_hPS2", "HABER_GUT"]
    d_human = (piv.loc["FOLI_mPS3", HUMAN_EPI] - piv.loc["FOLI_hPS2", HUMAN_EPI]
               if HUMAN_EPI in piv.columns
               else piv.loc["FOLI_mPS3", HUMAN_GUT].mean() - piv.loc["FOLI_hPS2", HUMAN_GUT].mean())
    print(f"\n  mPS3 minus hPS2 on MOUSE gut: {d_mouse:+.1f} pp")
    print(f"  mPS3 minus hPS2 on HUMAN gut epithelium (compartment-matched): {d_human:+.1f} pp")
    print(f"  -> species-match swing: {d_mouse - d_human:+.1f} pp"
          f"{'  (each panel favours its own species)' if d_mouse > 0 > d_human else ''}")

# ---- figure ----
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
BLUE, RED, ORANGE = "#2a78d6", "#e34948", "#e38f1f"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 200, "savefig.bbox": "tight",
                     "figure.facecolor": "white"})
fig, axs = plt.subplots(1, 2, figsize=(12.4, 5.0))
ordr = piv["SPECIES_gut"].sort_values().index
# the two Foli-seq panels must be distinguishable: they are the comparison the figure exists for
FOLI_STYLE = {"FOLI_mPS3": (RED, "o", "-"), "FOLI_hPS2": (ORANGE, "s", "--")}

ax = axs[0]
xs = ["MM_SLN", "HABER_GUT"] + ([HUMAN_EPI] if HUMAN_EPI in piv.columns else []) + [c for c in HUMAN_GUT if c in piv.columns]
lab = ["mouse\nimmune", "mouse gut\nepithelium", "human gut\nepithelium",
       "human gut\nwhole (healthy)", "human gut\nwhole (IBD)", "human gut\nwhole (CRC)"][:len(xs)]
for p in ordr:
    if p in FOLI_STYLE:
        c, mk, ls = FOLI_STYLE[p]
        ax.plot(range(len(xs)), piv.loc[p, xs], marker=mk, ls=ls, ms=5.5, lw=2.0,
                color=c, alpha=.98, label=f3.short(p), zorder=5)
    else:
        ax.plot(range(len(xs)), piv.loc[p, xs], "o-", ms=3.5, lw=.8, color=BLUE, alpha=.40)
for xv, txt in [(.5, " tissue "), (1.5, " species "), (2.5, " compartment ")]:
    if xv < len(xs) - .5:
        ax.axvline(xv, color=MUTED, lw=.8, ls=":")
        ax.text(xv, ax.get_ylim()[1], txt, fontsize=7, color=MUTED, va="top", ha="center")
ax.set_xticks(range(len(xs))); ax.set_xticklabels(lab, fontsize=8)
ax.set_ylabel("lift over size-matched random (pp)", fontsize=8.5)
ax.set_title("a  Each step changes exactly one thing", fontsize=9.5, loc="left", color=INK)
ax.legend(fontsize=8, frameon=False, loc="lower right"); ax.grid(axis="y", color=GRID); ax.set_axisbelow(True)

ax = axs[1]
y = np.arange(len(ordr))
ax.barh(y - .2, piv.loc[ordr, "TISSUE_mouse"], height=.38, color="#6a9fd8", label="tissue (mouse immune → mouse gut)")
key = "SPECIES_epi_matched" if "SPECIES_epi_matched" in piv.columns else "SPECIES_gut"
ax.barh(y + .2, piv.loc[ordr, key], height=.38, color="#d88a6a",
        label="species (mouse → human gut epithelium)" if key == "SPECIES_epi_matched"
              else "species (mouse gut → human gut)")
ax.axvline(0, color=MUTED, lw=.9)
ax.set_yticks(y); ax.set_yticklabels([f3.short(p) for p in ordr], fontsize=7.5)
ax.set_xlabel("change in lift (pp)", fontsize=8.5)
ax.set_title("b  Decomposed, with the other factor held fixed", fontsize=9.5, loc="left", color=INK)
ax.legend(fontsize=7.5, frameon=False, loc="upper left", bbox_to_anchor=(0, 1.02))
ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
ax.set_ylim(-0.7, len(ordr) - 0.3)
fig.suptitle("Separating species from tissue — the 2×2 the Haber mouse-gut reference completes",
             fontsize=10.5, color=INK, x=.01, ha="left")
fig.tight_layout(rect=[0, 0, 1, .95])
fig.savefig(F + "v3_fig8_species_vs_tissue.png"); plt.close(fig)
print("\nwrote figures/v3/v3_fig8_species_vs_tissue.png")
