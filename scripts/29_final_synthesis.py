"""Final synthesis: the complete species x tissue grid, with error bars where they exist.

WHAT THIS PULLS TOGETHER.

1. ERROR BARS. Every accuracy in the repo had been a single classifier run, so a 2-point difference
   could not be told from noise. 17_multiref_eval.py now passes n_seeds=5 and keeps the standard
   deviation matrix that spapros returns as out[1] -- which 12_spapros_eval.py discarded, and which is
   the reason nothing had an interval before. HABER_GUT and GUT_HEALTHY_EPI have been re-run at
   5 seeds; MM_SLN and HUMAN_IMMUNE remain at 1 seed, matched to each other.

2. THE COMPLETED GRID. With the Human Immune Health Atlas added, all four cells are filled:

                      immune                        gut
      mouse           MM_SLN (26 types)             HABER_GUT (15 types)
      human           HUMAN_IMMUNE (37 types)       GUT_HEALTHY_EPI (12 types)

   Each effect can now be estimated TWICE, holding the other factor fixed, and checked for
   replication:
      tissue effect  in mouse   = HABER_GUT  - MM_SLN
      tissue effect  in human   = GUT_EPI    - HUMAN_IMMUNE
      species effect on immune  = HUMAN_IMMUNE - MM_SLN
      species effect on gut     = GUT_EPI      - HABER_GUT

3. WHY REPLICATION MATTERS MORE THAN EITHER ESTIMATE. Lift still correlates with the number of cell
   types (r = +0.78), so each single estimate carries that artifact. But the artifact's DIRECTION
   depends on whether class count rises or falls:
      tissue in mouse   26 -> 15 types  (falls)
      tissue in human   37 -> 12 types  (falls)      -- same direction, artifact could drive both
      species on immune 26 -> 37 types  (RISES)
      species on gut    15 -> 12 types  (falls)      -- OPPOSITE directions
   So a species effect with the same sign in both rows cannot be produced by the class-count
   artifact, because the artifact would push the two estimates apart. That is the strongest
   available check short of matching class counts everywhere.
"""
import os, sys, glob, importlib.util
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ

O = PROJ + "results/v2/"
spec = importlib.util.spec_from_file_location("f3", os.path.join(os.path.dirname(os.path.abspath(__file__)), "18_figures_v3.py"))
f3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f3)

PANELS = ["FOLI_mPS3", "FOLI_hPS2", "COSMX_HIO", "MER_M_BRAINTRAIN", "MER_M_PAN", "MER_M_PANNEURO",
          "COSMX_MNEURO", "COSMX_MUCC", "COSMX_HUCC"]
GRID = {"MM_SLN": ("Mouse", "Immune", 1), "HUMAN_IMMUNE": ("Human", "Immune", 1),
        "HABER_GUT": ("Mouse", "Gut epithelium", 5), "GUT_HEALTHY_EPI": ("Human", "Gut epithelium", 5)}


def table(ref, seeds):
    """Accuracy table, preferring the 5-seed version where it exists."""
    if seeds == 5:
        cols = {}
        for f in sorted(glob.glob(O + f"q4_confusion_s5_{ref}__*.csv")):
            p = os.path.basename(f).split("__", 1)[1][:-4]
            cm = pd.read_csv(f, index_col=0)
            cols[p] = pd.Series(np.diag(cm.values), index=cm.index)
        if cols:
            return pd.DataFrame(cols) * 100
    A, _ = f3.accuracy_table(ref)
    return A


def sd_table(ref):
    cols = {}
    for f in sorted(glob.glob(O + f"q4_confusion_sd_s5_{ref}__*.csv")):
        p = os.path.basename(f).split("__", 1)[1][:-4]
        cm = pd.read_csv(f, index_col=0)
        cols[p] = pd.Series(np.diag(cm.values), index=cm.index)
    return pd.DataFrame(cols) * 100 if cols else None


rows = []
for ref, (sp_, ti, seeds) in GRID.items():
    A = table(ref, seeds)
    if A is None:
        print(f"{ref}: missing"); continue
    SD = sd_table(ref) if seeds == 5 else None
    base_cols = [c for c in A.columns if c.startswith("RANDOM_535")]
    floor = A[base_cols].values.mean()
    print(f"{ref:18s} {A.shape[0]:3d} types  {len(A.columns):2d} panels  "
          f"random-535 floor {floor:5.1f}%  seeds={seeds}")
    for p in PANELS:
        if p not in A.columns:
            continue
        acc = A[p].mean()
        # spread of the per-cell-type mean across seeds+folds, propagated as the mean of per-type SDs
        # divided by sqrt(n_types): an approximation, and stated as such
        sd = (SD[p].mean() / np.sqrt(A.shape[0])) if SD is not None and p in SD.columns else np.nan
        rows.append(dict(reference=ref, species=sp_, tissue=ti, seeds=seeds, panel_id=p,
                         n_celltypes=A.shape[0], acc=acc, acc_se=sd, floor=floor, lift=acc - floor))
L = pd.DataFrame(rows)
L.to_csv(O + "q14_grid_lift.csv", index=False)

piv = L.pivot_table(index="panel_id", columns="reference", values="lift")
se = L.pivot_table(index="panel_id", columns="reference", values="acc_se")
piv["TISSUE_in_mouse"] = piv["HABER_GUT"] - piv["MM_SLN"]
piv["TISSUE_in_human"] = piv["GUT_HEALTHY_EPI"] - piv["HUMAN_IMMUNE"]
piv["SPECIES_on_immune"] = piv["HUMAN_IMMUNE"] - piv["MM_SLN"]
piv["SPECIES_on_gut"] = piv["GUT_HEALTHY_EPI"] - piv["HABER_GUT"]
piv.round(1).to_csv(O + "q14_grid_effects.csv")

print("\n=== LIFT OVER RANDOM-535, ALL FOUR CELLS (pp) ===")
cols = ["MM_SLN", "HUMAN_IMMUNE", "HABER_GUT", "GUT_HEALTHY_EPI"]
print(piv[cols].round(1).to_string())

print("\n=== EACH EFFECT ESTIMATED TWICE ===")
eff = ["TISSUE_in_mouse", "TISSUE_in_human", "SPECIES_on_immune", "SPECIES_on_gut"]
print(piv[eff].round(1).to_string())

print("\nreplication check (does the sign agree across the two independent estimates?)")
for a, b, name in [("TISSUE_in_mouse", "TISSUE_in_human", "tissue"),
                   ("SPECIES_on_immune", "SPECIES_on_gut", "species")]:
    agree = int((np.sign(piv[a]) == np.sign(piv[b])).sum())
    r = piv[[a, b]].corr().iloc[0, 1]
    print(f"  {name:8s} sign agrees for {agree}/{len(piv)} panels   corr of the two estimates r={r:+.2f}")

print("\n=== 5-seed error bars on the gut references (approx SE of the per-panel mean, pp) ===")
for ref in ["HABER_GUT", "GUT_HEALTHY_EPI"]:
    s = L[L.reference == ref].set_index("panel_id")
    print(f"\n{ref}")
    for p in PANELS:
        if p in s.index and np.isfinite(s.loc[p, "acc_se"]):
            print(f"   {p:18s} acc {s.loc[p,'acc']:5.1f} +/- {s.loc[p,'acc_se']:4.1f}   lift {s.loc[p,'lift']:+6.1f}")

if {"FOLI_mPS3", "FOLI_hPS2"} <= set(piv.index):
    dm = piv.loc["FOLI_mPS3", "HABER_GUT"] - piv.loc["FOLI_hPS2", "HABER_GUT"]
    dh = piv.loc["FOLI_mPS3", "GUT_HEALTHY_EPI"] - piv.loc["FOLI_hPS2", "GUT_HEALTHY_EPI"]
    sem = np.sqrt(se.loc["FOLI_mPS3", "HABER_GUT"]**2 + se.loc["FOLI_hPS2", "HABER_GUT"]**2)
    seh = np.sqrt(se.loc["FOLI_mPS3", "GUT_HEALTHY_EPI"]**2 + se.loc["FOLI_hPS2", "GUT_HEALTHY_EPI"]**2)
    print(f"\n=== Foli-seq species match, now with intervals ===")
    print(f"  mPS3 - hPS2 on MOUSE gut epithelium: {dm:+.1f} +/- {sem:.1f} pp")
    print(f"  mPS3 - hPS2 on HUMAN gut epithelium: {dh:+.1f} +/- {seh:.1f} pp")
    print(f"  swing: {dm-dh:+.1f} pp   (each panel favouring its own species: {dm > 0 > dh})")
    di = piv.loc["FOLI_mPS3", "HUMAN_IMMUNE"] - piv.loc["FOLI_hPS2", "HUMAN_IMMUNE"]
    dmi = piv.loc["FOLI_mPS3", "MM_SLN"] - piv.loc["FOLI_hPS2", "MM_SLN"]
    print(f"  for comparison on immune tissue: mouse {dmi:+.1f}, human {di:+.1f}")


def figure():
    """The completed grid, and the replication check that matters most."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    INK, MUTED, GRID_C = "#0b0b0b", "#52514e", "#e6e5e1"
    RED, ORANGE, BLUE = "#e34948", "#e38f1f", "#2a78d6"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 200, "savefig.bbox": "tight",
                         "figure.facecolor": "white"})
    L = pd.read_csv(O + "q14_grid_lift.csv")
    E = pd.read_csv(O + "q14_grid_effects.csv", index_col=0)
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 5.0))

    # (a) the four cells, with error bars where 5 seeds were run
    ax = axs[0]
    xs = ["MM_SLN", "HUMAN_IMMUNE", "HABER_GUT", "GUT_HEALTHY_EPI"]
    lab = ["mouse\nimmune\n(26)", "human\nimmune\n(37)", "mouse gut\nepithelium\n(15)",
           "human gut\nepithelium\n(12)"]
    st = {"FOLI_mPS3": (RED, "o", "-"), "FOLI_hPS2": (ORANGE, "s", "--")}
    for p in sorted(set(L.panel_id)):
        sub = L[L.panel_id == p].set_index("reference").reindex(xs)
        if p in st:
            c, mk, ls = st[p]
            ax.errorbar(range(len(xs)), sub.lift, yerr=sub.acc_se, color=c, marker=mk, ls=ls,
                        ms=5.5, lw=2.0, capsize=3, label=f3.short(p), zorder=5)
        else:
            ax.plot(range(len(xs)), sub.lift, "o-", color=BLUE, ms=3.5, lw=.8, alpha=.40)
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xticks(range(len(xs))); ax.set_xticklabels(lab, fontsize=7.5)
    ax.set_ylabel("lift over size-matched random (pp)", fontsize=8.5)
    ax.set_title("a  All four cells of the grid\nerror bars = 5 seeds (gut refs only); 0 = no better than random",
                 fontsize=9, loc="left", color=INK)
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    ax.grid(axis="y", color=GRID_C); ax.set_axisbelow(True)

    # (b) does each effect replicate across its two independent estimates?
    ax = axs[1]
    ax.scatter(E.TISSUE_in_mouse, E.TISSUE_in_human, s=46, color=BLUE, label="tissue effect", zorder=4)
    ax.scatter(E.SPECIES_on_immune, E.SPECIES_on_gut, s=46, color="#d88a6a", marker="^",
               label="species effect", zorder=4)
    for p in ["FOLI_mPS3", "FOLI_hPS2", "COSMX_HIO"]:
        if p in E.index:
            ax.annotate(f3.short(p), (E.loc[p, "TISSUE_in_mouse"], E.loc[p, "TISSUE_in_human"]),
                        fontsize=6.5, color=MUTED, xytext=(4, -8), textcoords="offset points")
    lim = [-40, 18]
    ax.plot(lim, lim, color=MUTED, lw=.8, ls=":")
    ax.axhline(0, color=MUTED, lw=.8); ax.axvline(0, color=MUTED, lw=.8)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("estimate 1  (tissue: in mouse  /  species: on immune)", fontsize=8.5)
    ax.set_ylabel("estimate 2  (tissue: in human  /  species: on gut)", fontsize=8.5)
    rt = E[["TISSUE_in_mouse", "TISSUE_in_human"]].corr().iloc[0, 1]
    rs = E[["SPECIES_on_immune", "SPECIES_on_gut"]].corr().iloc[0, 1]
    ax.set_title(f"b  Replication: tissue r={rt:+.2f}, species r={rs:+.2f}\n"
                 f"points on the dotted line replicate; tissue does, species does not",
                 fontsize=9, loc="left", color=INK)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.grid(color=GRID_C); ax.set_axisbelow(True)
    fig.suptitle("The completed species x tissue grid", fontsize=10.5, color=INK, x=.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, .95])
    fig.savefig(PROJ + "figures/v3/v3_fig11_complete_grid.png"); plt.close(fig)
    print("wrote figures/v3/v3_fig11_complete_grid.png")


if __name__ == "__main__":
    figure()
