"""Control for the number of cell types, which the lift metric does not fully remove.

THE PROBLEM. Lift subtracts a size-matched random baseline measured on the same reference, which is
meant to cancel out how hard that reference's classification problem is. It does not fully succeed:
across the six references, corr(number of cell types, mean lift) = +0.78. References with more cell
types show systematically higher lift, because a harder problem leaves more headroom for curation.

WHY THAT THREATENS THE TISSUE RESULT. The tissue effect is measured MM_SLN (26 cell types) ->
HABER_GUT (15). Cell-type count drops, and the artifact predicts lift should drop with it. The
observed tissue effect is negative for every panel, which is exactly the direction the artifact would
produce. So "panels are immune-weighted" and "Haber simply has fewer cell types" cannot be told apart
from the existing runs.

(The species step is not threatened: HABER_GUT 15 types -> GUT_HEALTHY_EPI 12 types DECREASES the
count while lift RISES, against the artifact's direction.)

THE CONTROL. Re-run MM_SLN restricted to a random 15 of its 26 cell types, repeated over independent
draws, so the immune reference poses a 15-class problem exactly like Haber. Then:

    tissue effect, cell-type count matched = lift(HABER_GUT) - lift(MM_SLN restricted to 15 types)

If that stays negative for every panel, the tissue result is not an artifact of class count.
Everything else is held identical: same panels, same classifier settings, same baselines.
"""
import os, sys, time, warnings
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ, EXT
from common import to_human
from spapros.evaluation import metrics as spm
import scipy.sparse as sp
warnings.filterwarnings("ignore")

D, O = PROJ + "data/", PROJ + "results/v2/"
N_TYPES, N_DRAWS = 15, 3          # 15 = HABER_GUT's cell-type count
PANELS = ["FOLI_mPS3", "FOLI_hPS2", "COSMX_HIO", "MER_M_BRAINTRAIN", "MER_M_PAN", "MER_M_PANNEURO",
          "COSMX_MNEURO", "COSMX_MUCC", "COSMX_HUCC"]
N_RANDOM, RANDOM_SIZE = 3, 535
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)


def to_symbols(a, species):
    hs = [to_human(g, species) for g in a.var_names.astype(str)]
    keep = np.array([isinstance(h, str) and len(h) > 0 for h in hs])
    a = a[:, keep].copy()
    h = np.array(hs, dtype=object)[keep]
    codes, uniq = pd.factorize(h)
    Mx = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
    return ad.AnnData(sp.csr_matrix(a.X) @ Mx, obs=a.obs.copy(), var=pd.DataFrame(index=uniq))


def sample_per_type(a, key, n=200, min_cells=40, seed=0):
    rng = np.random.default_rng(seed)
    vc = a.obs[key].value_counts()
    idx = []
    for t in vc.index[vc >= min_cells]:
        ii = np.where(a.obs[key].values == t)[0]
        idx += list(rng.choice(ii, min(n, len(ii)), replace=False))
    return a[sorted(idx)].copy()


def load_mm_sln():
    """Identical preparation to 12_spapros_eval.py, so the restricted run differs only in class count."""
    a1 = ad.read_h5ad(EXT + "expr/sln_111.h5ad"); a2 = ad.read_h5ad(EXT + "expr/sln_208.h5ad")
    a = ad.concat([a1, a2], join="inner"); a.obs_names_make_unique()
    a = a[~a.obs.cell_types.str.contains("doublet|Low quality|Erythrocytes", case=False)].copy()
    a = to_symbols(a, "Mouse")
    a.obs["celltype"] = a.obs.cell_types.astype(str)
    sc.pp.filter_genes(a, min_cells=10)
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    return a


def main():
    a = load_mm_sln()
    types = sorted(a.obs.celltype.value_counts().pipe(lambda s: s[s >= 40]).index)
    print(f"MM_SLN: {a.shape}, {len(types)} usable cell types", flush=True)
    measured = set(a.var_names)
    rows = []
    for draw in range(N_DRAWS):
        rng = np.random.default_rng(100 + draw)
        chosen = sorted(rng.choice(types, N_TYPES, replace=False))
        sub = a[a.obs.celltype.isin(chosen)].copy()
        clf_ad = sample_per_type(sub, "celltype", seed=draw)
        print(f"\n-- draw {draw}: {len(chosen)} types, {clf_ad.shape[0]} cells", flush=True)
        sets = {p: [g for g in M.index[M[p] == 1] if g in measured] for p in PANELS}
        r2 = np.random.default_rng(draw)
        for k in range(N_RANDOM):
            sets[f"RANDOM_{RANDOM_SIZE}_{k}"] = list(r2.choice(a.var_names, RANDOM_SIZE, replace=False))
        for p, genes in sets.items():
            cache = O + f"q9_confusion_MM_SLN{N_TYPES}_d{draw}__{p}.csv"
            t0 = time.time()
            if os.path.exists(cache):
                cm = pd.read_csv(cache, index_col=0)
            else:
                cm = spm.xgboost_forest_classification(clf_ad, genes, celltypes="all", ct_key="celltype",
                                                       n_cells_min=40, n_seeds=1, n_jobs=2, verbosity=0)[0]
                cm.to_csv(cache)
            acc = float(np.diag(cm.values).mean()) * 100
            rows.append(dict(draw=draw, panel_id=p, n_genes=len(genes), n_celltypes=len(cm), mean_acc=acc))
            print(f"   {p:18s} {len(genes):5d} genes  acc {acc:5.1f}%  {time.time()-t0:4.0f}s", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(O + f"q9_mm_sln_{N_TYPES}types_raw.csv", index=False)

    base = R[R.panel_id.str.startswith("RANDOM")].groupby("draw").mean_acc.mean()
    R["lift"] = R.apply(lambda r: r.mean_acc - base[r.draw], axis=1)
    pan = R[~R.panel_id.str.startswith("RANDOM")]
    summ = pan.groupby("panel_id").agg(mm_sln15_lift=("lift", "mean"), sd=("lift", "std")).round(2)
    print(f"\nrandom-{RANDOM_SIZE} baseline per draw: " + ", ".join(f"{v:.1f}%" for v in base))

    hb = pd.read_csv(O + "q8_lift_by_reference.csv")
    hb = hb[hb.reference == "HABER_GUT"].set_index("panel_id").lift
    old = pd.read_csv(O + "q8_lift_by_reference.csv")
    old = old[old.reference == "MM_SLN"].set_index("panel_id").lift
    summ["haber_lift"] = hb
    summ["mm_sln26_lift"] = old
    summ["TISSUE_matched"] = (summ.haber_lift - summ.mm_sln15_lift).round(2)
    summ["TISSUE_unmatched"] = (summ.haber_lift - summ.mm_sln26_lift).round(2)
    summ.to_csv(O + "q9_celltype_count_control.csv")
    print("\n=== TISSUE EFFECT, CELL-TYPE COUNT MATCHED (both 15 classes) ===")
    print(summ.sort_values("TISSUE_matched").to_string())
    t = summ.TISSUE_matched.dropna()
    print(f"\nnegative for every panel? {(t < 0).all()}   range {t.min():.1f} to {t.max():.1f}")
    print(f"unmatched comparison ranged {summ.TISSUE_unmatched.min():.1f} to {summ.TISSUE_unmatched.max():.1f}")




def figure():
    """Matched vs unmatched tissue effect, so the size of the artifact is visible per panel."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    S = pd.read_csv(O + "q9_celltype_count_control.csv", index_col=0).dropna(subset=["TISSUE_matched"])
    S = S.sort_values("TISSUE_matched")
    SHORT = {"FOLI_mPS3": "Foli-seq mPS3", "FOLI_hPS2": "Foli-seq hPS2", "COSMX_HIO": "CosMx IO (100)",
             "MER_M_BRAINTRAIN": "MERSCOPE Brain (128)", "MER_M_PAN": "MERSCOPE Pan-Mouse (815)",
             "MER_M_PANNEURO": "MERSCOPE PanNeuro (815)", "COSMX_MNEURO": "CosMx Mouse Neuro (1K)",
             "COSMX_MUCC": "CosMx Mouse UCC (1K)", "COSMX_HUCC": "CosMx Human UCC (1K)"}
    INK, MUTED, GRID, RED, BLUE = "#0b0b0b", "#52514e", "#e6e5e1", "#e34948", "#2a78d6"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 200, "savefig.bbox": "tight",
                         "figure.facecolor": "white"})
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    y = np.arange(len(S))
    ax.barh(y - .2, S.TISSUE_unmatched, height=.38, color="#b9b8b4",
            label="as originally measured (26 immune types vs 15 gut types)")
    ax.barh(y + .2, S.TISSUE_matched, height=.38, color=BLUE,
            label="cell-type count matched (15 vs 15)")
    ax.errorbar(S.TISSUE_matched, y + .2, xerr=S["sd"], fmt="none", ecolor=INK, elinewidth=.9, capsize=2)
    for i, p in enumerate(S.index):
        if p.startswith("FOLI"):
            ax.get_yticklabels()
    ax.axvline(0, color=MUTED, lw=1)
    ax.set_yticks(y)
    labs = ax.set_yticklabels([SHORT.get(p, p) for p in S.index], fontsize=8)
    for lab, p in zip(labs, S.index):
        if p.startswith("FOLI"):
            lab.set_color(RED); lab.set_fontweight("bold")
    ax.set_xlabel("tissue effect: lift on mouse gut epithelium minus lift on mouse immune (pp)", fontsize=8.5)
    ax.set_title("Controlling for the number of cell types\n"
                 "Direction survives for every panel; magnitude shrinks, and the Foli-seq panels fall to ~0",
                 fontsize=9.5, loc="left", color=INK)
    ax.legend(fontsize=7.8, frameon=False, loc="upper left")
    ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(PROJ + "figures/v3/v3_fig9_celltype_count_control.png"); plt.close(fig)
    print("wrote figures/v3/v3_fig9_celltype_count_control.png")


if __name__ == "__main__":
    if "--figure-only" not in sys.argv:
        main()
    figure()
