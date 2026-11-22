"""Q3 + Q4 recomputed on the gut references built by 16_build_gut_refs.py.

Same Spapros metrics as 12_spapros_eval.py, run against annotated human intestinal tissue instead of
immune tissue, so that the panel comparison can be tested for robustness to the choice of reference.

THREE INGREDIENTS, held apart deliberately:
  (1) panels      the gene lists being compared        -- fixed across runs
  (2) metrics     accuracy / redundancy definitions    -- fixed across runs
  (3) reference   the annotated single-cell data       -- VARIED, one run per reference
Varying only (3) is what makes the comparison a robustness test: absolute scores are expected to move
(a reference with more, finer cell types is harder, so the whole distribution shifts left), but the
ORDERING of panels should survive if it reflects the panels rather than the reference.

NOTE ON A BUG IN 12_spapros_eval.py: its confusion-matrix cache key is `q4_confusion_{panel}.csv`
with no dataset in the name, so running a second reference through that script would silently reload
the first reference's results. Cache keys here are `q4_confusion_{REF}__{panel}.csv`.

Outputs per reference REF:
  results/v2/q4_accuracy_by_celltype_{REF}.csv    cell type x panel accuracy (the histogram input)
  results/v2/q4_celltype_identification_{REF}.csv per-panel summary
  results/v2/q3_redundancy_summary_{REF}.csv      redundancy, with size-normalised rates
  results/v2/q3_correlated_pairs_{REF}.csv        pairs with |r| >= 0.8
"""
import os, sys, time, warnings
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ, EXT
from spapros.evaluation import metrics as spm
warnings.filterwarnings("ignore")

D, O = PROJ + "data/", PROJ + "results/v2/"
os.makedirs(O, exist_ok=True)
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
PANELS = ["FOLI_mPS3", "FOLI_hPS2", "COSMX_HIO", "MER_M_BRAINTRAIN", "MER_M_PAN", "MER_M_PANNEURO",
          "COSMX_MNEURO", "COSMX_MUCC", "COSMX_HUCC", "GEOMX_H_CTA", "GEOMX_DOG_CCA",
          "XEN_M5K", "XEN_H5K", "COSMX_H6K"]
RANDOM_SIZES = [330, 535]      # 330 = Foli-seq mPS3 measured size; 535 = baseline used in v2
N_RANDOM = 3


def lognorm(a):
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    return a


def sample_per_type(a, key, n=200, min_cells=40, seed=0):
    rng = np.random.default_rng(seed)
    vc = a.obs[key].value_counts()
    idx = []
    for t in vc.index[vc >= min_cells]:
        ii = np.where(a.obs[key].values == t)[0]
        idx += list(rng.choice(ii, min(n, len(ii)), replace=False))
    return a[sorted(idx)].copy()


def corr_summary(adata, genes, tag, panel):
    cm = spm.correlation_matrix(adata, var_names=genes, verbosity=0)
    v = np.array(cm.fillna(0).values, dtype=float, copy=True)
    np.fill_diagonal(v, 0)
    absv, names = np.abs(v), list(cm.index)
    pairs = [(names[i], names[j], float(v[i, j])) for i, j in zip(*np.where(np.triu(absv, 1) >= 0.8))]
    n = len(genes)
    n_possible = n * (n - 1) // 2
    return dict(
        dataset=tag, panel_id=panel, n_genes_measured=n,
        pct_genes_nonredundant=float(np.mean(absv.max(axis=0) < 0.8)) * 100,
        # the size-normalised quantities the raw count cannot provide
        pct_genes_redundant=float(np.mean(absv.max(axis=0) >= 0.8)) * 100,
        mean_abs_corr=float(absv[np.triu_indices(n, 1)].mean()) if n > 1 else np.nan,
        n_pairs_r_ge_0_8=len(pairs),
        n_possible_pairs=n_possible,
        redundant_pair_rate_per_1k=1000.0 * len(pairs) / n_possible if n_possible else np.nan,
    ), pairs


def run(ref_tag, h5ad_name):
    a = ad.read_h5ad(EXT + f"expr/{h5ad_name}")
    title = a.uns.get("reference_title", ref_tag)
    sc.pp.filter_genes(a, min_cells=10)
    a = lognorm(a)
    clf_ad = sample_per_type(a, "celltype")
    n_corr = min(8000, a.n_obs)
    corr_ad = a[np.random.default_rng(0).choice(a.n_obs, n_corr, replace=False)].copy()
    measured = set(a.var_names)
    print(f"\n=== {ref_tag} ({title}) {a.shape}, clf {clf_ad.shape}, "
          f"{clf_ad.obs.celltype.nunique()} cell types", flush=True)

    sets = {p: [g for g in M.index[M[p] == 1] if g in measured] for p in PANELS}
    rng = np.random.default_rng(0)
    for size in RANDOM_SIZES:
        for k in range(N_RANDOM):
            sets[f"RANDOM_{size}_{k}"] = list(rng.choice(a.var_names, size, replace=False))

    rows, per_ct, red, allpairs = [], [], [], []
    for p, genes in sets.items():
        if len(genes) < 10:
            print(f"  {p}: only {len(genes)} genes measured, skipped", flush=True); continue
        t0 = time.time()
        cached = O + f"q4_confusion_{ref_tag}__{p}.csv"       # dataset-keyed cache
        if os.path.exists(cached):
            cm = pd.read_csv(cached, index_col=0)
        else:
            cm = spm.xgboost_forest_classification(clf_ad, genes, celltypes="all", ct_key="celltype",
                                                   n_cells_min=40, n_seeds=1, n_jobs=2, verbosity=0)[0]
            cm.to_csv(cached)
        diag = pd.Series(np.diag(cm.values), index=cm.index)
        rows.append(dict(reference=ref_tag, panel_id=p, n_genes_measured=len(genes),
                         n_celltypes=len(diag),
                         mean_celltype_accuracy=float(diag.mean()),
                         median_celltype_accuracy=float(diag.median()),
                         sd_celltype_accuracy=float(diag.std()),
                         min_celltype_accuracy=float(diag.min()),
                         pct_celltypes_acc_gt_0_8=float(np.mean(diag > 0.8)) * 100))
        per_ct += [dict(panel_id=p, celltype=c, accuracy=float(v)) for c, v in diag.items()]
        if not p.startswith("RANDOM"):
            s, pairs = corr_summary(corr_ad, genes, ref_tag, p)
            red.append(s)
            allpairs += [dict(dataset=ref_tag, panel_id=p, gene_a=x, gene_b=y, r=r) for x, y, r in pairs]
        print(f"  {p}: {len(genes)} genes, mean acc {diag.mean():.3f}, {time.time()-t0:.0f}s", flush=True)

    pd.DataFrame(rows).to_csv(O + f"q4_celltype_identification_{ref_tag}.csv", index=False)
    (pd.DataFrame(per_ct).pivot(index="celltype", columns="panel_id", values="accuracy")
       .to_csv(O + f"q4_accuracy_by_celltype_{ref_tag}.csv"))
    if red:
        pd.DataFrame(red).to_csv(O + f"q3_redundancy_summary_{ref_tag}.csv", index=False)
        pd.DataFrame(allpairs).to_csv(O + f"q3_correlated_pairs_{ref_tag}.csv", index=False)
    print(f"=== {ref_tag} done", flush=True)


REFS = {"GUT_HEALTHY": "gut_healthy.h5ad", "GUT_IBD": "gut_ibd.h5ad", "GUT_CRC": "gut_crc.h5ad"}

if __name__ == "__main__":
    # usage: 17_multiref_eval.py REF [REF ...] [--panels A,B,C]
    # --panels restricts the run to a subset, so a shared panel set can be completed across every
    # reference before the expensive 5K/6K panels are attempted on any one of them. Results are
    # cached per (reference, panel), so later runs extend rather than repeat.
    args = sys.argv[1:]
    subset = None
    if "--panels" in args:
        i = args.index("--panels")
        subset = [x for x in args[i + 1].split(",") if x]
        args = args[:i]
    if subset:
        PANELS = [p for p in PANELS if p in subset]
        print(f"panel subset: {len(PANELS)} panels", flush=True)
    for tag in (args or list(REFS)):
        run(tag, REFS[tag])
