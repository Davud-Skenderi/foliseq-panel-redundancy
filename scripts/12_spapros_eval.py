"""Q3 + Q4. Panel evaluation with the published Spapros metrics (Kuemmerle et al. 2024, Nature Methods).

Q4  Cell-type identification: can the panel's genes alone tell the cell types apart?
    Spapros `xgboost_forest_classification`: gradient-boosted tree classifiers (XGBoost) trained on the panel genes,
    5-fold cross-validation, confusion matrix normalised per true cell type (balanced). Per-cell-type accuracy is the
    diagonal; Spapros summarises it as mean accuracy and % of cell types identified with accuracy > 0.8.
Q3  Redundancy: Spapros `correlation_matrix` = pairwise Pearson correlation of the panel genes on log-normalised
    single-cell data; summarised as % of panel genes whose strongest partner has |r| < 0.8 (non-redundant), plus the
    explicit list of gene pairs with |r| >= 0.8.

Reference data available now (gut references are added the same way once provided — see README):
  * MM_SLN: mouse spleen + lymph nodes, CITE-seq (Gayoso et al. 2021, Nature Methods), 26 annotated immune cell types.
    All panels are evaluated in human-ortholog gene space, so human panels are scored on mouse cells via orthologs.
  * HS_PBMC: human PBMC 10k (10x Genomics) — unannotated, so redundancy (Q3) only.
Deviations from Spapros defaults, for compute: 1 seed instead of 5 (5 CV folds kept); at most 200 cells per cell type;
correlations on 8,000 sampled cells. Finished panels are cached in results/v2/q4_confusion_*.csv.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT  # noqa: E402
import os, json, time, warnings
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
from common import to_human
from spapros.evaluation import metrics as spm
warnings.filterwarnings("ignore")
D = PROJ + "data/"; O = PROJ + "results/v2/"; os.makedirs(O, exist_ok=True)
rng = np.random.default_rng(0)
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
PANELS = ["FOLI_mPS3", "FOLI_hPS2", "COSMX_HIO", "MER_M_BRAINTRAIN", "MER_M_PAN", "MER_M_PANNEURO", "COSMX_MNEURO",
          "COSMX_MUCC", "COSMX_HUCC", "GEOMX_H_CTA", "GEOMX_DOG_CCA", "XEN_M5K", "XEN_H5K", "COSMX_H6K"]


def to_symbols(a, species):
    hs = [to_human(g, species) for g in a.var_names.astype(str)]
    keep = np.array([isinstance(h, str) and len(h) > 0 for h in hs])
    a = a[:, keep].copy()
    h = np.array(hs, dtype=object)[keep]
    codes, uniq = pd.factorize(h)
    Mx = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
    return ad.AnnData(sp.csr_matrix(a.X) @ Mx, obs=a.obs.copy(), var=pd.DataFrame(index=uniq))


def lognorm(a):
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    return a


def sample_per_type(a, key, n=200, min_cells=40):
    vc = a.obs[key].value_counts()
    keep_types = vc.index[vc >= min_cells]
    idx = []
    for t in keep_types:
        ii = np.where(a.obs[key].values == t)[0]
        idx += list(rng.choice(ii, min(n, len(ii)), replace=False))
    return a[sorted(idx)].copy()


def corr_summary(adata, genes, tag, panel):
    cm = spm.correlation_matrix(adata, var_names=genes, verbosity=0)
    v = np.array(cm.fillna(0).values, dtype=float, copy=True)
    np.fill_diagonal(v, 0)
    a = np.abs(v)
    names = list(cm.index)
    pairs = [(names[i], names[j], float(v[i, j])) for i, j in zip(*np.where(np.triu(a, 1) >= 0.8))]
    return dict(dataset=tag, panel_id=panel, n_genes_measured=len(genes),
                pct_genes_nonredundant=float(np.mean(a.max(axis=0) < 0.8)) * 100,
                mean_abs_corr=float(a[np.triu_indices(len(genes), 1)].mean()) if len(genes) > 1 else np.nan,
                n_pairs_r_ge_0_8=len(pairs)), pairs


# ---------------------------------------------------------------- mouse spleen / lymph node (annotated)
a1 = ad.read_h5ad(EXT + "expr/sln_111.h5ad"); a2 = ad.read_h5ad(EXT + "expr/sln_208.h5ad")
a = ad.concat([a1, a2], join="inner"); a.obs_names_make_unique()
a = a[~a.obs.cell_types.str.contains("doublet|Low quality|Erythrocytes", case=False)].copy()
a = to_symbols(a, "Mouse")
a.obs["celltype"] = a.obs.cell_types.astype(str)
sc.pp.filter_genes(a, min_cells=10)
a = lognorm(a)
clf_ad = sample_per_type(a, "celltype")
corr_ad = a[rng.choice(a.n_obs, 8000, replace=False)].copy()
print("MM_SLN:", a.shape, "classification set", clf_ad.shape, clf_ad.obs.celltype.nunique(), "cell types", flush=True)

sets = {p: [g for g in M.index[M[p] == 1] if g in set(a.var_names)] for p in PANELS}
for k in range(3):  # size-matched random baseline (Spapros benchmarks against random selections)
    sets[f"RANDOM_535_{k}"] = list(rng.choice(a.var_names, 535, replace=False))

rows, per_ct, red, allpairs = [], [], [], []
for p, genes in sets.items():
    t = time.time()
    cached = O + f"q4_confusion_{p}.csv"
    if os.path.exists(cached):  # resume: reuse a finished panel
        cm = pd.read_csv(cached, index_col=0)
    else:
        cm, cm_sd = spm.xgboost_forest_classification(clf_ad, genes, celltypes="all", ct_key="celltype", n_cells_min=40,
                                                      n_seeds=1, n_jobs=2, verbosity=0)[:2]
    diag = pd.Series(np.diag(cm.values), index=cm.index)
    rows.append(dict(panel_id=p, n_genes_measured=len(genes), mean_celltype_accuracy=float(diag.mean()),
                     pct_celltypes_acc_gt_0_8=float(np.mean(diag > 0.8)) * 100,
                     spapros_perct_acc_gt_0_8=float(spm.summary_metric_diagonal_confusion_percentage(cm, threshold=0.8))))
    per_ct += [dict(panel_id=p, celltype=c, accuracy=float(v)) for c, v in diag.items()]
    if not p.startswith(("RANDOM", "ALL")):
        s, pairs = corr_summary(corr_ad, genes, "MM_SLN", p)
        red.append(s); allpairs += [dict(dataset="MM_SLN", panel_id=p, gene_a=x, gene_b=y, r=r) for x, y, r in pairs]
    print(f"{p}: {len(genes)} genes, acc {diag.mean():.3f}, {time.time()-t:.0f}s", flush=True)
    cm.to_csv(O + f"q4_confusion_{p}.csv")

pd.DataFrame(rows).to_csv(O + "q4_celltype_identification.csv", index=False)
pd.DataFrame(per_ct).pivot(index="celltype", columns="panel_id", values="accuracy").to_csv(O + "q4_accuracy_by_celltype.csv")

# ---------------------------------------------------------------- human PBMC (unannotated): redundancy only
b = ad.read_h5ad(EXT + "expr/pbmc_10k_protein_v3.h5ad")
b = to_symbols(b, "Human"); sc.pp.filter_genes(b, min_cells=10); b = lognorm(b)
for p in PANELS:
    genes = [g for g in M.index[M[p] == 1] if g in set(b.var_names)]
    s, pairs = corr_summary(b, genes, "HS_PBMC", p)
    red.append(s); allpairs += [dict(dataset="HS_PBMC", panel_id=p, gene_a=x, gene_b=y, r=r) for x, y, r in pairs]
pd.DataFrame(red).to_csv(O + "q3_redundancy_summary.csv", index=False)
pd.DataFrame(allpairs).sort_values(["panel_id", "r"], ascending=[True, False]).to_csv(O + "q3_correlated_pairs.csv", index=False)
print(pd.DataFrame(rows).round(3).to_string()); print(pd.DataFrame(red).round(2).to_string())
