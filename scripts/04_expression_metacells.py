"""Step 4: load public expression data, harmonise to human symbols, build metacells for covariance estimation.

Public datasets reachable from this workspace (GitHub-hosted):
  HS_PBMC   10x PBMC 10k (CITE-seq v3), human, 6,855 cells          [YosefLab/scVI-data]
  MM_SLN    mouse spleen + lymph node CITE-seq (Gayoso 2021), 32,648 cells, 2 runs  [YosefLab/scVI-data]
  HS_HEART  Human Heart Cell Atlas subsample (Litvinukova 2020), 20,000 cells/nuclei [YosefLab/scVI-data]
  HS_TISSUE HPA (35 tissues) + GTEx (29 tissues) median TPM, 64 tissue profiles incl. colon/SI/rectum/appendix [bioc/TissueEnrich]

Single-cell data are aggregated into metacells (k-means in PCA space, ~20 cells each) so that gene-gene
correlations are not dominated by dropout. Output: one log1p(CPM) gene x metacell matrix per dataset.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import sys, warnings
sys.path.insert(0, PROJ + "scripts")
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp, rdata
from sklearn.cluster import MiniBatchKMeans
from common import canon_h, to_human
warnings.filterwarnings("ignore")
X = EXT + "expr/"
D = PROJ + "data/"
rng = np.random.default_rng(0)


def to_symbols(a, species):
    names = a.var_names.astype(str)
    hs = [to_human(g, species) for g in names]
    keep = np.array([isinstance(h, str) and len(h) > 0 for h in hs])
    a = a[:, keep].copy()
    a.var["human"] = np.array(hs, dtype=object)[keep]
    # collapse duplicates by summing counts
    if a.var.human.duplicated().any():
        codes, uniq = pd.factorize(a.var.human)
        Mx = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
        Xs = sp.csr_matrix(a.X) @ Mx
        a = ad.AnnData(Xs, obs=a.obs, var=pd.DataFrame(index=uniq))
    else:
        a.var_names = a.var.human.values
    return a


def metacells(a, name, cells_per_mc=20, label_col=None):
    a.X = sp.csr_matrix(a.X)
    sc.pp.filter_cells(a, min_genes=200)
    b = a.copy()
    sc.pp.normalize_total(b, 1e4); sc.pp.log1p(b)
    sc.pp.highly_variable_genes(b, n_top_genes=3000)
    b = b[:, b.var.highly_variable].copy()
    sc.pp.scale(b, max_value=10); sc.tl.pca(b, n_comps=40)
    k = max(50, a.n_obs // cells_per_mc)
    lab = MiniBatchKMeans(k, random_state=0, n_init=3, batch_size=4096).fit_predict(b.obsm["X_pca"])
    G = sp.csr_matrix((np.ones(a.n_obs), (lab, np.arange(a.n_obs))), shape=(k, a.n_obs))
    S = np.asarray((G @ a.X).todense())  # summed raw counts per metacell
    S = S[S.sum(1) > 0]
    cpm = S / S.sum(1, keepdims=True) * 1e6
    E = pd.DataFrame(np.log1p(cpm / 10).T, index=a.var_names, columns=[f"{name}_mc{i}" for i in range(S.shape[0])])
    # keep genes detected in >=1% of cells (metacell correlations of near-absent genes are noise)
    det = np.asarray((a.X > 0).mean(0)).ravel()
    E = E.loc[det >= 0.01]
    mc_label = None
    if label_col:
        mc_label = pd.Series(a.obs[label_col].values).groupby(lab).agg(lambda s: s.value_counts().index[0])
    print(f"{name}: {a.n_obs} cells -> {E.shape[1]} metacells, {E.shape[0]} genes detected")
    return E, mc_label


out = {}
labels = {}
a = ad.read_h5ad(X + "pbmc_10k_protein_v3.h5ad"); a = to_symbols(a, "Human")
out["HS_PBMC"], _ = metacells(a, "HS_PBMC")

a1 = ad.read_h5ad(X + "sln_111.h5ad"); a2 = ad.read_h5ad(X + "sln_208.h5ad")
a = ad.concat([a1, a2], join="inner"); a.obs_names_make_unique()
a = a[~a.obs.cell_types.str.contains("doublet|Low quality", case=False)].copy()
a = to_symbols(a, "Mouse")
out["MM_SLN"], labels["MM_SLN"] = metacells(a, "MM_SLN", label_col="cell_types")

a = ad.read_h5ad(X + "hca_subsampled_20k.h5ad")
a = a[~a.obs.cell_type.isin(["doublets", "NotAssigned"])].copy()
a = to_symbols(a, "Human")
out["HS_HEART"], labels["HS_HEART"] = metacells(a, "HS_HEART", label_col="cell_type")

te = rdata.conversion.convert(rdata.parser.parse_file(EXT + "chk_bioc_TissueEnrich/inst/extdata/combine-expression.rda"))["dataset"]
hm = te["humanGeneMapping"]
ens2sym = dict(zip(hm.iloc[:, 0], hm.iloc[:, 1])) if hm.shape[1] >= 2 else {}
frames = []
for key, tag in [("Protein-Atlas", "HPA"), ("GTEx-Combine", "GTEx")]:
    e = te[key]["expressionData"].copy()
    e.index = [canon_h(ens2sym.get(i)) for i in e.index]
    e = e[e.index.notna()].groupby(level=0).mean()
    e.columns = [f"{tag}_{c}" for c in e.columns]
    frames.append(e)
T = pd.concat(frames, axis=1, join="inner")
T = T.loc[(T > 1).sum(1) >= 3]  # expressed (>1 TPM) in at least 3 tissues
T = np.log1p(T)
# HPA and GTEx have different scales: centre each gene within source before pooling the 64 tissue profiles
for tag in ("HPA", "GTEx"):
    cols = [c for c in T.columns if c.startswith(tag)]
    T[cols] = T[cols].sub(T[cols].mean(1), axis=0)
out["HS_TISSUE"] = T
print("HS_TISSUE:", T.shape)

for k, v in out.items():
    v.astype("float32").to_parquet(D + f"expr_{k}.parquet")
pd.to_pickle(labels, D + "metacell_labels.pkl")
