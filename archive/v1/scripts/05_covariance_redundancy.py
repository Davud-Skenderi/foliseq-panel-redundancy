"""Step 5: gene covariance + redundancy of each panel in public expression data.

For each dataset (metacells / tissue profiles) and each expression panel:
  n_measured            panel genes detected in the dataset
  frac_pair_redundant   fraction of panel genes with a partner on the same panel at |r| >= 0.8
  eff_dim_frac          participation ratio of the panel's correlation spectrum / n_measured
                        (1 = every gene carries independent signal; ->0 = genes collapse onto few axes)
  self_pred_R2_median   median held-out R^2 of predicting each panel gene from the *other* panel genes
  frac_self_predictable fraction of panel genes with held-out R^2 >= 0.7 (i.e., recoverable from the rest)
  tx_recon_R2           mean held-out R^2 for 1,000 off-panel genes predicted from the panel (information captured)
  ...and the same three for size-matched random gene sets (baseline).
Pathway level (Hallmark + gut/IBD lens): fidelity of panel-restricted pathway score vs full score,
minimum member genes needed to reproduce the full score (r >= 0.9) -> pathway redundancy index.
Also writes consensus gene-gene correlations (Fisher-z averaged over datasets) for the targeted-panel gene union.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import sys
sys.path.insert(0, PROJ + "scripts")
from common import top_pcs
import numpy as np, pandas as pd, warnings
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.utils.extmath import randomized_svd
import functools
print = functools.partial(print, flush=True)
warnings.filterwarnings("ignore")
D = PROJ + "data/"
O = PROJ + "results/"
rng = np.random.default_rng(1)
DATASETS = ["HS_PBMC", "MM_SLN", "HS_HEART", "HS_TISSUE"]
E = {k: pd.read_parquet(D + f"expr_{k}.parquet") for k in DATASETS}
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
meta = pd.read_csv(D + "panel_meta.csv")
gs = pd.read_pickle(D + "genesets.pkl")
EXPR_PANELS = meta[(meta.modality == "RNA") & (meta.n_targets > 0)].panel_id.tolist()
PANELS = {p: set(M.index[M[p] == 1]) for p in EXPR_PANELS}
FOLI = [p for p in PANELS if p.startswith("FOLI_")]


def zscore(df):
    v = df.values.astype(np.float64)
    v = v - v.mean(1, keepdims=True)
    sd = v.std(1, keepdims=True); sd[sd == 0] = 1
    return pd.DataFrame(v / sd, index=df.index, columns=df.columns)


Z = {k: zscore(e.loc[e.var(1) > 1e-6]) for k, e in E.items()}


def cv_predict_r2(Xp, Y, n_pc=40, n_splits=5, alpha=1.0):
    """held-out R^2 of each column of Y (samples x targets) predicted from PCs of Xp (samples x predictors)."""
    n = Xp.shape[0]
    pred = np.zeros_like(Y)
    for tr, te in KFold(n_splits, shuffle=True, random_state=0).split(Xp):
        mu = Xp[tr].mean(0); Xc = Xp[tr] - mu
        k = int(min(n_pc, Xc.shape[1], len(tr) - 2))
        Vt = top_pcs(Xc, k)
        P_tr = Xc @ Vt[:k].T; P_te = (Xp[te] - mu) @ Vt[:k].T
        m = Ridge(alpha=alpha).fit(P_tr, Y[tr]); pred[te] = m.predict(P_te)
    ss_res = ((Y - pred) ** 2).sum(0); ss_tot = ((Y - Y.mean(0)) ** 2).sum(0)
    return 1 - ss_res / np.where(ss_tot == 0, 1, ss_tot)


def panel_metrics(genes, z, n_pc, n_splits, off_targets):
    g = [x for x in genes if x in z.index]
    n = len(g)
    if n < 5:
        return dict(n_measured=n)
    Xg = z.loc[g].values.T  # samples x genes (rows already z-scored per gene)
    Xs = (Xg - Xg.mean(0)) / np.where(Xg.std(0) == 0, 1, Xg.std(0))
    ns = Xs.shape[0]
    maxr = np.zeros(n)
    for a in range(0, n, 2000):  # chunked max |r| to a partner gene
        c = (Xs[:, a:a + 2000].T @ Xs) / ns
        for i in range(c.shape[0]):
            c[i, a + i] = 0
        maxr[a:a + 2000] = np.abs(c).max(1)
    sv = np.linalg.svd(Xs, compute_uv=False)
    lam = sv ** 2 / ns  # non-zero eigenvalues of the correlation matrix
    pr = lam.sum() ** 2 / (lam ** 2).sum()
    # self-predictability with gene folds: genes in fold f predicted from PCs of genes outside fold f
    # cap self-prediction targets at 500 genes (random subset) to keep whole-transcriptome panels tractable
    tgt = rng.permutation(n)[:500]
    folds = np.array_split(tgt, min(5, len(tgt)))
    r2 = np.full(n, np.nan)
    for f in folds:
        rest = np.setdiff1d(np.arange(n), f)
        r2[f] = cv_predict_r2(Xg[:, rest], Xg[:, f], n_pc=n_pc, n_splits=n_splits)
    gset = set(g)
    off = [x for x in off_targets if x not in gset][:1000]
    tx = cv_predict_r2(Xg, z.loc[off].values.T, n_pc=n_pc, n_splits=n_splits)
    return dict(n_measured=n, frac_pair_redundant=(maxr >= 0.8).mean(), mean_max_abs_r=maxr.mean(),
                eff_dim=pr, eff_dim_frac=pr / n, self_pred_R2_median=np.nanmedian(r2),
                frac_self_predictable=np.mean(r2[~np.isnan(r2)] >= 0.7), tx_recon_R2=np.mean(np.clip(tx, -1, 1)),
                _gene_r2=dict(zip(g, r2)), _gene_maxr=dict(zip(g, maxr)))


rows, gene_rows = [], []
import os
_prev = O + "panel_redundancy_by_dataset.csv"
DONE = set()
if os.path.exists(_prev):  # reuse vendor-panel results; recompute only panels not yet scored (e.g. new Foli-seq lists)
    _p = pd.read_csv(_prev); _p = _p[_p.panel_id.isin(PANELS) & ~_p.panel_id.str.startswith("FOLI_")]
    rows = _p.to_dict("records"); DONE = set(_p.panel_id)
    _g = pd.read_csv(O + "gene_redundancy_by_panel.csv.gz"); _g = _g[_g.panel_id.isin(DONE)]
    gene_rows = list(_g.itertuples(index=False, name=None))
for ds, z in Z.items():
    n_pc, n_splits = (40, 5) if z.shape[1] > 100 else (15, 8)
    off_targets = list(rng.permutation(z.index))
    for p, genes in PANELS.items():
        if p in DONE:
            continue
        m = panel_metrics(genes, z, n_pc, n_splits, off_targets)
        gr, gm = m.pop("_gene_r2", {}), m.pop("_gene_maxr", {})
        for g_ in gr:
            gene_rows.append((ds, p, g_, gr[g_], gm[g_]))
        # size-matched random baseline (3 draws)
        if m.get("n_measured", 0) >= 5:
            base = [panel_metrics(list(rng.choice(z.index, m["n_measured"], replace=False)), z, n_pc, n_splits, off_targets)
                    for _ in range(2)]
            for key in ["frac_pair_redundant", "eff_dim_frac", "self_pred_R2_median", "frac_self_predictable", "tx_recon_R2"]:
                m["random_" + key] = np.mean([b[key] for b in base])
        rows.append(dict(dataset=ds, panel_id=p, **m))
        print(ds, p, {k: round(v, 3) if isinstance(v, float) else v for k, v in m.items() if k in ("n_measured", "eff_dim_frac", "frac_self_predictable", "tx_recon_R2", "random_tx_recon_R2")})
pd.DataFrame(rows).to_csv(O + "panel_redundancy_by_dataset.csv", index=False)
pd.DataFrame(gene_rows, columns=["dataset", "panel_id", "gene", "self_pred_R2", "max_abs_r_on_panel"]).to_csv(
    O + "gene_redundancy_by_panel.csv.gz", index=False)

# ------------------------------------------------------------------ pathway-level redundancy
lens = gs[(gs.collection == "Hallmark") | gs.gut_lens]
prow = []
pscores = {}
for ds, z in Z.items():
    S = {}
    for name, genes in zip(lens.gs_name, lens.genes):
        mem = [g for g in genes if g in z.index]
        if len(mem) < 10:
            continue
        X = z.loc[mem].values
        full = X.mean(0)
        S[name] = full
        # greedy: fewest member genes whose mean reproduces the full score at r >= 0.9
        chosen, cur, best_r = [], np.zeros_like(full), 0
        while best_r < 0.9 and len(chosen) < len(mem):
            cand = [i for i in range(len(mem)) if i not in chosen]
            rs = [np.corrcoef((cur * len(chosen) + X[i]) / (len(chosen) + 1), full)[0, 1] for i in cand]
            j = cand[int(np.nanargmax(rs))]; best_r = np.nanmax(rs)
            cur = (cur * len(chosen) + X[j]) / (len(chosen) + 1); chosen.append(j)
        rec = dict(dataset=ds, gene_set=name, n_members_measured=len(mem), min_genes_r90=len(chosen),
                   redundancy_index=1 - len(chosen) / len(mem), mean_intra_r=(np.corrcoef(X)[np.triu_indices(len(mem), 1)]).mean(),
                   core_genes=";".join(mem[i] for i in chosen))
        for p, pg in PANELS.items():
            on = [i for i, g in enumerate(mem) if g in pg]
            rec[f"fid_{p}"] = np.corrcoef(X[on].mean(0), full)[0, 1] if len(on) >= 2 else np.nan
            rec[f"n_{p}"] = len(on)
        prow.append(rec)
    pscores[ds] = pd.DataFrame(S)
pw = pd.DataFrame(prow)
pw.to_csv(O + "pathway_redundancy_fidelity.csv", index=False)
# pathway-pathway score correlation (biological redundancy between pathways), averaged over datasets
pc = []
for ds, S in pscores.items():
    c = S.corr(); c["dataset"] = ds; pc.append(c)
pcc = pd.concat(pc)
pcc_mean = pcc.drop(columns="dataset").groupby(level=0).mean()
pcc_mean.to_csv(O + "pathway_score_correlation_mean.csv")

# ------------------------------------------------------------------ consensus gene-gene correlation (targeted union)
targeted = [c for c in M.columns if c not in ("GEOMX_H_WTA", "GEOMX_M_WTA", "AMPLI_H_TX", "AMPLI_COMP_V3", "IDT_XGEN_PANCANCER")]
U = sorted(set(M.index[M[targeted].sum(axis=1) > 0]) | set().union(*[PANELS[p] for p in FOLI]))
idx = {g: i for i, g in enumerate(U)}
Ssum = np.zeros((len(U), len(U)), np.float32); Cnt = np.zeros((len(U), len(U)), np.uint8)
for ds, z in Z.items():
    g = [x for x in U if x in z.index]
    ii = np.array([idx[x] for x in g])
    c = np.corrcoef(z.loc[g].values).astype(np.float32)
    np.clip(c, -0.999, 0.999, out=c)
    Ssum[np.ix_(ii, ii)] += np.arctanh(c); Cnt[np.ix_(ii, ii)] += 1
    print("consensus add", ds, len(g))
with np.errstate(invalid="ignore", divide="ignore"):
    R = np.tanh(Ssum / Cnt)
np.save(D + "consensus_corr.npy", R.astype(np.float16)); np.save(D + "consensus_cnt.npy", Cnt)
pd.Series(U).to_csv(D + "consensus_genes.csv", index=False)
print("done")
