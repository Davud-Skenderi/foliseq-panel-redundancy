"""Step 6: redundancy-aware panel optimisation.

Objective: choose genes that best reconstruct (held-out) (a) the transcriptome's major axes of variation and
(b) Hallmark + gut/IBD-lens pathway scores, pooled across all public datasets (each dataset weighted equally).
Greedy forward selection with orthogonal deflation (column-subset selection), trained on 80% of samples,
evaluated on the held-out 20% with PC-ridge regression.

Also:
  * evaluates every vendor panel (and size-matched random sets) under the same held-out metric;
  * Foli-seq: backward 'unique contribution' ranking of the 65 recovered genes (which are redundant),
    and forward additions that most improve reconstruction (what Foli-seq should add).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import sys
sys.path.insert(0, PROJ + "scripts")
from common import top_pcs
import numpy as np, pandas as pd, warnings, json
from sklearn.linear_model import Ridge
from sklearn.utils.extmath import randomized_svd
import functools
print = functools.partial(print, flush=True)
warnings.filterwarnings("ignore")
D = PROJ + "data/"; O = PROJ + "results/"
rng = np.random.default_rng(7)
DATASETS = ["HS_PBMC", "MM_SLN", "HS_HEART", "HS_TISSUE"]
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
meta = pd.read_csv(D + "panel_meta.csv")
gs = pd.read_pickle(D + "genesets.pkl")
lens = gs[(gs.collection == "Hallmark") | gs.gut_lens]
EXPR = meta[(meta.modality == "RNA") & (meta.n_targets > 0)].panel_id.tolist()
PANELS = {p: set(M.index[M[p] == 1]) for p in EXPR}


def z(df):
    v = df.values.astype(np.float64); v -= v.mean(1, keepdims=True)
    sd = v.std(1, keepdims=True); sd[sd == 0] = 1
    return pd.DataFrame(v / sd, index=df.index, columns=df.columns)


Z = {k: z(pd.read_parquet(D + f"expr_{k}.parquet")) for k in DATASETS}
# candidate pool: genes measured in >= 3 of 4 datasets (missing -> 0 after z-scoring, i.e. uninformative)
cnt = pd.Series(0, index=sorted(set().union(*[set(v.index) for v in Z.values()])))
for v in Z.values():
    cnt[v.index] += 1
POOL = cnt.index[cnt >= 3].tolist()
print("candidate pool", len(POOL))

blocks = []
for ds, zz in Z.items():
    n = zz.shape[1]
    te = rng.random(n) < 0.2
    X = zz.reindex(POOL).fillna(0).values.T  # samples x pool
    # targets: top PCs of the full transcriptome (fit on train) + pathway scores
    Xall = zz.values.T
    mu = Xall[~te].mean(0)
    U, S, Vt = np.linalg.svd(Xall[~te] - mu, full_matrices=False)
    k = min(30, (~te).sum() - 2)
    pcs = (Xall - mu) @ Vt[:k].T
    pcs = pcs / np.sqrt((pcs[~te] ** 2).sum(0).sum())  # weight PCs by their variance; block total = 1
    path = []
    for genes in lens.genes:
        mem = [g for g in genes if g in zz.index]
        if len(mem) >= 10:
            path.append(zz.loc[mem].values.mean(0))
    path = np.array(path).T
    path = (path - path[~te].mean(0)) / path[~te].std(0)
    path = path / np.sqrt((path[~te] ** 2).sum())  # block total = 1
    Y = np.hstack([pcs, path])
    blocks.append(dict(ds=ds, X=X, Y=Y, te=te, n_pc=k, n_path=path.shape[1]))


def heldout_R2(genes_idx, blocks, n_pc=40):
    """mean over datasets of held-out variance-weighted R^2 of targets from selected genes (PC-ridge)."""
    out = []
    for b in blocks:
        X, Y, te = b["X"][:, genes_idx], b["Y"], b["te"]
        keep = X[~te].std(0) > 0
        X = X[:, keep]
        if X.shape[1] == 0:
            out.append(0.0); continue
        mu = X[~te].mean(0); Xc = X[~te] - mu
        kk = int(min(n_pc, Xc.shape[1], Xc.shape[0] - 2))
        Vt = top_pcs(Xc, kk)
        m = Ridge(alpha=1.0).fit(Xc @ Vt[:kk].T, Y[~te])
        pred = m.predict((X[te] - mu) @ Vt[:kk].T)
        out.append(1 - ((Y[te] - pred) ** 2).sum() / ((Y[te] - Y[~te].mean(0)) ** 2).sum())
    return float(np.mean(out)), out


def greedy(blocks, n_select, start=(), allowed=None):
    """forward selection on training rows; returns ordered list of pool indices."""
    sel = list(start)
    Xr = [b["X"][~b["te"]].copy() for b in blocks]
    Yr = [b["Y"][~b["te"]].copy() for b in blocks]
    for b in range(len(blocks)):  # centre
        Xr[b] -= Xr[b].mean(0); Yr[b] -= Yr[b].mean(0)

    def deflate(j):
        for b in range(len(blocks)):
            q = Xr[b][:, j].copy(); nq = np.linalg.norm(q)
            if nq < 1e-8:
                continue
            q /= nq
            Xr[b] -= np.outer(q, q @ Xr[b]); Yr[b] -= np.outer(q, q @ Yr[b])
    for j in sel:
        deflate(j)
    mask = np.ones(len(POOL), bool) if allowed is None else np.isin(np.arange(len(POOL)), allowed)
    mask[sel] = False
    while len(sel) < len(start) + n_select:
        score = np.zeros(len(POOL))
        for b in range(len(blocks)):
            num = ((Yr[b].T @ Xr[b]) ** 2).sum(0)
            den = (Xr[b] ** 2).sum(0)
            score += np.where(den > 1e-8, num / np.maximum(den, 1e-8), 0)
        score[~mask] = -1
        j = int(np.argmax(score)); sel.append(j); mask[j] = False; deflate(j)
    return sel


gi = {g: i for i, g in enumerate(POOL)}
# 1) de novo optimal panel, up to 500 genes
opt = greedy(blocks, 500)
sizes = [10, 25, 50, 75, 100, 150, 200, 300, 400, 500]
curve = []
for s in sizes:
    r, per = heldout_R2(opt[:s], blocks)
    rr = [heldout_R2(list(rng.choice(len(POOL), s, replace=False)), blocks)[0] for _ in range(3)]
    curve.append(dict(size=s, optimized_R2=r, random_R2=np.mean(rr), **{f"opt_{b['ds']}": v for b, v in zip(blocks, per)}))
    print("size", s, round(r, 3), round(np.mean(rr), 3))
pd.DataFrame(curve).to_csv(O + "optimization_curve.csv", index=False)
pd.DataFrame(dict(rank=np.arange(1, 501), gene=[POOL[i] for i in opt])).to_csv(O + "optimized_panel_500.csv", index=False)

# 2) vendor panels under the same metric (genes in pool)
FINE = list(range(5, 101, 5)) + list(range(110, 501, 10))
FINE_R2 = {s_: heldout_R2(opt[:s_], blocks)[0] for s_ in FINE}
pd.Series(FINE_R2, name="optimized_R2").rename_axis("size").to_csv(O + "optimization_curve_fine.csv")
vp = []
for p, genes in PANELS.items():
    idx = [gi[g] for g in genes if g in gi]
    if len(idx) < 5:
        continue
    r, per = heldout_R2(idx, blocks)
    rr = np.mean([heldout_R2(list(rng.choice(len(POOL), len(idx), replace=False)), blocks)[0] for _ in range(3)])
    # size needed by the optimised ordering to match this panel
    match = next((s_ for s_ in FINE if FINE_R2[s_] >= r), None)
    vp.append(dict(panel_id=p, n_in_pool=len(idx), heldout_R2=r, random_same_size_R2=rr,
                   optimized_genes_to_match=match, **{f"R2_{b['ds']}": v for b, v in zip(blocks, per)}))
    print(p, len(idx), round(r, 3), round(rr, 3), match)
pd.DataFrame(vp).to_csv(O + "vendor_panel_reconstruction.csv", index=False)

# 3) Foli-seq: backward unique contribution + forward additions
import os
FOLI_ID = os.environ.get("FOLI_PANEL", "FOLI_mPS3")
foli = [gi[g] for g in PANELS[FOLI_ID] if g in gi]
base, _ = heldout_R2(foli, blocks)
uniq = []
for j in foli:
    r, _ = heldout_R2([x for x in foli if x != j], blocks)
    uniq.append(dict(gene=POOL[j], R2_without=r, unique_contribution=base - r))
pd.DataFrame(uniq).sort_values("unique_contribution").to_csv(O + "foli_gene_unique_contribution.csv", index=False)
add = greedy(blocks, 100, start=foli)[len(foli):]
fc = []
for s in [0, 10, 25, 50, 100]:
    fc.append(dict(n_added=s, R2=heldout_R2(foli + add[:s], blocks)[0]))
pd.DataFrame(fc).to_csv(O + "foli_augmentation_curve.csv", index=False)
pd.DataFrame(dict(rank=np.arange(1, 101), gene=[POOL[i] for i in add])).to_csv(O + "foli_recommended_additions.csv", index=False)
json.dump(dict(pool=len(POOL), foli_in_pool=len(foli), foli_R2=base), open(O + "opt_summary.json", "w"))
print("foli base", base, fc)
