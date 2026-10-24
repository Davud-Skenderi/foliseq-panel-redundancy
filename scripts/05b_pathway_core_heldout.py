"""Step 5b: pathway redundancy with split-half validation.
Core genes are chosen greedily on a random half of samples (until the core-gene mean reaches r>=0.9 with the
full-membership score), then the same core is scored on the other half. Repeated for 3 splits.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import numpy as np, pandas as pd
D = PROJ + "data/"; O = PROJ + "results/"
gs = pd.read_pickle(D + "genesets.pkl")
lens = gs[(gs.collection == "Hallmark") | gs.gut_lens]
rng = np.random.default_rng(5)
rows = []
for ds in ["HS_PBMC", "MM_SLN", "HS_HEART", "HS_TISSUE"]:
    e = pd.read_parquet(D + f"expr_{ds}.parquet")
    v = e.values.astype(float); v -= v.mean(1, keepdims=True); sd = v.std(1, keepdims=True); sd[sd == 0] = 1
    z = pd.DataFrame(v / sd, index=e.index)
    n = z.shape[1]
    for name, coll, genes in zip(lens.gs_name, lens.collection, lens.genes):
        mem = [g for g in genes if g in z.index]
        if len(mem) < 10:
            continue
        X = z.loc[mem].values
        for rep in range(3):
            tr = rng.random(n) < 0.5; te = ~tr
            full = X.mean(0)
            chosen, s, best = [], np.zeros(n), 0.0
            while best < 0.9 and len(chosen) < len(mem):
                cand = [i for i in range(len(mem)) if i not in chosen]
                k = len(chosen)
                rs = np.array([np.corrcoef(((s * k + X[i]) / (k + 1))[tr], full[tr])[0, 1] for i in cand])
                j = cand[int(np.nanargmax(rs))]; best = np.nanmax(rs)
                s = (s * k + X[j]) / (k + 1); chosen.append(j)
            r_te = np.corrcoef(s[te], full[te])[0, 1]
            rows.append(dict(dataset=ds, gene_set=name, collection=coll, n_members_measured=len(mem), rep=rep,
                             core_size=len(chosen), heldout_r=r_te, core_frac=len(chosen) / len(mem),
                             core_genes=";".join(mem[i] for i in chosen)))
    print(ds, "done", flush=True)
R = pd.DataFrame(rows)
R.to_csv(O + "pathway_core_heldout.csv", index=False)
S = R.groupby(["gene_set", "collection"]).agg(n_members=("n_members_measured", "median"), core_size=("core_size", "median"),
                                             heldout_r=("heldout_r", "median"), core_frac=("core_frac", "median")).reset_index()
S.to_csv(O + "pathway_core_summary.csv", index=False)
print(S.describe().round(3)); print(S.groupby("collection")[["core_size", "heldout_r", "core_frac"]].median().round(3))
