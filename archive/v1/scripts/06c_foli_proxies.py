"""Step 6c: can a spatial panel stand in for each Foli-seq gene? Direct target, or best co-expressed proxy
(consensus correlation, Fisher-z mean over public datasets, pairs seen in >=2 datasets)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import numpy as np, pandas as pd
D = PROJ + "data/"; O = PROJ + "results/"
R = np.load(D + "consensus_corr.npy").astype(np.float32); C = np.load(D + "consensus_cnt.npy")
U = pd.read_csv(D + "consensus_genes.csv").iloc[:, 0].tolist(); ix = {g: i for i, g in enumerate(U)}
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
panels = ["XEN_H5K", "XEN_M5K", "COSMX_H6K", "COSMX_HUCC", "COSMX_MUCC", "COSMX_HIO", "MER_M_PAN", "GEOMX_H_CTA"]
rows = []
FOLI = [c for c in M.columns if c.startswith("FOLI_")]
for fp in FOLI:
  foli = sorted(M.index[M[fp] == 1])
  for p in panels:
    on = [g for g in M.index[M[p] == 1] if g in ix]
    oi = np.array([ix[g] for g in on])
    for g in foli:
        rec = dict(foli_panel=fp, panel_id=p, foli_gene=g, direct=bool(M.loc[g, p]) if g in M.index else False)
        if g in ix and not rec["direct"]:
            r = R[ix[g], oi].copy(); r[C[ix[g], oi] < 2] = np.nan
            if np.isfinite(r).any():
                j = int(np.nanargmax(np.abs(r))); rec.update(best_proxy=on[j], proxy_r=float(r[j]))
        rows.append(rec)
P = pd.DataFrame(rows); P.to_csv(O + "foli_spatial_proxies.csv", index=False)
S = P.groupby(["foli_panel", "panel_id"]).agg(direct=("direct", "sum"), proxy_r06=("proxy_r", lambda s: int((s.abs() >= 0.6).sum())),
                              proxy_r04=("proxy_r", lambda s: int((s.abs() >= 0.4).sum())))
S["n_foli_genes"] = P.groupby(["foli_panel", "panel_id"]).size()
S["no_coverage"] = S.n_foli_genes - S.direct - S.proxy_r04
S.to_csv(O + "foli_spatial_proxy_summary.csv"); print(S)
