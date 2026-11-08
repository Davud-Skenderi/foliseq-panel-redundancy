"""Step 6b: Foli-seq-focused optimisation.

Same machinery as step 6, but
  * targets = gut/IBD-lens pathway & intestinal cell-type scores only (no generic transcriptome axes),
  * candidates = members of those gene sets, minus ribosomal / mitochondrial / housekeeping-like transcripts,
so the selected genes are interpretable biology rather than generic variance carriers.
Outputs: foli_focus_* tables (unique contribution, additions with pathway membership, curves) and a
de novo gut/IBD-focused panel ranking.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import os, re, json
import numpy as np, pandas as pd
src = open(PROJ + "scripts/06_optimize_panel.py").read()
exec(src.split("gi = {g: i for i, g in enumerate(POOL)}")[0])  # loads data, POOL, heldout_R2, greedy

GL = gs[gs.gut_lens]
EXCL = re.compile(r"^(RPL|RPS|MRPL|MRPS|MT-|EEF1|TMSB|HIST|H[1-4]C)|^(MALAT1|NEAT1|ACTB|ACTG1|GAPDH|B2M|FTL|FTH1|TPT1|UBA52|UBB|UBC)$")
lens_genes = set(g for gl in GL.genes for g in gl)
ALLOWED = [i for i, g in enumerate(POOL) if g in lens_genes and not EXCL.search(g)]
print("lens-restricted candidates", len(ALLOWED))

rng = np.random.default_rng(11)
fblocks = []
for b, (ds, zz) in zip(blocks, Z.items()):
    te = b["te"]
    path, names = [], []
    for name, genes in zip(GL.gs_name, GL.genes):
        mem = [g for g in genes if g in zz.index]
        if len(mem) >= 10:
            path.append(zz.loc[mem].values.mean(0)); names.append(name)
    path = np.array(path).T
    path = (path - path[~te].mean(0)) / path[~te].std(0)
    fblocks.append(dict(ds=ds, X=b["X"], Y=path, te=te))
gi = {g: i for i, g in enumerate(POOL)}
FOLI_ID = os.environ.get("FOLI_PANEL", "FOLI_mPS3")
SUF = "_" + FOLI_ID
foli = [gi[g] for g in PANELS[FOLI_ID] if g in gi]

base, per = heldout_R2(foli, fblocks)
uniq = []
for j in foli:
    r, _ = heldout_R2([x for x in foli if x != j], fblocks)
    uniq.append(dict(gene=POOL[j], R2_without=r, unique_contribution=base - r))
U = pd.DataFrame(uniq).sort_values("unique_contribution")
U.to_csv(O + "foli_focus_unique_contribution" + SUF + ".csv", index=False)

add = greedy(fblocks, 100, start=foli, allowed=ALLOWED)[len(foli):]
memb = {}
for name, genes in zip(GL.gs_name, GL.genes):
    for g in genes:
        memb.setdefault(g, []).append(name)
A = pd.DataFrame(dict(rank=np.arange(1, 101), gene=[POOL[i] for i in add]))
A["n_gut_lens_sets"] = A.gene.map(lambda g: len(memb.get(g, [])))
A["example_sets"] = A.gene.map(lambda g: "; ".join(sorted(memb.get(g, []), key=len)[:4]))
A.to_csv(O + "foli_focus_additions" + SUF + ".csv", index=False)
curve = []
for s in [0, 5, 10, 25, 50, 100]:
    r, per = heldout_R2(foli + add[:s], fblocks)
    curve.append(dict(n_added=s, total_genes=len(foli) + s, R2=r, **{f"R2_{b['ds']}": v for b, v in zip(fblocks, per)}))
# random lens-gene additions as a baseline
for s in [10, 25, 50, 100]:
    rr = [heldout_R2(foli + list(rng.choice(ALLOWED, s, replace=False)), fblocks)[0] for _ in range(3)]
    curve.append(dict(n_added=s, total_genes=len(foli) + s, R2=np.mean(rr), baseline="random lens genes"))
pd.DataFrame(curve).to_csv(O + "foli_focus_augmentation_curve" + SUF + ".csv", index=False)

if os.environ.get("SKIP_DENOVO"):
    raise SystemExit(0)
# de novo gut/IBD-focused panel and its curve vs random lens genes; plus each vendor panel on this objective
opt = greedy(fblocks, 300, allowed=ALLOWED)
pd.DataFrame(dict(rank=np.arange(1, 301), gene=[POOL[i] for i in opt],
                  n_gut_lens_sets=[len(memb.get(POOL[i], [])) for i in opt])).to_csv(O + "gut_focused_panel_300.csv", index=False)
rows = []
for s in [10, 25, 50, 75, 100, 150, 200, 300]:
    rows.append(dict(size=s, optimized_R2=heldout_R2(opt[:s], fblocks)[0],
                     random_lens_R2=np.mean([heldout_R2(list(rng.choice(ALLOWED, s, replace=False)), fblocks)[0] for _ in range(3)])))
pd.DataFrame(rows).to_csv(O + "gut_focused_curve.csv", index=False)
vp = []
for p, genes in PANELS.items():
    idx = [gi[g] for g in genes if g in gi]
    if len(idx) >= 5:
        vp.append(dict(panel_id=p, n_in_pool=len(idx), gut_lens_R2=heldout_R2(idx, fblocks)[0]))
pd.DataFrame(vp).to_csv(O + "vendor_gut_lens_R2.csv", index=False)
json.dump(dict(foli_panel=FOLI_ID, foli_in_pool=len(foli), foli_R2=base, candidates=len(ALLOWED)), open(O + "foli_focus_summary" + SUF + ".json", "w"))
print(pd.DataFrame(curve).round(3).to_string()); print(pd.DataFrame(rows).round(3).to_string()); print(pd.DataFrame(vp).round(3).to_string())
print(U.round(4).to_string()); print(A.head(30).to_string())
