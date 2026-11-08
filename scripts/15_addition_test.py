"""Q4b. Does adding the vendor-consensus genes Foli-seq lacks improve immune cell-type identification?

The consensus genes (Q1: on all 4 vendors' general panels, not on mPS3) are chosen without looking at the reference
data, so this is a fair test. Control: mPS3 plus the same number of random genes (3 draws). Same Spapros XGBoost
classifier and cell sample as 12_spapros_eval.py.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT  # noqa: E402
import os, warnings
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
from spapros.evaluation import metrics as spm
import importlib.util
warnings.filterwarnings("ignore")
O = PROJ + "results/v2/"
spec = importlib.util.spec_from_file_location("ev", _os.path.join(_os.path.dirname(__file__), "12_spapros_eval.py"))
# reuse the data preparation of 12_spapros_eval.py without rerunning its loop
src = open(spec.origin).read().split("sets = {p:")[0]
ns = {"__file__": spec.origin}
exec(compile(src, spec.origin, "exec"), ns)
clf_ad, M, a = ns["clf_ad"], ns["M"], ns["a"]
rng = np.random.default_rng(1)
V = pd.read_csv(O + "q1_gene_vendor_matrix.csv", index_col=0)
var = set(clf_ad.var_names)
base = [g for g in M.index[M["FOLI_mPS3"] == 1] if g in var]
cons = [g for g in V.index[(V.n_vendors == 4) & (V.on_mPS3 == 0)] if g in var]
others = sorted(var - set(base))
sets = {"mPS3": base, "mPS3+consensus": base + cons}
for k in range(3):
    sets[f"mPS3+random_{k}"] = base + list(rng.choice(others, len(cons), replace=False))
rows, per = [], []
for name, genes in sets.items():
    cf = O + f"q4b_confusion_{name}.csv"
    if os.path.exists(cf):
        cm = pd.read_csv(cf, index_col=0)
    else:
        cm = spm.xgboost_forest_classification(clf_ad, genes, celltypes="all", ct_key="celltype", n_cells_min=40,
                                               n_seeds=1, n_jobs=2, verbosity=0)[0]
        cm.to_csv(cf)
    d = pd.Series(np.diag(cm.values), index=cm.index)
    rows.append(dict(set=name, n_genes=len(genes), n_added=len(genes) - len(base), mean_acc=d.mean(), pct_gt_0_8=100 * (d > 0.8).mean()))
    per += [dict(set=name, celltype=c, accuracy=v) for c, v in d.items()]
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(O + "q4b_addition_test.csv", index=False)
pd.DataFrame(per).pivot(index="celltype", columns="set", values="accuracy").to_csv(O + "q4b_addition_by_celltype.csv")
print("consensus genes added:", " ".join(cons))
