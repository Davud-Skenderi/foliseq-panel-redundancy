"""Q2. Gaps: which cell types and pathways does each panel cover?

Three readouts, each a plain count:
 (a) Intestinal epithelial cell types: published cell-type signatures from single-cell atlases, as distributed in
     MSigDB collection C8 — Busslinger et al. 2021 (Cell Reports; adult human duodenum, 14 epithelial types) and
     Gao et al. 2018 (Nature Cell Biology; human colon). For each signature: how many of its genes are on the panel.
 (b) Immune cell types: the top 20 marker genes of each annotated cell type in the mouse spleen/lymph-node reference
     (Gayoso et al. 2021), found with the Wilcoxon rank-sum test in scanpy (Wolf et al. 2018; the Wilcoxon test is
     among the best-performing single-cell differential-expression tests, Soneson & Robinson 2018, Nature Methods).
     Markers must have log2 fold-change > 1 and adjusted p < 0.05.
 (c) Hallmark pathways (Liberzon et al. 2015, Cell Systems): coverage (% of the pathway's genes on the panel) and
     over-representation — a one-sided hypergeometric test of whether the panel holds more of the pathway's genes
     than a random gene list of the same size would, with Benjamini–Hochberg false discovery rate correction.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT  # noqa: E402
import os, re, warnings
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from common import to_human
warnings.filterwarnings("ignore")
D = PROJ + "data/"; O = PROJ + "results/v2/"; os.makedirs(O, exist_ok=True)
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
gs = pd.read_pickle(D + "genesets.pkl")
PANELS = ["FOLI_mPS3", "FOLI_hPS2", "COSMX_HIO", "MER_M_PAN", "COSMX_MUCC", "COSMX_HUCC", "GEOMX_H_CTA",
          "XEN_M5K", "XEN_H5K", "COSMX_H6K"]
P = {p: set(M.index[M[p] == 1]) for p in PANELS}

# (a) epithelial signatures
EPI = gs[gs.gs_name.str.match(r"BUSSLINGER_DUODENAL_|GAO_LARGE_INTESTINE_ADULT_C[AEFGH]_|GAO_LARGE_INTESTINE_24W_C10_ENTEROCYTE")]
rows = []
for name, genes in zip(EPI.gs_name, EPI.genes):
    label = re.sub(r"^(BUSSLINGER_DUODENAL_|GAO_LARGE_INTESTINE_ADULT_C[A-Z]_|GAO_LARGE_INTESTINE_24W_C10_)", "", name).replace("_", " ").lower()
    src = "Busslinger 2021 (duodenum)" if name.startswith("BUSSLINGER") else "Gao 2018 (colon)"
    r = dict(source=src, cell_type=label, msigdb_set=name, n_signature_genes=len(genes))
    for p in PANELS:
        on = sorted(set(genes) & P[p]); r[f"{p}_n"] = len(on); r[f"{p}_genes"] = ", ".join(on[:25])
    rows.append(r)
pd.DataFrame(rows).to_csv(O + "q2a_epithelial_signature_coverage.csv", index=False)

# (b) immune markers from the annotated mouse reference
a1 = ad.read_h5ad(EXT + "expr/sln_111.h5ad"); a2 = ad.read_h5ad(EXT + "expr/sln_208.h5ad")
a = ad.concat([a1, a2], join="inner"); a.obs_names_make_unique()
a = a[~a.obs.cell_types.str.contains("doublet|Low quality|Erythrocytes", case=False)].copy()
vc = a.obs.cell_types.value_counts(); a = a[a.obs.cell_types.isin(vc.index[vc >= 40])].copy()
hs = [to_human(g, "Mouse") for g in a.var_names.astype(str)]
keep = np.array([isinstance(h, str) and h != "" for h in hs]); a = a[:, keep].copy()
a.var["human"] = np.array(hs, dtype=object)[keep]
sc.pp.filter_genes(a, min_cells=10); sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
a.obs["celltype"] = a.obs.cell_types.astype(str)
sc.tl.rank_genes_groups(a, "celltype", method="wilcoxon", n_genes=200)
rg = sc.get.rank_genes_groups_df(a, None)
rg = rg[(rg.logfoldchanges > 1) & (rg.pvals_adj < 0.05)]
rg["human"] = rg.names.map(dict(zip(a.var_names, a.var.human)))
top = rg.groupby("group").head(20)
top[["group", "names", "human", "logfoldchanges", "pvals_adj"]].rename(columns={"group": "celltype", "names": "mouse_gene"}).to_csv(O + "q2b_immune_markers_top20.csv", index=False)
rows = []
for ct, sub in top.groupby("group"):
    hm = set(sub.human)
    r = dict(cell_type=ct, n_markers=len(hm))
    for p in PANELS:
        on = sorted(hm & P[p]); r[f"{p}_n"] = len(on); r[f"{p}_genes"] = ", ".join(on)
    rows.append(r)
pd.DataFrame(rows).to_csv(O + "q2b_immune_marker_coverage.csv", index=False)

# (c) Hallmark coverage + over-representation
H = gs[gs.collection == "Hallmark"]
U = set(g for gl in gs.genes for g in gl)  # background: genes annotated in MSigDB collections used here
N = len(U)
rows = []
for p in PANELS:
    Pg = P[p] & U; n = len(Pg)
    for name, genes in zip(H.gs_name, H.genes):
        S = set(genes) & U; K = len(S); k = len(S & Pg)
        rows.append(dict(panel_id=p, pathway=name.replace("HALLMARK_", ""), pathway_size=K, genes_on_panel=k,
                         coverage_pct=100 * k / K, expected_by_chance=K * n / N,
                         p_value=hypergeom.sf(k - 1, N, K, n) if k else 1.0))
Hc = pd.DataFrame(rows)
Hc["q_value_BH"] = Hc.groupby("panel_id").p_value.transform(lambda s: multipletests(s, method="fdr_bh")[1])
Hc.to_csv(O + "q2c_hallmark_coverage.csv", index=False)
print(pd.DataFrame(rows).groupby("panel_id").coverage_pct.median().round(1))
