"""Step 3: what does each panel look at? Pathway coverage, enrichment, and panel-panel redundancy.

Outputs in results/:
  coverage_long.csv        panel x gene set: overlap, coverage, odds ratio, hypergeometric p, BH q
  panel_scope_summary.csv  per-panel breadth metrics
  panel_jaccard.csv        gene-level Jaccard between panels
  panel_pathway_similarity.csv  similarity of pathway-coverage profiles (biological-scope redundancy)
  foli_gene_platform_presence.csv  which platforms measure each Foli-seq gene
  top_enriched_per_panel.csv
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import numpy as np, pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

D = PROJ + "data/"
O = PROJ + "results/"
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
meta = pd.read_csv(D + "panel_meta.csv")
gs = pd.read_pickle(D + "genesets.pkl")

# universe = genes annotated in the kept collections (protein-coding dominated), plus panel genes
U = set(g for gl in gs.genes for g in gl)
M = M.loc[M.index.isin(U)]
N = len(U)
panels = [p for p in M.columns]
P = {p: set(M.index[M[p] == 1]) for p in panels}

rows = []
for p in panels:
    Pg = P[p]
    n = len(Pg)
    for name, coll, size, lens, genes in zip(gs.gs_name, gs.collection, gs["size"], gs.gut_lens, gs.genes):
        S = set(genes) & U
        k = len(S & Pg)
        K = len(S)
        # enrichment: is the panel over-represented for this pathway relative to the genome?
        pval = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
        exp = K * n / N
        orr = ((k + .5) / (n - k + .5)) / ((K - k + .5) / (N - K - n + k + .5))
        rows.append((p, name, coll, K, k, k / K, exp, orr, pval, lens))
cov = pd.DataFrame(rows, columns=["panel_id", "gene_set", "collection", "set_size", "overlap", "coverage",
                                  "expected_overlap", "odds_ratio", "p_enrich", "gut_lens"])
cov["q_enrich"] = np.nan
for (p, c), idx in cov.groupby(["panel_id", "collection"]).groups.items():
    cov.loc[idx, "q_enrich"] = multipletests(cov.loc[idx, "p_enrich"], method="fdr_bh")[1]
cov.to_csv(O + "coverage_long.csv.gz", index=False)

# per-panel scope summary
summ = []
for p in panels:
    c = cov[cov.panel_id == p]
    h = c[c.collection == "Hallmark"]
    r = c[c.collection.isin(["Reactome", "KEGG"])]
    gl = c[c.gut_lens]
    summ.append(dict(
        panel_id=p, n_genes_in_universe=len(P[p]),
        hallmark_mean_coverage=h.coverage.mean(),
        hallmarks_scorable=(h.overlap >= 5).sum(),  # >=5 genes -> a usable pathway score
        reactome_kegg_scorable_frac=(r.overlap >= 5).mean(),
        gut_lens_scorable_frac=(gl.overlap >= 5).mean(),
        gut_lens_mean_coverage=gl.coverage.mean(),
        n_enriched_pathways_q05=((c.q_enrich < .05) & c.collection.isin(["Hallmark", "KEGG", "Reactome"])).sum(),
        n_enriched_celltypes_q05=((c.q_enrich < .05) & (c.collection == "C8 cell types")).sum(),
        # coverage per gene: how many pathways each panel gene participates in (pathway 'reuse' of genes)
    ))
summ = pd.DataFrame(summ).merge(meta[["panel_id", "platform", "panel_name", "species", "modality"]], on="panel_id")
# genes-per-pathway efficiency: fraction of genome-wide hallmark coverage achieved per 1k genes
summ["hallmark_cov_per_1k_genes"] = summ.hallmark_mean_coverage / (summ.n_genes_in_universe / 1000)
summ.to_csv(O + "panel_scope_summary.csv", index=False)

# panel-panel gene Jaccard and overlap coefficient
J = pd.DataFrame(index=panels, columns=panels, dtype=float)
OC = J.copy()
for a in panels:
    for b in panels:
        i = len(P[a] & P[b])
        J.loc[a, b] = i / len(P[a] | P[b])
        OC.loc[a, b] = i / min(len(P[a]), len(P[b]))
J.to_csv(O + "panel_jaccard.csv")
OC.to_csv(O + "panel_overlap_coefficient.csv")

# biological-scope similarity: correlate log-odds enrichment profiles over Hallmark+KEGG+Reactome
# enrichment strength vector per panel: -log10 p (only over-representation counts), cosine similarity
prof = cov[cov.collection.isin(["Hallmark", "KEGG", "Reactome"])].pivot(index="gene_set", columns="panel_id", values="p_enrich")
prof = np.clip(-np.log10(prof), 0, 30)
nv = prof / np.sqrt((prof ** 2).sum())
S = nv.T @ nv
S.to_csv(O + "panel_pathway_similarity.csv")

# top enriched pathways per panel (what is the panel 'about')
top = (cov[cov.collection.isin(["Hallmark", "KEGG", "Reactome", "C8 cell types"]) & (cov.overlap >= 5)]
       .sort_values("p_enrich").groupby(["panel_id", "collection"]).head(8))
top.to_csv(O + "top_enriched_per_panel.csv", index=False)

# Foli-seq genes: which platforms can measure them?
long = pd.read_csv(D + "panels_long.csv")
FOLI = [p for p in M.columns if p.startswith("FOLI_")]
spatial = [c for c in M.columns if c.startswith(("XEN", "COSMX", "MER", "GEOMX")) and c not in ("GEOMX_H_WTA", "GEOMX_M_WTA", "COSMX_MNEURO_PROT")]
allpres = []
for fp in FOLI:
    foli = sorted(M.index[M[fp] == 1])
    pres = M.loc[foli, spatial].astype(int)
    pres["n_spatial_panels"] = pres.sum(axis=1)
    pres.insert(0, "foli_panel", fp)
    allpres.append(pres)
pd.concat(allpres).rename_axis("human_symbol").to_csv(O + "foli_gene_platform_presence.csv")
print(summ.sort_values("n_genes_in_universe").to_string())
print(pres[["XEN_H5K", "XEN_M5K", "COSMX_H6K", "COSMX_HUCC", "MER_M_PAN", "GEOMX_H_CTA", "n_spatial_panels"]].sum())
