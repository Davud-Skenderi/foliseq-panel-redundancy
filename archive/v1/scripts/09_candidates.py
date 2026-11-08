"""Step 9: per-gene evidence table and keep / review / add recommendations for future Foli-seq panels.

Evidence per gene (human symbol space; mouse symbol kept for Foli-seq mouse panels):
  on_mPS1/2/3, on_hPS2      current Foli-seq panels (Huang et al. 2025 Tables S2, S8)
  n_spatial_panels, vendors  how many targeted spatial panels / vendors include it (consistency)
  stool_evidence             DEG in Foli-seq stool data in a colitis model (Tables S3-S7, mouse)
  gut_enriched               HPA or GTEx tissue-enriched / group-enriched in intestine (TissueEnrich)
  epithelial_programs        intestinal epithelial cell-type signatures it belongs to (MSigDB C8)
  n_pathway_cores            how many gut/IBD programs list it among their 5-10 core genes (step 5b)
  unique_contribution        drop-one held-out R2 loss on the gut/IBD objective (step 6b), per Foli panel
  best_partner, partner_r    most correlated other gene on the same Foli panel (consensus r, step 5)
  addition_rank              rank among covariance-chosen additions (step 6b)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import sys, re, glob, os
sys.path.insert(0, PROJ + "scripts")
import numpy as np, pandas as pd, rdata
from common import to_human, canon_h, m2h
D = PROJ + "data/"; O = PROJ + "results/"

M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
gs = pd.read_pickle(D + "genesets.pkl")
SPATIAL = ["XEN_H5K", "XEN_M5K", "COSMX_H6K", "COSMX_HUCC", "COSMX_MUCC", "COSMX_HIO", "COSMX_MNEURO",
           "MER_M_PAN", "MER_M_PANNEURO", "MER_M_BRAINTRAIN", "GEOMX_H_CTA", "GEOMX_DOG_CCA"]
VENDOR = {"XEN": "10x", "COSMX": "CosMx", "MER": "Vizgen", "GEOMX": "GeoMx"}

# Foli panels incl. mPS1/2 (not analysed as separate panels, but reported)
foli = {}
for fp in glob.glob(D + "foli_panels/*.csv") + glob.glob(D + "foli_panels_extra/*.csv"):
    pid, sp = os.path.basename(fp)[:-4].split("__")
    if pid == "FOLI_hPS1":
        continue
    g = pd.read_csv(fp).gene.astype(str)
    foli[pid] = {to_human(x, sp) or x.upper(): x for x in g}
h2m = {}
for pid in ("FOLI_mPS1", "FOLI_mPS2", "FOLI_mPS3"):
    for h, m in foli[pid].items():
        h2m.setdefault(h, m)
for m_, h_ in m2h.items():  # remaining genes: mouse ortholog (MGI/AGR/babelgene)
    if isinstance(h_, str):
        h2m.setdefault(canon_h(h_) or h_, m_)

# stool evidence (mouse) -> human
E = pd.read_csv(D + "foli_stool_degs.csv")
E["human"] = [to_human(g, "Mouse") or g.upper() for g in E.mouse_gene]
stool = E.groupby("human").agg(n_stool_evidence=("n_evidence", "max"), stool_evidence=("evidence", lambda s: ";".join(sorted(set(";".join(s).split(";"))))))

# gut-enriched genes (HPA + GTEx)
te = rdata.conversion.convert(rdata.parser.parse_file(EXT + "chk_bioc_TissueEnrich/inst/extdata/combine-expression.rda"))["dataset"]
hm = te["humanGeneMapping"]; e2s = dict(zip(hm.iloc[:, 0], hm.iloc[:, 1]))
GUT = re.compile(r"Colon|Small.Intestine|Duodenum|Rectum|Appendix|Stomach")
gut_en = {}
for key in ("Protein-Atlas", "GTEx-Combine"):
    t = te[key]["tissueSpecificGenes"]
    t = t[t.Tissue.astype(str).str.contains(GUT) & t.Group.isin(["Tissue-Enriched", "Group-Enriched", "Tissue-Enhanced"])]
    for ens, tis, grp in zip(t.Gene, t.Tissue, t.Group):
        s = canon_h(e2s.get(ens))
        if s:
            gut_en.setdefault(s, set()).add(f"{grp}:{tis}")

# epithelial cell-type programs
EPI = gs[(gs.collection == "C8 cell types") & gs.gs_name.str.contains(r"GAO_LARGE_INTESTINE|GAO_SMALL_INTESTINE|BUSSLINGER_DUODENAL|FETAL_INTESTINE_INTESTINAL_EPI|MAIN_FETAL_INTESTINAL_EPI")
         & ~gs.gs_name.str.contains("IMMUNE|MESENCHYMAL")]
epi = {}
for n, genes in zip(EPI.gs_name, EPI.genes):
    short = re.sub(r"^(GAO_LARGE_INTESTINE_|GAO_SMALL_INTESTINE_|BUSSLINGER_DUODENAL_|DESCARTES_)", "", n)
    for g in genes:
        epi.setdefault(g, []).append(short)
lens_n = {}
for genes in gs[gs.gut_lens].genes:
    for g in genes:
        lens_n[g] = lens_n.get(g, 0) + 1

# pathway core membership (split-half validated, step 5b)
core = pd.read_csv(O + "pathway_core_heldout.csv")
core_n = core.assign(g=core.core_genes.str.split(";")).explode("g").groupby("g").gene_set.nunique()

# covariance: consensus correlation among genes on the same Foli panel
R = np.load(D + "consensus_corr.npy").astype(np.float32); Cn = np.load(D + "consensus_cnt.npy")
U = pd.read_csv(D + "consensus_genes.csv").iloc[:, 0].tolist(); ix = {g: i for i, g in enumerate(U)}

rows = []
cand = set(M.index[M[[c for c in M.columns if c.startswith("FOLI_")]].sum(axis=1) > 0])
cand |= set(M.index[M[SPATIAL].sum(axis=1) >= 6]) & set(lens_n)          # consistent spatial picks in gut/IBD biology
cand |= set(g for g in gut_en if g in epi)                                  # gut-enriched epithelial genes
for fp in ("FOLI_mPS3", "FOLI_hPS2"):
    f = O + f"foli_focus_additions_{fp}.csv"
    if os.path.exists(f):
        cand |= set(pd.read_csv(f).gene.head(50))
uc = {fp: dict(zip(*pd.read_csv(O + f"foli_focus_unique_contribution_{fp}.csv")[["gene", "unique_contribution"]].T.values))
      for fp in ("FOLI_mPS3", "FOLI_hPS2") if os.path.exists(O + f"foli_focus_unique_contribution_{fp}.csv")}
addr = {fp: {g: i + 1 for i, g in enumerate(pd.read_csv(O + f"foli_focus_additions_{fp}.csv").gene)}
        for fp in ("FOLI_mPS3", "FOLI_hPS2") if os.path.exists(O + f"foli_focus_additions_{fp}.csv")}

for g in sorted(cand):
    on = {fp: g in foli[fp] for fp in ("FOLI_mPS1", "FOLI_mPS2", "FOLI_mPS3", "FOLI_hPS2")}
    sp = [p for p in SPATIAL if g in M.index and M.loc[g, p] == 1]
    rec = dict(human_symbol=g, mouse_symbol=h2m.get(g, ""), on_mPS1=on["FOLI_mPS1"], on_mPS2=on["FOLI_mPS2"],
               on_mPS3=on["FOLI_mPS3"], on_hPS2=on["FOLI_hPS2"], n_spatial_panels=len(sp),
               n_vendors=len({VENDOR[p.split("_")[0]] for p in sp}),
               n_stool_evidence=int(stool.n_stool_evidence.get(g, 0)), stool_evidence=stool.stool_evidence.get(g, ""),
               gut_enriched="; ".join(sorted(gut_en.get(g, []))), n_epithelial_programs=len(epi.get(g, [])),
               epithelial_programs="; ".join(sorted(set(epi.get(g, [])))[:5]), n_gut_lens_sets=lens_n.get(g, 0),
               n_pathway_cores=int(core_n.get(g, 0)), measured_in_public_data=g in ix)
    for fp, key in (("FOLI_mPS3", "mPS3"), ("FOLI_hPS2", "hPS2")):
        rec[f"unique_contribution_{key}"] = uc.get(fp, {}).get(g, np.nan)
        rec[f"addition_rank_{key}"] = addr.get(fp, {}).get(g, np.nan)
        if on[fp] and g in ix:
            others = [x for x in foli[fp] if x != g and x in ix]
            oi = np.array([ix[x] for x in others])
            r = R[ix[g], oi].copy(); r[Cn[ix[g], oi] < 2] = np.nan
            if np.isfinite(r).any():
                j = int(np.nanargmax(r)); rec[f"partner_{key}"] = others[j]; rec[f"partner_r_{key}"] = float(r[j])
    rows.append(rec)
T = pd.DataFrame(rows)


def call(r):
    on = r.on_mPS3 or r.on_hPS2
    if on:
        if r.n_stool_evidence > 0:
            return "Keep: responds in Foli-seq stool data"
        ucs = [v for v in (r.unique_contribution_mPS3, r.unique_contribution_hPS2) if pd.notna(v)]
        pr = max([v for v in (r.get("partner_r_mPS3"), r.get("partner_r_hPS2")) if pd.notna(v)], default=np.nan)
        if not r.measured_in_public_data or not ucs:
            return "Keep: no evidence either way"
        if max(ucs) <= 0 and pd.notna(pr) and pr >= 0.6:
            return "Review: redundant with a panel partner"
        if max(ucs) > 0.002 or r.n_pathway_cores >= 3:
            return "Keep: unique or core signal"
        return "Keep: supporting"
    epi_marker = r.n_epithelial_programs > 0 and r.gut_enriched != ""
    covar = any(pd.notna(v) and v <= 25 for v in (r.addition_rank_mPS3, r.addition_rank_hPS2))
    if epi_marker and r.n_epithelial_programs >= 2:
        return "Add: gut epithelial marker"
    if covar:
        return "Add: covariance-ranked"
    if r.n_spatial_panels >= 8 and r.n_gut_lens_sets >= 3:
        return "Add: spatial consensus"
    if epi_marker:
        return "Consider: gut epithelial (single program)"
    return "Consider"


T["recommendation"] = T.apply(call, axis=1)
order = {"Keep: responds in Foli-seq stool data": 0, "Keep: unique or core signal": 1, "Keep: supporting": 2,
         "Keep: no evidence either way": 3, "Review: redundant with a panel partner": 4,
         "Add: gut epithelial marker": 5, "Add: covariance-ranked": 6, "Add: spatial consensus": 7,
         "Consider: gut epithelial (single program)": 8, "Consider": 9}
T["_o"] = T.recommendation.map(order)
T["_score"] = (T.n_stool_evidence * 3 + T.n_pathway_cores + T.n_epithelial_programs + T.n_spatial_panels / 4
               - T[["addition_rank_mPS3", "addition_rank_hPS2"]].min(axis=1).fillna(100) / 50)
T = T.sort_values(["_o", "_score"], ascending=[True, False]).drop(columns=["_o"]).rename(columns={"_score": "priority_score"})
T.to_csv(O + "foli_candidate_genes.csv", index=False)
print(T.recommendation.value_counts().to_string())
for k in ["Review: redundant with a panel partner", "Add: gut epithelial marker", "Add: covariance-ranked", "Add: spatial consensus"]:
    sub = T[T.recommendation == k]
    print("\n==", k, len(sub))
    print(sub.head(30)[["human_symbol", "mouse_symbol", "n_spatial_panels", "n_epithelial_programs", "gut_enriched", "n_pathway_cores",
                        "addition_rank_mPS3", "partner_mPS3" if "partner_mPS3" in sub else "human_symbol"]].to_string()[:4000])

# ---- epithelial cell-type group for each gene (majority vote over its epithelial programs)
GROUPS = [("goblet", r"GOBLET"), ("enteroendocrine", r"ENTEROENDOCRINE|EC_CELLS|K_CELLS|I_CELLS|MX_CELLS"),
          ("tuft", r"TUFT"), ("Paneth", r"PANETH"), ("stem / progenitor", r"STEM|LGR5|OLFM4|TRANSIT|PROGENITOR|MKI67"),
          ("enterocyte", r"ENTEROCYTE"), ("pan-epithelial", r"EPITHELIAL|BCHE")]


def group_of(progs):
    votes = {}
    for p in str(progs).split("; "):
        for name, pat in GROUPS:
            if re.search(pat, p):
                votes[name] = votes.get(name, 0) + 1
                break
    return max(votes, key=votes.get) if votes else ""


T = pd.read_csv(O + "foli_candidate_genes.csv")
T["epithelial_group"] = [group_of(epi_progs) if n else "" for epi_progs, n in
                         zip(["; ".join(sorted(set(epi.get(g, [])))) for g in T.human_symbol], T.n_epithelial_programs)]
T.to_csv(O + "foli_candidate_genes.csv", index=False)
add = T[T.recommendation == "Add: gut epithelial marker"]
short = (add.groupby("epithelial_group").apply(lambda d: pd.Series(dict(n=len(d), top=", ".join(d.human_symbol.head(8)))))
         .sort_values("n", ascending=False))
short.to_csv(O + "foli_epithelial_additions_by_group.csv")
print(short.to_string())
# current Foli-seq coverage of each epithelial group (genes on mPS3 / hPS2 belonging mostly to that group)
cur = T[T.on_mPS3 | T.on_hPS2]
print(cur[cur.epithelial_group != ""].groupby("epithelial_group").agg(on_mPS3=("on_mPS3", "sum"), on_hPS2=("on_hPS2", "sum"),
      examples=("human_symbol", lambda s: ", ".join(s.head(8)))).to_string())
