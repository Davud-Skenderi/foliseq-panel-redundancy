"""Q5 + Q6. Stool relevance and explicit, rule-based recommendations for future Foli-seq panels.

Every call is a fixed rule over evidence columns that are each directly readable; there is no weighted score.

Evidence (per gene, human symbol space; mouse symbol kept):
  stool_evidence        Foli-seq stool DEG tables in Huang et al. 2025 (Tables S3–S7; mouse colitis models)
  n_vendors             number of vendors (of 4) whose general-purpose panels include the gene (Q1)
  immune_marker_of      immune cell types for which it is a top-20 marker (Q2b)
  epithelial_signatures intestinal epithelial signatures containing it (Q2a)
  sole_marker_for       cell types for which it is one of at most 2 top-20 markers on that Foli-seq panel (Q2b)
  partners_r>=0.8       other genes on the same Foli-seq panel correlated at |r| >= 0.8 in a reference (Q3)

Rules for genes on mPS3 / hPS2 (applied in order):
  1 Keep — stool evidence          in any stool DEG table
  2 Keep — sole marker             one of <= 2 markers of a cell type on that panel
  3 Review — redundant             has a partner at |r| >= 0.8 on the same panel (Spapros' redundancy threshold),
                                   and rules 1–2 do not apply
  4 Keep — default                 no evidence against it
Rules for genes on neither panel:
  A Add — cell type poorly identified   top-20 marker of an immune cell type the panel identifies with accuracy < 0.8 (Q4)
  B Add — field consensus               on all 4 vendors' general panels (Q1)
  C Consider — epithelial               in >= 2 intestinal epithelial signatures and on >= 1 vendor panel
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ  # noqa: E402
import os, json
import numpy as np, pandas as pd
from common import to_human, canon_h, m2h
D = PROJ + "data/"; O = PROJ + "results/v2/"
FOLI = {"FOLI_mPS3": "mPS3", "FOLI_hPS2": "hPS2"}

V = pd.read_csv(O + "q1_gene_vendor_matrix.csv", index_col=0)
epi = pd.read_csv(O + "q2a_epithelial_signature_coverage.csv")
gs = pd.read_pickle(D + "genesets.pkl").set_index("gs_name")
imm = pd.read_csv(O + "q2b_immune_markers_top20.csv")
acc = pd.read_csv(O + "q4_accuracy_by_celltype.csv", index_col=0)
pairs = pd.read_csv(O + "q3_correlated_pairs.csv")
E = pd.read_csv(D + "foli_stool_degs.csv")
E["human"] = [to_human(g, "Mouse") or g.upper() for g in E.mouse_gene]
stool = E.groupby("human").evidence.apply(lambda s: ";".join(sorted(set(";".join(s).split(";")))))

h2m = {}
for f in ["foli_panels/FOLI_mPS3__Mouse.csv", "foli_panels_extra/FOLI_mPS2__Mouse.csv", "foli_panels_extra/FOLI_mPS1__Mouse.csv"]:
    for g in pd.read_csv(D + f).gene:
        h = to_human(g, "Mouse")
        if h:
            h2m.setdefault(h, g)
for m_, h_ in m2h.items():
    if isinstance(h_, str):
        h2m.setdefault(canon_h(h_) or h_, m_)

epi_of = {}
for name, label in zip(epi.msigdb_set, epi.cell_type + " (" + epi.source.str.split().str[0] + ")"):
    for g in gs.loc[name, "genes"]:
        epi_of.setdefault(g, []).append(label)
imm_of = imm.groupby("human").celltype.apply(lambda s: sorted(set(s))).to_dict()

on = {k: set(V.index[V[f"on_{k}"] == 1]) for k in FOLI.values()}
sole = {k: {} for k in FOLI.values()}
for k in FOLI.values():
    for ct, sub in imm.groupby("celltype"):
        present = sorted(set(sub.human) & on[k])
        if 1 <= len(present) <= 2:
            for g in present:
                sole[k].setdefault(g, []).append(ct)
partners = {k: {} for k in FOLI.values()}
for pid, k in FOLI.items():
    for _, r in pairs[pairs.panel_id == pid].iterrows():
        for x, y in ((r.gene_a, r.gene_b), (r.gene_b, r.gene_a)):
            partners[k].setdefault(x, []).append(f"{y} ({r.r:.2f}, {r.dataset})")

poor = {pid: list(acc.index[acc[pid] < 0.8]) for pid in FOLI if pid in acc}
add_A = {}
for pid, k in FOLI.items():
    for ct in poor.get(pid, []):
        for g in imm[imm.celltype == ct].human:
            if g not in on[k]:
                add_A.setdefault(g, set()).add(f"{k}: {ct}")

genes = set(on["mPS3"]) | set(on["hPS2"]) | set(add_A) | set(V.index[V.n_vendors == 4]) | \
        set(g for g, l in epi_of.items() if len(l) >= 2 and g in V.index and V.loc[g, "n_vendors"] >= 1)
rows = []
for g in sorted(genes):
    nv = int(V.loc[g, "n_vendors"]) if g in V.index else 0
    rec = dict(gene=g, mouse_gene=h2m.get(g, ""), on_mPS3=g in on["mPS3"], on_hPS2=g in on["hPS2"], n_vendors=nv,
               stool_evidence=stool.get(g, ""), immune_marker_of="; ".join(imm_of.get(g, [])),
               epithelial_signatures="; ".join(sorted(set(epi_of.get(g, [])))),
               sole_marker_for="; ".join(f"{k}: {ct}" for k in FOLI.values() for ct in sole[k].get(g, [])),
               partners_r_ge_0_8="; ".join(f"{k}: {p}" for k in FOLI.values() for p in partners[k].get(g, [])[:4]),
               add_reason_poorly_identified="; ".join(sorted(add_A.get(g, []))))
    if rec["on_mPS3"] or rec["on_hPS2"]:
        if rec["stool_evidence"]:
            call = "Keep: responds in stool"
        elif rec["sole_marker_for"]:
            call = "Keep: sole marker of a cell type"
        elif rec["partners_r_ge_0_8"]:
            call = "Review: redundant (|r| >= 0.8 partner)"
        else:
            call = "Keep: default"
    elif rec["add_reason_poorly_identified"]:
        call = "Add: marker for a poorly identified cell type"
    elif nv == 4:
        call = "Add: field consensus (all 4 vendors)"
    else:
        call = "Consider: epithelial marker"
    rec["call"] = call
    rows.append(rec)
T = pd.DataFrame(rows)
order = ["Keep: responds in stool", "Keep: sole marker of a cell type", "Keep: default", "Review: redundant (|r| >= 0.8 partner)",
         "Add: marker for a poorly identified cell type", "Add: field consensus (all 4 vendors)", "Consider: epithelial marker"]
T["call"] = pd.Categorical(T.call, order, ordered=True)
T = T.sort_values(["call", "n_vendors", "gene"], ascending=[True, False, True])
T.to_csv(O + "q6_recommendations.csv", index=False)
summ = T.groupby("call", observed=True).size().to_dict()
summ = {str(k): int(v) for k, v in summ.items()}
summ["poorly_identified_celltypes"] = poor
json.dump(summ, open(O + "q6_summary.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
