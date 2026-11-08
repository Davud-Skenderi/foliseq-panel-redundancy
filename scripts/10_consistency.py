"""Q1. Consistency: which genes do the commercial vendors agree on, and how does Foli-seq relate to that consensus?

Method
  * All panels are compared in human gene symbols (step 01: HGNC-approved symbols; mouse/dog genes mapped to their
    human orthologs with MGI / Alliance of Genome Resources ortholog tables).
  * Only general-purpose RNA panels relevant to gut tissue are used (pan-tissue, universal cell characterisation,
    immuno-oncology, cancer transcriptome). Brain-specific, canine, whole-transcriptome and DNA panels are excluded,
    because they would distort a "what the field measures in tissue" consensus.
  * Panels are grouped by vendor, and a gene counts once per vendor, so a vendor with several panels does not dominate.
  * Set overlaps are shown as an UpSet plot (Lex et al. 2014, IEEE TVCG), the standard display for >3 sets.

Outputs (results/v2/)
  q1_gene_vendor_matrix.csv   gene x vendor presence (0/1), n_vendors, on_mPS3, on_hPS2
  q1_upset_counts.csv         size of every vendor/Foli-seq intersection
  q1_consensus_genes.csv      genes on all four vendors' general panels, annotated
  q1_summary.json
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ  # noqa: E402
import json, os
import pandas as pd

D = PROJ + "data/"
O = PROJ + "results/v2/"
os.makedirs(O, exist_ok=True)
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
gs = pd.read_pickle(D + "genesets.pkl")

VENDOR_PANELS = {
    "10x Xenium": ["XEN_H5K", "XEN_M5K"],
    "NanoString CosMx": ["COSMX_H6K", "COSMX_HUCC", "COSMX_MUCC", "COSMX_HIO"],
    "Vizgen MERSCOPE": ["MER_M_PAN"],
    "NanoString GeoMx": ["GEOMX_H_CTA"],
}
V = pd.DataFrame({v: M[p].max(axis=1) for v, p in VENDOR_PANELS.items()}).astype(int)
V["n_vendors"] = V.sum(axis=1)
V["on_mPS3"] = M["FOLI_mPS3"].astype(int)
V["on_hPS2"] = M["FOLI_hPS2"].astype(int)
V["on_FoliSeq"] = ((V.on_mPS3 + V.on_hPS2) > 0).astype(int)
V = V[(V.n_vendors > 0) | (V.on_FoliSeq > 0)]

# Hallmark annotation (Liberzon et al. 2015) so consensus genes can be read by function
hall = gs[gs.collection == "Hallmark"]
hmap = {}
for n, genes in zip(hall.gs_name, hall.genes):
    for g in genes:
        hmap.setdefault(g, []).append(n.replace("HALLMARK_", "").lower())
V["hallmark_sets"] = ["; ".join(sorted(hmap.get(g, []))) for g in V.index]
V.rename_axis("gene").to_csv(O + "q1_gene_vendor_matrix.csv")

# UpSet counts over the four vendors + Foli-seq
sets = list(VENDOR_PANELS) + ["on_FoliSeq"]
up = V[sets].astype(bool).value_counts().rename("n_genes").reset_index()
up["sets"] = up[sets].apply(lambda r: " & ".join(s.replace("on_FoliSeq", "Foli-seq") for s in sets if r[s]), axis=1)
up[["sets", "n_genes"]].to_csv(O + "q1_upset_counts.csv", index=False)

cons = V[V.n_vendors == 4].copy()
cons.rename_axis("gene").to_csv(O + "q1_consensus_genes.csv")
dist = V[V.n_vendors > 0].n_vendors.value_counts().sort_index().to_dict()
foli = V[V.on_FoliSeq == 1]
summary = dict(
    genes_on_any_vendor_general_panel=int((V.n_vendors > 0).sum()),
    genes_by_number_of_vendors={int(k): int(v) for k, v in dist.items()},
    consensus_all_four_vendors=int(len(cons)),
    consensus_on_FoliSeq=int(cons.on_FoliSeq.sum()),
    consensus_on_mPS3=int(cons.on_mPS3.sum()),
    consensus_on_hPS2=int(cons.on_hPS2.sum()),
    FoliSeq_genes=int(len(foli)),
    FoliSeq_genes_on_no_vendor_panel=int((foli.n_vendors == 0).sum()),
    FoliSeq_genes_on_all_four=int((foli.n_vendors == 4).sum()),
)
json.dump(summary, open(O + "q1_summary.json", "w"), indent=1)
print(json.dumps(summary, indent=1))
print(up.head(12)[["sets", "n_genes"]].to_string())
