"""How much host-microbiome interface biology does each panel actually carry?

The recurring finding across this project is that commercial panels are curated for other tissues and
lose power on gut epithelium. That says what they are bad at. This says what is MISSING, concretely,
which is what a panel designer can act on.

The interface set is the biology that makes the gut lining a barrier and an interface with the
microbiota: antimicrobial peptides, mucins, pattern-recognition receptors, tight junctions, IgA
transport, and the trefoil/goblet secretory programme. If Foli-seq reads cells shed from that lining,
this is the content it exists to measure.

Coverage is reported against the CosMx Human WTx panel as the denominator -- the set of interface
genes a whole-transcriptome imaging assay can actually see -- rather than against an abstract genome,
so "missing" means missing relative to what the platform could have included.
"""
import os, sys
import pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ

D, O = PROJ + "data/", PROJ + "results/v2/"
INTERFACE = {
    "antimicrobial peptides":   lambda g: g.startswith(("DEFA", "DEFB")) or g in
                                {"LYZ", "REG3A", "REG3G", "PLA2G2A", "LCN2", "CAMP", "ITLN1", "ITLN2", "PRSS2"},
    "mucins":                   lambda g: g.startswith("MUC"),
    "pattern-recognition":      lambda g: g.startswith(("TLR", "NLR", "NOD")) or g in {"CARD9", "RIPK2", "MYD88", "TICAM1"},
    "tight junction / barrier": lambda g: g.startswith(("CLDN", "TJP")) or g in {"OCLN", "CGN", "MARVELD2", "F11R", "CRB3"},
    "IgA transport":            lambda g: g in {"PIGR", "CCL28", "JCHAIN", "FCGRT", "TNFSF13B"},
    "trefoil / goblet":         lambda g: g.startswith("TFF") or g in {"FCGBP", "CLCA1", "ZG16", "AGR2", "SPDEF", "ATOH1"},
    "Paneth / secretory":       lambda g: g in {"LYZ1", "MMP7", "ANG4", "PLA2G2A", "DLL1", "DLL4", "SOX9"},
}
PANELS = ["FOLI_hPS2", "FOLI_mPS3", "COSMX_HIO", "MER_M_BRAINTRAIN", "MER_M_PAN", "MER_M_PANNEURO",
          "COSMX_MNEURO", "COSMX_MUCC", "COSMX_HUCC", "GEOMX_H_CTA", "XEN_H5K", "XEN_M5K", "COSMX_H6K"]
CEILING = "COSMX_H_WTX"


def main():
    M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
    assert CEILING in M.columns, f"{CEILING} missing; run 23_add_cosmx_wtx.py first"
    wtx = set(M.index[M[CEILING] == 1])
    cat_genes = {c: sorted(g for g in wtx if f(g)) for c, f in INTERFACE.items()}
    universe = sorted(set().union(*cat_genes.values()))
    print(f"interface set measurable on CosMx Human WTx: {len(universe)} genes")
    for c, gs in cat_genes.items():
        print(f"   {c:26s} {len(gs):4d}")

    rows = []
    for p in [c for c in PANELS if c in M.columns]:
        got = set(M.index[M[p] == 1])
        d = dict(panel_id=p, panel_size=int(M[p].sum()),
                 interface_genes=len(got & set(universe)),
                 interface_pct_of_set=round(100 * len(got & set(universe)) / len(universe), 1))
        # share of the panel's own budget spent on interface biology
        d["pct_of_panel_budget"] = round(100 * len(got & set(universe)) / max(int(M[p].sum()), 1), 2)
        for c, gs in cat_genes.items():
            d[c] = len(got & set(gs))
        rows.append(d)
    R = pd.DataFrame(rows).sort_values("interface_pct_of_set", ascending=False)
    R.to_csv(O + "q11_interface_coverage.csv", index=False)

    print(f"\n=== coverage of the {len(universe)}-gene interface set ===")
    cols = ["panel_id", "panel_size", "interface_genes", "interface_pct_of_set", "pct_of_panel_budget"]
    print(R[cols].to_string(index=False))

    print("\n=== by category (count of genes carried) ===")
    print(R[["panel_id"] + list(INTERFACE)].to_string(index=False))

    # what no targeted panel carries
    targeted = [c for c in PANELS if c in M.columns and M[c].sum() < 8000]
    carried = set()
    for p in targeted:
        carried |= set(M.index[M[p] == 1])
    orphan = [g for g in universe if g not in carried]
    print(f"\n=== interface genes carried by NO targeted panel: {len(orphan)} of {len(universe)} ===")
    for c, gs in cat_genes.items():
        miss = [g for g in gs if g in orphan]
        if miss:
            print(f"  {c:26s} {len(miss):3d}/{len(gs):3d}  e.g. {', '.join(miss[:8])}")
    pd.DataFrame({"gene_symbol": orphan}).to_csv(O + "q11_interface_orphans.csv", index=False)


if __name__ == "__main__":
    main()
