"""Audit the many-to-one collapses introduced by mouse-to-human ortholog mapping.

WHY THIS MATTERS HERE. Every cross-species comparison in this repo runs mouse data through
common.to_human. When several mouse genes map to one human symbol, to_symbols SUMS their counts into
a single column. That is a silent distortion, and it is not uniform across biology: it concentrates
on gene families that expanded in mouse relative to human.

The gut is one of the places where that expansion is largest. Mouse has a greatly expanded
alpha-defensin (cryptdin) family expressed by Paneth cells; human has far fewer. Defensins are
antimicrobial peptides -- core host-microbiome interface biology, and precisely the category this
project keeps identifying as the Foli-seq opportunity. If 19 mouse defensin genes collapse into one
human DEFA1 column, then any mouse-vs-human claim touching antimicrobial peptides is suspect.

This script quantifies the collapse on the Haber reference, flags the families affected, and checks
how much of the host-microbiome interface gene set is involved.
"""
import os, sys, collections
import pandas as pd, numpy as np, anndata as ad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ, EXT
from common import to_human

# host-microbiome interface categories, by human symbol prefix or explicit membership
INTERFACE = {
    "alpha/beta defensins":      lambda g: g.startswith(("DEFA", "DEFB")),
    "other antimicrobial":       lambda g: g in {"LYZ", "REG3A", "REG3G", "PLA2G2A", "LCN2", "CAMP", "ITLN1", "ITLN2"},
    "mucins":                    lambda g: g.startswith("MUC"),
    "pattern-recognition rcpt":  lambda g: g.startswith(("TLR", "NLR", "NOD")) or g in {"CARD9", "RIPK2"},
    "tight junction":            lambda g: g.startswith(("CLDN", "TJP")) or g in {"OCLN", "CGN", "MARVELD2"},
    "IgA transport / secretion": lambda g: g in {"PIGR", "CCL28", "JCHAIN", "FCGRT"},
    "trefoil / goblet":          lambda g: g.startswith("TFF") or g in {"FCGBP", "CLCA1", "ZG16", "AGR2"},
}

def main():
    a = ad.read_h5ad(EXT + "expr/haber_gut.h5ad", backed="r")
    mouse_genes = [str(g) for g in a.var_names]
    pairs = [(g, to_human(g, "Mouse")) for g in mouse_genes]
    pairs = [(m, h) for m, h in pairs if isinstance(h, str) and h]
    fwd = collections.defaultdict(list)
    for m, h in pairs:
        fwd[h].append(m)
    collapsed = {h: ms for h, ms in fwd.items() if len(ms) > 1}
    n_lost = sum(len(ms) - 1 for ms in collapsed.values())
    print(f"Haber: {len(mouse_genes)} mouse genes -> {len(fwd)} human symbols")
    print(f"  {len(collapsed)} human symbols receive 2+ mouse genes, absorbing {n_lost} extra genes "
          f"({100*n_lost/len(pairs):.1f}% of mapped genes)\n")

    rows = []
    for h, ms in collapsed.items():
        cat = next((c for c, f in INTERFACE.items() if f(h)), "")
        rows.append(dict(human_symbol=h, n_mouse_genes=len(ms), interface_category=cat,
                         mouse_genes="; ".join(sorted(ms))))
    C = pd.DataFrame(rows).sort_values("n_mouse_genes", ascending=False)
    C.to_csv(PROJ + "results/v2/q10_ortholog_collapses.csv", index=False)

    print("worst collapses overall:")
    for _, r in C.head(10).iterrows():
        tag = f"   [{r.interface_category}]" if r.interface_category else ""
        print(f"  {r.human_symbol:12s} <- {r.n_mouse_genes:3d} mouse genes{tag}")

    print("\n=== host-microbiome interface genes ===")
    inter = C[C.interface_category != ""]
    if len(inter):
        print(f"{len(inter)} interface symbols are collapse targets, absorbing "
              f"{int((inter.n_mouse_genes-1).sum())} extra mouse genes\n")
        for cat in INTERFACE:
            sub = inter[inter.interface_category == cat]
            if len(sub):
                worst = "; ".join(f"{r.human_symbol}x{r.n_mouse_genes}" for _, r in sub.head(4).iterrows())
                print(f"  {cat:28s} {len(sub):3d} collapsed symbols   worst: {worst}")
    else:
        print("none")

    # how much of the interface set is even measurable, and how much is collapsed
    print("\n=== coverage of the interface set in Haber (mapped human space) ===")
    present = set(fwd)
    for cat, f in INTERFACE.items():
        hits = sorted(g for g in present if f(g))
        coll = [g for g in hits if g in collapsed]
        print(f"  {cat:28s} {len(hits):3d} present, {len(coll):3d} of them collapsed")

if __name__ == "__main__":
    main()
