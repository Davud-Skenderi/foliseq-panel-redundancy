"""Add the CosMx Whole Transcriptome (WTx) panels, human and mouse, to the panel matrix.

These were the sales-gated lists this project could not previously obtain; they were supplied
directly. They matter less as "another panel" than as a CEILING: whole-transcriptome imaging measures
essentially everything a CosMx instrument can, so the accuracy it reaches bounds what ANY targeted
panel on that platform could achieve. A targeted panel is then readable as a fraction of that ceiling
rather than as an absolute number with nothing above it.

Each workbook also carries a "Known Missing Targets" sheet listing genes deliberately or
unavoidably excluded, with a rationale. That is a direct, vendor-stated gap list and is saved
alongside: even the whole transcriptome has holes, and knowing which ones are intentional
(housekeeping genes dropped to avoid saturating detection) versus technical matters for any claim
about what the platform cannot see.

Mouse symbols are mapped to human with common.to_human, the same mapping used everywhere else in
this repo, so the new columns join data/gene_by_panel.csv in human-symbol space.
"""
import os, sys, re
import pandas as pd, numpy as np, openpyxl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ
from common import to_human, canon_h

SRC = {
    "COSMX_H_WTX": ("/mnt/user-data/uploads/3c077ad5-LBL-11220-01_CosMx_SMI_Human_Whole_Transcriptome_WTX_GeneList.xlsx", "Human"),
    "COSMX_M_WTX": ("/mnt/user-data/uploads/dfc0ff13-LBL-11447-01_CosMx_SMI_Mouse_WTX_Gene_List.xlsx", "Mouse"),
}
D = PROJ + "data/"


def sheet_rows(path, sheet, header_row=1):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    hdr, out = None, []
    for i, r in enumerate(ws.iter_rows(values_only=True)):
        if i < header_row:
            continue
        if i == header_row:
            hdr = [str(c).strip() if c is not None else "" for c in r]
            continue
        out.append(r)
    wb.close()
    return hdr, out


def primary(sym):
    """Gene Symbol(s) separates alternatives with ' | '; the first is the primary symbol."""
    parts = [p.strip() for p in str(sym).split("|") if p.strip()]
    return parts[0] if parts else ""


def main():
    M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
    print(f"gene_by_panel.csv before: {M.shape[0]} genes x {M.shape[1]} panels")
    added = {}
    for tag, (path, species) in SRC.items():
        hdr, rows = sheet_rows(path, "Gene and Probe Details")
        ci = [h.lower() for h in hdr].index("gene symbol(s)")
        raw = [primary(r[ci]) for r in rows if ci < len(r) and r[ci] is not None]
        raw = [g for g in raw if g]
        if species == "Mouse":
            mapped = [to_human(g, "Mouse") for g in raw]
        else:
            mapped = [canon_h(g) for g in raw]
        hs = sorted({h for h in mapped if isinstance(h, str) and h})
        added[tag] = hs
        print(f"{tag}: {len(raw)} source genes -> {len(hs)} human symbols "
              f"({len(raw)-len(hs)} lost to dropping or many-to-one collapse)")

        # vendor-stated gaps
        h2, r2 = sheet_rows(path, "Known Missing Targets")
        gi = [h.lower() for h in h2].index("gene symbol")
        ri = [h.lower() for h in h2].index("rationale")
        miss = pd.DataFrame([{"panel_id": tag, "gene_symbol": str(r[gi]).strip(),
                              "rationale": str(r[ri]).strip() if ri < len(r) and r[ri] else ""}
                             for r in r2 if gi < len(r) and r[gi] is not None])
        miss.to_csv(D + f"known_missing_{tag}.csv", index=False)
        cats = miss.rationale.str.slice(0, 46).value_counts()
        print(f"   known missing targets: {len(miss)}")
        for k, v in cats.head(4).items():
            print(f"      {v:5d}  {k}")

    all_genes = sorted(set(M.index) | set().union(*[set(v) for v in added.values()]))
    M2 = M.reindex(all_genes).fillna(0).astype(int)
    for tag, hs in added.items():
        M2[tag] = 0
        M2.loc[[g for g in hs if g in M2.index], tag] = 1
    M2.to_csv(D + "gene_by_panel.csv")
    print(f"\ngene_by_panel.csv after: {M2.shape[0]} genes x {M2.shape[1]} panels")
    print("panel sizes:")
    for c in sorted(M2.columns, key=lambda c: -M2[c].sum()):
        print(f"   {c:20s} {M2[c].sum():6d}")


if __name__ == "__main__":
    main()
