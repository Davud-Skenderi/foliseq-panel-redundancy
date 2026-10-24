"""Step 0: Foli-seq panels and stool-responsive genes from the paper's supplement (Huang et al. 2025, Tables S1-S8).

Writes data/foli_panels/<ID>__<Species>.csv (gene column, ERCC spike-ins removed) and
data/foli_stool_degs.csv: one row per mouse gene with the colitis evidence found in Foli-seq stool data.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import os, re
import pandas as pd
SUP = SUP
D = PROJ + "data/"
os.makedirs(D + "foli_panels", exist_ok=True); os.makedirs(D + "foli_panels_extra", exist_ok=True)
x = pd.ExcelFile(SUP)


def blocks(sheet):
    d = pd.read_excel(x, sheet, header=None)
    out = {}
    for c in range(d.shape[1]):
        if str(d.iloc[1, c]).strip() == "Gene Symbol":
            title = next(str(d.iloc[0, k]) for k in (c - 1, c) if pd.notna(d.iloc[0, k]))
            out[title] = d.iloc[2:, c].dropna().astype(str).str.strip()
    return out


panels = {}
for t, g in blocks("S2").items():
    panels["FOLI_" + t.split()[0] + "__Mouse"] = g
for t, g in blocks("S8").items():
    panels["FOLI_" + t.split()[0] + "__Human"] = g
for k, g in panels.items():
    g = sorted(set(g[~g.str.startswith("ERCC")]))
    sub = "foli_panels" if k.split("__")[0] in ("FOLI_mPS3", "FOLI_hPS2") else "foli_panels_extra"  # analysed vs reported
    pd.DataFrame({"gene": g}).to_csv(D + f"{sub}/{k}.csv", index=False)
    print(k, len(g))

# ---- stool-responsive genes (all mouse)
ev = {}


def add(genes, tag):
    for g in set(genes):
        ev.setdefault(g, set()).add(tag)


s3 = pd.read_excel(x, "S3")
add(s3["Gene Symbol"].dropna().astype(str), "DSS_temporal")
s4 = pd.read_excel(x, "S4", header=None)
add(s4.iloc[2:, 1].dropna().astype(str), "shared_3_models")
s5 = pd.read_excel(x, "S5", header=None).iloc[2:]
s5[0] = s5[0].ffill()
for m, sub in s5.groupby(0):
    add(sub[2].dropna().astype(str), "unique_" + re.sub(r"\W+", "_", str(m)).strip("_"))
s6 = pd.read_excel(x, "S6")
add(s6["Gene Symbol"].dropna().astype(str), "antiP19_treatment_response")
s7 = pd.read_excel(x, "S7", header=None)
sec = None
for _, r in s7.iterrows():
    a, b = r[0], r[1]
    if isinstance(a, str) and not a.startswith("ENSMUSG") and pd.isna(b) and a != "Gene ID":
        sec = "S7_" + re.sub(r"\W+", "_", a).strip("_")[:40]
    elif isinstance(a, str) and a.startswith("ENSMUSG") and isinstance(b, str):
        add([b.strip()], sec)
E = pd.DataFrame([(g, len(v), ";".join(sorted(v))) for g, v in ev.items()], columns=["mouse_gene", "n_evidence", "evidence"])
E.to_csv(D + "foli_stool_degs.csv", index=False)
print("stool-responsive genes:", len(E)); print(E.evidence.str.split(";").explode().value_counts())
