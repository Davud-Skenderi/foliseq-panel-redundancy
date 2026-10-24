"""Step 1: parse every vendor panel in the attached zip, map to human symbols, build gene x panel matrix.

Outputs (in data/):
  panel_meta.csv        one row per panel (platform, species, modality, n genes, mapping rate)
  panels_long.csv       one row per (panel, original gene) with human ortholog
  gene_by_panel.csv     binary human-gene x panel membership
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import re, warnings
import sys
sys.path.insert(0, PROJ + "scripts")
import pandas as pd, numpy as np, rdata

warnings.filterwarnings("ignore")
RAW = RAW + ""
EXT = EXT + ""
OUT = PROJ + "data/"

from common import *  # noqa  (canon_h, to_human, orthology tables)
import sys

# ---------------------------------------------------------------- parsers
def nano_sheet(path, sheet, col="Gene Symbol(s)", hcol=None):
    d = pd.read_excel(RAW + path, sheet, header=1)
    d = d[d[col].notna()]
    rows = []
    for _, r in d.iterrows():
        # multi-gene probes are separated by comma / slash; keep each target
        for g in re.split(r"\s*[,/;|]\s*", str(r[col])):
            if g.strip():
                rows.append((g.strip(), r[hcol] if hcol else None))
    return rows


def simple_list(path, sheet, skip):
    d = pd.read_excel(RAW + path, sheet, header=None)
    v = d.iloc[skip:, 0].dropna().astype(str).str.strip()
    return [(g, None) for g in v if g and g.upper() != "COUNT"]


P = []  # (panel_id, platform, panel_name, species, modality, technology_class, rows)
xh = pd.read_csv(RAW + "10X Genomics Xenium/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")
xm = pd.read_csv(RAW + "10X Genomics Xenium/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv")
P += [
    ("XEN_H5K", "10x Xenium", "Xenium Prime 5K Pan-Tissue & Pathways (Human)", "Human", "RNA", "Imaging, in situ", [(g, None) for g in xh.gene_name]),
    ("XEN_M5K", "10x Xenium", "Xenium Prime 5K Pan-Tissue & Pathways (Mouse)", "Mouse", "RNA", "Imaging, in situ", [(g, None) for g in xm.gene_name]),
    ("COSMX_H6K", "NanoString CosMx", "CosMx Human 6K Discovery", "Human", "RNA", "Imaging, in situ",
     nano_sheet("NanoString CosMx SMI/LBL-11190-04-Human-6k-Discovery-Gene-List.xlsx", "Gene and Probe Details")),
    ("COSMX_HUCC", "NanoString CosMx", "CosMx Human Universal Cell Characterization (1K)", "Human", "RNA", "Imaging, in situ",
     nano_sheet("NanoString CosMx SMI/LBL-11178-04-Human-Universal-Cell-Characterization-Panel-Gene-Target-List.xlsx", "Gene and Probe Details")),
    ("COSMX_MUCC", "NanoString CosMx", "CosMx Mouse Universal Cell Characterization (1K)", "Mouse", "RNA", "Imaging, in situ",
     nano_sheet("NanoString CosMx SMI/LBL-11176-04-Mouse-Universal-Cell-Characterization-Gene-List.xlsx", "Gene and Probe Details", hcol="Human Gene")),
    ("COSMX_HIO", "NanoString CosMx", "CosMx Human Immuno-Oncology add-on (100)", "Human", "RNA", "Imaging, in situ",
     nano_sheet("NanoString CosMx SMI/LBL-11180-04-Human-IO-RNA-Panel-Gene-Target-List.xlsx", "Gene and Probe Details")),
    ("COSMX_MNEURO", "NanoString CosMx", "CosMx Mouse Neuroscience (1K)", "Mouse", "RNA", "Imaging, in situ",
     nano_sheet("NanoString CosMx SMI/LBL-11179-05-Mouse-Neuroscience-Panel-Gene-List.xlsx", "Gene and Probe Details")),
    ("COSMX_MNEURO_PROT", "NanoString CosMx", "CosMx Mouse Neuroscience Protein panel", "Mouse", "Protein", "Imaging, in situ",
     nano_sheet("NanoString CosMx SMI/CosMx#U2122-Mouse-Neuroscience-Protein-Panel-Protein-Target-List.xlsx", "Protein Details", col="Gene")),
]
for pid, f, name in [
    ("MER_M_PANNEURO", "FinalGeneList_MERSCOPE815mousepanneuron_D2M1763.csv", "MERSCOPE 815 Mouse PanNeuro"),
    ("MER_M_PAN", "FinalGeneList_MERSCOPE815panmouse_D2M1791.csv", "MERSCOPE 815 Pan-Mouse"),
    ("MER_M_BRAINTRAIN", "FinalGeneList_mousebraintrainingpanel_A2M1620.csv", "MERSCOPE Mouse Brain Training (128)"),
]:
    d = pd.read_csv(RAW + "Vizgen MERSCOPE/" + f)
    P.append((pid, "Vizgen MERSCOPE", name, "Mouse", "RNA", "Imaging, in situ (MERFISH)", [(g, None) for g in d["Vizgen Gene"]]))
P += [
    ("GEOMX_H_WTA", "NanoString GeoMx", "GeoMx Human Whole Transcriptome Atlas", "Human", "RNA", "ROI profiling (DSP)",
     nano_sheet("Nanostring GeoMx/GeoMx-Human-Whole-Transcriptome-Atlas-Probe-List_v0.1.xlsx", "Gene Details")),
    ("GEOMX_M_WTA", "NanoString GeoMx", "GeoMx Mouse Whole Transcriptome Atlas v2", "Mouse", "RNA", "ROI profiling (DSP)",
     nano_sheet("Nanostring GeoMx/lbl-10861-02_geomx_mouse_whole_transcriptome_atlas_probe_list.xlsx", "Gene Details")),
    ("GEOMX_H_CTA", "NanoString GeoMx", "GeoMx Cancer Transcriptome Atlas", "Human", "RNA", "ROI profiling (DSP)",
     nano_sheet("Nanostring GeoMx/LBL-10809-03_GeoMx_Cancer_Transcriptome_Atlas_Probe_List.xlsx", "Gene Details")),
    ("GEOMX_DOG_CCA", "NanoString GeoMx", "GeoMx Canine Cancer Atlas", "Dog", "RNA", "ROI profiling (DSP)",
     nano_sheet("Nanostring GeoMx/LBL-11121-01-Canine-Cancer-Atlas-Gene-List.xlsx", "Gene Details")),
    ("AMPLI_H_TX", "Illumina AmpliSeq", "AmpliSeq Transcriptome Human Gene Expression", "Human", "RNA", "Bulk targeted sequencing",
     simple_list("Illumina AmpliSeq/ampliseq-for-illumina-transcriptome-gene-list.xlsx", "Sheet1", 4)),
    ("AMPLI_COMP_V3", "Illumina AmpliSeq", "AmpliSeq Comprehensive Panel v3 (DNA)", "Human", "DNA", "Bulk targeted sequencing",
     simple_list("Illumina AmpliSeq/ampliseq-for-illumina-comprehensive-v3-gene-list.xlsx", "All Genes", 4)),
    ("IDT_XGEN_PANCANCER", "IDT xGen", "xGen Pan-Cancer hybridization panel (DNA)", "Human", "DNA", "Bulk targeted sequencing",
     simple_list("IDT xGen/xgen-pan-cancer-gene-list.xlsx", "GeneList", 1)),
]
# canine: use NanoString's own "Human Gene" mapping from the Annotations sheet when available
can_ann = pd.read_excel(RAW + "Nanostring GeoMx/LBL-11121-01-Canine-Cancer-Atlas-Gene-List.xlsx", "Annotations", header=1)
CAN_H = dict(zip(can_ann["Gene"], can_ann["Human Gene"]))
P = [(a, b, c, d, e, f, [(g, CAN_H.get(g)) for g, _ in rows] if a == "GEOMX_DOG_CCA" else rows) for a, b, c, d, e, f, rows in P]

# Foli-seq panels: one CSV per panel in data/foli_panels/ (written by 00_foli_inputs.py), named <PANEL_ID>__<Species>.csv
# (e.g. FOLI_mPS3__Mouse.csv, FOLI_hPS2__Human.csv) with a column 'gene'. They are picked up automatically.
import glob, os
for fp in sorted(glob.glob(OUT + "foli_panels/*.csv")):
    pid, sp = os.path.basename(fp)[:-4].split("__")
    g = pd.read_csv(fp)["gene"].dropna().astype(str).str.strip().unique()
    P.append((pid, "Foli-seq (Wang Lab)", f"Foli-seq {pid}", sp, "RNA", "Fecal exfoliome targeted multiplex PCR", [(x, None) for x in g]))
# the public 70-gene excerpt is a subset of hPS2: drop it once the full panels are present
if any(p[0].startswith("FOLI_") and p[0] != "FOLI_hIBD_70" for p in P):
    P = [p for p in P if p[0] != "FOLI_hIBD_70"]

# ---------------------------------------------------------------- assemble
long = []
for pid, plat, name, sp, mod, tech, rows in P:
    for g, hint in rows:
        long.append(dict(panel_id=pid, platform=plat, panel_name=name, species=sp, modality=mod, technology=tech,
                         original_symbol=g, human_symbol=to_human(g, sp, hint)))
long = pd.DataFrame(long).drop_duplicates(["panel_id", "original_symbol"])

# Foli excerpt: a few primers resolved to transcript IDs (ENST...) — resolve via Ensembl transcript table if possible
long.loc[long.original_symbol.str.startswith("ENST", na=False), "human_symbol"] = None

meta = long.groupby(["panel_id", "platform", "panel_name", "species", "modality", "technology"]).agg(
    n_targets=("original_symbol", "nunique"), n_human_mapped=("human_symbol", "nunique")).reset_index()
meta["mapping_rate"] = (long.groupby("panel_id").human_symbol.apply(lambda s: s.notna().mean()).reindex(meta.panel_id).values).round(3)
# missing panels listed in the master spreadsheet but absent/empty in the zip
missing = pd.DataFrame([
    dict(panel_id="MER_H_IO", platform="Vizgen MERSCOPE", panel_name="MERSCOPE 815 Immuno-Oncology (Human)", species="Human",
         modality="RNA", technology="Imaging, in situ (MERFISH)", n_targets=0, n_human_mapped=0, mapping_rate=np.nan,
         note="files in zip are empty (0-byte CSV, blank XLSX)"),
])
meta["note"] = ""
meta.loc[meta.panel_id.str.startswith("FOLI_m"), "note"] = "Foli-seq mouse panel, Huang et al. 2025 Table S2 (ERCC spike-ins removed)"
meta.loc[meta.panel_id == "FOLI_hPS2", "note"] = "Foli-seq human panel, Huang et al. 2025 Table S8"
meta.loc[meta.panel_id == "FOLI_hIBD_70", "note"] = "alphabetical excerpt of the real human panel recovered from public folitools test primers; full mPS/hPS lists needed"
meta.loc[meta.panel_id.isin(["AMPLI_COMP_V3", "IDT_XGEN_PANCANCER"]), "note"] = "DNA (mutation/CNV) panel - not an expression panel"
meta.loc[meta.panel_id == "COSMX_MNEURO_PROT", "note"] = "protein targets mapped to their genes"
meta = pd.concat([meta, missing], ignore_index=True)

M = (long.dropna(subset=["human_symbol"]).assign(v=1)
     .pivot_table(index="human_symbol", columns="panel_id", values="v", aggfunc="max", fill_value=0))
long.to_csv(OUT + "panels_long.csv", index=False)
meta.to_csv(OUT + "panel_meta.csv", index=False)
M.to_csv(OUT + "gene_by_panel.csv")
print(meta[["panel_id", "species", "modality", "n_targets", "n_human_mapped", "mapping_rate"]].to_string())
print("union human genes:", M.shape)
print("unmapped examples:", long[long.human_symbol.isna()].groupby("panel_id").original_symbol.apply(lambda s: list(s)[:6]).to_string())
