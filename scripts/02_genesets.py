"""Step 2: build gene-set library (MSigDB 7.5.1 via msigdbr) with HGNC-canonical symbols.

Collections kept: Hallmark, KEGG, Reactome, BioCarta, PID, WikiPathways, GO:BP, C8 cell-type signatures.
Also defines a curated 'gut/IBD lens' subset used for the Foli-seq-focused analysis.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import sys, re, warnings
import pandas as pd, rdata
warnings.filterwarnings("ignore")
sys.path.insert(0, PROJ + "scripts")
from common import canon_h  # noqa

EXT = EXT + ""
OUT = PROJ + "data/"
d = rdata.conversion.convert(rdata.parser.parse_file(EXT + "msigdbr/R/sysdata.rda"))
gs, gg, genes = d["msigdbr_genesets"], d["msigdbr_geneset_genes"], d["msigdbr_genes"]
keep = {("H", ""): "Hallmark", ("C2", "CP:KEGG"): "KEGG", ("C2", "CP:REACTOME"): "Reactome",
        ("C2", "CP:BIOCARTA"): "BioCarta", ("C2", "CP:PID"): "PID", ("C2", "CP:WIKIPATHWAYS"): "WikiPathways",
        ("C5", "GO:BP"): "GO:BP", ("C8", ""): "C8 cell types"}
gs["collection"] = [keep.get((c, s)) for c, s in zip(gs.gs_cat, gs.gs_subcat)]
gs = gs[gs.collection.notna()]
sym = dict(zip(genes.gene_id, genes.human_gene_symbol))
gg = gg[gg.gs_id.isin(gs.gs_id)].copy()
gg["sym"] = gg.gene_id.map(sym)
cmap = {s: canon_h(s) or s for s in gg.sym.dropna().unique()}
gg["sym"] = gg.sym.map(cmap)
members = gg.dropna().groupby("gs_id").sym.apply(lambda s: sorted(set(s)))
gs["genes"] = gs.gs_id.map(members)
gs = gs[gs.genes.notna()]
gs["size"] = gs.genes.str.len()
gs = gs[(gs["size"] >= 5) & (gs["size"] <= 500)]

# gut / IBD lens: intestinal cell-type signatures + immune/barrier pathways most relevant to Foli-seq
gut_pat = re.compile(r"INTESTIN|COLON|DUODEN|ILEUM|GAO_LARGE|GAO_SMALL|BUSSLINGER_DUODENAL")
ibd_pat = re.compile(r"INTERLEUKIN_(1|2|4|6|10|12|13|17|20|22|23)|IL17|IL23|IL22|TNF|INTERFERON|NFKB|JAK_STAT|TOLL|"
                     r"INFLAMMATORY|CHEMOKINE|ANTIGEN_PROCESSING|COMPLEMENT|TIGHT_JUNCTION|MUCIN|O_LINKED_GLYCOSYLATION_OF_MUCINS|"
                     r"T_CELL_RECEPTOR|B_CELL_RECEPTOR|NATURAL_KILLER|TH1_AND_TH2|TH17|APOPTOSIS|EPITHELIAL_MESENCHYMAL|"
                     r"WNT|NOTCH|HYPOXIA|REACTIVE_OXYGEN|XENOBIOTIC|BILE_ACID|FATTY_ACID|CELL_CYCLE|G2M|E2F|"
                     r"INTESTINAL_IMMUNE_NETWORK_FOR_IGA|CYTOKINE_CYTOKINE_RECEPTOR")
gs["gut_lens"] = [
    (c == "C8 cell types" and bool(gut_pat.search(n))) or
    (c in ("Hallmark", "KEGG", "Reactome") and bool(ibd_pat.search(n)))
    for n, c in zip(gs.gs_name, gs.collection)]
gs[["gs_name", "collection", "size", "gut_lens", "gs_description", "genes"]].to_pickle(OUT + "genesets.pkl")
print(gs.groupby("collection").size().to_string())
print("gut lens sets:", gs.gut_lens.sum())
print(gs[gs.gut_lens].groupby("collection").size())
