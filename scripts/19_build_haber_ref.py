"""Build the Haber et al. 2017 mouse small-intestinal epithelium reference (GEO GSE92332).

Why this reference specifically. With MM_SLN (mouse immune) and the three human gut references, the
contrast "immune vs gut" is perfectly confounded with "mouse vs human", so neither effect can be
estimated. Haber is mouse AND gut, which completes a 2x2:

                 immune            gut
    mouse        MM_SLN            HABER_GUT      <- tissue effect, species held fixed
    human        (HS_PBMC, unann.) GUT_HEALTHY/IBD/CRC
                                   ^ species effect, tissue held fixed (mouse vs human gut)

Cell types are encoded in the count matrix's column names as barcode fields:
    B1_AAACATTGTTTGGG_Enterocyte.Immature.Distal  ->  celltype "Enterocyte.Immature.Distal"

Gene symbols are left in MOUSE space; 17_multiref_eval.py maps them with common.to_human, the same
mapping 12_spapros_eval.py applied to MM_SLN, so the two mouse references are directly comparable.
"""
import os, sys, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import EXT

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_atlas_UMIcounts.txt.gz"
SRC = os.environ.get("HABER_SRC", "/home/claude/haber_raw/atlas.txt.gz")
OUT = EXT + "expr/haber_gut.h5ad"
MIN_CELLS_PER_TYPE = 40


def main():
    if not os.path.exists(SRC):
        raise SystemExit(f"download first:\n  curl -sSL -o {SRC} {URL}")
    df = pd.read_csv(SRC, sep="\t", index_col=0)          # genes x cells, raw UMI counts
    print(f"loaded {df.shape[0]} genes x {df.shape[1]} cells", flush=True)

    celltype = pd.Series([c.split("_")[2] if len(c.split("_")) > 2 else "unknown" for c in df.columns],
                         index=df.columns)
    vc = celltype.value_counts()
    keep = celltype.isin(vc.index[vc >= MIN_CELLS_PER_TYPE])
    df, celltype = df.loc[:, keep.values], celltype[keep.values]
    print(f"kept {keep.sum()} cells across {celltype.nunique()} cell types "
          f"(>= {MIN_CELLS_PER_TYPE} cells each)", flush=True)

    X = sp.csr_matrix(df.T.values.astype(np.float32))      # cells x genes
    a = ad.AnnData(X=X,
                   obs=pd.DataFrame({"celltype": celltype.values,
                                     "batch": [c.split("_")[0] for c in df.columns]},
                                    index=df.columns.astype(str)),
                   var=pd.DataFrame(index=df.index.astype(str)))
    a.var_names_make_unique()
    keep_g = np.asarray((a.X > 0).sum(axis=0)).ravel() >= 10
    a = a[:, keep_g].copy()
    a.uns["reference_title"] = "Haber et al. 2017 mouse small intestinal epithelium (GSE92332)"
    a.uns["species"] = "Mouse"
    a.uns["source_url"] = URL
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    a.write_h5ad(OUT, compression="gzip")
    print(f"-> {OUT}  {a.shape}  {a.obs.celltype.nunique()} cell types", flush=True)
    print("   types:", ", ".join(sorted(a.obs.celltype.unique())), flush=True)


if __name__ == "__main__":
    main()
