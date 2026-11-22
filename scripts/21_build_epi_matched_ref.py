"""Derive an epithelium-only human gut reference, to match Haber's compartment.

HABER_GUT is mouse small-intestinal EPITHELIUM (15 types). The human gut references are whole tissue:
epithelial, immune and stromal cells together (72/54/45 types). So the contrast HABER_GUT -> human gut
is not species alone — it is also epithelium -> mixed tissue, and a different annotation granularity.

Restricting the Gut Cell Atlas to its epithelial cell types gives a human reference matched to Haber
on compartment (12 vs 15 epithelial types), so that contrast isolates species far more cleanly.

The epithelial type list is matched on the reference's own ontology labels and is printed on every run
so the selection stays auditable rather than hidden in a keyword regex.
"""
import os, sys, anndata as ad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import EXT

EPITHELIAL = [
    "M cell of gut", "type EC enteroendocrine cell", "transit amplifying cell", "stem cell",
    "colon epithelial cell", "colonocyte", "epithelial cell", "intestine goblet cell",
    "enteroendocrine cell", "intestinal tuft cell", "type I enteroendocrine cell",
    "type D enteroendocrine cell",
]
MIN_CELLS = 40

def main():
    src = EXT + "expr/gut_healthy.h5ad"
    a = ad.read_h5ad(src)
    before = a.obs.celltype.nunique()
    a = a[a.obs.celltype.isin(EPITHELIAL)].copy()
    vc = a.obs.celltype.value_counts()
    a = a[a.obs.celltype.isin(vc.index[vc >= MIN_CELLS])].copy()
    a.uns["reference_title"] = "Gut Cell Atlas, epithelial compartment only (Haber-matched)"
    a.uns["derived_from"] = "gut_healthy.h5ad"
    out = EXT + "expr/gut_healthy_epi.h5ad"
    a.write_h5ad(out, compression="gzip")
    print(f"{before} cell types -> {a.obs.celltype.nunique()} epithelial types, {a.n_obs} cells")
    for t, n in a.obs.celltype.value_counts().items():
        print(f"   {n:6d}  {t}")
    print(f"-> {out}  {a.shape}")

if __name__ == "__main__":
    main()
