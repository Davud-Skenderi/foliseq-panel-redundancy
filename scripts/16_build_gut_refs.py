"""Build annotated gut single-cell references for Q3/Q4 from the CELLxGENE Census.

The README's open item: Q3/Q4 ran on immune references only (mouse spleen/lymph node, human PBMC),
so epithelial conclusions rested on gene lists alone. These three human intestinal references close
that gap and let the same metrics be recomputed on gut tissue across three disease contexts:

  GUT_HEALTHY  Gut Cell Atlas (Elmentaite et al. 2021)      healthy intestine
  GUT_IBD      metaplasia in inflammatory gut disease        inflamed
  GUT_CRC      colorectal cancer immune hubs                 malignant

Census genes already carry HGNC symbols (var feature_name), so no ortholog mapping is needed and
common.py is not imported. Raw counts are kept: 12/17 normalise themselves.
The FULL measured transcriptome is kept, not only panel genes, so that random baselines are drawn
from the transcriptome rather than from the union of commercial panels.
"""
import os, sys, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
import cellxgene_census

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import EXT

CENSUS_VERSION = "2025-11-08"
PER_TYPE  = 250     # cells sampled per annotated cell type
MIN_CELLS = 40      # Spapros n_cells_min; drop rarer types here
TISSUE    = "tissue_general in ['large intestine','small intestine']"

REFS = {
    "GUT_HEALTHY": ("80a2c5b6-02e7-4fc0-9f12-179f5247c1bc", "Gut Cell Atlas (Elmentaite et al. 2021)"),
    "GUT_IBD":     ("19053a82-9c89-4fb8-bd19-d7b1800b0b7b", "Metaplasia in inflammatory gut disease"),
    "GUT_CRC":     ("16023185-de21-4c0d-a9c8-73abdd52d142", "CRC immune hubs"),
}

def build(tag, dataset_id, title, out_dir):
    rng = np.random.default_rng(0)
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as c:
        obs = c["census_data"]["homo_sapiens"].obs.read(
            value_filter=f"dataset_id == '{dataset_id}' and is_primary_data == True and {TISSUE}",
            column_names=["soma_joinid", "cell_type", "tissue", "disease", "donor_id", "assay"],
        ).concat().to_pandas()
        vc = obs.cell_type.value_counts()
        obs = obs[obs.cell_type.isin(vc.index[vc >= MIN_CELLS])]
        idx = (obs.groupby("cell_type", observed=True)["soma_joinid"]
                  .apply(lambda s: pd.Series(rng.choice(s.values, min(len(s), PER_TYPE), replace=False))))
        coords = np.sort(idx.values.astype(np.int64))
        print(f"{tag}: {len(coords)} cells, {obs.cell_type.nunique()} cell types", flush=True)
        a = cellxgene_census.get_anndata(
            census=c, organism="homo_sapiens", obs_coords=coords,
            column_names={"obs": ["cell_type", "tissue", "disease", "donor_id", "assay"],
                          "var": ["feature_id", "feature_name"]},
        )
    a.var_names = a.var["feature_name"].astype(str).values
    a.var_names_make_unique()
    a.obs["celltype"] = a.obs["cell_type"].astype(str)
    # drop genes seen in <10 cells (same filter 12_spapros_eval applies)
    keep = np.asarray((a.X > 0).sum(axis=0)).ravel() >= 10
    a = a[:, keep].copy()
    a.uns["reference_title"] = title
    a.uns["census_version"] = CENSUS_VERSION
    a.uns["dataset_id"] = dataset_id
    out = os.path.join(out_dir, f"{tag.lower()}.h5ad")
    a.write_h5ad(out, compression="gzip")
    print(f"  -> {out}  {a.shape}  {a.obs.celltype.nunique()} types", flush=True)

if __name__ == "__main__":
    out_dir = EXT + "expr"
    os.makedirs(out_dir, exist_ok=True)
    only = sys.argv[1:] or list(REFS)
    for tag in only:
        did, title = REFS[tag]
        build(tag, did, title, out_dir)
