"""Build a human immune reference, to fill the empty corner of the species-by-tissue design.

THE DESIGN PROBLEM THIS SOLVES. Comparisons in this repo are made across references that differ in
two ways at once unless care is taken. The grid as it stood:

                  immune                    gut
    mouse         MM_SLN (26 types)         HABER_GUT (15 types)
    human         -- missing --             GUT_HEALTHY_EPI (12), GUT_HEALTHY/IBD/CRC

With the human-immune cell empty, the species effect could only be estimated on gut tissue, so there
was no way to check whether it replicates elsewhere. Filling it makes the grid complete and allows
the species effect to be estimated twice, independently, on different tissue.

Source: Human Immune Health Atlas, via the CELLxGENE Census (release 2025-11-08), restricted to
normal blood and primary data. Prepared exactly as the gut references in 16_build_gut_refs.py:
250 cells per annotated cell type, types under 40 cells dropped, raw counts kept, full measured
transcriptome retained so random baselines are drawn from the transcriptome.
"""
import os, sys
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import EXT
import cellxgene_census

CENSUS_VERSION = "2025-11-08"
# Matched on dataset_id prefix, not title: the Census holds several datasets whose titles all begin
# "Human Immune Health Atlas", and one of the siblings ("Other cells") has no cells under this filter.
DATASET_PREFIX = "e522d2cd"
PER_TYPE, MIN_CELLS = 250, 40


def main():
    rng = np.random.default_rng(0)
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as c:
        ds = c["census_info"]["datasets"].read().concat().to_pandas()
        hit = ds[ds.dataset_id.str.startswith(DATASET_PREFIX)]
        if hit.empty:
            raise SystemExit(f"no census dataset with id starting {DATASET_PREFIX!r}")
        did = hit.iloc[0].dataset_id
        print(f"dataset {did}  {hit.iloc[0].dataset_title}", flush=True)
        obs = c["census_data"]["homo_sapiens"].obs.read(
            value_filter=f"dataset_id == '{did}' and is_primary_data == True "
                         f"and tissue_general == 'blood' and disease == 'normal'",
            column_names=["soma_joinid", "cell_type", "donor_id", "assay"],
        ).concat().to_pandas()
        vc = obs.cell_type.value_counts()
        obs = obs[obs.cell_type.isin(vc.index[vc >= MIN_CELLS])]
        idx = (obs.groupby("cell_type", observed=True)["soma_joinid"]
                  .apply(lambda s: pd.Series(rng.choice(s.values, min(len(s), PER_TYPE), replace=False))))
        coords = np.sort(idx.values.astype(np.int64))
        print(f"{len(coords)} cells across {obs.cell_type.nunique()} cell types", flush=True)
        a = cellxgene_census.get_anndata(
            census=c, organism="homo_sapiens", obs_coords=coords,
            column_names={"obs": ["cell_type", "donor_id", "assay"],
                          "var": ["feature_id", "feature_name"]})
    a.var_names = a.var["feature_name"].astype(str).values
    a.var_names_make_unique()
    a.obs["celltype"] = a.obs["cell_type"].astype(str)
    a.obs["batch"] = a.obs["donor_id"].astype(str)
    keep = np.asarray((a.X > 0).sum(axis=0)).ravel() >= 10
    a = a[:, keep].copy()
    a.uns["reference_title"] = "Human Immune Health Atlas (blood, normal)"
    a.uns["census_version"] = CENSUS_VERSION
    out = EXT + "expr/human_immune.h5ad"
    a.write_h5ad(out, compression="gzip")
    print(f"-> {out}  {a.shape}  {a.obs.celltype.nunique()} cell types", flush=True)


if __name__ == "__main__":
    main()
