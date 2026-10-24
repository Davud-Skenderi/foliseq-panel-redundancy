"""Step 8: collect every result table into one workbook."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import pandas as pd, os
D = PROJ + "data/"; O = PROJ + "results/"
out = PROJ + "results/Panel_redundancy_results.xlsx"
os.makedirs(os.path.dirname(out), exist_ok=True)
gs = pd.read_pickle(D + "genesets.pkl")
cov = pd.read_csv(O + "coverage_long.csv.gz")
lens = cov[cov.gut_lens | (cov.collection == "Hallmark")]
wide = lens.pivot_table(index=["collection", "gene_set", "set_size"], columns="panel_id", values="coverage").round(3).reset_index()
_names = [
    "panels", "panel_scope", "gene_jaccard", "pathway_scope_similarity", "gut_pathway_coverage", "top_enriched",
    "foli_on_platforms", "panel_redundancy", "panel_redundancy_by_dataset", "pathway_redundancy", "pathway_score_corr",
    "vendor_reconstruction", "optimization_curve", "optimized_panel_500", "foli_candidate_genes", "foli_stool_degs", "foli_epithelial_add_groups", "foli_mPS3_unique", "foli_mPS3_additions", "foli_mPS3_curve", "foli_hPS2_unique", "foli_hPS2_additions", "foli_hPS2_curve",
    "foli_spatial_proxies", "foli_proxy_summary", "gut_focused_panel_300", "gut_focused_curve", "vendor_gut_lens_R2", "pathway_core_summary"]
_desc = [
    "Every panel parsed from the zip: platform, species, targets, human-ortholog mapping rate (MERSCOPE human IO files were empty)",
    "Per-panel biological breadth: Hallmark coverage, % of gut/IBD-lens pathways scorable (>=5 genes), # enriched pathways",
    "Gene-level Jaccard overlap between panels (human-ortholog space)",
    "Cosine similarity of pathway-enrichment profiles (-log10 p over Hallmark+KEGG+Reactome): do two panels probe the same biology?",
    "Fraction of each Hallmark / gut-IBD-lens gene set on each panel",
    "Most enriched Hallmark / KEGG / Reactome / cell-type sets per panel (what the panel is 'about')",
    "Which spatial/targeted platforms carry each recovered Foli-seq gene",
    "Within-panel redundancy averaged over datasets vs size-matched random genes",
    "Same, per public dataset",
    "Per pathway & dataset: genes needed to reproduce the full pathway score (r>=0.9), redundancy index, core genes, panel fidelity",
    "Pathway-pathway score correlation (mean over datasets): pathways that move together",
    "Held-out reconstruction R2 of transcriptome axes + pathway scores for each vendor panel vs random and vs optimised gene count",
    "R2 vs panel size for the covariance-optimised panel and random genes",
    "Greedy covariance-optimised gene ranking (top 500)",
    "MAIN OUTPUT: every gene on mPS3/hPS2 plus candidate additions, with all evidence columns and a Keep / Review / Add / Consider call",
    "Mouse genes differentially expressed in Foli-seq stool data across colitis models (paper Tables S3-S7)",
    "Recommended gut epithelial additions grouped by epithelial cell type",
    "mPS3: genes ranked by unique contribution to gut/IBD pathway reconstruction (<=0 = redundant given the others)",
    "mPS3: top 100 covariance-ranked additions with the gut/IBD gene sets each belongs to",
    "mPS3: held-out R2 as genes are added, covariance-ranked vs random gut/IBD-lens genes",
    "hPS2: genes ranked by unique contribution",
    "hPS2: top 100 covariance-ranked additions",
    "hPS2: held-out R2 as genes are added",
    "For each spatial panel and Foli-seq gene: on panel, or best co-expressed proxy on that panel (consensus r)",
    "Per spatial panel: Foli-seq genes direct / proxy |r|>=0.6 / 0.4-0.6 / uncovered",
    "De novo covariance-optimised gut/IBD panel ranking (top 300, ribosomal/mito/housekeeping excluded)",
    "Gut/IBD objective: optimised vs random lens genes by panel size",
    "Gut/IBD objective: held-out R2 of every vendor panel",
    "Split-half validated pathway core size (genes needed to reproduce the full pathway score) and held-out r"]
print(len(_names), len(_desc))
readme = pd.DataFrame({"sheet": _names, "what it holds": _desc})
sheets = {"README": readme,
          "foli_candidate_genes": None,
          "panels": pd.read_csv(D + "panel_meta.csv"),
          "panel_scope": pd.read_csv(O + "panel_scope_summary.csv"),
          "gene_jaccard": pd.read_csv(O + "panel_jaccard.csv", index_col=0).round(3).reset_index(),
          "pathway_scope_similarity": pd.read_csv(O + "panel_pathway_similarity.csv", index_col=0).round(3).reset_index(),
          "gut_pathway_coverage": wide,
          "top_enriched": pd.read_csv(O + "top_enriched_per_panel.csv"),
          "foli_on_platforms": pd.read_csv(O + "foli_gene_platform_presence.csv"),
          "panel_redundancy": pd.read_csv(O + "panel_redundancy_summary.csv"),
          "panel_redundancy_by_dataset": pd.read_csv(O + "panel_redundancy_by_dataset.csv"),
          "pathway_redundancy": pd.read_csv(O + "pathway_redundancy_fidelity.csv"),
          "pathway_score_corr": pd.read_csv(O + "pathway_score_correlation_mean.csv", index_col=0).round(3).reset_index(),
          "vendor_reconstruction": pd.read_csv(O + "vendor_panel_reconstruction.csv"),
          "optimization_curve": pd.read_csv(O + "optimization_curve.csv"),
          "optimized_panel_500": pd.read_csv(O + "optimized_panel_500.csv"),
          "foli_candidate_genes": pd.read_csv(O + "foli_candidate_genes.csv"),
          "foli_stool_degs": pd.read_csv(D + "foli_stool_degs.csv"),
          "foli_epithelial_add_groups": pd.read_csv(O + "foli_epithelial_additions_by_group.csv"),
          "foli_mPS3_unique": pd.read_csv(O + "foli_focus_unique_contribution_FOLI_mPS3.csv"),
          "foli_mPS3_additions": pd.read_csv(O + "foli_focus_additions_FOLI_mPS3.csv"),
          "foli_mPS3_curve": pd.read_csv(O + "foli_focus_augmentation_curve_FOLI_mPS3.csv"),
          "foli_hPS2_unique": pd.read_csv(O + "foli_focus_unique_contribution_FOLI_hPS2.csv"),
          "foli_hPS2_additions": pd.read_csv(O + "foli_focus_additions_FOLI_hPS2.csv"),
          "foli_hPS2_curve": pd.read_csv(O + "foli_focus_augmentation_curve_FOLI_hPS2.csv"),
          "foli_spatial_proxies": pd.read_csv(O + "foli_spatial_proxies.csv"),
          "foli_proxy_summary": pd.read_csv(O + "foli_spatial_proxy_summary.csv"),
          "gut_focused_panel_300": pd.read_csv(O + "gut_focused_panel_300.csv"),
          "gut_focused_curve": pd.read_csv(O + "gut_focused_curve.csv"),
          "vendor_gut_lens_R2": pd.read_csv(O + "vendor_gut_lens_R2.csv"),
          "pathway_core_summary": pd.read_csv(O + "pathway_core_summary.csv")}
sheets = {k: v for k, v in sheets.items()}
with pd.ExcelWriter(out, engine="openpyxl") as w:
    for name, df in sheets.items():
        df.to_excel(w, sheet_name=name[:31], index=False)
        ws = w.sheets[name[:31]]
        ws.freeze_panes = "B2"
        for col in ws.columns:
            width = min(60, max(10, max(len(str(c.value)) if c.value is not None else 0 for c in col[:200]) + 2))
            ws.column_dimensions[col[0].column_letter].width = width
print("wrote", out)
