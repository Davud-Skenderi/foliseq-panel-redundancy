"""Q3 + Q4 recomputed on the gut references built by 16_build_gut_refs.py.

Same Spapros metrics as 12_spapros_eval.py, run against annotated human intestinal tissue instead of
immune tissue, so that the panel comparison can be tested for robustness to the choice of reference.

THREE INGREDIENTS, held apart deliberately:
  (1) panels      the gene lists being compared        -- fixed across runs
  (2) metrics     accuracy / redundancy definitions    -- fixed across runs
  (3) reference   the annotated single-cell data       -- VARIED, one run per reference
Varying only (3) is what makes the comparison a robustness test: absolute scores are expected to move
(a reference with more, finer cell types is harder, so the whole distribution shifts left), but the
ORDERING of panels should survive if it reflects the panels rather than the reference.

NOTE ON A BUG IN 12_spapros_eval.py: its confusion-matrix cache key is `q4_confusion_{panel}.csv`
with no dataset in the name, so running a second reference through that script would silently reload
the first reference's results. Cache keys here are `q4_confusion_{REF}__{panel}.csv`.

Outputs per reference REF:
  results/v2/q4_accuracy_by_celltype_{REF}.csv    cell type x panel accuracy (the histogram input)
  results/v2/q4_celltype_identification_{REF}.csv per-panel summary
  results/v2/q3_redundancy_summary_{REF}.csv      redundancy, with size-normalised rates
  results/v2/q3_correlated_pairs_{REF}.csv        pairs with |r| >= 0.8
"""
import os, sys, time, warnings
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PROJ, EXT
from common import to_human
from spapros.evaluation import metrics as spm
warnings.filterwarnings("ignore")

D, O = PROJ + "data/", PROJ + "results/v2/"
os.makedirs(O, exist_ok=True)
M = pd.read_csv(D + "gene_by_panel.csv", index_col=0)
PANELS = ["FOLI_mPS3", "FOLI_hPS2", "COSMX_HIO", "MER_M_BRAINTRAIN", "MER_M_PAN", "MER_M_PANNEURO",
          "COSMX_MNEURO", "COSMX_MUCC", "COSMX_HUCC", "GEOMX_H_CTA", "GEOMX_DOG_CCA",
          "XEN_M5K", "XEN_H5K", "COSMX_H6K",
          # whole-transcriptome panels: not competitors but a CEILING -- what the platform achieves
          # measuring essentially everything, so targeted panels can be read as a fraction of it
          "COSMX_H_WTX", "COSMX_M_WTX"]
RANDOM_SIZES = [330, 535]      # 330 = Foli-seq mPS3 measured size; 535 = baseline used in v2
CORR_MAX_GENES = 8000          # above this the n x n correlation matrix exhausts memory; see run()
N_SEEDS = int(os.environ.get("N_SEEDS", "1"))   # Spapros default is 5; v2 used 1 for compute.
                                                # With 1 seed no accuracy carries a spread, so a 2-point
                                                # difference cannot be told from run-to-run noise.
N_RANDOM = 3


def to_symbols(a, species):
    """Collapse a reference onto human symbols. Identical to 12_spapros_eval.to_symbols, so that a
    mouse reference added here is mapped exactly as MM_SLN was and the two stay comparable."""
    hs = [to_human(g, species) for g in a.var_names.astype(str)]
    keep = np.array([isinstance(h, str) and len(h) > 0 for h in hs])
    a = a[:, keep].copy()
    h = np.array(hs, dtype=object)[keep]
    codes, uniq = pd.factorize(h)
    Mx = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
    return ad.AnnData(sp.csr_matrix(a.X) @ Mx, obs=a.obs.copy(), var=pd.DataFrame(index=uniq))


def lognorm(a):
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    return a


def sample_per_type(a, key, n=200, min_cells=40, seed=0):
    rng = np.random.default_rng(seed)
    vc = a.obs[key].value_counts()
    idx = []
    for t in vc.index[vc >= min_cells]:
        ii = np.where(a.obs[key].values == t)[0]
        idx += list(rng.choice(ii, min(n, len(ii)), replace=False))
    return a[sorted(idx)].copy()


def corr_summary(adata, genes, tag, panel):
    cm = spm.correlation_matrix(adata, var_names=genes, verbosity=0)
    v = np.array(cm.fillna(0).values, dtype=float, copy=True)
    np.fill_diagonal(v, 0)
    absv, names = np.abs(v), list(cm.index)
    pairs = [(names[i], names[j], float(v[i, j])) for i, j in zip(*np.where(np.triu(absv, 1) >= 0.8))]
    n = len(genes)
    n_possible = n * (n - 1) // 2
    return dict(
        dataset=tag, panel_id=panel, n_genes_measured=n,
        pct_genes_nonredundant=float(np.mean(absv.max(axis=0) < 0.8)) * 100,
        # the size-normalised quantities the raw count cannot provide
        pct_genes_redundant=float(np.mean(absv.max(axis=0) >= 0.8)) * 100,
        mean_abs_corr=float(absv[np.triu_indices(n, 1)].mean()) if n > 1 else np.nan,
        n_pairs_r_ge_0_8=len(pairs),
        n_possible_pairs=n_possible,
        redundant_pair_rate_per_1k=1000.0 * len(pairs) / n_possible if n_possible else np.nan,
    ), pairs


def load_mm_sln():
    """Same preparation as 12_spapros_eval.py: concat the two donors, drop doublets and low quality."""
    a1 = ad.read_h5ad(EXT + "expr/sln_111.h5ad"); a2 = ad.read_h5ad(EXT + "expr/sln_208.h5ad")
    a = ad.concat([a1, a2], join="inner"); a.obs_names_make_unique()
    a = a[~a.obs.cell_types.str.contains("doublet|Low quality|Erythrocytes", case=False)].copy()
    a.obs["celltype"] = a.obs.cell_types.astype(str)
    a.uns["reference_title"] = "Mouse spleen / lymph node CITE-seq"
    return a


def run(ref_tag, spec):
    h5ad_name, species = spec
    a = load_mm_sln() if h5ad_name == "__MM_SLN__" else ad.read_h5ad(EXT + f"expr/{h5ad_name}")
    title = a.uns.get("reference_title", ref_tag)
    if species == "Mouse":
        n0 = a.n_vars
        a = to_symbols(a, "Mouse")
        print(f"  {ref_tag}: {n0} mouse genes -> {a.n_vars} human symbols", flush=True)
    sc.pp.filter_genes(a, min_cells=10)
    a = lognorm(a)
    clf_ad = sample_per_type(a, "celltype")
    n_corr = min(8000, a.n_obs)
    corr_ad = a[np.random.default_rng(0).choice(a.n_obs, n_corr, replace=False)].copy()
    measured = set(a.var_names)
    print(f"\n=== {ref_tag} ({title}) {a.shape}, clf {clf_ad.shape}, "
          f"{clf_ad.obs.celltype.nunique()} cell types", flush=True)

    sets = {p: [g for g in M.index[M[p] == 1] if g in measured] for p in PANELS}
    rng = np.random.default_rng(0)
    for size in RANDOM_SIZES:
        for k in range(N_RANDOM):
            sets[f"RANDOM_{size}_{k}"] = list(rng.choice(a.var_names, size, replace=False))

    rows, per_ct, red, allpairs = [], [], [], []
    for p, genes in sets.items():
        if len(genes) < 10:
            print(f"  {p}: only {len(genes)} genes measured, skipped", flush=True); continue
        t0 = time.time()
        # cache key carries the seed count as well as the dataset: a 5-seed result is a different
        # quantity from a 1-seed result and must not silently overwrite or reload it
        suffix = "" if N_SEEDS == 1 else f"_s{N_SEEDS}"
        cached = O + f"q4_confusion{suffix}_{ref_tag}__{p}.csv"
        if os.path.exists(cached):
            cm = pd.read_csv(cached, index_col=0)
        else:
            out = spm.xgboost_forest_classification(clf_ad, genes, celltypes="all", ct_key="celltype",
                                                    n_cells_min=40, n_seeds=N_SEEDS, n_jobs=2,
                                                    verbosity=0)
            # out[0] is the confusion matrix averaged over CV folds and seeds; out[1] is its standard
            # deviation over the same. 12_spapros_eval.py discarded out[1], which is why no number in
            # this repo had an error bar. Keep it.
            cm = out[0]
            cm.to_csv(cached)
            if len(out) > 1 and out[1] is not None:
                out[1].to_csv(cached.replace("q4_confusion", "q4_confusion_sd"))
        diag = pd.Series(np.diag(cm.values), index=cm.index)
        rows.append(dict(reference=ref_tag, panel_id=p, n_genes_measured=len(genes),
                         n_celltypes=len(diag),
                         mean_celltype_accuracy=float(diag.mean()),
                         median_celltype_accuracy=float(diag.median()),
                         sd_celltype_accuracy=float(diag.std()),
                         min_celltype_accuracy=float(diag.min()),
                         pct_celltypes_acc_gt_0_8=float(np.mean(diag > 0.8)) * 100))
        per_ct += [dict(panel_id=p, celltype=c, accuracy=float(v)) for c, v in diag.items()]
        # The redundancy step builds an n x n correlation matrix. At whole-transcriptome size that is
        # tens of millions of float64 entries plus copies, which exhausts memory and kills the process
        # AFTER the classification result has been cached -- silently, with no error in the log.
        # Redundancy is not a meaningful question for a whole-transcriptome panel anyway: it is not a
        # curated selection, so "what fraction of its genes are redundant" has no design interpretation.
        if not p.startswith("RANDOM") and len(genes) <= CORR_MAX_GENES:
            s, pairs = corr_summary(corr_ad, genes, ref_tag, p)
            red.append(s)
            allpairs += [dict(dataset=ref_tag, panel_id=p, gene_a=x, gene_b=y, r=r) for x, y, r in pairs]
        print(f"  {p}: {len(genes)} genes, mean acc {diag.mean():.3f}, {time.time()-t0:.0f}s", flush=True)

    sfx = "" if N_SEEDS == 1 else f"_s{N_SEEDS}"
    pd.DataFrame(rows).to_csv(O + f"q4_celltype_identification{sfx}_{ref_tag}.csv", index=False)
    (pd.DataFrame(per_ct).pivot(index="celltype", columns="panel_id", values="accuracy")
       .to_csv(O + f"q4_accuracy_by_celltype{sfx}_{ref_tag}.csv"))
    if red:
        pd.DataFrame(red).to_csv(O + f"q3_redundancy_summary_{ref_tag}.csv", index=False)
        pd.DataFrame(allpairs).to_csv(O + f"q3_correlated_pairs_{ref_tag}.csv", index=False)
    print(f"=== {ref_tag} done", flush=True)


REFS = {"GUT_HEALTHY": ("gut_healthy.h5ad", "Human"),
        "GUT_IBD":     ("gut_ibd.h5ad", "Human"),
        "GUT_CRC":     ("gut_crc.h5ad", "Human"),
        "HABER_GUT":   ("haber_gut.h5ad", "Mouse"),
        # epithelium-only human reference, compartment-matched to HABER_GUT (see 21_)
        "GUT_HEALTHY_EPI": ("gut_healthy_epi.h5ad", "Human"),
        # fills the empty human-immune corner of the species x tissue grid (see 28_)
        "HUMAN_IMMUNE": ("human_immune.h5ad", "Human"),
        # MM_SLN is two files that must be concatenated and filtered, so it is loaded specially
        # below rather than from a single h5ad. Registered here so it can be re-run at 5 seeds and
        # the whole grid share one seed setting instead of mixing 1-seed and 5-seed numbers.
        "MM_SLN": ("__MM_SLN__", "Mouse")}

if __name__ == "__main__":
    # usage: 17_multiref_eval.py REF [REF ...] [--panels A,B,C]
    # --panels restricts the run to a subset, so a shared panel set can be completed across every
    # reference before the expensive 5K/6K panels are attempted on any one of them. Results are
    # cached per (reference, panel), so later runs extend rather than repeat.
    args = sys.argv[1:]
    subset = None
    if "--panels" in args:
        i = args.index("--panels")
        subset = [x for x in args[i + 1].split(",") if x]
        args = args[:i]
    if subset:
        PANELS = [p for p in PANELS if p in subset]
        print(f"panel subset: {len(PANELS)} panels", flush=True)
    for tag in (args or list(REFS)):
        run(tag, REFS[tag])
