#!/usr/bin/env bash
# Downloads the public resources the pipeline uses into external/ (about 400 MB).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p external/expr && cd external
git clone -q --depth 1 --branch 7.5.1 https://github.com/cran/msigdbr.git msigdbr            # MSigDB 7.5.1 gene sets
git clone -q --depth 1 https://github.com/cran/babelgene.git babelgene                        # orthologs
git clone -q --depth 1 https://github.com/cran/HGNChelper.git HGNChelper                      # HGNC symbol updates
git clone -q --depth 1 --filter=blob:none --no-checkout https://github.com/bioc/TissueEnrich.git chk_bioc_TissueEnrich
git -C chk_bioc_TissueEnrich checkout HEAD -- inst/extdata/combine-expression.rda           # HPA + GTEx tissue expression
for f in pbmc_10k_protein_v3.h5ad sln_111.h5ad sln_208.h5ad hca_subsampled_20k.h5ad; do      # public scRNA-seq
  curl -sSL --retry 3 -o "expr/$f" "https://raw.githubusercontent.com/YosefLab/scVI-data/master/$f"
done
echo "done"
