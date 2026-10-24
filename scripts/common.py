"""Shared helpers: HGNC-canonical human symbols and mouse/dog -> human orthology."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from paths import PROJ, EXT, RAW, SUP  # noqa: E402

import warnings
import numpy as np, pandas as pd, rdata
warnings.filterwarnings("ignore")
EXT = EXT + ""
# ---- reference: symbols + orthologs (cached)
import os, pickle
_CACHE = PROJ + "data/_symbol_cache.pkl"
if os.path.exists(_CACHE):
    HUMAN_SYMBOLS, HGNC_APPROVED, H_ALIAS, M_ALIAS, m2h, d2h = pickle.load(open(_CACHE, "rb"))
else:
    msig = rdata.conversion.convert(rdata.parser.parse_file(EXT + "msigdbr/R/sysdata.rda"))
    HUMAN_SYMBOLS = set(msig["msigdbr_genes"]["human_gene_symbol"].dropna())
    bg = rdata.conversion.convert(rdata.parser.parse_file(EXT + "babelgene/R/sysdata.rda"))
    orth = bg["orthologs_df"]
    # mouse -> human: prefer MGI table, then AGR, then babelgene consensus (highest support)
    m2h = {}
    for a, b in zip(bg["mgi_orthologs_df"].mouse_symbol, bg["mgi_orthologs_df"].human_symbol):
        m2h.setdefault(a, b)
    agr = bg["agr_orthologs_df"]
    _a = agr[agr.species_name == "Mus musculus"]
    for a, b in zip(_a.species_symbol, _a.human_symbol):
        m2h.setdefault(a, b)
    mo = orth[orth.taxon_id == 10090].sort_values("support_n", ascending=False)
    for s, h in zip(mo.symbol, mo.human_symbol):
        if isinstance(s, str):
            m2h.setdefault(s, h)
    dog = orth[orth.taxon_id == 9615].sort_values("support_n", ascending=False)
    d2h = {}
    for s, h in zip(dog.symbol, dog.human_symbol):
        if isinstance(s, str):
            d2h.setdefault(s, h)
    HUMAN_UPPER = {s.upper(): s for s in HUMAN_SYMBOLS}
    # HGNChelper alias/previous-symbol tables -> current approved symbols
    _h = rdata.conversion.convert(rdata.parser.parse_file(EXT + "HGNChelper/data/hgnc.table.rda"))["hgnc.table"]
    HGNC_APPROVED = set(_h["Approved.Symbol"].dropna())
    _h = _h[_h.Symbol.notna() & _h["Approved.Symbol"].notna()]
    # ambiguous aliases (one alias -> several approved) are dropped unless the alias is itself approved
    _amb = _h.groupby("Symbol")["Approved.Symbol"].nunique().to_dict()
    H_ALIAS = {a: b for a, b in zip(_h.Symbol, _h["Approved.Symbol"]) if _amb[a] == 1 or a in HGNC_APPROVED}
    _m = rdata.conversion.convert(rdata.parser.parse_file(EXT + "HGNChelper/data/mouse.table.rda"))["mouse.table"]
    _m = _m[_m.Symbol.notna() & _m["Approved.Symbol"].notna()]
    _mamb = _m.groupby("Symbol")["Approved.Symbol"].nunique().to_dict()
    M_ALIAS = {a: b for a, b in zip(_m.Symbol, _m["Approved.Symbol"]) if _mamb[a] == 1}



    pickle.dump((HUMAN_SYMBOLS, HGNC_APPROVED, H_ALIAS, M_ALIAS, m2h, d2h), open(_CACHE, "wb"))

def canon_h(sym):
    """Return the current HGNC-approved human symbol (or None)."""
    if not isinstance(sym, str):
        return None
    sym = sym.strip()
    if sym in HGNC_APPROVED:
        return sym
    if sym in H_ALIAS:
        return H_ALIAS[sym]
    u = sym.upper()
    if u in HGNC_APPROVED:
        return u
    return H_ALIAS.get(u)


def to_human(sym, species, provided_human=None):
    if sym is None or (isinstance(sym, float) and np.isnan(sym)):
        return None
    sym = str(sym).strip()
    if provided_human and isinstance(provided_human, str) and provided_human.strip():
        h = canon_h(provided_human)
        if h:
            return h
    if species == "Human":
        return canon_h(sym)
    if species == "Mouse":
        for cand in (sym, M_ALIAS.get(sym)):
            if cand and cand in m2h:
                return canon_h(m2h[cand]) or m2h[cand]
        return canon_h(sym)  # many mouse symbols == human symbol capitalised
    if species == "Dog":
        if sym in d2h:
            return canon_h(d2h[sym]) or d2h[sym]
        return canon_h(sym)
    return None





def top_pcs(Xc, k):
    """Top-k right singular vectors of a centred (samples x features) matrix via the samples Gram matrix.
    Fast when samples << features (metacells x genes)."""
    import numpy as _np
    G = Xc @ Xc.T
    w, U = _np.linalg.eigh(G)
    order = _np.argsort(w)[::-1][:k]
    w, U = _np.clip(w[order], 1e-12, None), U[:, order]
    return (Xc.T @ U / _np.sqrt(w)).T  # k x features, rows orthonormal
