"""Repository-relative paths. Override any of them with environment variables."""
import os
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
PROJ = os.environ.get("PANEL_PROJ", str(_ROOT)) + "/"                     # repo root (data/, results/ live here)
EXT = os.environ.get("PANEL_EXT", str(_ROOT / "external")) + "/"          # public downloads (scripts/fetch_public_data.sh)
RAW = os.environ.get("PANEL_RAW", str(_ROOT / "inputs" / "vendor_panels")) + "/"
SUP = os.environ.get("PANEL_SUP", str(_ROOT / "inputs" / "foli_supplement" / "NIHMS2146370-supplement-Suppl_Table_S1-S9.xlsx"))
