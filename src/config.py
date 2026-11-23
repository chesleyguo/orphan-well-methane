"""Shared paths. Run every script from the repository root: python src/<script>.py"""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data' / 'raw'
PROC = ROOT / 'data' / 'processed'
OUT = ROOT / 'outputs'
SNAPSHOT = '2026-11-22'
DOI_WELLS = RAW / 'doi' / f'doi_owp_wells_{SNAPSHOT}.json'
USGS_2022 = RAW / 'usgs' / 'usgs_dow_2022.csv'
USGS_2026 = RAW / 'usgs' / 'usgs_dowdb_2026_subset.json'
STATES = RAW / 'states'
DRISCOLL = RAW / 'research' / 'driscoll2025_la_chamber.xlsx'
for p in (PROC, OUT):
    p.mkdir(parents=True, exist_ok=True)
