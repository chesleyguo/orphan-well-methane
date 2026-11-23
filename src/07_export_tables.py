"""Export analysis-ready CSVs from the processed pickles.
Run from the repo root after 01–02: python src/07_export_tables.py
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
import pandas as pd

e = pd.read_pickle(PROC / 'doi_enriched.pkl')
for c in [c for c in e.columns if 'date' in c and e[c].dtype != 'datetime64[ns]']:
    if pd.api.types.is_numeric_dtype(e[c]):
        e[c] = pd.to_datetime(e[c], unit='ms', errors='coerce').dt.date
keep_first = ['OBJECTID', 'us_well_id', 'api10', 'st', 'county', 'name', 'well_type', 'lat', 'lon', 'witness_date',
              'methane_detection_method', 'method_class', 'method_reported', 'measured_real',
              'methane_pre_plugging_emissions', 'methane_post_plugging_emissions', 'total_cost_complete_pa',
              'uid', 'match', 'spud', 'depth', 'last_prod', 'age_at_plug', 'src', 'Program', 'Recipient']
cols = [c for c in keep_first if c in e.columns] + [c for c in e.columns if c not in keep_first]
e[cols].to_csv(PROC / 'doi_wells_analysis.csv', index=False)

u = pd.read_pickle(PROC / 'inv.pkl')
u.drop(columns=[c for c in ['Latitude', 'Longitude'] if c in u.columns]).to_csv(PROC / 'usgs2022_inventory_linked.csv', index=False)
print('wrote', PROC / 'doi_wells_analysis.csv', e.shape, '|', PROC / 'usgs2022_inventory_linked.csv', u.shape)
