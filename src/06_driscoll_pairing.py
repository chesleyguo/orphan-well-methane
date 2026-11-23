"""Pair Driscoll et al. (2025) chamber measurements (northern Louisiana, 132 unplugged wells) with
later DOI program (Hi-Flow) pre-plugging records for the same wells.

Driscoll IDs are SONRIS well serial numbers; they are converted to API numbers via the SONRIS
GIS service (cached to data/raw/states/LA_serial_to_api.csv). The placeholder API '0000000000'
is treated as missing.
Run from the repo root: python src/06_driscoll_pairing.py
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *

import json, re, urllib.request, urllib.parse
import pandas as pd

PLACEHOLDER = '0000000000'
c = pd.read_excel(DRISCOLL).iloc[1:].copy()          # row 0 holds units
c['ID'] = c['ID'].astype(int).astype(str)
c['chamber_g_h'] = pd.to_numeric(c['CH4 Flow Rate'], errors='coerce')
print('Driscoll wells:', len(c), '| above control (0.0001 g/h):', int((c.chamber_g_h > 0.0001).sum()))

cache = STATES / 'LA_serial_to_api.csv'
if not cache.exists():
    rows = []
    ids = c.ID.tolist()
    for i in range(0, len(ids), 40):
        q = urllib.parse.urlencode({'where': 'WELL_SERIAL_NUM IN (' + ','.join(ids[i:i + 40]) + ')',
                                    'outFields': 'WELL_SERIAL_NUM,API_NUM', 'returnGeometry': 'false', 'f': 'json'})
        url = 'https://sonris-gis.dnr.la.gov/arcgis/rest/services/DNRSvc/OC/MapServer/0/query?' + q
        r = json.load(urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=60))
        rows += [f['attributes'] for f in r['features']]
    pd.DataFrame(rows).to_csv(cache, index=False)
s = pd.read_csv(cache, dtype=str)
s['ID'] = s.WELL_SERIAL_NUM.astype(float).astype(int).astype(str)
s['api10'] = s.API_NUM.str[:10]
c = c.merge(s[['ID', 'api10']].drop_duplicates('ID'), on='ID', how='left')
c.loc[c.api10 == PLACEHOLDER, 'api10'] = None
print('Driscoll wells with a valid API:', int(c.api10.notna().sum()))

d = pd.read_pickle(PROC / 'doi.pkl')
d['api10'] = d.us_well_id.map(lambda x: (lambda y: y[:10] if len(y) >= 10 else None)(re.sub(r'\D', '', str(x or ''))))
d.loc[d.api10 == PLACEHOLDER, 'api10'] = None
dl = d[(d.st == 'LA') & d.api10.notna()]
dl = dl[~dl.api10.duplicated(keep=False)]             # only API numbers that identify one DOI record

m = c.dropna(subset=['api10']).merge(
    dl[['api10', 'methane_pre_plugging_emissions', 'methane_detection_method', 'witness_date']], on='api10', how='inner')
m['plug_date'] = pd.to_datetime(m.witness_date, unit='ms', errors='coerce').dt.date
m['days_between'] = (pd.to_datetime(m.plug_date) - pd.to_datetime(m.Date)).dt.days
m = m.rename(columns={'Date': 'chamber_date', 'methane_pre_plugging_emissions': 'doi_hiflow_g_h'})
cols = ['ID', 'api10', 'Product Type', 'chamber_date', 'chamber_g_h', 'plug_date', 'days_between', 'doi_hiflow_g_h', 'methane_detection_method']
m[cols].sort_values('chamber_g_h', ascending=False).to_csv(OUT / 'driscoll_doi_pairs.csv', index=False)
print('paired wells:', len(m))
print(m[cols].sort_values('chamber_g_h', ascending=False).to_string(index=False))
print('unmatched for lack of API (candidates for coordinate matching):', int(c.api10.isna().sum()))
