"""Pair Driscoll et al. (2025) chamber measurements (northern Louisiana, 132 unplugged wells) with
later DOI program (Hi-Flow) pre-plugging records for the same wells.

Driscoll IDs are SONRIS well serial numbers; they are converted to API numbers via the SONRIS
GIS service (cached to data/raw/states/LA_serial_to_api.csv). The placeholder API '0000000000'
is treated as missing.

Three linkage routes, in order: (1) API number; (2) SONRIS serial number, which a few DOI LA
records use as us_well_id; (3) coordinates, accepted only within 50 m (D17 in docs/decision_log.md).
SONRIS coordinates are cached in data/raw/states/LA_serial_coords.csv. The 4 API-linked wells
sit 0-36 m from their DOI record, while sister wells on the same lease sit about 60-110 m apart,
so a wider radius would pair different wells.
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

keep = ['methane_pre_plugging_emissions', 'methane_detection_method', 'witness_date', 'name']
m1 = c.dropna(subset=['api10']).merge(dl[['api10'] + keep], on='api10', how='inner').assign(link='api')

ds = d[(d.st == 'LA') & d.us_well_id.astype(str).str.strip().str.fullmatch(r'\d{1,6}')].copy()
ds['ID'] = ds.us_well_id.astype(str).str.strip()
m2 = c[~c.ID.isin(m1.ID)].merge(ds[['ID'] + keep], on='ID', how='inner').assign(link='serial')

import numpy as np
from sklearn.neighbors import BallTree
xcache = STATES / 'LA_serial_coords.csv'
if not xcache.exists():
    rows = []
    ids = c.ID.tolist()
    for i in range(0, len(ids), 40):
        q = urllib.parse.urlencode({'where': 'WELL_SERIAL_NUM IN (' + ','.join(ids[i:i + 40]) + ')',
                                    'outFields': 'WELL_SERIAL_NUM,API_NUM,WELL_NAME,WELL_STATUS_CODE',
                                    'returnGeometry': 'true', 'outSR': '4326', 'f': 'json'})
        url = 'https://sonris-gis.dnr.la.gov/arcgis/rest/services/DNRSvc/OC/MapServer/0/query?' + q
        r = json.load(urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=60))
        rows += [dict(f['attributes'], lon=f['geometry']['x'], lat=f['geometry']['y']) for f in r['features']]
    pd.DataFrame(rows).to_csv(xcache, index=False)
xy = pd.read_csv(xcache, dtype={'WELL_SERIAL_NUM': str})
xy['ID'] = xy.WELL_SERIAL_NUM.astype(float).astype(int).astype(str)
c = c.merge(xy[['ID', 'lat', 'lon', 'WELL_NAME']].drop_duplicates('ID'), on='ID', how='left')
lag = d[(d.st == 'LA') & d.lat.notna()].reset_index(drop=True)
dist, idx = BallTree(np.radians(lag[['lat', 'lon']].values), metric='haversine').query(np.radians(c[['lat', 'lon']].values), k=1)
c['nn_m'] = dist[:, 0] * 6371000
c['nn_name'] = lag.name.values[idx[:, 0]]
rest = c[~c.ID.isin(set(m1.ID) | set(m2.ID)) & (c.nn_m <= 50)]
m3 = rest.join(lag[keep].iloc[idx[rest.index, 0]].set_index(rest.index)).assign(link='coord')
m = pd.concat([m1, m2, m3], ignore_index=True)
print('pairs by route:', m.link.value_counts().to_dict())
nl = c[c.api10.isna() & ~c.ID.isin(m.ID)]
print('placeholder-API wells unmatched: %d | nearest DOI record (m): min %.0f, median %.0f'
      % (len(nl), nl.nn_m.min(), nl.nn_m.median()))
c[['ID', 'api10', 'WELL_NAME', 'nn_m', 'nn_name']].round(1).to_csv(OUT / 'driscoll_nearest_doi.csv', index=False)

m['plug_date'] = pd.to_datetime(m.witness_date, unit='ms', errors='coerce').dt.date
m['days_between'] = (pd.to_datetime(m.plug_date) - pd.to_datetime(m.Date)).dt.days
m = m.rename(columns={'Date': 'chamber_date', 'methane_pre_plugging_emissions': 'doi_hiflow_g_h'})
cols = ['ID', 'api10', 'link', 'name', 'Product Type', 'chamber_date', 'chamber_g_h', 'plug_date', 'days_between', 'doi_hiflow_g_h', 'methane_detection_method']
m[cols].sort_values('chamber_g_h', ascending=False).to_csv(OUT / 'driscoll_doi_pairs.csv', index=False)
print('paired wells:', len(m))
print(m[cols].sort_values('chamber_g_h', ascending=False).to_string(index=False))
