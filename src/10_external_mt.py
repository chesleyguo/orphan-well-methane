"""USGS Kevin-Sunburst Dome (Toole County, MT) orphan-well methane measurements, linked to the
2022 USGS orphan inventory and to the nearest DOI plugged-well record.

Source: USGS data release doi:10.5066/P13A7RXN (ScienceBase item 677bfa33d34e711d10b7c891),
cached in data/raw/usgs_mt_kevin_sunburst/. Note: in KD_OW_Geochemistry_Data.csv the columns
Decimal_Longitude and Decimal_Latitude are swapped (the metadata definitions are swapped too);
latitude ~48.8 is in the column named Decimal_Longitude.
Run from the repo root: python src/10_external_mt.py
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *

import re
import numpy as np, pandas as pd
from sklearn.neighbors import BallTree

src = RAW / 'usgs_mt_kevin_sunburst'
e = pd.read_csv(src / 'KD_OW_Emissions_Data.csv', encoding='utf-8-sig')
g = pd.read_csv(src / 'KD_OW_Geochemistry_Data.csv', encoding='utf-8-sig')
g = g.groupby('API').agg(lat=('Decimal_Longitude', 'first'), lon=('Decimal_Latitude', 'first')).reset_index()

w = e.groupby(['API', 'Well_Name']).agg(
    sessions=('Average_g_h', 'size'), total_hours=('Duration_h', 'sum'),
    mean_of_sessions_g_h=('Average_g_h', 'mean'), min_session_g_h=('Average_g_h', 'min'),
    max_session_g_h=('Average_g_h', 'max')).reset_index()
w['cv_between_sessions'] = e.groupby('API').Average_g_h.std().reindex(w.API).values / w.mean_of_sessions_g_h
w = w.merge(g, on='API', how='left')
w['api10'] = w.API.str.replace('-', '')

inv = pd.read_pickle(PROC / 'inv.pkl')
inv_api = set(inv['Well identifier'].astype(str).str.replace(r'\D', '', regex=True).str[:10])
w['in_usgs_2022_inventory'] = w.api10.isin(inv_api)

d = pd.read_pickle(PROC / 'doi.pkl')
d['a'] = d.us_well_id.map(lambda x: re.sub(r'\D', '', str(x or ''))[:10])
w['in_doi_plugged'] = w.api10.isin(set(d.a))
dl = d[d.lat.notna()].reset_index(drop=True)
dist, idx = BallTree(np.radians(dl[['lat', 'lon']].values), metric='haversine').query(np.radians(w[['lat', 'lon']].values), k=1)
w['nearest_doi_m'] = (dist[:, 0] * 6371000).round(0)
w['nearest_doi_value_g_h'] = dl.methane_pre_plugging_emissions.values[idx[:, 0]]
w['nearest_doi_method'] = dl.methane_detection_method.values[idx[:, 0]]

w = w.drop(columns='api10').sort_values('mean_of_sessions_g_h', ascending=False)
w.round(3).to_csv(OUT / 'usgs_mt_wells.csv', index=False)
print('USGS MT: sessions %d, wells with emission data %d, wells with coordinates %d' % (len(e), e.API.nunique(), len(g)))
print(w[['API', 'Well_Name', 'sessions', 'mean_of_sessions_g_h', 'min_session_g_h', 'max_session_g_h',
         'in_usgs_2022_inventory', 'in_doi_plugged', 'nearest_doi_m']].round(1).to_string(index=False))
print('in 2022 inventory: %d/%d | in DOI plugged records: %d | nearest DOI record: min %.0f m'
      % (w.in_usgs_2022_inventory.sum(), len(w), w.in_doi_plugged.sum(), w.nearest_doi_m.min()))
print('median between-session CV (wells with >=3 sessions): %.2f'
      % w.loc[w.sessions >= 3, 'cv_between_sessions'].median())
