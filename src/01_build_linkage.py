"""Build analysis dataset: USGS 2022 inventory (pre-program) linked to DOI plugged wells.
Match by API10 first, then by nearest coordinates within 100 m (same state)."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *

import json, re, numpy as np, pandas as pd
from sklearn.neighbors import BallTree

ST = {'Alabama':'AL','Alaska':'AK','Arizona':'AZ','Arkansas':'AR','California':'CA','Colorado':'CO','Florida':'FL',
      'Illinois':'IL','Indiana':'IN','Kansas':'KS','Kentucky':'KY','Louisiana':'LA','Michigan':'MI','Mississippi':'MS',
      'Missouri':'MO','Montana':'MT','Nebraska':'NE','New Mexico':'NM','New York':'NY','North Dakota':'ND','Ohio':'OH',
      'Oklahoma':'OK','Pennsylvania':'PA','Tennessee':'TN','Texas':'TX','Utah':'UT','Virginia':'VA','West Virginia':'WV',
      'Wyoming':'WY','South Dakota':'SD','Nevada':'NV','Oregon':'OR','Idaho':'ID','Maryland':'MD','New Jersey':'NJ'}
NA = {'N/A', 'NA', 'NR', 'NP', 'ND', ''}

def api10(s):
    x = re.sub(r'\D', '', str(s or ''))
    return x[:10] if len(x) >= 10 else None

u = pd.read_csv(USGS_2022, dtype=str, encoding='utf-8-sig')
u['st'] = u.State.map(ST).fillna(u.State)
u['api10'] = u['Well identifier'].map(lambda s: api10(s) if str(s).upper().startswith('API') else None)
u.loc[u.api10 == '0000000000', 'api10'] = None
u['lat'] = pd.to_numeric(u.Latitude, errors='coerce'); u['lon'] = pd.to_numeric(u.Longitude, errors='coerce')
u = u.reset_index(drop=True); u['uid'] = u.index

d = pd.DataFrame(json.load(open(DOI_WELLS)))
d['st'] = d.state.str.upper().str.strip()
d['api10'] = d.us_well_id.map(api10)
d.loc[d.api10 == '0000000000', 'api10'] = None  # placeholder API used by some LA records
d['lat'] = pd.to_numeric(d.lat, errors='coerce'); d['lon'] = pd.to_numeric(d.lon, errors='coerce')

# 1) API match
k = u.dropna(subset=['api10']).drop_duplicates('api10').set_index('api10')['uid']
d['uid'] = d.api10.map(k); d['match'] = np.where(d.uid.notna(), 'api', None)

# 2) spatial match within state, 100 m, one-to-one greedy by distance
R = 6371000.0
for st, g in d[d.uid.isna() & d.lat.notna()].groupby('st'):
    pool = u[(u.st == st) & u.lat.notna() & ~u.uid.isin(d.uid.dropna())]
    if pool.empty: continue
    tree = BallTree(np.radians(pool[['lat', 'lon']].values), metric='haversine')
    dist, idx = tree.query(np.radians(g[['lat', 'lon']].values), k=1)
    cand = pd.DataFrame({'i': g.index, 'uid': pool.uid.values[idx[:, 0]], 'm': dist[:, 0] * R})
    cand = cand[cand.m <= 100].sort_values('m').drop_duplicates('uid')
    d.loc[cand.i, 'uid'] = cand.uid.values; d.loc[cand.i, 'match'] = 'spatial'

pre = d.methane_pre_plugging_emissions
meth = d.methane_detection_method.fillna('').str.strip()
d['method_reported'] = ~meth.isin(NA)
d['measured_real'] = pre.notna() & d.method_reported
d['detect'] = pre > 0

def method_class(m):
    s = m.lower()
    if s in {x.lower() for x in NA}: return 'none'
    if any(w in s for w in ['flowmeter', 'flow meter', 'chromatograph', 'gc']): return 'flow+GC'
    if any(w in s for w in ['hi-flow', 'hiflow', 'high flow', 'hi flow', 'semtech', 'bacharach']): return 'hi-flow'
    if any(w in s for w in ['federal methane measurement', 'guideline']): return 'federal-guidelines'
    if any(w in s for w in ['flow through', 'chamber', 'flux']): return 'chamber/flow-through'
    if any(w in s for w in ['laser', 'lgr', 'mgga', 'picarro', 'tdlas', 'rmld', 'absorption']): return 'laser'
    if any(w in s for w in ['ogi', 'flir', 'gx320', 'optical gas']): return 'OGI camera'
    if any(w in s for w in ['sensit', 'bascom', 'altair', 'rki', 'gas detector', 'monitor', '4 gas', 'four gas', 'lel']): return 'handheld detector'
    if 'calc' in s or 'estimat' in s: return 'calculated'
    return 'other'
d['method_class'] = meth.map(method_class)

u['plugged'] = u.uid.isin(d.uid.dropna())
lk = d.dropna(subset=['uid']).drop_duplicates('uid').set_index('uid')
for c in ['methane_pre_plugging_emissions', 'method_class', 'measured_real', 'total_cost_complete_pa', 'Program']:
    u[c] = u.uid.map(lk[c])

d.to_pickle(PROC / 'doi.pkl'); u.to_pickle(PROC / 'inv.pkl')
print('DOI matched:', d.match.value_counts(dropna=False).to_dict())
print(d.groupby('st').match.apply(lambda x: x.notna().mean()).round(2).sort_values().to_dict())
print('inventory plugged:', u.plugged.sum(), 'of', len(u))
print('method classes among measured:', d[pre.notna()].method_class.value_counts().to_dict())
