"""Link DOI plugged orphaned wells to state well databases by API10; report attribute coverage.

LA, NM, CA, WV are queried live from state ArcGIS services the first time and cached to
data/raw/states/<ST>_query.csv; later runs read the cache, so results stay tied to the snapshot.
CO and NY come from bulk downloads in data/raw/states/. KS (49 MB, ks_kgs_wells.zip) is optional:
download https://www.kgs.ku.edu/PRS/Ora_Archive/ks_wells.zip into data/raw/states/ to include it.
Run from the repo root: python src/02_state_attributes.py
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
import json, re, zipfile, urllib.request, urllib.parse, numpy as np, pandas as pd
from dbfread import DBF

S = str(STATES) + '/'
d = pd.read_pickle(PROC / 'doi.pkl')
d['api10'] = d.us_well_id.map(lambda s: (lambda x: x[:10] if len(x) >= 10 else None)(re.sub(r'\D', '', str(s or ''))))

def arcgis(base, field, vals, fmt, outf, cache):
    cache = STATES / cache
    if cache.exists():
        return pd.read_csv(cache, dtype=str)
    rows = []
    vals = [v for v in vals if v]
    for i in range(0, len(vals), 40):
        chunk = [fmt(v) for v in vals[i:i+40]]
        where = f"{field} IN (" + ",".join("'" + c + "'" for c in chunk) + ")"
        q = urllib.parse.urlencode({'where': where, 'outFields': outf, 'returnGeometry': 'false', 'f': 'json'})
        r = json.load(urllib.request.urlopen(urllib.request.Request(base + '/query?' + q, headers={'User-Agent': 'Mozilla/5.0 (research data check)'}), timeout=90))
        rows += [f['attributes'] for f in r.get('features', [])]
    df = pd.DataFrame(rows); df.to_csv(cache, index=False)  # epoch-ms dates are kept as numbers
    return pd.read_csv(cache, dtype=str)

def ms2date(x):
    return pd.to_datetime(pd.to_numeric(x, errors='coerce'), unit='ms', errors='coerce')

out = {}
def pick(st):
    g = d[(d.st == st) & d.api10.notna()]
    return g, sorted(set(g.api10))

# Louisiana (SONRIS)
g, apis = pick('LA')
la = arcgis('https://sonris-gis.dnr.la.gov/arcgis/rest/services/DNRSvc/OC/MapServer/0', 'API_NUM', apis,
            lambda a: a + '0000', 'API_NUM,SPUD_DATE,MEASURED_DEPTH,TRUE_VERTICAL_DEPTH,PRODUCT_TYPE_CODE,ORIGINAL_COMPLETION_DATE,LAST_TEST_DATE,FIELD_NAME', 'LA_query.csv')
la['api10'] = la.API_NUM.astype(str).str[:10]; la['spud'] = ms2date(la.SPUD_DATE)
la['depth'] = pd.to_numeric(la.TRUE_VERTICAL_DEPTH, errors='coerce').where(lambda x: x > 0, pd.to_numeric(la.MEASURED_DEPTH, errors='coerce'))
la['last_prod'] = ms2date(la.LAST_TEST_DATE)
out['LA'] = la.drop_duplicates('api10')[['api10', 'spud', 'depth', 'last_prod']]

# New Mexico (OCD)
g, apis = pick('NM')
nm = arcgis('https://gis.emnrd.nm.gov/arcgis/rest/services/OCDView/Wells_Public/FeatureServer/0', 'id', apis,
            lambda a: f'{a[:2]}-{a[2:5]}-{a[5:10]}', 'id,spud_date,true_vertical_depth,measured_vertical_depth,last_production_date', 'NM_query.csv')
nm['api10'] = nm['id'].str.replace('-', '').str[:10]; nm['spud'] = ms2date(nm.spud_date)
nm['depth'] = pd.to_numeric(nm.true_vertical_depth, errors='coerce').where(lambda x: x > 0, pd.to_numeric(nm.measured_vertical_depth, errors='coerce'))
nm['last_prod'] = ms2date(nm.last_production_date)
out['NM'] = nm.drop_duplicates('api10')[['api10', 'spud', 'depth', 'last_prod']]

# California (WellSTAR) — spud only
g, apis = pick('CA')
ca = arcgis('https://gis.conservation.ca.gov/server/rest/services/WellSTAR/Wells/MapServer/0', 'API', apis, lambda a: a, 'API,SpudDate', 'CA_query.csv')
ca['api10'] = ca.API.astype(str).str.zfill(10).str[:10]; ca['spud'] = pd.to_datetime(ca.SpudDate, errors='coerce'); ca['depth'] = np.nan; ca['last_prod'] = pd.NaT
out['CA'] = ca.drop_duplicates('api10')[['api10', 'spud', 'depth', 'last_prod']]

# West Virginia (TAGIS)
g, apis = pick('WV')
wv = arcgis('https://tagis.dep.wv.gov/arcgis/rest/services/WVDEP_enterprise/oil_gas/MapServer/7', 'api', apis, lambda a: a, 'api,welldepth,compdate,issuedate', 'WV_query.csv')
wv['api10'] = wv.api.astype(str).str[:10]; wv['depth'] = pd.to_numeric(wv.welldepth, errors='coerce')
# compdate/issuedate in this layer are dates of the current (plugging) permit, not spud dates -> no age for WV
wv['spud'] = pd.NaT; wv['last_prod'] = pd.NaT
out['WV'] = wv.drop_duplicates('api10')[['api10', 'spud', 'depth', 'last_prod']]

# Colorado (ECMC shapefile dbf)
z = zipfile.ZipFile(S + 'co_ecmc_wells_shp.zip'); z.extract('Wells.dbf', PROC)
co = pd.DataFrame(iter(DBF(PROC / 'Wells.dbf', encoding='latin1', char_decode_errors='ignore')))
co['api10'] = co.API_Label.astype(str).str.replace('-', '').str[:10]
co['spud'] = pd.to_datetime(co.Spud_Date, errors='coerce')
co['depth'] = pd.to_numeric(co.Max_TVD, errors='coerce').where(lambda x: x > 0, pd.to_numeric(co.Max_MD, errors='coerce'))
co['last_prod'] = pd.NaT
out['CO'] = co.drop_duplicates('api10')[['api10', 'spud', 'depth', 'last_prod']]

# New York (DEC, data.ny.gov)
ny = pd.read_csv(S + 'ny_dec_wells.csv', dtype=str)
ny['api10'] = ny['API Well Number'].str[:10]; ny['spud'] = pd.to_datetime(ny['Date Spudded'], errors='coerce')
ny['depth'] = pd.to_numeric(ny['True Vertical Depth, ft'], errors='coerce').where(lambda x: x > 0, pd.to_numeric(ny['Measured Depth, ft'], errors='coerce'))
ny['last_prod'] = pd.NaT
out['NY'] = ny.drop_duplicates('api10')[['api10', 'spud', 'depth', 'last_prod']]

# Kansas (KGS) — optional; KS has no real methane measurements, so it only matters for selection models
if (STATES / 'ks_kgs_wells.zip').exists():
    ks = pd.read_csv(zipfile.ZipFile(S + 'ks_kgs_wells.zip').open('ks_wells.txt'), dtype=str, encoding='latin1')
    ks['api10'] = ks.API_NUM_NODASH.str[:10]; ks['spud'] = pd.to_datetime(ks.SPUD, format='%d-%b-%Y', errors='coerce')
    ks['depth'] = pd.to_numeric(ks.DEPTH, errors='coerce'); ks['last_prod'] = pd.NaT
    out['KS'] = ks.dropna(subset=['api10']).drop_duplicates('api10')[['api10', 'spud', 'depth', 'last_prod']]

def clean(v):
    v = v.copy()
    for c in ['spud', 'last_prod']:
        v[c] = pd.to_datetime(v[c], errors='coerce')
        v[c] = v[c].where((v[c] > pd.Timestamp('1850-01-01')) & (v[c] < pd.Timestamp('2026-12-31')))
        v[c] = v[c].astype('datetime64[ns]')
    return v
allattr = pd.concat([clean(v).assign(src=k) for k, v in out.items()])
allattr = allattr[allattr.api10.notna() & (allattr.api10 != '0000000000')].drop_duplicates('api10')
d.loc[d.api10 == '0000000000', 'api10'] = None  # placeholder API used by some LA records
e = d.merge(allattr, on='api10', how='left')
e.loc[e.depth <= 0, 'depth'] = np.nan
e['age_at_plug'] = (pd.to_datetime(e.witness_date, unit='ms', errors='coerce') - e.spud).dt.days / 365.25
e.loc[(e.age_at_plug < 0) | (e.age_at_plug > 170), 'age_at_plug'] = np.nan
e.to_pickle(PROC / 'doi_enriched.pkl')
e.to_csv(PROC / 'doi_wells_enriched.csv', index=False)

rows = []
for st in ['LA', 'NM', 'CO', 'CA', 'WV', 'NY', 'KS']:
    g = e[e.st == st]; m = g[g.measured_real & (g.method_class != 'handheld detector')]
    rows.append(dict(state=st, doi_wells=len(g), matched=g.src.notna().sum(), spud=g.spud.notna().sum(), depth=g.depth.notna().sum(),
                     last_prod=g.last_prod.notna().sum(), measured_real=len(m), measured_with_age=m.age_at_plug.notna().sum(),
                     measured_with_depth=m.depth.notna().sum(), detects=int((m.methane_pre_plugging_emissions > 0).sum())))
print(pd.DataFrame(rows).to_string(index=False))
