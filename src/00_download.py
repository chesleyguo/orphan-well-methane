"""Re-download the raw data. NOTE: the DOI layer changes over time, so a fresh pull will NOT match
the frozen 2026-11-22 snapshot. Use this only to create a new, separately dated snapshot.

Usage (from repo root): python src/00_download.py YYYY-MM-DD
Writes data/raw_<date>/... and never overwrites data/raw/.
"""
import sys, json, urllib.request, urllib.parse, pathlib

DATE = sys.argv[1] if len(sys.argv) > 1 else 'new'
ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / f'data/raw_{DATE}'
for sub in ('doi', 'usgs', 'states', 'research'):
    (OUT / sub).mkdir(parents=True, exist_ok=True)
UA = {'User-Agent': 'Mozilla/5.0 (research data download)'}

def get(url, path):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r, open(path, 'wb') as f:
        f.write(r.read())
    print('saved', path)

def arcgis_all(base, path, fields='*', geometry=True):
    rows, off = [], 0
    while True:
        q = urllib.parse.urlencode({'where': '1=1', 'outFields': fields, 'returnGeometry': str(geometry).lower(), 'outSR': 4326,
                                    'f': 'json', 'resultOffset': off, 'resultRecordCount': 2000, 'orderByFields': 'OBJECTID'})
        d = json.load(urllib.request.urlopen(urllib.request.Request(base + '/query?' + q, headers=UA), timeout=120))
        fs = d.get('features', [])
        for f in fs:
            a = f['attributes']
            if geometry:
                g = f.get('geometry') or {}; a['lon'] = g.get('x'); a['lat'] = g.get('y')
            rows.append(a)
        off += len(fs)
        if not fs or not d.get('exceededTransferLimit'):
            break
    json.dump(rows, open(path, 'w')); print('saved', path, len(rows), 'records')

# DOI Orphaned Wells Program, well-level layer behind the program dashboard (methane, method, cost)
arcgis_all('https://services.arcgis.com/4OV0eRKiLAYkbH2J/arcgis/rest/services/DOI_BIL_Orphaned_Wells_Web_Map_Joined/FeatureServer/0',
           OUT / f'doi/doi_owp_wells_{DATE}.json')
# USGS Documented Orphaned Well Database 2026 (subset of fields)
arcgis_all('https://services.arcgis.com/v01gqwM5QqNysAAi/arcgis/rest/services/US_Documented_Unplugged_Orphaned_Oil_and_Gas_Well_Dataset/FeatureServer/0',
           OUT / 'usgs/usgs_dowdb_2026_subset.json',
           fields='OBJECTID,primary_identifier_type,primary_identifier,state,county,type,status,latitude,longitude,data_file_date', geometry=False)
# USGS 2022 DOW dataset (pre-program inventory; static)
get('https://www.sciencebase.gov/catalog/file/get/62ebd67bd34eacf539724c56?f=__disk__11%2Fe9%2F27%2F11e927c652d995f46129f282b400063b5d262369',
    OUT / 'usgs/usgs_dow_2022.csv')
# State bulk files
get('https://ecmc.state.co.us/documents/data/downloads/gis/WELLS_SHP.ZIP', OUT / 'states/co_ecmc_wells_shp.zip')
get('https://data.ny.gov/api/views/szye-wmt3/rows.csv?accessType=DOWNLOAD', OUT / 'states/ny_dec_wells.csv')
# Optional, 49 MB: get('https://www.kgs.ku.edu/PRS/Ora_Archive/ks_wells.zip', OUT / 'states/ks_kgs_wells.zip')
# Driscoll et al. (2025) chamber data, Zenodo record 14827052 (concept DOI 10.5281/zenodo.13646591)
get('https://zenodo.org/api/records/14827052/files/Final_Data.xlsx/content', OUT / 'research/driscoll2025_la_chamber.xlsx')
