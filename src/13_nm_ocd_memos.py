"""New Mexico OCD well files: contractor methane memos for wells plugged under the program.

For each NM DOI record with an API number: list the well's documents on the OCD imaging site,
download PDFs filed from 2022-06 on, keep only their extracted text (PDFs are deleted), and parse
whether the pre-plugging value was measured or a field average, the reported rate, and the
instrument. Resumable: texts are cached in data/raw/states/nm_ocd_memos.jsonl (one line per PDF).
Polite: one request at a time with a pause.

Usage (from repo root): python src/13_nm_ocd_memos.py [--limit N] [--minutes M]
"""
import sys, pathlib, time, json, re, subprocess, tempfile, argparse, urllib.request
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--limit', type=int, default=0)
ap.add_argument('--minutes', type=float, default=9.0)
ap.add_argument('--parse-only', action='store_true')
args = ap.parse_args()

CACHE = STATES / 'nm_ocd_memos.jsonl'
LISTS = STATES / 'nm_ocd_listings.json'
UA = {'User-Agent': 'Mozilla/5.0 (research; orphan-well-methane)'}
PAUSE = 0.5
t_end = time.time() + args.minutes * 60


def get(url, timeout=60):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout).read()


d = pd.read_pickle(PROC / 'doi.pkl')
nm = d[(d.st == 'NM')].copy()
nm['api14'] = nm.us_well_id.astype(str).str.replace(r'\D', '', regex=True).str.ljust(14, '0').str[:14]
nm = nm[nm.api14.str.startswith('30') & (nm.api14 != '30000000000000')]
apis = sorted(nm.api14.unique())
if args.limit:
    apis = apis[:args.limit]

listings = json.load(open(LISTS)) if LISTS.exists() else {}
done = set()
if CACHE.exists():
    done = {json.loads(l)['url'] for l in open(CACHE)}

if not args.parse_only:
    for api in apis:
        if time.time() > t_end:
            break
        if api not in listings:
            try:
                h = get('https://ocdimage.emnrd.nm.gov/imaging/WellFileView.aspx?RefType=WF&RefID=' + api).decode('utf-8', 'ignore')
            except Exception as ex:
                print('listing failed', api, ex); continue
            urls = sorted(set(re.findall(r'https://ocdimage\.emnrd\.nm\.gov/Imaging/FileStore/santafe/wf/(20\d{6})/(\d{14}_[\d_]+\.pdf)', h, re.I)))
            listings[api] = ['https://ocdimage.emnrd.nm.gov/Imaging/FileStore/santafe/wf/%s/%s' % u for u in urls if u[0] >= '20220601']
            json.dump(listings, open(LISTS, 'w'))
            time.sleep(PAUSE)
        for url in listings[api]:
            if url in done or time.time() > t_end:
                continue
            try:
                pdf = get(url, 120)
            except Exception as ex:
                print('pdf failed', url, ex); continue
            with tempfile.TemporaryDirectory() as td:
                pth = pathlib.Path(td) / 'f.pdf'
                pth.write_bytes(pdf)
                txt = subprocess.run(['pdftotext', '-layout', str(pth), '-'], capture_output=True, text=True, timeout=120).stdout
            with open(CACHE, 'a') as f:
                f.write(json.dumps({'api': api, 'url': url, 'chars': len(txt), 'text': txt}) + '\n')
            done.add(url)
            time.sleep(PAUSE)
    print('wells listed: %d of %d | PDFs cached: %d' % (sum(a in listings for a in apis), len(apis), len(done)))

# ---- parse -------------------------------------------------------------------------------
recs = [json.loads(l) for l in open(CACHE)] if CACHE.exists() else []
num = r'([0-9]+(?:\.[0-9]+)?)'
rows = []
for r in recs:
    t = ' '.join(r['text'].split())
    if not re.search(r'methane', t, re.I):
        continue
    pre = bool(re.search(r'pre[- ]plug', t, re.I))
    avg = bool(re.search(r'averaged and applied|field averaging|were not measured|average[d]? .{0,40}applied to', t, re.I))
    m = (re.search(r'emission rate in grams per hour \(g/hr\)\s*' + num, t, re.I)
         or re.search(r'Average CH4 Mass\s*' + num + r'\s*g/hr', t, re.I)
         or re.search(r'calculated at\s*' + num + r'\s*grams per hour', t, re.I)
         or re.search(r'Calculated Average Wellhead Methane Gas Emissions:\s*' + num, t, re.I))
    dur = re.search(r'Flow Duration\s*([0-9]+)\s*hrs?', t, re.I)
    rows.append({'api': r['api'], 'url': r['url'], 'pre_plugging_memo': pre, 'field_average': avg,
                 'reported_g_h': float(m.group(1)) if m else None,
                 'flow_duration_h': float(dur.group(1)) if dur else None,
                 'ventbuster': bool(re.search(r'Ventbuster', t, re.I)),
                 'contractor_wdf': bool(re.search(r'Well Done', t, re.I))})
memo = pd.DataFrame(rows)
if len(memo):
    per = memo.groupby('api').agg(memos=('url', 'size'), any_pre=('pre_plugging_memo', 'any'),
                                  field_average=('field_average', 'any'), ventbuster=('ventbuster', 'any'),
                                  memo_g_h=('reported_g_h', lambda s: s.dropna().iloc[0] if s.notna().any() else None),
                                  flow_duration_h=('flow_duration_h', 'max')).reset_index()
    nm1 = nm.drop_duplicates('api14')[['api14', 'name', 'methane_pre_plugging_emissions', 'methane_detection_method']]
    out = nm1.merge(per, left_on='api14', right_on='api', how='left')
    out.to_csv(OUT / 'nm_memo_audit.csv', index=False)
    has = out.memos.notna()
    print('NM wells: %d | with a methane memo: %d | field average flagged: %d | memo value found: %d'
          % (len(out), has.sum(), (out.field_average == True).sum(), out.memo_g_h.notna().sum()))
    both = out.dropna(subset=['memo_g_h', 'methane_pre_plugging_emissions'])
    if len(both):
        agree = (both.memo_g_h - both.methane_pre_plugging_emissions).abs() <= 0.01 + 0.05 * both.methane_pre_plugging_emissions.abs()
        print('memo vs DOI value agree (within 5%%): %d of %d' % (agree.sum(), len(both)))
