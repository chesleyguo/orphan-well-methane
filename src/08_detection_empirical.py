"""Empirical view of detection by method class in the DOI program data: detection rate, smallest
non-zero values, share of detects below 1 g/h, and recurring (rounded) values.
These describe what each method class actually reported; they are inputs for choosing
detection-limit priors (see docs/detection_limits.md), not estimates of the limits themselves.
Run from the repo root after 01: python src/08_detection_empirical.py
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
import re
import pandas as pd, numpy as np

def method_group(s):
    """Refined instrument grouping (proposed, see decision log D14). Separates quantitative rate
    methods from screening devices that read concentration or show a plume."""
    s0 = str(s).strip()
    s = s0.lower()
    if s0 in ('', 'nan', 'N/A', 'NA', 'NR', 'NP'):
        return 'none'
    if s0 == 'ND':
        return 'ND (ambiguous: no data or not detected)'
    if 'calc' in s:
        return 'calculated'
    if any(w in s for w in ['hi-flow', 'hiflow', 'high flow', 'semtech']):
        return 'hi-flow sampler'
    if 'flow through' in s or 'chamber' in s:
        return 'chamber / flow-through cell'
    if 'chromatograph' in s:
        return 'flow meter + GC'
    if 'mgga' in s or s == 'laser absorption spectroscopy':
        return 'laser analyzer (quantitative)'
    if any(w in s for w in ['rmld', 'tdlas', 'laser methane']):
        return 'screening: standoff laser (RMLD-type)'
    if any(w in s for w in ['ogi', 'flir', 'gx320']):
        return 'screening: OGI camera'
    if any(w in s for w in ['sensit', 'bascom', 'altair', 'rki', 'eagle', '4-gas', 'gas meter', 'sniffer', 'arearae', 'detector', 'monitor']):
        return 'screening: handheld concentration detector'
    if 'federal methane' in s:
        return 'federal guidelines (instrument unspecified)'
    if s == 'fid':
        return 'FID (concentration)'
    return 'other / unspecified'

d = pd.read_pickle(PROC / 'doi.pkl')
m = d[d.methane_pre_plugging_emissions.notna()].copy()
m['y'] = m.methane_pre_plugging_emissions
m['method_group'] = m.methane_detection_method.map(method_group)
m[['OBJECTID', 'method_class', 'method_group']].to_csv(PROC / 'method_group.csv', index=False)
rows = []
for cls, g in m.groupby('method_group'):
    pos = g.y[g.y > 0]
    top = pos.round(4).value_counts()
    rows.append(dict(
        method_group=cls, n=len(g), states=', '.join(sorted(g.st.dropna().unique()))[:60],
        n_detect=len(pos), detect_rate=round(len(pos) / len(g), 3),
        min_pos=pos.min() if len(pos) else np.nan,
        p05_pos=round(pos.quantile(.05), 4) if len(pos) else np.nan,
        p25_pos=round(pos.quantile(.25), 3) if len(pos) else np.nan,
        median_pos=round(pos.median(), 2) if len(pos) else np.nan,
        share_pos_below_1=round((pos < 1).mean(), 2) if len(pos) else np.nan,
        top_values='; '.join(f'{v:g}×{c}' for v, c in top.head(3).items() if c > 1)))
t = pd.DataFrame(rows).sort_values('n', ascending=False)
t.to_csv(OUT / 'detection_by_method_empirical.csv', index=False)
pd.set_option('display.width', 250)
print(t.to_string(index=False))

# raw method strings within each class (to audit the keyword mapping)
aud = (m.groupby(['method_group', 'methane_detection_method']).agg(n=('y', 'size'), detect=('y', lambda x: int((x > 0).sum())))
         .reset_index().sort_values(['method_group', 'n'], ascending=[True, False]))
aud.to_csv(OUT / 'method_string_audit.csv', index=False)
print('\nmethod strings:', len(aud), '-> outputs/method_string_audit.csv')


# --- Hi-Flow maker/model named in the method string, by state (D18) ---
mm = d.methane_detection_method.fillna('').astype(str)
is_hf = mm.str.contains(r'hi.?flow|high.?flow|semtec|bacharach', case=False, regex=True)
def hf_model(t):
    t = t.lower()
    if 'semtec' in t: return 'Semtech Hi-Flow 2 (named)'
    if 'bacharach' in t: return 'Bacharach (named)'
    if re.search(r'2\.0|flow ?2|flow ii', t): return 'Hi-Flow 2/II, maker not named'
    return 'model not named'
hf = pd.crosstab(d.loc[is_hf, 'st'], mm[is_hf].map(hf_model), margins=True)
hf.to_csv(OUT / 'hiflow_model_by_state.csv')
print('Hi-Flow records by named model:', hf.loc['All'].to_dict())
