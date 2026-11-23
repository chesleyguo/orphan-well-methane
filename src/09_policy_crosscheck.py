"""Cross-check state measurement policy against what the DOI records show.

For each state: records, share with a pre-plugging value before/after the July 2025 federal guidance
change (measurement no longer required, estimates allowed), share of values that are exactly zero,
share with no method recorded, and values that repeat three or more times among positive values
(candidates for field averaging, as in New Mexico, or for instrument resolution steps, as in Louisiana).
Run from the repo root: python src/09_policy_crosscheck.py
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *

import pandas as pd

CHANGE = pd.Timestamp('2025-07-11')   # revised State Formula Grant guidance and FAQ (D19)

d = pd.read_pickle(PROC / 'doi.pkl')
g = pd.read_csv(PROC / 'method_group.csv')
d = d.merge(g[['OBJECTID', 'method_group']], on='OBJECTID', how='left')
d['date'] = pd.to_datetime(d.witness_date, unit='ms', errors='coerce')
d['v'] = d.methane_pre_plugging_emissions
d['has'] = d.v.notna()
d['post'] = d.date >= CHANGE


def repeats(s):
    vc = s[s > 0].round(4).value_counts()
    vc = vc[vc >= 3]
    return '; '.join(f'{k:g} x{n}' for k, n in vc.head(4).items())


rows = []
for st, x in d.groupby('st'):
    pre, post = x[~x.post & x.date.notna()], x[x.post]
    vals = x.v.dropna()
    rows.append({
        'state': st, 'records': len(x),
        'with_value': int(x.has.sum()),
        'share_with_value_before_2025_07': round(pre.has.mean(), 2) if len(pre) else None,
        'share_with_value_after_2025_07': round(post.has.mean(), 2) if len(post) else None,
        'records_after_2025_07': len(post),
        'share_zero_of_values': round((vals == 0).mean(), 2) if len(vals) else None,
        'share_no_method_of_values': round((x.loc[x.has, 'method_group'] == 'none').mean(), 2) if len(vals) else None,
        'share_above_1g_h': round((vals > 1).mean(), 2) if len(vals) else None,
        'repeated_positive_values': repeats(vals),
    })
t = pd.DataFrame(rows).sort_values('records', ascending=False)
t.to_csv(OUT / 'state_measurement_crosscheck.csv', index=False)
print(t.to_string(index=False))

allzero = t[(t.share_zero_of_values == 1) & (t.with_value >= 10)]
print('\nstates where every reported value is 0 (>=10 values):', ', '.join(allzero.state),
      '| records:', int(allzero.with_value.sum()))
print('records with a value: before change %.2f, after %.2f' % (d[~d.post & d.date.notna()].has.mean(), d[d.post].has.mean()))
