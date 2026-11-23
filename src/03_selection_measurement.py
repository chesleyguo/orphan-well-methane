import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold, StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
import statsmodels.formula.api as smf

u = pd.read_pickle(PROC / 'inv.pkl'); d = pd.read_pickle(PROC / 'doi.pkl')
out = {}

# ---------- Selection: which inventory wells got plugged ----------
u['type_c'] = u.Type.str.upper().str.strip().fillna('UNK')
top = u.type_c.value_counts().index[:8]; u['type_c'] = u.type_c.where(u.type_c.isin(top), 'OTHER')
u['status_c'] = u.Status.str.upper().str.strip().fillna('UNK')
tops = u.status_c.value_counts().index[:10]; u['status_c'] = u.status_c.where(u.status_c.isin(tops), 'OTHER')
u['has_api'] = u.api10.notna()
# local density of inventory wells (wells within ~5 km) as a crude "cluster" covariate
from sklearn.neighbors import BallTree
ok = u.lat.notna()
tree = BallTree(np.radians(u.loc[ok, ['lat', 'lon']].values), metric='haversine')
u.loc[ok, 'dens5km'] = tree.query_radius(np.radians(u.loc[ok, ['lat', 'lon']].values), r=5000/6371000, count_only=True)
u['plug'] = u.plugged.astype(int)
v = u[ok & u.st.isin(u[u.plugged].st.unique())].copy()
print('selection universe (states with any plugging):', len(v), 'plugged', v.plug.sum())
print('plug rate by state:\n', v.groupby('st').plug.agg(['size', 'mean']).sort_values('size', ascending=False).round(3).head(15).T.to_string())

# within-state predictability: GBM on lat/lon/type/status/density, AUC computed within each state then averaged
X = pd.get_dummies(v[['lat', 'lon', 'dens5km', 'type_c', 'status_c', 'has_api']], columns=['type_c', 'status_c']).astype(float)
p = cross_val_predict(HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, random_state=0),
                      X, v.plug, cv=StratifiedKFold(5, shuffle=True, random_state=0), method='predict_proba')[:, 1]
v['p_sel'] = p
aucs = {s: roc_auc_score(g.plug, g.p_sel) for s, g in v.groupby('st') if g.plug.sum() >= 30 and (1 - g.plug).sum() >= 30}
print('within-state AUC of selection model:', {k: round(x, 2) for k, x in sorted(aucs.items())})
print('median within-state AUC:', round(np.median(list(aucs.values())), 2))
out['sel_auc'] = aucs
v[['uid', 'st', 'plug', 'p_sel', 'dens5km', 'type_c']].to_pickle(PROC / 'sel.pkl')

# ---------- Measurement: does detection depend on method? ----------
m = d[d.methane_pre_plugging_emissions.notna()].copy()
m['det'] = (m.methane_pre_plugging_emissions > 0).astype(int)
t = m.groupby('method_class').agg(n=('det', 'size'), det_rate=('det', 'mean'),
        min_pos=('methane_pre_plugging_emissions', lambda x: x[x > 0].min()),
        med_pos=('methane_pre_plugging_emissions', lambda x: x[x > 0].median())).sort_values('n', ascending=False)
print('\nDetection by method class:\n', t.round(3).to_string())
mr = m[m.method_class != 'none'].copy()
mr['wt'] = mr.well_type.fillna('Unknown')
fit = smf.logit('det ~ C(method_class, Treatment("hi-flow")) + C(st) + C(wt)', data=mr).fit(disp=0, method='bfgs', maxiter=2000)
coef = fit.params.filter(like='method_class').round(2)
print('\nlogit coefs (vs hi-flow), controlling for state + well type:\n', coef.to_string())
# within-state variation in method: states using >1 method class
mix = mr.groupby('st').method_class.nunique(); print('states using >1 method class:', mix[mix > 1].to_dict())
for s in mix[mix > 1].index:
    g = mr[mr.st == s]; tt = g.groupby('method_class').det.agg(['size', 'mean']); tt = tt[tt['size'] >= 15]
    if len(tt) > 1: print(' ', s, tt.round(2).to_dict('index'))
