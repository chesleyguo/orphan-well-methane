import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score

d = pd.read_pickle(PROC / 'doi.pkl'); u = pd.read_pickle(PROC / 'inv.pkl')
m = d[d.measured_real & d.lat.notna() & (d.method_class != 'handheld detector')].copy()   # drop methods that never detect
m['y'] = m.methane_pre_plugging_emissions
m['det'] = (m.y > 0).astype(int)
m['wt'] = m.well_type.fillna('Unknown')
m['inv_status'] = m.uid.map(u.set_index('uid').Status.str.upper().str.strip()).fillna('NOT_IN_2022')
m['in2022'] = m.uid.notna().astype(int)
m['year'] = pd.to_datetime(m.witness_date, unit='ms', errors='coerce').dt.year
print('modeling sample:', len(m), 'detects', m.det.sum(), 'states', m.st.nunique())

X = pd.get_dummies(m[['lat', 'lon', 'wt', 'method_class', 'in2022']], columns=['wt', 'method_class']).astype(float)
# spatial blocks ~0.5 degree to avoid neighbour leakage
grp = (np.floor(m.lat / 0.5)).astype(int).astype(str) + '_' + (np.floor(m.lon / 0.5)).astype(int).astype(str)
def cvpred(X, y, groups, clf):
    p = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, y, groups):
        clf.fit(X.iloc[tr], y.iloc[tr]); p[te] = clf.predict_proba(X.iloc[te])[:, 1]
    return p
clf = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=15, random_state=0)
m['p_det'] = cvpred(X, m.det, grp, clf)
print('detection AUC (spatial-block CV):', round(roc_auc_score(m.det, m.p_det), 3),
      ' AP:', round(average_precision_score(m.det, m.p_det), 3), ' base rate:', round(m.det.mean(), 3))
# without location (attributes + method only) -> how much is just geography?
Xn = X.drop(columns=['lat', 'lon'])
pn = cvpred(Xn, m.det, grp, clf); print('AUC without lat/lon:', round(roc_auc_score(m.det, pn), 3))
# high emitter (>=10 g/h) prediction
m['hi'] = (m.y >= 10).astype(int); m['p_hi'] = cvpred(X, m.hi, grp, clf)
print('high-emitter (>=10 g/h) n=', m.hi.sum(), ' AUC:', round(roc_auc_score(m.hi, m.p_hi), 3))

# ---- prioritisation: plug top-k by model vs random vs actual order; emissions captured ----
tot = m.y.sum(); n = len(m); rng = np.random.default_rng(0)
res = []
for frac in [0.05, 0.1, 0.2, 0.3, 0.5]:
    k = int(n * frac)
    mod = m.nlargest(k, 'p_hi').y.sum() / tot
    rnd = np.mean([m.sample(k, random_state=s).y.sum() / tot for s in range(200)])
    orc = m.nlargest(k, 'y').y.sum() / tot
    res.append((frac, round(rnd, 3), round(mod, 3), round(orc, 3)))
print('share of emissions captured plugging top fraction [frac, random, model, oracle]:'); [print(' ', r) for r in res]

# ---- cost-effectiveness ----
c = m.total_cost_complete_pa
cc = m[(c > 0) & (m.y > 0)]
t_yr = cc.y * 8760 / 1e6
print('cost per detected well median $', round(cc.total_cost_complete_pa.median()),
      ' corr(log cost, log emission)=', round(np.corrcoef(np.log(cc.total_cost_complete_pa), np.log(cc.y))[0, 1], 2))
print('$ per t CH4/yr avoided (detected wells): median', round((cc.total_cost_complete_pa / t_yr).median()),
      ' program-wide (all measured):', round(m.total_cost_complete_pa.sum() / (m.y.sum() * 8760 / 1e6)))

# ---- active measurement simulation: estimate statewide total with budget B, random vs uncertainty-guided ----
# Here: pool = modeling sample; truth = sum(y). Strategy A random sample + expansion; B stratified by model score (oversample high p_hi) + Horvitz-Thompson.
truth = tot; B = 200; sims = 500
errA, errB = [], []
w = (m.p_hi.clip(0.02, None)) ; w = w / w.sum()
for s in range(sims):
    r = np.random.default_rng(s)
    ia = r.choice(n, B, replace=False); errA.append(m.y.values[ia].mean() * n / truth - 1)
    ib = r.choice(n, B, replace=True, p=w.values); errB.append(np.mean(m.y.values[ib] / (w.values[ib] * n)) * n / truth - 1)
for nm, e in [('random', errA), ('model-guided PPS', errB)]:
    e = np.array(e); print(f'B={B} {nm}: RMSE of total = {np.sqrt((e**2).mean()):.2f}, 90% interval of rel. error [{np.quantile(e,.05):.2f}, {np.quantile(e,.95):.2f}]')
m.to_pickle(PROC / 'model_sample.pkl')

# ---- robustness: tail sensitivity and within-state ranking (the decision states actually face) ----
print('top emitters g/h:', m.y.nlargest(5).round(0).tolist(), ' top-3 share:', round(m.y.nlargest(3).sum() / m.y.sum(), 3))
for drop in [0, 1, 3, 10]:
    mm = m.drop(m.y.nlargest(drop).index) if drop else m
    print('drop top', drop, '-> model capture of emissions at 10/20/30%:',
          [round(mm.nlargest(int(len(mm) * f), 'p_hi').y.sum() / mm.y.sum(), 2) for f in (.1, .2, .3)])
for f in [.1, .2, .3]:
    k = int(n * f)
    sel = m.groupby('st', group_keys=False).apply(lambda g: g.nlargest(max(1, int(round(len(g) * f))), 'p_hi'))
    print(f'top {int(f*100)}%: pooled ranking captures {m.nlargest(k, "p_hi").hi.sum() / m.hi.sum():.2f} of high emitters; '
          f'within-state ranking captures {sel.hi.sum() / m.hi.sum():.2f} (random = {f:.2f})')
wa = [roc_auc_score(g.hi, g.p_hi) for s_, g in m.groupby('st') if g.hi.sum() >= 10 and (1 - g.hi).sum() >= 10]
print('within-state AUC (high emitter):', [round(x, 2) for x in wa], ' mean', round(float(np.mean(wa)), 2))
