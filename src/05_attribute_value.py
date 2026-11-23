import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
import zlib, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
e = pd.read_pickle(PROC / 'doi_enriched.pkl')
m = e[e.measured_real & (e.method_class != 'handheld detector') & e.st.isin(['LA','NM','CO','CA','NY','WV']) & e.lat.notna()].copy()
m['y'] = m.methane_pre_plugging_emissions; m['det'] = (m.y > 0).astype(int); m['hi'] = (m.y >= 10).astype(int)
m['wt'] = m.well_type.fillna('Unknown')
m['yrs_since_prod'] = (pd.to_datetime(m.witness_date, unit='ms', errors='coerce') - m.last_prod).dt.days / 365.25
m['spud_year'] = m.spud.dt.year
print('n', len(m), 'det', m.det.sum(), 'hi', m.hi.sum(), ' age avail', m.age_at_plug.notna().mean().round(2), ' depth avail', m.depth.notna().mean().round(2))
base = pd.get_dummies(m[['lat','lon','wt','method_class','st']], columns=['wt','method_class','st']).astype(float)
plus = pd.concat([base, m[['age_at_plug','depth','spud_year','yrs_since_prod']].astype(float)], axis=1)
grp = (np.floor(m.lat/0.5)).astype(int).astype(str)+'_'+(np.floor(m.lon/0.5)).astype(int).astype(str)
def cv(X, y, reps=5):
    aucs=[]
    for r in range(reps):
        p=np.zeros(len(y)); g=grp.map(lambda s: zlib.crc32(f'{s}|{r}'.encode())%1000)
        for tr,te in GroupKFold(5).split(X,y,g):
            c=HistGradientBoostingClassifier(max_iter=200,learning_rate=0.05,max_leaf_nodes=15,random_state=r).fit(X.iloc[tr],y.iloc[tr]); p[te]=c.predict_proba(X.iloc[te])[:,1]
        aucs.append(roc_auc_score(y,p))
    return np.mean(aucs), np.std(aucs)
for t in ['det','hi']:
    a=cv(base,m[t]); b=cv(plus,m[t]); print(f'{t}: base AUC {a[0]:.3f}±{a[1]:.3f}  +age/depth {b[0]:.3f}±{b[1]:.3f}')
# simple univariate looks
m['age_bin']=pd.cut(m.age_at_plug,[0,30,50,70,200]); print(m.groupby('age_bin',observed=True).agg(n=('det','size'),det=('det','mean'),hi=('hi','mean')).round(2).to_string())
m['dep_bin']=pd.qcut(m.depth,4); print(m.groupby('dep_bin',observed=True).agg(n=('det','size'),det=('det','mean'),hi=('hi','mean')).round(2).to_string())

# ---- single-state logistic models (LA, NM): do attributes help within a state? ----
from sklearn.linear_model import LogisticRegression
for st in ['LA', 'NM']:
    g = e[(e.st == st) & e.measured_real & e.lat.notna()].copy(); g['det'] = (g.methane_pre_plugging_emissions > 0).astype(int)
    g = g[g.depth.notna() & g.age_at_plug.notna()]
    grp2 = (np.floor(g.lat / 0.25)).astype(int).astype(str) + '_' + (np.floor(g.lon / 0.25)).astype(int).astype(str)
    g['ld'] = np.log(g.depth.clip(1)); g['gas'] = g.well_type.fillna('').str.contains('Gas').astype(int)
    def cvl(cols):
        X = g[cols].astype(float); X = (X - X.mean()) / X.std(); p = np.zeros(len(g))
        for tr, te in GroupKFold(5).split(X, g.det, grp2):
            p[te] = LogisticRegression(max_iter=1000).fit(X.iloc[tr], g.det.iloc[tr]).predict_proba(X.iloc[te])[:, 1]
        return round(roc_auc_score(g.det, p), 3)
    print(st, 'n', len(g), 'det', g.det.sum(), ' loc', cvl(['lat', 'lon']), ' loc+gas', cvl(['lat', 'lon', 'gas']),
          ' +depth+age', cvl(['lat', 'lon', 'gas', 'ld', 'age_at_plug']), ' depth+age+gas', cvl(['ld', 'age_at_plug', 'gas']))
