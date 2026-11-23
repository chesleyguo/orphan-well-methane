"""PILOT (preliminary, for discussion only): how much do non-detects and method-specific detection
limits change the picture, and what is screening worth?

Simplifications relative to the planned analysis (do not cite these numbers as results):
- one pooled hurdle model for all states: a share pi of wells does not leak; leak rates of
  leaking wells are lognormal(mu, sigma). A non-detect is either a non-leaker or a leak below the
  method's limit; limits spanning five orders of magnitude make pi and the lognormal separable;
- leak rates are capped at 50 kg/h, just above the largest DOI value (46 kg/h, CO); without a cap
  the fitted lognormal (sigma ~4.5) puts mass far beyond anything observed;
- no selection model (plugged wells are treated as the candidate population);
- each method is a single left-censoring limit (no two-stage leak-finding for Hi-Flow);
- placeholder zeros (no method) and NM field averages are dropped (D20, D21).

Outputs: outputs/pilot/pilot_fit.csv, pilot_breakeven.csv, pilot_costs.csv and two figures.
Run from the repo root: python src/11_pilot_censored.py
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *

import numpy as np, pandas as pd
from scipy import optimize, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker

P = OUT / 'pilot'; P.mkdir(exist_ok=True)
rng = np.random.default_rng(0)

# ---- censoring limit per method group (g/h); midpoints of data/reference/detection_limits.csv
LIMIT = {
    'chamber / flow-through cell': 0.001, 'flow meter + GC': 0.01, 'hi-flow sampler': 1.0,
    'laser analyzer (quantitative)': 1.0, 'federal guidelines (instrument unspecified)': 1.0,
    'screening: OGI camera': 50.0, 'screening: handheld concentration detector': 10.0,
    'screening: standoff laser (RMLD-type)': 10.0, 'FID (concentration)': 10.0,
}

d = pd.read_pickle(PROC / 'doi.pkl')
g = pd.read_csv(PROC / 'method_group.csv')
d = d.merge(g[['OBJECTID', 'method_group']], on='OBJECTID', how='left')
d['v'] = d.methane_pre_plugging_emissions
x = d[d.v.notna() & d.method_group.isin(LIMIT)].copy()
nm_rep = x[(x.st == 'NM') & (x.v > 0)].v.round(4).value_counts()
x = x[~((x.st == 'NM') & x.v.round(4).isin(nm_rep[nm_rep >= 3].index))]       # D21
x['L'] = x.method_group.map(LIMIT)
x['det'] = x.v > 0
print('pilot sample: %d wells (%d detects, %d non-detects), %d states'
      % (len(x), x.det.sum(), (~x.det).sum(), x.st.nunique()))


def negll(theta, ly, lL):
    mu, ls, lp = theta; sd = np.exp(ls); pi = 1 / (1 + np.exp(-lp))
    return -(len(ly) * np.log1p(-pi) + stats.norm.logpdf(ly, mu, sd).sum()
             + np.log(pi + (1 - pi) * stats.norm.cdf((lL - mu) / sd)).sum())


def fit(df):
    ly = np.log(df.loc[df.det, 'v'].values); lL = np.log(df.loc[~df.det, 'L'].values)
    best = None
    for p0 in (-1.0, 0.0, 1.0):
        r = optimize.minimize(negll, [np.median(ly), np.log(ly.std()), p0], args=(ly, lL), method='Nelder-Mead',
                              options={'xatol': 1e-7, 'fatol': 1e-7, 'maxiter': 8000})
        if best is None or r.fun < best.fun: best = r
    mu, ls, lp = best.x
    return mu, np.exp(ls), 1 / (1 + np.exp(-lp))


mu, s, pi = fit(x)
boot = np.array([fit(x.sample(len(x), replace=True, random_state=int(i))) for i in range(200)])
print('hurdle fit: share not leaking %.2f | leakers: median %.2f g/h, sigma %.2f' % (pi, np.exp(mu), s))


CAP = 50000.0   # g/h


def capped_mean(m, sd):
    lc = np.log(CAP)
    return np.exp(m + sd ** 2 / 2) * stats.norm.cdf((lc - m - sd ** 2) / sd) + CAP * stats.norm.sf((lc - m) / sd)


def summ(m, sd, p):
    q = 1 - p
    return {'share_leaking': q, 'median_leaker_g_h': np.exp(m), 'mean_g_h': q * capped_mean(m, sd),
            'P_gt_1': q * stats.norm.sf(-m / sd), 'P_gt_10': q * stats.norm.sf((np.log(10) - m) / sd),
            'P_gt_100': q * stats.norm.sf((np.log(100) - m) / sd)}


pos = x.loc[x.det, 'v']
rows = [
    {'estimate': 'quantified wells only (as in RFF 2026)', 'n': len(pos), 'share_leaking': 1.0, 'median_leaker_g_h': pos.median(),
     'mean_g_h': pos.mean(), 'P_gt_1': (pos > 1).mean(), 'P_gt_10': (pos > 10).mean(), 'P_gt_100': (pos > 100).mean()},
    {'estimate': 'non-detects as zero (as in DOI averages)', 'n': len(x), 'share_leaking': x.det.mean(), 'median_leaker_g_h': pos.median(),
     'mean_g_h': x.v.mean(), 'P_gt_1': (x.v > 1).mean(), 'P_gt_10': (x.v > 10).mean(), 'P_gt_100': (x.v > 100).mean()},
    dict(estimate='censoring-corrected hurdle model (pilot)', n=len(x), mu=mu, sigma=s, pi_not_leaking=pi, **summ(mu, s, pi)),
]
bs = pd.DataFrame([summ(*b) for b in boot])
for k in ['share_leaking', 'median_leaker_g_h', 'mean_g_h', 'P_gt_1', 'P_gt_10', 'P_gt_100']:
    rows[2][k + '_lo'], rows[2][k + '_hi'] = bs[k].quantile([.05, .95])
fitdf = pd.DataFrame(rows)
fitdf.to_csv(P / 'pilot_fit.csv', index=False, float_format='%.6g')
print(fitdf[['estimate', 'n', 'share_leaking', 'median_leaker_g_h', 'mean_g_h', 'P_gt_1', 'P_gt_10', 'P_gt_100']].to_string(index=False, float_format='%.4g'))

# ---- break-even screening cost vs detection limit
SCC = 1900.0                     # $/t CH4, EPA 2023 social cost of methane, 2% discount, ~2025-2030
DECLINE, DISC, YEARS = 0.064, 0.02, 25
PVF = sum(((1 - DECLINE) / (1 + DISC)) ** t for t in range(YEARS))
VALUE = 8.766e-3 * SCC * PVF     # $ present value per g/h abated
SIG_M = np.sqrt(np.log(1 + 0.65 ** 2))   # snapshot noise from USGS MT repeat visits (D24)
M = 20000
limits = np.logspace(-3, 4, 50)


def gain_per_screened(m, sd, p, L, frac, seed, p_find=1.0):
    r = np.random.default_rng(seed)
    z, u, n, tie = r.normal(size=M), r.random(M), r.normal(size=M), r.random(M)   # common random numbers
    e = np.minimum(np.exp(m + sd * z), CAP) * (u > p)
    obs = e * np.exp(-SIG_M ** 2 / 2 + SIG_M * n)
    found = np.random.default_rng(seed + 7).random(M) < p_find
    detected = (obs > L) & found
    order = np.lexsort((tie, -np.where(detected, obs, 0)))
    k = int(frac * M)
    return (e[order[:k]].sum() - frac * e.sum()) * VALUE / M


res = []
for frac in (0.1, 0.3):
    for i, L in enumerate(limits):
        central = gain_per_screened(mu, s, pi, L, frac, 1000)
        band = [gain_per_screened(bm, bsd, bp, L, frac, 2000 + j) for j, (bm, bsd, bp) in enumerate(boot[:50])]
        res.append({'plug_fraction': frac, 'detection_limit_g_h': L, 'breakeven_cost_usd': central,
                    'lo': np.quantile(band, .05), 'hi': np.quantile(band, .95)})
be = pd.DataFrame(res)

# sensitivity: leak point found with probability p_find regardless of size (Hi-Flow finding stage;
# LSU: Hi-Flow detected 23% of wells vs 96% by chamber -> p_find ~0.25 is a pessimistic bound)
sens = []
for pf in (1.0, 0.75, 0.5, 0.25):
    central = gain_per_screened(mu, s, pi, 1.0, 0.1, 1000, pf)
    band = [gain_per_screened(bm, bsd, bp, 1.0, 0.1, 2000 + j, pf) for j, (bm, bsd, bp) in enumerate(boot[:50])]
    sens.append({'p_find': pf, 'detection_limit_g_h': 1.0, 'plug_fraction': 0.1, 'breakeven_cost_usd': central,
                 'lo': np.quantile(band, .05), 'hi': np.quantile(band, .95)})
sens = pd.DataFrame(sens)
sens.to_csv(P / 'pilot_breakeven_findprob.csv', index=False, float_format='%.6g')
print('break-even (L=1 g/h, plug 10%) by probability the leak is found:')
print(sens[['p_find', 'breakeven_cost_usd', 'lo', 'hi']].round({'breakeven_cost_usd': 0, 'lo': 0, 'hi': 0}).to_string(index=False))

# state heterogeneity: same hurdle model fitted separately where data allow
st_rows = []
for st, xs_ in x.groupby('st'):
    if xs_.det.sum() >= 20 and (~xs_.det).sum() >= 20:
        m_, s_, p_ = fit(xs_)
        st_rows.append(dict(state=st, n=len(xs_), detects=int(xs_.det.sum()),
                            main_method=xs_.method_group.value_counts().index[0], **summ(m_, s_, p_), sigma=s_))
stdf = pd.DataFrame(st_rows).sort_values('n', ascending=False)
stdf.to_csv(P / 'pilot_by_state.csv', index=False, float_format='%.6g')
print(stdf[['state', 'n', 'detects', 'main_method', 'share_leaking', 'median_leaker_g_h', 'P_gt_10', 'sigma']].to_string(index=False, float_format='%.3g'))
be.to_csv(P / 'pilot_breakeven.csv', index=False, float_format='%.6g')
print('value per g/h abated: $%.0f (PV, %d yr, %.1f%%/yr decline, %.0f%% discount)' % (VALUE, YEARS, DECLINE * 100, DISC * 100))
for frac in (0.1, 0.3):
    b = be[be.plug_fraction == frac].set_index('detection_limit_g_h').breakeven_cost_usd
    print('plug %d%%: break-even $/screened well at L=0.01: %.0f | L=1: %.0f | L=50: %.0f'
          % (frac * 100, np.interp(-2, np.log10(b.index), b), np.interp(0, np.log10(b.index), b),
             np.interp(np.log10(50), np.log10(b.index), b)))

# ---- plugging cost context (state median; many states report contract averages)
c = d.assign(c=pd.to_numeric(d.total_cost_complete_pa, errors='coerce'))
c = c[c.c > 0].groupby('st').c.agg(n='size', median='median', share_repeated=lambda s: s.duplicated(keep=False).mean())
c.round(2).to_csv(P / 'pilot_costs.csv')
print('median reported plugging cost per well, all states: $%.0f' % d.total_cost_complete_pa[d.total_cost_complete_pa > 0].median())

# ---- figures (reference palette, light mode)
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.facecolor': SURF, 'figure.facecolor': SURF, 'font.family': 'DejaVu Sans'})


def style(ax):
    ax.grid(True, which='major', color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)


xs = np.logspace(-4, 4, 400)
fig, ax = plt.subplots(figsize=(7.2, 4.4))
emp = lambda v: np.array([(v > t).mean() for t in xs])
ax.plot(xs, emp(pos.values), color=C2, lw=2, label='Quantified wells only (as in RFF 2026)')
ax.plot(xs, emp(x.v.values), color=C1, lw=2, label='Non-detects counted as zero (as in DOI averages)')
sf = np.array([(1 - bp) * stats.norm.sf((np.log(xs) - bm) / bsd) * (xs < CAP) for bm, bsd, bp in boot])
ax.fill_between(xs, np.quantile(sf, .05, 0), np.quantile(sf, .95, 0), color=C3, alpha=0.2, lw=0)
ax.plot(xs, (1 - pi) * stats.norm.sf((np.log(xs) - mu) / s) * (xs < CAP), color=C3, lw=2, label='Censoring-corrected (pilot, 90% band)')
ax.set_xscale('log'); ax.set_xlim(1e-4, 1e4); ax.set_ylim(0, 1)
ax.set_xlabel('Pre-plugging methane emission rate (g/h)'); ax.set_ylabel('Share of wells above this rate')
ax.set_title('Same wells, three readings of the data (PRELIMINARY)', loc='left', color=INK, fontsize=11)
style(ax); ax.legend(frameon=False, fontsize=9, loc='upper right')
fig.text(0.01, 0.01, 'DOI plugged-well records with a named method (n=%d, %d states). Pooled hurdle model capped at 50 kg/h; no selection model.'
         % (len(x), x.st.nunique()), fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(P / 'fig_pilot_distribution.png', dpi=200); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.2, 4.4))
for frac, col in ((0.1, C1), (0.3, C2)):
    b = be[be.plug_fraction == frac]
    ax.fill_between(b.detection_limit_g_h, b.lo, b.hi, color=col, alpha=0.18, lw=0)
    ax.plot(b.detection_limit_g_h, b.breakeven_cost_usd, color=col, lw=2,
            label='Budget plugs %d%% of candidate wells' % (frac * 100))
ax.set_ylim(0, be.hi.max() * 1.15)
for L, name in ((0.01, 'chamber /\nflow meter'), (1, 'Hi-Flow'), (50, 'OGI camera')):
    ax.axvline(L, color=GRID, lw=1, zorder=0)
    ax.text(L, ax.get_ylim()[1] * 0.98, name, fontsize=8, color=INK2, ha='center', va='top',
            bbox=dict(facecolor=SURF, edgecolor='none', pad=1))
ax.set_xscale('log')
ax.set_xlabel('Detection limit of the screening measurement (g/h)')
ax.set_ylabel('Break-even cost per screened well (US$)')
ax.set_title('What screening before plugging is worth (PRELIMINARY)', loc='left', color=INK, fontsize=11)
style(ax); ax.legend(frameon=False, fontsize=9, loc='lower left')
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: '${:,.0f}'.format(v)))
fig.text(0.01, 0.01, 'Plug top-ranked wells after one noisy visit (CV 0.65) vs random choice; SC-CH4 $%d/t, %d-yr horizon, %.1f%%/yr decline.'
         % (SCC, YEARS, DECLINE * 100), fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(P / 'fig_pilot_breakeven.png', dpi=200); plt.close(fig)
print('figures written to', P.relative_to(ROOT))
