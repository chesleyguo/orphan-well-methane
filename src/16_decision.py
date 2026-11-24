"""Decision analysis (analysis_plan.md, Sec. 5): what is measuring before plugging worth?

For each posterior draw, leak status and leak rate are drawn for every unplugged 2022-inventory well,
conditional on that well not having been plugged (the selection model with rho makes leaking wells
more likely to be plugged in leak-rule states, so the remaining wells leak a little less). Each state
then plugs a fraction f of its candidates under four policies:
  P0 the observed selection model (state, spatial cell, and rho * leaking);
  P1 random;
  P2 rank by posterior mean emission given attributes (state is fixed by the budget, so this is
     gas wells first when the posterior says gas wells emit more);
  P3 screen every candidate once with a method of detection limit L and leak-finding probability
     p_find, plug in order of the reading, fill the rest by P2.
Outputs (outputs/decision/):
  policy_value.csv     avoided emission (g/h, t/yr) and present value ($) per policy and budget
  breakeven.csv        E5: c*(L, p_find) = (value P3 - value P2) / wells screened
  e6_rule_change.csv   E6: value of screening with the formula-grant-era method mix minus P0,
                       per plugged well and for the wells plugged after 2025-07-11
  by_state.csv         P3 (Hi-Flow-like: L = 1 g/h, posterior p_find) gain over P2, per state
  fig_breakeven.png
Until the NPE ensemble passes SBC and coverage (src/15), outputs are labelled TRIAL.
Run from the repo root: python src/16_decision.py [--draws N] [--prior]
"""
import sys, pathlib, argparse; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
from sim.model import build_design, sample_prior, NAMES, METHODS, LA_STEP

import numpy as np, pandas as pd
from scipy.special import expit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker

ap = argparse.ArgumentParser()
ap.add_argument('--draws', type=int, default=200)
ap.add_argument('--prior', action='store_true', help='use prior draws (pipeline test only)')
ap.add_argument('--label', default='TRIAL')
args = ap.parse_args()
P = OUT / 'decision'; P.mkdir(exist_ok=True)
rng = np.random.default_rng(16)

# ---- value of 1 g/h avoided (same constants as the pilot and the analysis plan)
SCC = 1900.0                      # $/t CH4 (EPA 2023, 2% discount)
DECLINE, DISC, YEARS = 0.064, 0.02, 25
PVF = sum(((1 - DECLINE) / (1 + DISC)) ** t for t in range(YEARS))
T_PER_GH = 8.766e-3               # t/yr per g/h
VALUE = T_PER_GH * SCC * PVF      # $ per g/h avoided
FRACS = (0.1, 0.3)
LGRID = np.logspace(-3, 4, 15)
PGRID = (0.25, 0.5, 0.75, 1.0)
N_AFTER_RULE = 2259               # DOI wells with witness date on or after 2025-07-11 (snapshot 2026-11-22)

d = build_design()
S = len(d['states'])
cand = ~d['inv_plugged']
st, gas, cell = d['inv_state'][cand], d['inv_gas'][cand], d['inv_cell'][cand]
rule = d['leak_rule'][st]
N = len(st)
n_s = np.bincount(st, minlength=S)

if args.prior:
    th = sample_prior(rng, args.draws); label = 'PRIOR'
else:
    z = np.load(PROC / 'npe_ensemble' / 'real_posterior_draws.npz')
    assert list(z['names']) == NAMES
    th = z['theta'][rng.choice(len(z['theta']), args.draws, replace=False)]; label = args.label
print('decision analysis on %d %s draws, %d candidate wells in %d states' % (len(th), label, N, S))


def draw_wells(p):
    a_s = rng.normal(p['a0'], p['tau_a'], S); m_s = rng.normal(p['m0'], p['tau_m'], S)
    d_s = rng.normal(p['d0'], p['tau_d'], S); u_c = rng.normal(0, p['tau_u'], d['n_cells'])
    sigma = np.exp(p['log_sigma0'])
    q = expit(a_s[st] + p['b_gas'] * gas)
    lp0 = d_s[st] + u_c[cell]
    s1, s0 = expit(lp0 + p['rho'] * rule), expit(lp0)
    q = q * (1 - s1) / (q * (1 - s1) + (1 - q) * (1 - s0))            # condition on not yet plugged
    zz = rng.random(N) < q
    e = np.where(zz, np.exp(rng.normal(m_s[st] + p['g_gas'] * gas, sigma)), 0.0)
    w0 = np.where(zz, s1, s0)                                           # P0 selection weights
    mean_e = expit(a_s + p['b_gas'] * np.array([[0.], [1.]])) * np.exp(m_s + p['g_gas'] * np.array([[0.], [1.]]) + sigma ** 2 / 2)
    return e, w0, mean_e                                                # mean_e: [gas, state]


def plug_top(score, frac):
    """Plug the top ceil(frac * n_s) wells of each state by score; returns a boolean mask."""
    order = np.lexsort((-score, st))
    start = np.concatenate([[0], np.cumsum(n_s)[:-1]])
    pos = np.arange(N) - start[st[order]]
    k = np.ceil(frac * n_s).astype(int)
    m = np.zeros(N, bool); m[order[pos < k[st[order]]]] = True
    return m


def screen(e, L, pf, sv, noise=True):
    """Reading of a one-visit screen; 0 if not found or below the limit."""
    r = e * (np.exp(rng.normal(-sv ** 2 / 2, sv, N)) if noise else 1.0)
    return np.where((rng.random(N) < pf) & (r > L), r, 0.0)


# P2 needs the posterior mean emission by (state, gas) before the loop
draws = [dict(zip(NAMES, t)) for t in th]
lat = [draw_wells(p) for p in draws]
mean_e = np.mean([m for _, _, m in lat], axis=0)
p2_score = mean_e[gas.astype(int), st]

# formula-grant-era (period 1) method mix per state, for E6
mix1 = d['method_mix'][:, 1]
Lm = np.array([METHODS[k][1] or np.inf for k in range(len(METHODS))])
two = np.array([METHODS[k][2] for k in range(len(METHODS))])

rows, be, e6, bys = [], [], [], []
for i, (p, (e, w0, _)) in enumerate(zip(draws, lat)):
    tie = rng.random(N) * 1e-9
    tot = e.sum()
    for f in FRACS:
        k = np.ceil(f * n_s).astype(int)
        # P0: weighted sampling without replacement within state (Efraimidis-Spirakis keys)
        m0 = plug_top(np.log(rng.random(N)) / w0, f)
        m1 = plug_top(rng.random(N), f)
        m2 = plug_top(p2_score + tie, f)
        v = {'P0': e[m0].sum(), 'P1': e[m1].sum(), 'P2': e[m2].sum(), 'oracle': e[plug_top(e + tie, f)].sum()}
        # P3 at the Hi-Flow-like setting with this draw's p_find
        r = screen(e, 1.0, p['p_find'], p['s_v'])
        m3 = plug_top(np.where(r > 0, 1e6 + r, p2_score + tie), f)
        v['P3_hiflow'] = e[m3].sum()
        for pol, val in v.items():
            rows.append(dict(draw=i, frac=f, policy=pol, avoided_g_h=val, share_of_total=val / tot))
        if f == 0.1:
            for s in range(S):
                ins = st == s
                bys.append(dict(draw=i, state=d['states'][s], n=int(n_s[s]),
                                gain_per_screened_usd=(e[m3 & ins].sum() - e[m2 & ins].sum()) * VALUE / max(n_s[s], 1),
                                p3_over_p2=e[m3 & ins].sum() / max(e[m2 & ins].sum(), 1e-9)))
        for L in LGRID:
            for pf in PGRID:
                r = screen(e, L, pf, p['s_v'])
                m3 = plug_top(np.where(r > 0, 1e6 + r, p2_score + tie), f)
                be.append(dict(draw=i, frac=f, L=L, p_find=pf, c_star=(e[m3].sum() - v['P2']) * VALUE / N))
        # E6: screening with each state's 2023-07 to 2025-07 method mix vs P0, per plugged well
        cum = np.cumsum(mix1[st], axis=1)
        code = np.minimum((rng.random(N)[:, None] > cum).sum(1), len(METHODS) - 1)
        found = ~two[code] | (rng.random(N) < p['p_find'])
        r = e * np.exp(rng.normal(-p['s_v'] ** 2 / 2, p['s_v'], N))
        r = np.where((code > 0) & found & (r > Lm[code] * np.exp(rng.normal(0, .5, N))), r, 0.0)
        m3 = plug_top(np.where(r > 0, 1e6 + r, p2_score + tie), f)
        e6.append(dict(draw=i, frac=f, gain_g_h_per_plugged=(e[m3].sum() - v['P0']) / k.sum(),
                       gain_vs_P2_g_h_per_plugged=(e[m3].sum() - v['P2']) / k.sum()))
    if (i + 1) % 50 == 0:
        print('  %d draws done' % (i + 1))

q = lambda s: pd.Series({'median': s.median(), 'q05': s.quantile(.05), 'q95': s.quantile(.95)})
R = pd.DataFrame(rows)
R['t_per_yr'] = R.avoided_g_h * T_PER_GH
R['pv_usd_millions'] = R.avoided_g_h * VALUE / 1e6
pv = R.groupby(['frac', 'policy'])[['avoided_g_h', 't_per_yr', 'pv_usd_millions', 'share_of_total']].median()
pv.insert(0, 'label', label)
pv.to_csv(P / 'policy_value.csv', float_format='%.4g')
B = pd.DataFrame(be).groupby(['frac', 'L', 'p_find']).c_star.apply(q).unstack().reset_index()
B.insert(0, 'label', label); B.to_csv(P / 'breakeven.csv', index=False, float_format='%.4g')
E = pd.DataFrame(e6)
E['pv_usd_per_plugged'] = E.gain_g_h_per_plugged * VALUE
E['pv_usd_after_rule_millions'] = E.pv_usd_per_plugged * N_AFTER_RULE / 1e6
E6 = E.groupby('frac')[['gain_g_h_per_plugged', 'gain_vs_P2_g_h_per_plugged', 'pv_usd_per_plugged',
                        'pv_usd_after_rule_millions']].agg(['median', lambda s: s.quantile(.05), lambda s: s.quantile(.95)])
E6.columns = ['%s_%s' % (a, b.replace('<lambda_0>', 'q05').replace('<lambda_1>', 'q95')) for a, b in E6.columns]
E6.insert(0, 'label', label); E6.to_csv(P / 'e6_rule_change.csv', float_format='%.4g')
BS = pd.DataFrame(bys).groupby(['state', 'n'])[['gain_per_screened_usd', 'p3_over_p2']].median().reset_index()
BS.insert(0, 'label', label); BS.sort_values('gain_per_screened_usd', ascending=False).to_csv(P / 'by_state.csv', index=False, float_format='%.4g')

print('value per g/h avoided: $%.0f (PVF %.2f)' % (VALUE, PVF))
print(pv.round(3).to_string())
b1 = B[(B.frac == 0.1)]
for pf in PGRID:
    x = b1[b1.p_find == pf].set_index('L')['median']
    print('f=0.1 p_find=%.2f  c* at L=0.01: $%.0f  L=1: $%.0f  L=10: $%.0f  L=100: $%.0f'
          % (pf, np.interp(-2, np.log10(x.index), x), np.interp(0, np.log10(x.index), x),
             np.interp(1, np.log10(x.index), x), np.interp(2, np.log10(x.index), x)))
print(E6.round(1).to_string())

# ---- figure: break-even screening cost (reference palette, light mode; sequential ramp for p_find)
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
RAMP = {0.25: '#a9c8f0', 0.5: '#6aa1e6', 0.75: '#2a78d6', 1.0: '#1a4f94'}
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': GRID, 'axes.labelcolor': INK2, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.facecolor': SURF, 'figure.facecolor': SURF, 'font.family': 'DejaVu Sans'})
fig, ax = plt.subplots(figsize=(7.2, 4.4))
for pf in PGRID:
    x = b1[b1.p_find == pf]
    if pf == 1.0:
        ax.fill_between(x.L, x.q05, x.q95, color=RAMP[pf], alpha=0.15, lw=0)
    ax.plot(x.L, x['median'], color=RAMP[pf], lw=2, label='leak found with probability %.2f' % pf)
for L, name in ((0.001, 'chamber'), (1, 'Hi-Flow'), (10, 'handheld'), (50, 'OGI')):
    ax.axvline(L, color=GRID, lw=1, zorder=0)
    ax.text(L, ax.get_ylim()[1] * 0.98, name, fontsize=8, color=INK2, ha='center', va='top',
            bbox=dict(facecolor=SURF, edgecolor='none', pad=1))
ax.set_xscale('log'); ax.set_ylim(bottom=0)
ax.set_xlabel('Detection limit of the screening measurement (g/h)')
ax.set_ylabel('Break-even cost per screened well (US$)')
ax.set_title('What one screening visit is worth, budget = 10%% of candidates (%s)' % label, loc='left', color=INK, fontsize=11)
ax.grid(True, which='major', color=GRID, lw=0.8); ax.set_axisbelow(True)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: '${:,.0f}'.format(v)))
ax.legend(frameon=False, fontsize=9, loc='upper right')
fig.text(0.01, 0.01, 'Posterior median over %d draws; band: 90%% interval at p=1. SC-CH4 $%.0f/t, %d yr, decline %.1f%%/yr. '
         'Compared with ranking by state and well type.' % (len(th), SCC, YEARS, DECLINE * 100), fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(P / 'fig_breakeven.png', dpi=200); plt.close(fig)
