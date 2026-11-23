"""Forward simulator for the selection / emission / measurement model (docs/analysis_plan.md, Sec. 3).

SKELETON (v0.1): global parameters only; state effects are drawn inside the simulator as nuisance
latents. Attributes are reduced to a gas-well indicator and a 0.5-degree spatial cell. Measurement
regimes are represented by each state's empirical method mix among plugged wells (period split not
yet used). The same summary function is applied to simulated and real data, so the pair
(simulate, real_summaries) is what neural posterior estimation will be trained on.

Usage:
    design = build_design()          # once; cached to data/processed/sim_design.pkl
    theta = sample_prior(rng, n)
    s = simulate(theta[0], design, rng)
    s_obs = real_summaries(design)
"""
import pathlib, pickle, re, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from config import PROC  # noqa: E402

import numpy as np
import pandas as pd
from scipy.special import expit

# ---- parameters (name, prior) ------------------------------------------------------------
# Priors are weakly informative; see analysis_plan.md. Order fixes the theta vector layout.
PARAMS = [
    ('a0', 'normal', 0.0, 1.5),        # logit share leaking, average state
    ('tau_a', 'halfnormal', 1.0),      # between-state SD of a_s
    ('b_gas', 'normal', 0.0, 1.0),     # gas-well effect on leaking (logit)
    ('m0', 'normal', 0.0, 2.0),        # mean log leak rate (log g/h), average state
    ('tau_m', 'halfnormal', 1.0),      # between-state SD of m_s
    ('g_gas', 'normal', 0.0, 1.0),     # gas-well effect on log leak rate
    ('log_sigma0', 'normal', np.log(2.0), 0.4),   # log of within-state SD of log rate
    ('d0', 'normal', -2.5, 1.0),       # logit plugging probability, average state
    ('tau_d', 'halfnormal', 1.0),
    ('tau_u', 'halfnormal', 1.0),      # spatial-cell SD (selection clustering)
    ('rho', 'normal', 0.0, 1.0),       # selection on leaking, in states with leak-based rules
    ('p_find', 'beta', 4.0, 2.0),      # leak found by two-stage methods
    ('s_v', 'halfnormal_loc', 0.59, 0.2),  # single-visit log-noise SD
]
NAMES = [p[0] for p in PARAMS]

LEAK_RULE_STATES = {'NM', 'LA', 'NY', 'KY', 'TX', 'CO', 'OH', 'KS', 'MT', 'AK', 'IN', 'FL', 'MS'}

# measurement codes: (name, detection-limit midpoint g/h, two-stage leak finding?, quantitative?)
METHODS = {
    0: ('not measured', None, False, False),
    1: ('chamber / flow-through cell', 0.001, False, True),
    2: ('flow meter + GC', 0.01, False, True),
    3: ('hi-flow sampler', 1.0, True, True),
    4: ('laser analyzer (quantitative)', 1.0, False, True),
    5: ('federal guidelines (instrument unspecified)', 1.0, True, True),
    6: ('screening: OGI camera', 50.0, True, False),
    7: ('screening: handheld concentration detector', 10.0, False, False),
    8: ('screening: standoff laser (RMLD-type)', 10.0, False, False),
    9: ('FID (concentration)', 10.0, False, False),
}
GROUP_TO_CODE = {v[0]: k for k, v in METHODS.items() if k > 0}
GROUP_TO_CODE['ND (ambiguous: no data or not detected)'] = 3     # analysis plan: censored at Hi-Flow limit
QUANT = [k for k, v in METHODS.items() if v[3]]
SCREEN = [k for k, v in METHODS.items() if k > 0 and not v[3]]


def is_gas(s):
    t = str(s).upper()
    return bool(re.search(r'GAS|CBM|^G$', t)) and 'OIL' not in t


def sample_prior(rng, n):
    out = np.empty((n, len(PARAMS)))
    for j, p in enumerate(PARAMS):
        kind = p[1]
        if kind == 'normal':
            out[:, j] = rng.normal(p[2], p[3], n)
        elif kind == 'halfnormal':
            out[:, j] = np.abs(rng.normal(0, p[2], n))
        elif kind == 'halfnormal_loc':
            out[:, j] = np.abs(rng.normal(p[2], p[3], n))
        elif kind == 'beta':
            out[:, j] = rng.beta(p[2], p[3], n)
    return out


# ---- design: everything about the real data that the simulator conditions on --------------
def build_design(cache=True):
    f = PROC / 'sim_design.pkl'
    if cache and f.exists():
        return pickle.load(open(f, 'rb'))
    inv = pd.read_pickle(PROC / 'inv.pkl')
    inv = inv[inv.lat.notna() & inv.st.notna()].copy()
    a = pd.read_csv(PROC / 'doi_wells_analysis.csv', low_memory=False)
    g = pd.read_csv(PROC / 'method_group.csv')
    a = a.merge(g[['OBJECTID', 'method_group']], on='OBJECTID', how='left')
    a['code'] = a.method_group.map(GROUP_TO_CODE).fillna(0).astype(int)
    a.loc[a.methane_pre_plugging_emissions.isna(), 'code'] = 0
    # NM field averages (D21) are not measurements
    v = a.methane_pre_plugging_emissions
    nm_rep = v[(a.st == 'NM') & (v > 0)].round(4).value_counts()
    a.loc[(a.st == 'NM') & v.round(4).isin(nm_rep[nm_rep >= 3].index), 'code'] = 0

    states = sorted(set(inv.st) & set(a.st))
    sidx = {s: i for i, s in enumerate(states)}
    inv = inv[inv.st.isin(sidx)]
    a = a[a.st.isin(sidx)]
    cell = (np.floor(inv.lat / 0.5).astype(int).astype(str) + '_' + np.floor(inv.lon / 0.5).astype(int).astype(str)
            + '_' + inv.st)
    cell_codes, cell_id = np.unique(cell.values, return_inverse=True)
    mix = np.zeros((len(states), len(METHODS)))
    for s, grp in a.groupby('st'):
        mix[sidx[s]] = np.bincount(grp.code, minlength=len(METHODS)) / len(grp)
    design = dict(
        states=states,
        inv_state=inv.st.map(sidx).values.astype(int),
        inv_gas=inv.Type.map(is_gas).values.astype(float),
        inv_cell=cell_id.astype(int), n_cells=len(cell_codes),
        inv_cell_state=np.array([sidx[c.rsplit('_', 1)[1]] for c in cell_codes]),
        inv_plugged=inv.plugged.values.astype(bool),
        method_mix=mix,
        leak_rule=np.array([s in LEAK_RULE_STATES for s in states]),
        doi_state=a.st.map(sidx).values.astype(int),
        doi_gas=a.well_type.map(is_gas).values.astype(float),
        doi_code=a.code.values.astype(int),
        doi_value=a.methane_pre_plugging_emissions.fillna(0).values,
    )
    if cache:
        pickle.dump(design, open(f, 'wb'))
    return design


# ---- summaries ---------------------------------------------------------------------------
Q = (0.1, 0.5, 0.9)
EMPTY = -12.0     # sentinel for an undefined quantile (log g/h); paired with a count summary


def _summaries(n_states, inv_state, inv_cell, inv_cell_state, n_cells, plugged,
               w_state, w_gas, w_code, w_value):
    """Same function for real and simulated data. Plugging summaries use the inventory; measurement
    summaries use plugged wells (w_*: one entry per plugged well)."""
    out = []
    n_inv = np.bincount(inv_state, minlength=n_states)
    out.append(np.bincount(inv_state, weights=plugged, minlength=n_states) / np.maximum(n_inv, 1))
    # spatial clustering: SD of plugged share across cells with >= 20 wells, per state
    cn = np.bincount(inv_cell, minlength=n_cells)
    cp = np.bincount(inv_cell, weights=plugged, minlength=n_cells) / np.maximum(cn, 1)
    sd = np.zeros(n_states)
    for s in range(n_states):
        m = (inv_cell_state == s) & (cn >= 20)
        sd[s] = cp[m].std() if m.sum() >= 2 else 0.0
    out.append(sd)
    n_w = np.bincount(w_state, minlength=n_states)
    measured = w_code > 0
    out.append(np.bincount(w_state, weights=measured, minlength=n_states) / np.maximum(n_w, 1))
    q = np.isin(w_code, QUANT)
    sc = np.isin(w_code, SCREEN)
    det = w_value > 0
    nq = np.bincount(w_state, weights=q, minlength=n_states)
    ns = np.bincount(w_state, weights=sc, minlength=n_states)
    out.append(np.bincount(w_state, weights=q & det, minlength=n_states) / np.maximum(nq, 1))
    out.append(np.bincount(w_state, weights=sc & det, minlength=n_states) / np.maximum(ns, 1))
    lq = np.full((n_states, len(Q)), EMPTY)
    lv = np.log(np.where(det, w_value, 1.0))
    for s in range(n_states):
        m = (w_state == s) & q & det
        if m.sum() >= 3:
            lq[s] = np.quantile(lv[m], Q)
    out.append(lq.ravel())
    # pooled: detection among quantitative measurements, gas vs other wells; pooled positive quantiles
    pooled = []
    for gmask in (w_gas == 1, w_gas == 0):
        m = q & gmask
        pooled.append((m & det).sum() / max(m.sum(), 1))
    m = q & det
    pooled += list(np.quantile(lv[m], (0.1, 0.5, 0.9, 0.99))) if m.sum() >= 5 else [EMPTY] * 4
    out.append(np.array(pooled))
    # between-state dispersion (informs the tau parameters): SD across states with enough wells
    lg = lambda x: np.log((x + 0.01) / (1.01 - x))
    ok_inv = n_inv >= 200
    ok_q = nq >= 20
    med = lq[:, 1]
    ok_m = ok_q & (med > EMPTY)
    plug_share = out[0]
    det_q = out[3]
    out.append(np.array([lg(plug_share[ok_inv]).std() if ok_inv.sum() >= 3 else 0.0,
                         lg(det_q[ok_q]).std() if ok_q.sum() >= 3 else 0.0,
                         med[ok_m].std() if ok_m.sum() >= 3 else 0.0]))
    return np.concatenate(out)


def summary_names(design):
    st = design['states']
    names = [f'plug_share_{s}' for s in st] + [f'cell_sd_{s}' for s in st] + [f'measured_{s}' for s in st]
    names += [f'detect_quant_{s}' for s in st] + [f'detect_screen_{s}' for s in st]
    names += [f'logq{int(q * 100)}_{s}' for s in st for q in Q]
    names += ['detect_gas', 'detect_other', 'pooled_logq10', 'pooled_logq50', 'pooled_logq90', 'pooled_logq99']
    names += ['between_state_sd_plug', 'between_state_sd_detect', 'between_state_sd_logq50']
    return names


def real_summaries(design):
    d = design
    return _summaries(len(d['states']), d['inv_state'], d['inv_cell'], d['inv_cell_state'], d['n_cells'],
                      d['inv_plugged'].astype(float), d['doi_state'], d['doi_gas'], d['doi_code'], d['doi_value'])


# ---- simulator ---------------------------------------------------------------------------
def simulate(theta, design, rng, return_latent=False):
    p = dict(zip(NAMES, theta))
    d = design
    S = len(d['states'])
    a_s = rng.normal(p['a0'], p['tau_a'], S)
    m_s = rng.normal(p['m0'], p['tau_m'], S)
    d_s = rng.normal(p['d0'], p['tau_d'], S)
    u_c = rng.normal(0, p['tau_u'], d['n_cells'])
    sigma = np.exp(p['log_sigma0'])

    st, gas = d['inv_state'], d['inv_gas']
    n = len(st)
    z = rng.random(n) < expit(a_s[st] + p['b_gas'] * gas)                       # leaking
    e = np.where(z, np.exp(rng.normal(m_s[st] + p['g_gas'] * gas, sigma)), 0.0)  # true g/h
    lp = d_s[st] + u_c[d['inv_cell']] + p['rho'] * z * d['leak_rule'][st]
    plugged = rng.random(n) < expit(lp)

    # measurement of plugged wells: method from the state's empirical mix
    w = np.flatnonzero(plugged)
    ws = st[w]
    cum = np.cumsum(d['method_mix'][ws], axis=1)
    code = (rng.random(len(w))[:, None] > cum).sum(axis=1)
    code = np.minimum(code, len(METHODS) - 1)
    L = np.array([METHODS[k][1] or np.inf for k in range(len(METHODS))])
    L = L * np.exp(rng.normal(0, 0.5, len(L)))                                 # limit uncertainty
    two_stage = np.array([METHODS[k][2] for k in range(len(METHODS))])
    found = ~two_stage[code] | (rng.random(len(w)) < p['p_find'])
    reading = e[w] * np.exp(rng.normal(-p['s_v'] ** 2 / 2, p['s_v'], len(w)))
    value = np.where((code > 0) & found & (reading > L[code]), reading, 0.0)

    s = _summaries(S, st, d['inv_cell'], d['inv_cell_state'], d['n_cells'], plugged.astype(float),
                   ws, gas[w], code, value)
    if return_latent:
        return s, dict(z=z, e=e, plugged=plugged)
    return s
