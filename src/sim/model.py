"""Forward simulator for the selection / emission / measurement model (docs/analysis_plan.md, Sec. 3).

v0.2: global parameters only; state effects are drawn inside the simulator as nuisance latents.
Attributes are reduced to a gas-well indicator and a 0.5-degree spatial cell. Measurement regimes are
each state's empirical method mix among plugged wells in each of three periods (split at 2023-07-07 and
2025-07-11). Includes the reporting-error term (D28), LA Hi-Flow quantization (D22), and the external
data: Driscoll chamber sample and the 4 usable Driscoll-DOI pairs, the Montana repeat visits, and the
New Mexico memo-record pairs. The same summary function is applied to simulated and real data, so the pair
(simulate, real_summaries) is what neural posterior estimation will be trained on.

Usage:
    design = build_design()          # once; cached to data/processed/sim_design.pkl
    theta = sample_prior(rng, n)
    s = simulate(theta[0], design, rng)
    s_obs = real_summaries(design)
"""
import pathlib, pickle, re, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from config import PROC, RAW, OUT, DRISCOLL  # noqa: E402

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
    ('omega', 'beta', 1.0, 30.0),      # reporting error: measured leak recorded as zero
    ('kappa', 'beta', 1.0, 9.0),       # reporting error: recorded value off by a random factor
    ('s_r', 'halfnormal', 0.7),        # log-SD of that factor
]
PERIOD_BREAKS = ('2023-07-07', '2025-07-11')
DESIGN_VERSION = 2   # bump whenever build_design changes; a stale cache is rebuilt
LA_STEP = 1.1504   # g/h, LA Hi-Flow reporting resolution (D22)
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
        d = pickle.load(open(f, 'rb'))
        if d.get('version') == DESIGN_VERSION:
            return d
    inv = pd.read_pickle(PROC / 'inv.pkl')
    inv = inv[inv.lat.notna() & inv.st.notna()].copy()
    a = pd.read_csv(PROC / 'doi_wells_analysis.csv', low_memory=False)
    g = pd.read_csv(PROC / 'method_group.csv')
    a = a.merge(g[['OBJECTID', 'method_group']], on='OBJECTID', how='left')
    a['code'] = a.method_group.map(GROUP_TO_CODE).fillna(0).astype(int)
    dt = pd.to_datetime(a.witness_date, unit='ms', errors='coerce')
    a['period'] = np.select([dt < PERIOD_BREAKS[0], dt < PERIOD_BREAKS[1]], [0, 1], 2)
    a.loc[dt.isna(), 'period'] = 1
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
    mix = np.zeros((len(states), 3, len(METHODS)))
    pshare = np.zeros((len(states), 3))
    for s, grp in a.groupby('st'):
        pshare[sidx[s]] = np.bincount(grp.period, minlength=3) / len(grp)
        for t, gt in grp.groupby('period'):
            mix[sidx[s], t] = np.bincount(gt.code, minlength=len(METHODS)) / len(gt)
    for i in range(len(states)):                    # a period with no wells borrows the state's overall mix
        tot = (mix[i] * pshare[i][:, None]).sum(0)
        for t in range(3):
            if mix[i, t].sum() == 0:
                mix[i, t] = tot
    design = dict(
        version=DESIGN_VERSION,
        states=states,
        inv_state=inv.st.map(sidx).values.astype(int),
        inv_gas=inv.Type.map(is_gas).values.astype(float),
        inv_cell=cell_id.astype(int), n_cells=len(cell_codes),
        inv_cell_state=np.array([sidx[c.rsplit('_', 1)[1]] for c in cell_codes]),
        inv_plugged=inv.plugged.values.astype(bool),
        method_mix=mix, period_share=pshare,
        doi_period=a.period.values.astype(int),
        leak_rule=np.array([s in LEAK_RULE_STATES for s in states]),
        doi_state=a.st.map(sidx).values.astype(int),
        doi_gas=a.well_type.map(is_gas).values.astype(float),
        doi_code=a.code.values.astype(int),
        doi_value=a.methane_pre_plugging_emissions.fillna(0).values,
    )
    # external data (analysis plan Sec. 3): Driscoll chamber sample of unplugged LA wells, and the
    # USGS Montana repeat visits (used only for the single-visit noise s_v)
    dr = pd.read_excel(DRISCOLL).iloc[1:]
    design['dr_state'] = sidx['LA']
    design['dr_gas'] = dr['Product Type'].astype(str).str.upper().str.contains('GAS').values.astype(float)
    design['dr_value'] = pd.to_numeric(dr['CH4 Flow Rate'], errors='coerce').fillna(0).values
    pr = pd.read_csv(OUT / 'driscoll_doi_pairs.csv').dropna(subset=['doi_hiflow_g_h'])
    design['pair_chamber'] = pr.chamber_g_h.values
    design['pair_doi'] = pr.doi_hiflow_g_h.values
    nmm = pd.read_csv(OUT / 'nm_memo_audit.csv')
    nmm = nmm.dropna(subset=['memo_g_h', 'methane_pre_plugging_emissions'])
    nmm = nmm[nmm.field_average != True]
    design['nm_memo'] = nmm.memo_g_h.values
    design['nm_record'] = nmm.methane_pre_plugging_emissions.values
    design['nm_state'] = sidx['NM']
    mt = pd.read_csv(RAW / 'usgs_mt_kevin_sunburst' / 'KD_OW_Emissions_Data.csv', encoding='utf-8-sig')
    design['mt_sessions'] = mt.groupby('API').size().values
    design['mt_logsd'] = float(np.median([np.log(g).std(ddof=1) for _, g in mt.groupby('API').Average_g_h if len(g) >= 3]))
    if cache:
        pickle.dump(design, open(f, 'wb'))
    return design


# ---- summaries ---------------------------------------------------------------------------
Q = (0.1, 0.5, 0.9)
EMPTY = -12.0     # sentinel for an undefined quantile (log g/h); paired with a count summary


def _summaries(n_states, inv_state, inv_cell, inv_cell_state, n_cells, plugged,
               w_state, w_gas, w_code, w_value, w_period):
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
    # by period (pooled over states): share measured, detection among quantitative measurements
    per = []
    for t in range(3):
        m = w_period == t
        per.append(measured[m].mean() if m.any() else 0.0)
        mq = m & q
        per.append((mq & det).sum() / max(mq.sum(), 1))
    out.append(np.array(per))
    return np.concatenate(out)


DR_CONTROL = 0.0001   # g/h; Driscoll et al. control-blank level


def _pair_nm_summaries(pair_chamber, pair_doi, nm_memo, nm_record):
    big = pair_chamber > 10
    s = [big.sum(), (big & (pair_doi == 0)).sum()]
    m1 = nm_memo > 1
    s.append((m1 & (nm_record == 0)).sum() / max(m1.sum(), 1))
    both = (nm_memo > 0.1) & (nm_record > 0.1)
    lr = np.log(nm_record[both] / nm_memo[both]) if both.any() else np.array([0.0])
    off = np.abs(lr) > 0.1
    s += [off.mean(), lr[off].std() if off.sum() >= 3 else 0.0]
    return np.array(s, dtype=float)


def _external_summaries(dr_value, mt_logsd, la_hiflow_detect):
    det = dr_value > DR_CONTROL
    lv = np.log(dr_value[det]) if det.sum() >= 3 else None
    q = list(np.quantile(lv, (0.1, 0.5, 0.9))) if lv is not None else [EMPTY] * 3
    # leak-finding contrast: LA Hi-Flow detection share relative to the chamber share above 1 g/h
    above1 = (dr_value > 1.0).mean()
    ratio = la_hiflow_detect / max(above1, 1 / len(dr_value))
    return np.array([det.mean()] + q + [mt_logsd, ratio])


def _la_hiflow_detect(state, code, value, la):
    m = (state == la) & (code == 3)
    return (value[m] > 0).mean() if m.sum() else 0.0


def summary_names(design):
    st = design['states']
    names = [f'plug_share_{s}' for s in st] + [f'cell_sd_{s}' for s in st] + [f'measured_{s}' for s in st]
    names += [f'detect_quant_{s}' for s in st] + [f'detect_screen_{s}' for s in st]
    names += [f'logq{int(q * 100)}_{s}' for s in st for q in Q]
    names += ['detect_gas', 'detect_other', 'pooled_logq10', 'pooled_logq50', 'pooled_logq90', 'pooled_logq99']
    names += ['between_state_sd_plug', 'between_state_sd_detect', 'between_state_sd_logq50']
    names += [f'{k}_period{t}' for t in range(3) for k in ('measured', 'detect_quant')]
    names += ['driscoll_detect', 'driscoll_logq10', 'driscoll_logq50', 'driscoll_logq90', 'mt_within_well_logsd',
              'la_hiflow_vs_chamber_ratio']
    names += ['pairs_chamber_gt10', 'pairs_chamber_gt10_doi_zero', 'nm_zeroed_share', 'nm_off_share', 'nm_off_logsd']
    return names


def real_summaries(design):
    d = design
    s = _summaries(len(d['states']), d['inv_state'], d['inv_cell'], d['inv_cell_state'], d['n_cells'],
                   d['inv_plugged'].astype(float), d['doi_state'], d['doi_gas'], d['doi_code'], d['doi_value'],
                   d['doi_period'])
    la = _la_hiflow_detect(d['doi_state'], d['doi_code'], d['doi_value'], d['dr_state'])
    return np.concatenate([s, _external_summaries(d['dr_value'], d['mt_logsd'], la),
                           _pair_nm_summaries(d['pair_chamber'], d['pair_doi'], d['nm_memo'], d['nm_record'])])


# ---- simulator ---------------------------------------------------------------------------
def report(value, p, rng):
    """Reporting error between the reading and the program record (D28): a measured leak is recorded as
    zero with probability omega; with probability kappa the value is off by a lognormal factor."""
    n = len(value)
    u = rng.random(n)
    zero = u < p['omega']
    off = (u >= p['omega']) & (u < p['omega'] + p['kappa'])
    out = np.where(off, value * np.exp(rng.normal(0, p['s_r'], n)), value)
    return np.where(zero, 0.0, out)


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
    pcum = np.cumsum(d['period_share'][ws], axis=1)
    wper = np.minimum((rng.random(len(w))[:, None] > pcum).sum(axis=1), 2)
    cum = np.cumsum(d['method_mix'][ws, wper], axis=1)
    code = (rng.random(len(w))[:, None] > cum).sum(axis=1)
    code = np.minimum(code, len(METHODS) - 1)
    L = np.array([METHODS[k][1] or np.inf for k in range(len(METHODS))])
    L = L * np.exp(rng.normal(0, 0.5, len(L)))                                 # limit uncertainty
    two_stage = np.array([METHODS[k][2] for k in range(len(METHODS))])
    found = ~two_stage[code] | (rng.random(len(w)) < p['p_find'])
    reading = e[w] * np.exp(rng.normal(-p['s_v'] ** 2 / 2, p['s_v'], len(w)))
    value = np.where((code > 0) & found & (reading > L[code]), reading, 0.0)
    value = report(value, p, rng)
    la_hf = (ws == d['dr_state']) & (code == 3)
    value[la_hf] = np.round(value[la_hf] / LA_STEP) * LA_STEP

    s = _summaries(S, st, d['inv_cell'], d['inv_cell_state'], d['n_cells'], plugged.astype(float),
                   ws, gas[w], code, value, wper)
    # Driscoll: chamber (no leak-finding stage, limit ~0.001 g/h) on a sample of unplugged LA wells
    k = d['dr_state']; dg = d['dr_gas']; nd = len(dg)
    dz = rng.random(nd) < expit(a_s[k] + p['b_gas'] * dg)
    de = np.where(dz, np.exp(rng.normal(m_s[k] + p['g_gas'] * dg, sigma)), 0.0)
    dread = de * np.exp(rng.normal(-p['s_v'] ** 2 / 2, p['s_v'], nd))
    dval = np.where(dread > 0.001 * np.exp(rng.normal(0, 0.5)), dread, 0.0)
    # Montana repeat visits: within-well SD of log readings depends only on s_v
    sds = [rng.normal(0, p['s_v'], n_k).std(ddof=1) for n_k in d['mt_sessions'] if n_k >= 3]
    la = _la_hiflow_detect(ws, code, value, k)
    noise = lambda n: np.exp(rng.normal(-p['s_v'] ** 2 / 2, p['s_v'], n))
    # Driscoll-DOI pairs: chamber reading, then a program Hi-Flow reading at a later visit
    npair = len(d['pair_chamber'])
    pz = rng.random(npair) < expit(a_s[k])
    pe = np.where(pz, np.exp(rng.normal(m_s[k], sigma, npair)), 0.0)
    pch = pe * noise(npair)
    pch = np.where(pch > 0.001, pch, 0.0)
    phf = pe * noise(npair)
    phf = np.where((rng.random(npair) < p['p_find']) & (phf > L[3]), phf, 0.0)
    phf = np.round(report(phf, p, rng) / LA_STEP) * LA_STEP
    # NM memo-record pairs: flow-meter reading (memo) and the program record after reporting error
    kn = d['nm_state']; nn = len(d['nm_memo'])
    nz = rng.random(nn) < expit(a_s[kn])
    ne = np.where(nz, np.exp(rng.normal(m_s[kn], sigma, nn)), 0.0) * noise(nn)
    memo = np.where(ne > L[2], ne, 0.0)
    rec = report(memo, p, rng)
    s = np.concatenate([s, _external_summaries(dval, float(np.median(sds)), la),
                        _pair_nm_summaries(pch, phf, memo, rec)])
    if return_latent:
        return s, dict(z=z, e=e, plugged=plugged)
    return s
