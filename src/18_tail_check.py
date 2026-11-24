"""Tail family check (analysis_plan.md Sec. 3 "Tail" and Sec. 6 item 5), on the reduced exact-likelihood
model of src/17 (no selection, no reporting error; s_v = 0.59 fixed).

Families for log leak rate of leaking wells (state shift m_s + g_gas * gas in all):
  lognormal  N(mu, sigma)                                  (primary)
  logt       Student-t(nu, mu, sigma), nu ~ Gamma(2, 0.1)  (pre-specified alternative)
  mix2       (1 - pi_h) N(mu, sigma1) + pi_h N(mu + delta, sigma2), pi_h ~ Beta(1, 9), delta ~ HalfNormal(4)
             (added: a low/high-emitter mixture, as in the super-emitter literature; D35)
For log-t the single-visit noise is folded into the scale (t convolved with a normal has no closed
form); for the normal families it is added exactly.

Checks:
  1. 5-fold held-out (random wells, seed 18): log predictive density of the held-out records, and for
     each fold the predicted vs observed number of held-out records above 1 kg/h and above 10 kg/h
     (the top ~1% and ~0.1%), with 90% predictive intervals.
  2. Driscoll et al. (132 unplugged LA wells, chamber, no selection), predicted from the full fit:
     shares above 10, 100, 1000 g/h and the maximum.
  3. Implied mean leak rate of a leaking well in an average state (what drives the decision values).
Outputs: outputs/tail_check/{heldout.csv, driscoll.csv, summary.csv, posterior.csv}
Run from the repo root: python src/18_tail_check.py
"""
import sys, pathlib, time; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
from sim.model import build_design, METHODS

import numpy as np, pandas as pd
import jax, jax.numpy as jnp
import numpyro, numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS, init_to_median
from jax.scipy.stats import norm as jnorm, t as jt
from jax.scipy.special import logsumexp
from scipy.special import expit
P = OUT / 'tail_check'; P.mkdir(exist_ok=True)
SV = 0.59
FAMS = ('lognormal', 'logt', 'mix2')

d = build_design()
keep = d['doi_code'] > 0
ST = d['doi_state'][keep]; GAS = d['doi_gas'][keep]; CODE = d['doi_code'][keep]; Y = d['doi_value'][keep]
S = len(d['states'])
_GX, _GW = np.polynomial.legendre.leggauss(96)
_W, _X = jnp.array(_GW / 2), jnp.array((_GX + 1) / 2)


def t_sf(x, nu):
    """Survival function of Student-t by Gauss-Legendre quadrature (JAX has no gradient of betainc in nu)."""
    a = jnp.pi * _X / 2
    u = x[..., None] + jnp.tan(a)
    return (jnp.exp(jt.logpdf(u, nu)) * (jnp.pi / 2) / jnp.cos(a) ** 2 * _W).sum(-1)


LIM = np.array([METHODS[k][1] or np.inf for k in range(len(METHODS))])
TWO = np.array([METHODS[k][2] for k in range(len(METHODS))])


def loglik(p, fam, st, gas, code, y):
    a_s, m_s = p['a0'] + p['tau_a'] * p['za'], p['m0'] + p['tau_m'] * p['zm']
    eta = a_s[st] + p['b_gas'] * gas
    lq, l1q = jax.nn.log_sigmoid(eta), jax.nn.log_sigmoid(-eta)
    mu = m_s[st] + p['g_gas'] * gas - SV ** 2 / 2
    pfi = jnp.where(TWO[code], p['p_find'], 1.0)
    pos = y > 0
    ly = jnp.log(jnp.where(pos, y, 1.0)); lL = jnp.log(LIM[code])
    if fam == 'lognormal':
        sd = jnp.sqrt(jnp.exp(2 * p['log_sigma0']) + SV ** 2)
        lpdf = jnorm.logpdf(ly, mu, sd); sf = jnorm.sf((lL - mu) / sd)
    elif fam == 'logt':
        sd = jnp.sqrt(jnp.exp(2 * p['log_sigma0']) + SV ** 2)
        lpdf = jt.logpdf((ly - mu) / sd, p['nu']) - jnp.log(sd); sf = t_sf((lL - mu) / sd, p['nu'])
    else:
        s1 = jnp.sqrt(jnp.exp(2 * p['log_sigma0']) + SV ** 2); s2 = jnp.sqrt(p['sigma2'] ** 2 + SV ** 2)
        lpdf = jnp.logaddexp(jnp.log1p(-p['pi_h']) + jnorm.logpdf(ly, mu, s1),
                             jnp.log(p['pi_h']) + jnorm.logpdf(ly, mu + p['delta'], s2))
        sf = (1 - p['pi_h']) * jnorm.sf((lL - mu) / s1) + p['pi_h'] * jnorm.sf((lL - mu - p['delta']) / s2)
    lpos = lq + jnp.log(pfi) + lpdf - ly
    lzero = jnp.logaddexp(l1q, lq + jnp.log1p(-jnp.clip(pfi * sf, 0, 1 - 1e-12)))
    return jnp.where(pos, lpos, lzero)


def model(fam, st, gas, code, y):
    p = dict(a0=numpyro.sample('a0', dist.Normal(0, 1.5)), tau_a=numpyro.sample('tau_a', dist.HalfNormal(1.0)),
             b_gas=numpyro.sample('b_gas', dist.Normal(0, 1)), m0=numpyro.sample('m0', dist.Normal(0, 2)),
             tau_m=numpyro.sample('tau_m', dist.HalfNormal(1.0)), g_gas=numpyro.sample('g_gas', dist.Normal(0, 1)),
             log_sigma0=numpyro.sample('log_sigma0', dist.Normal(np.log(2.0), 0.4)),
             p_find=numpyro.sample('p_find', dist.Beta(4, 2)))
    if fam == 'logt':
        p['nu'] = numpyro.sample('nu', dist.Gamma(2, 0.1))
    if fam == 'mix2':
        p['pi_h'] = numpyro.sample('pi_h', dist.Beta(1, 9)); p['delta'] = numpyro.sample('delta', dist.HalfNormal(4))
        p['sigma2'] = numpyro.sample('sigma2', dist.HalfNormal(2))
    with numpyro.plate('s', S):
        p['za'] = numpyro.sample('za', dist.Normal(0, 1)); p['zm'] = numpyro.sample('zm', dist.Normal(0, 1))
    numpyro.factor('lik', loglik(p, fam, st, gas, code, y).sum())


def fit(fam, idx, seed, n=600):
    m = MCMC(NUTS(model, target_accept_prob=0.9, init_strategy=init_to_median), num_warmup=n, num_samples=n,
             num_chains=2, chain_method='sequential', progress_bar=False)
    m.run(jax.random.PRNGKey(seed), fam=fam, st=jnp.array(ST[idx]), gas=jnp.array(GAS[idx]),
          code=jnp.array(CODE[idx]), y=jnp.array(Y[idx]))
    return {k: np.asarray(v) for k, v in m.get_samples().items()}, int(m.get_extra_fields()['diverging'].sum())


def draw_log_rate(p, i, fam, mu, rng):
    n = len(mu)
    if fam == 'lognormal':
        return rng.normal(mu, np.exp(p['log_sigma0'][i]), n)
    if fam == 'logt':
        return mu + np.exp(p['log_sigma0'][i]) * rng.standard_t(p['nu'][i], n)
    h = rng.random(n) < p['pi_h'][i]
    return np.where(h, rng.normal(mu + p['delta'][i], p['sigma2'][i], n), rng.normal(mu, np.exp(p['log_sigma0'][i]), n))


def predict(p, fam, st, gas, code, rng, nd=400, limit_override=None):
    out = []
    for i in rng.choice(len(p['a0']), nd, replace=False):
        a_s = p['a0'][i] + p['tau_a'][i] * p['za'][i]; m_s = p['m0'][i] + p['tau_m'][i] * p['zm'][i]
        n = len(st)
        z = rng.random(n) < expit(a_s[st] + p['b_gas'][i] * gas)
        le = draw_log_rate(p, i, fam, m_s[st] + p['g_gas'][i] * gas, rng)
        r = np.exp(np.minimum(le, 30) + rng.normal(-SV ** 2 / 2, SV, n))
        L = LIM[code] if limit_override is None else limit_override
        found = (~TWO[code] | (rng.random(n) < p['p_find'][i])) if limit_override is None else np.ones(n, bool)
        out.append(np.where(z & found & (r > L), r, 0.0))
    return np.array(out)


rng = np.random.default_rng(18)
t0 = time.time()
fold = rng.permutation(len(Y)) % 5
H, DR, SUM, POST = [], [], [], []
for fam in FAMS:
    elpd, nd_tot = 0.0, 0
    for k in range(5):
        tr, te = fold != k, fold == k
        s, div = fit(fam, tr, seed=100 + k)
        nd_tot += div
        ll = jax.vmap(lambda pp: loglik(pp, fam, jnp.array(ST[te]), jnp.array(GAS[te]), jnp.array(CODE[te]), jnp.array(Y[te])))(
            {kk: jnp.array(v) for kk, v in s.items()})
        elpd += float((logsumexp(ll, axis=0) - np.log(ll.shape[0])).sum())
        yp = predict(s, fam, ST[te], GAS[te], CODE[te], rng)
        for thr in (1000, 10000):
            c = (yp > thr).sum(1); lo, hi = np.quantile(c, [.05, .95])
            ob = int((Y[te] > thr).sum())
            H.append(dict(family=fam, fold=k, threshold_g_h=thr, observed=ob, pred_median=np.median(c), pred_q05=lo,
                          pred_q95=hi, inside_90=bool(lo <= ob <= hi), divergences=div))
        print('%s fold %d done (%d divergences, %.0f s)' % (fam, k, div, time.time() - t0), flush=True)
    s, div = fit(fam, np.ones(len(Y), bool), seed=7)
    for kk in ('a0', 'tau_a', 'm0', 'tau_m', 'log_sigma0', 'p_find', 'nu', 'pi_h', 'delta', 'sigma2'):
        if kk in s:
            POST.append(dict(family=fam, parameter=kk, median=np.median(s[kk]), q05=np.quantile(s[kk], .05), q95=np.quantile(s[kk], .95)))
    # Driscoll: LA, chamber (limit 0.001 g/h, no leak-finding stage)
    k_la = d['dr_state']; dg = d['dr_gas']; dv = d['dr_value']
    yp = predict(s, fam, np.full(len(dg), k_la), dg, np.full(len(dg), 1), rng, limit_override=0.001)
    for thr in (10, 100, 1000):
        sh = (yp > thr).mean(1); lo, hi = np.quantile(sh, [.05, .95])
        DR.append(dict(family=fam, statistic='share > %d g/h' % thr, observed=(dv > thr).mean(), pred_median=np.median(sh),
                       pred_q05=lo, pred_q95=hi, inside_90=bool(lo <= (dv > thr).mean() <= hi)))
    mx = yp.max(1); lo, hi = np.quantile(mx, [.05, .95])
    DR.append(dict(family=fam, statistic='max g/h', observed=dv.max(), pred_median=np.median(mx), pred_q05=lo, pred_q95=hi,
                   inside_90=bool(lo <= dv.max() <= hi)))
    # implied mean rate of a leaking well, average state, non-gas (Monte Carlo, capped at e^30 to avoid overflow)
    means = []
    for i in rng.choice(len(s['a0']), 400, replace=False):
        le = draw_log_rate(s, i, fam, np.full(20000, s['m0'][i]), rng)
        means.append(np.exp(np.minimum(le, 30)).mean())
    hh = pd.DataFrame([h for h in H if h['family'] == fam])
    SUM.append(dict(family=fam, elpd_heldout=elpd, divergences_folds=nd_tot, divergences_full=div,
                    tail_inside_90=int(hh.inside_90.sum()), tail_checks=len(hh),
                    driscoll_inside_90=int(sum(r['inside_90'] for r in DR if r['family'] == fam)),
                    mean_leaking_g_h_median=np.median(means), mean_leaking_g_h_q05=np.quantile(means, .05),
                    mean_leaking_g_h_q95=np.quantile(means, .95)))
    print(pd.DataFrame(SUM).round(2).to_string(index=False), flush=True)

pd.DataFrame(H).to_csv(P / 'heldout.csv', index=False, float_format='%.4g')
pd.DataFrame(DR).to_csv(P / 'driscoll.csv', index=False, float_format='%.4g')
pd.DataFrame(SUM).to_csv(P / 'summary.csv', index=False, float_format='%.4g')
pd.DataFrame(POST).to_csv(P / 'posterior.csv', index=False, float_format='%.4g')
print(pd.DataFrame(DR).round(3).to_string(index=False))
print(pd.DataFrame(H).groupby(['family', 'threshold_g_h'])[['observed', 'pred_median', 'inside_90']].sum().to_string())
print('total %.0f s' % (time.time() - t0))
