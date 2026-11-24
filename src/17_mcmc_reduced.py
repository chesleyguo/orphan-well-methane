"""Exact-likelihood check of the inference (analysis_plan.md, Sec. 4): NUTS vs NPE on a reduced model,
and leave-one-state-out predictions for LA, CA, NM, CO.

Reduced model (likelihood is tractable; drops selection, reporting error, quantization and
detection-limit uncertainty from simulator v0.2):
  wells with a pre-plugging measurement (DOI codes 1-9, NM field averages excluded);
  leaking: logit q = a_s + b_gas * gas, a_s ~ N(a0, tau_a);
  reading if leaking: log y ~ N(m_s + g_gas * gas - s_v^2 / 2, sqrt(sigma^2 + s_v^2)), m_s ~ N(m0, tau_m),
  s_v fixed at 0.59 (USGS Montana repeat visits, D24);
  recorded value: y if the leak is found (probability p_find for two-stage methods, 1 otherwise)
  and y > L_method, else 0.
Part A: 6 synthetic data sets from the prior, same wells and methods. Compare the NUTS posterior
  (exact) with an NPE trained on 20,000 simulations of this model and summary statistics of the same
  kind as the full simulator's. Pass if, for every hyperparameter, |mean difference| < 0.3 posterior
  SD and the SD ratio is within 0.7-1.4 on at least 5 of 6 data sets.
Part B: NUTS on the real records, all states, then leaving out LA, CA, NM, CO one at a time and
  predicting that state's detection share and positive-value quantiles from a new state drawn from the
  population. Reported: observed value, 90% predictive interval, inside or not.
Outputs: outputs/mcmc_reduced/{npe_vs_nuts.csv, real_posterior.csv, loso.csv}
Run from the repo root: python src/17_mcmc_reduced.py
"""
import sys, pathlib, time; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
from sim.model import build_design, METHODS

import numpy as np, pandas as pd
import jax, jax.numpy as jnp
import numpyro, numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS
from jax.scipy.stats import norm as jnorm
numpyro.set_host_device_count(1)
P = OUT / 'mcmc_reduced'; P.mkdir(exist_ok=True)
SV = 0.59
HYPER = ['a0', 'tau_a', 'b_gas', 'm0', 'tau_m', 'g_gas', 'log_sigma0', 'p_find']

d = build_design()
keep = d['doi_code'] > 0
ST_ALL = d['doi_state'][keep]; GAS = d['doi_gas'][keep]; CODE = d['doi_code'][keep]; Y_REAL = d['doi_value'][keep]
S = len(d['states'])
LIM = np.array([METHODS[k][1] or np.inf for k in range(len(METHODS))])
TWO = np.array([METHODS[k][2] for k in range(len(METHODS))])


def model(st, gas, code, y, n_states):
    a0 = numpyro.sample('a0', dist.Normal(0, 1.5)); tau_a = numpyro.sample('tau_a', dist.HalfNormal(1.0))
    b = numpyro.sample('b_gas', dist.Normal(0, 1)); m0 = numpyro.sample('m0', dist.Normal(0, 2))
    tau_m = numpyro.sample('tau_m', dist.HalfNormal(1.0)); g = numpyro.sample('g_gas', dist.Normal(0, 1))
    ls = numpyro.sample('log_sigma0', dist.Normal(np.log(2.0), 0.4)); pf = numpyro.sample('p_find', dist.Beta(4, 2))
    with numpyro.plate('s', n_states):
        za = numpyro.sample('za', dist.Normal(0, 1)); zm = numpyro.sample('zm', dist.Normal(0, 1))
    a_s, m_s = a0 + tau_a * za, m0 + tau_m * zm
    lq = jax.nn.log_sigmoid(a_s[st] + b * gas); l1q = jax.nn.log_sigmoid(-(a_s[st] + b * gas))
    mu = m_s[st] + g * gas - SV ** 2 / 2
    sd = jnp.sqrt(jnp.exp(2 * ls) + SV ** 2)
    pfi = jnp.where(TWO[code], pf, 1.0)
    pos = y > 0
    ly = jnp.log(jnp.where(pos, y, 1.0))
    lpos = lq + jnp.log(pfi) + jnorm.logpdf(ly, mu, sd) - ly
    pdet = pfi * jnorm.sf((jnp.log(LIM[code]) - mu) / sd)
    lzero = jnp.logaddexp(l1q, lq + jnp.log1p(-pdet))
    numpyro.factor('lik', jnp.where(pos, lpos, lzero).sum())


def nuts(st, gas, code, y, seed, n=1000):
    m = MCMC(NUTS(model, target_accept_prob=0.9), num_warmup=n, num_samples=n, num_chains=2,
             chain_method='sequential', progress_bar=False)
    m.run(jax.random.PRNGKey(seed), st=jnp.array(st), gas=jnp.array(gas), code=jnp.array(code), y=jnp.array(y), n_states=S)
    s = m.get_samples()
    div = int(m.get_extra_fields()['diverging'].sum())
    return {k: np.asarray(v) for k, v in s.items()}, div


def simulate(th, rng, st=ST_ALL, gas=GAS, code=CODE):
    a0, tau_a, b, m0, tau_m, g, ls, pf = th
    a_s, m_s = rng.normal(a0, tau_a, S), rng.normal(m0, tau_m, S)
    n = len(st)
    z = rng.random(n) < 1 / (1 + np.exp(-(a_s[st] + b * gas)))
    r = np.exp(rng.normal(m_s[st] + g * gas - SV ** 2 / 2, np.sqrt(np.exp(2 * ls) + SV ** 2), n))
    found = ~TWO[code] | (rng.random(n) < pf)
    return np.where(z & found & (r > LIM[code]), r, 0.0)


def summaries(y, st=ST_ALL, gas=GAS, code=CODE):
    out = []
    for grp, ng in ((st, S), (code, len(METHODS)), (gas.astype(int), 2)):
        for k in range(ng):
            m = grp == k
            if grp is st and k not in set(ST_ALL):
                continue
            if grp is code and k == 0:
                continue
            yy = y[m & (y > 0)]
            out += [(y[m] > 0).mean() if m.any() else 0.0]
            out += list(np.quantile(np.log(yy), (.1, .5, .9))) if len(yy) >= 3 else [-12.0, -12.0, -12.0]
    return np.array(out)


def prior_draw(rng, n):
    return np.column_stack([rng.normal(0, 1.5, n), abs(rng.normal(0, 1, n)), rng.normal(0, 1, n), rng.normal(0, 2, n),
                            abs(rng.normal(0, 1, n)), rng.normal(0, 1, n), rng.normal(np.log(2), .4, n), rng.beta(4, 2, n)])


rng = np.random.default_rng(17)
t0 = time.time()
# ---- Part A: NPE vs exact posterior on synthetic data
import torch
from sbi.inference import NPE
from sbi.neural_nets import posterior_nn
from scipy.special import logit, expit
NS = 20000
TH = prior_draw(rng, NS)
X = np.array([summaries(simulate(t, rng)) for t in TH])
U = TH.copy(); U[:, [1, 4]] = np.log(TH[:, [1, 4]]); U[:, 7] = logit(TH[:, 7])
mu, sd = X.mean(0), X.std(0) + 1e-9
zf = lambda x: torch.tensor(np.clip((x - mu) / sd, -10, 10), dtype=torch.float32)
um, us = U.mean(0), U.std(0)
prior = torch.distributions.Independent(torch.distributions.Normal(torch.tensor(um, dtype=torch.float32),
                                                                   torch.tensor(3 * us, dtype=torch.float32)), 1)
torch.manual_seed(17)
inf = NPE(prior=prior, density_estimator=posterior_nn('nsf', hidden_features=64, num_transforms=4), show_progress_bars=False)
est = inf.append_simulations(torch.tensor(U, dtype=torch.float32), zf(X)).train(max_num_epochs=300, training_batch_size=256)
post = inf.build_posterior(est)
print('reduced NPE trained on %d sims in %.0f s' % (NS, time.time() - t0))


def npe_draws(x, n=4000):
    u = post.sample((n,), x=zf(x), show_progress_bars=False).numpy().astype(float)
    u[:, [1, 4]] = np.exp(u[:, [1, 4]]); u[:, 7] = expit(u[:, 7])
    return u


rows = []
test = prior_draw(np.random.default_rng(170), 6)
for j, th in enumerate(test):
    y = simulate(th, rng)
    s, div = nuts(ST_ALL, GAS, CODE, y, seed=j)
    nd = npe_draws(summaries(y))
    for i, nme in enumerate(HYPER):
        a, b = s[nme], nd[:, i]
        rows.append(dict(dataset=j, parameter=nme, truth=th[i], nuts_mean=a.mean(), nuts_sd=a.std(), npe_mean=b.mean(),
                         npe_sd=b.std(), z_diff=(b.mean() - a.mean()) / a.std(), sd_ratio=b.std() / a.std(), divergences=div))
    print('  synthetic data set %d done (%.0f s)' % (j, time.time() - t0))
A = pd.DataFrame(rows)
A['ok'] = (A.z_diff.abs() < 0.3) & A.sd_ratio.between(0.7, 1.4)
A.to_csv(P / 'npe_vs_nuts.csv', index=False, float_format='%.4g')
okp = A.groupby('parameter').ok.sum()
print(A.groupby('parameter')[['z_diff', 'sd_ratio']].agg(['median', 'min', 'max']).round(2).to_string())
print('data sets passing per parameter (need >= 5 of 6):', okp.to_dict())

# ---- Part B: real records, all states and leave-one-state-out
s_all, div = nuts(ST_ALL, GAS, CODE, Y_REAL, seed=100)
R = pd.DataFrame([dict(parameter=k, median=np.median(s_all[k]), q05=np.quantile(s_all[k], .05), q95=np.quantile(s_all[k], .95))
                  for k in HYPER])
R['npe_reduced_median'] = np.median(npe_draws(summaries(Y_REAL)), 0)
R.to_csv(P / 'real_posterior.csv', index=False, float_format='%.4g')
print('real records, NUTS (%d divergences) vs reduced NPE:' % div); print(R.round(3).to_string(index=False))


def stats_of(y):
    p = y[y > 0]
    return [(y > 0).mean()] + (list(np.quantile(np.log(p), (.1, .5, .9))) if len(p) >= 3 else [np.nan] * 3)


L = []
for name in ('LA', 'CA', 'NM', 'CO'):
    k = d['states'].index(name)
    tr, te = ST_ALL != k, ST_ALL == k
    s, div = nuts(ST_ALL[tr], GAS[tr], CODE[tr], Y_REAL[tr], seed=200 + k)
    pr = []
    for i in rng.choice(len(s['a0']), 1000, replace=False):
        th = [s[h][i] for h in HYPER]
        a_new, m_new = rng.normal(th[0], th[1]), rng.normal(th[3], th[4])
        n = te.sum(); g_ = GAS[te]; c_ = CODE[te]
        z = rng.random(n) < expit(a_new + th[2] * g_)
        r = np.exp(rng.normal(m_new + th[5] * g_ - SV ** 2 / 2, np.sqrt(np.exp(2 * th[6]) + SV ** 2), n))
        pr.append(stats_of(np.where(z & (~TWO[c_] | (rng.random(n) < th[7])) & (r > LIM[c_]), r, 0.0)))
    pr = np.array(pr); ob = stats_of(Y_REAL[te])
    for j, stat in enumerate(['detect_share', 'log_q10', 'log_q50', 'log_q90']):
        lo, hi = np.nanquantile(pr[:, j], [.05, .95])
        L.append(dict(state=name, n=int(te.sum()), statistic=stat, observed=ob[j], pred_median=np.nanmedian(pr[:, j]),
                      pred_q05=lo, pred_q95=hi, inside_90=bool(lo <= ob[j] <= hi), divergences=div))
L = pd.DataFrame(L); L.to_csv(P / 'loso.csv', index=False, float_format='%.4g')
print(L.round(3).to_string(index=False))
print('LOSO: %d of %d observed statistics inside the 90%% predictive interval' % (L.inside_90.sum(), L.inside_90.notna().sum()))
print('total %.0f s' % (time.time() - t0))
