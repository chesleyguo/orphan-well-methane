"""Simulator skeleton check: timing, summary dimensions, and a prior predictive check (where the
real summaries fall within the distribution of summaries simulated from the prior).

A real summary far outside the prior predictive range (quantile < 0.005 or > 0.995) means the
simulator or the prior cannot produce data like ours yet; those are listed for revision.
Run from the repo root: python src/12_sim_prior_check.py
"""
import sys, pathlib, time; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
from sim.model import build_design, sample_prior, simulate, real_summaries, summary_names, NAMES

import numpy as np, pandas as pd

N = 1000
rng = np.random.default_rng(0)
design = build_design(cache=False)
names = summary_names(design)
s_obs = real_summaries(design)
print('states: %d | inventory wells: %d | DOI wells: %d | summaries: %d | parameters: %d'
      % (len(design['states']), len(design['inv_state']), len(design['doi_state']), len(s_obs), len(NAMES)))

theta = sample_prior(rng, N)
t0 = time.time()
sims = np.array([simulate(th, design, rng) for th in theta])
dt = (time.time() - t0) / N
print('seconds per simulation: %.3f (10,000 sims ~ %.0f min)' % (dt, dt * 1e4 / 60))

q = (sims < s_obs).mean(0) + 0.5 * (sims == s_obs).mean(0)
pp = pd.DataFrame({'summary': names, 'observed': s_obs, 'prior_q05': np.quantile(sims, .05, 0),
                   'prior_q50': np.quantile(sims, .5, 0), 'prior_q95': np.quantile(sims, .95, 0),
                   'quantile_of_observed': q})
P = OUT / 'sim'; P.mkdir(exist_ok=True)
pp.to_csv(P / 'prior_predictive.csv', index=False, float_format='%.4g')
bad = pp[(pp.quantile_of_observed < 0.005) | (pp.quantile_of_observed > 0.995)]
print('summaries outside the 99%% prior predictive range: %d of %d' % (len(bad), len(pp)))
print(bad.round(3).to_string(index=False))

# which parameters do the summaries carry information about? (max |rank correlation| over summaries)
from scipy.stats import spearmanr
ok = sims.std(0) > 0
rc = spearmanr(np.column_stack([theta, sims[:, ok]])).correlation[:len(NAMES), len(NAMES):]
ok_names = np.array(names)[ok]
info = pd.DataFrame({'parameter': NAMES, 'max_abs_rank_corr': np.abs(rc).max(1),
                     'most_informative_summary': ok_names[np.abs(rc).argmax(1)]})
info.to_csv(P / 'parameter_information.csv', index=False, float_format='%.3f')
print('\nparameter information (max |Spearman| between parameter and any summary):')
print(info.round(2).to_string(index=False))
