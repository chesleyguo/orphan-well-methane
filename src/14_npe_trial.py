"""TRIAL RUN of neural posterior estimation on the simulator skeleton (src/sim/model.py, D26).

Checks that the inference pipeline works end to end: simulate from the prior, train a neural
posterior estimator, run simulation-based calibration (SBC) on held-out simulations, then apply it to
the real summaries. Not a result: the simulator is still a skeleton (analysis_plan.md governs the
real fit). Outputs in outputs/npe_trial/.
Run from the repo root: python src/14_npe_trial.py [--n 20000]
"""
import sys, pathlib, time, argparse; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
from sim.model import build_design, sample_prior, simulate, real_summaries, PARAMS, NAMES

import numpy as np, pandas as pd, torch
from scipy.special import expit, logit
from scipy import stats
from sbi.inference import NPE
from sbi.neural_nets import posterior_nn

ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, default=100000); ap.add_argument('--sbc', type=int, default=200); ap.add_argument('--hidden', type=int, default=100)
args = ap.parse_args()
torch.manual_seed(0); rng = np.random.default_rng(1)
P = OUT / 'npe_trial'; P.mkdir(exist_ok=True)

# unconstrained transform of the parameters
KIND = [p[1] for p in PARAMS]
def to_u(th):
    u = th.copy()
    for j, k in enumerate(KIND):
        if k.startswith('halfnormal'): u[:, j] = np.log(np.maximum(th[:, j], 1e-6))
        if k == 'beta': u[:, j] = logit(np.clip(th[:, j], 1e-6, 1 - 1e-6))
    return u
def from_u(u):
    th = u.copy()
    for j, k in enumerate(KIND):
        if k.startswith('halfnormal'): th[:, j] = np.exp(u[:, j])
        if k == 'beta': th[:, j] = expit(u[:, j])
    return th

design = build_design()
cache = PROC / ('npe_sims_v%d_n%d.npz' % (design['version'], args.n + args.sbc))
if cache.exists():
    z0 = np.load(cache); theta, X = z0['theta'], z0['X']
    print('loaded %d cached simulations' % len(X))
else:
    t0 = time.time()
    theta = sample_prior(rng, args.n + args.sbc)
    X = np.array([simulate(t, design, rng) for t in theta])
    np.savez_compressed(cache, theta=theta, X=X)
    print('simulated %d datasets' % len(X))
x_obs = real_summaries(design)

keep = X[:args.n].std(0) > 1e-9
mu, sd = X[:args.n, keep].mean(0), X[:args.n, keep].std(0)
z = lambda x: np.clip((x[..., keep] - mu) / sd, -10, 10)
U = to_u(theta)
umu, usd = U[:args.n].mean(0), U[:args.n].std(0)

prior = torch.distributions.Independent(torch.distributions.Normal(torch.tensor(umu, dtype=torch.float32),
                                                                   torch.tensor(3 * usd, dtype=torch.float32)), 1)
inf = NPE(prior=prior, density_estimator=posterior_nn('nsf', hidden_features=args.hidden, num_transforms=5),
          show_progress_bars=False)
t0 = time.time()
inf.append_simulations(torch.tensor(U[:args.n], dtype=torch.float32), torch.tensor(z(X[:args.n]), dtype=torch.float32))
est = inf.train(max_num_epochs=300, training_batch_size=256)
post = inf.build_posterior(est)
print('trained in %.0f s, %d epochs' % (time.time() - t0, len(inf.summary['training_loss'])))

# ---- SBC on held-out simulations
ranks = []
for i in range(args.n, args.n + args.sbc):
    s = post.sample((500,), x=torch.tensor(z(X[i]), dtype=torch.float32), show_progress_bars=False).numpy()
    ranks.append((s < U[i]).sum(0))
ranks = np.array(ranks)
sbc = []
for j, nme in enumerate(NAMES):
    h = np.histogram(ranks[:, j], bins=10, range=(0, 500))[0]
    sbc.append({'parameter': nme, 'chi2_p': stats.chisquare(h).pvalue,
                'mean_rank': ranks[:, j].mean() / 500,
                'posterior_sd_over_prior_sd': None})
sbc = pd.DataFrame(sbc)

# ---- contraction (how much the data narrow each parameter) on held-out sims, and the real-data fit
contr = []
for i in range(args.n, args.n + min(args.sbc, 50)):
    s = post.sample((500,), x=torch.tensor(z(X[i]), dtype=torch.float32), show_progress_bars=False).numpy()
    contr.append(s.std(0) / usd)
sbc['posterior_sd_over_prior_sd'] = np.mean(contr, 0)
sbc.to_csv(P / 'sbc.csv', index=False, float_format='%.3f')
print(sbc.round(3).to_string(index=False))

s = post.sample((4000,), x=torch.tensor(z(x_obs), dtype=torch.float32), show_progress_bars=False).numpy()
th = from_u(s.astype(float))
res = pd.DataFrame({'parameter': NAMES, 'median': np.median(th, 0), 'q05': np.quantile(th, .05, 0), 'q95': np.quantile(th, .95, 0)})
# derived: average-state share leaking (non-gas), median leak rate, P(E > 10 g/h)
p = dict(zip(NAMES, th.T))
q = expit(p['a0']); sig = np.exp(p['log_sigma0'])
der = {'share_leaking_avg_state': q, 'median_leak_rate_g_h': np.exp(p['m0']),
       'P_gt_10_g_h_avg_state': q * stats.norm.sf((np.log(10) - p['m0']) / sig)}
for k, v in der.items():
    res.loc[len(res)] = [k, np.median(v), np.quantile(v, .05), np.quantile(v, .95)]
res.to_csv(P / 'real_data_posterior.csv', index=False, float_format='%.4g')
print('\nreal data (TRIAL, skeleton simulator):')
print(res.round(3).to_string(index=False))
n_out = (np.abs((x_obs[keep] - mu) / sd) > 10).sum()
print('real summaries clipped at |z|>10: %d of %d' % (n_out, keep.sum()))
