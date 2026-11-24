"""Ensemble of neural posterior estimators on simulator v0.2 (D30 follow-up).

Trains K flows with different seeds on the cached simulations, mixes their posteriors equally, and
runs SBC (20 rank bins, fail if chi-square p < 0.01, as in the analysis plan; 200 held-out
simulations instead of 300) plus expected-coverage checks. Each trained flow is saved so a
restarted environment resumes where it stopped. Real-data posterior draws are saved for the decision
analysis (src/16) but are not results until SBC and coverage pass (analysis_plan.md, Sec. 4).
Run from the repo root: python src/15_npe_ensemble.py [--k 5]
"""
import sys, pathlib, time, argparse; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from config import *
from sim.model import build_design, real_summaries, PARAMS, NAMES

import numpy as np, pandas as pd, torch
from scipy.special import expit, logit
from scipy import stats
from sbi.inference import NPE
from sbi.neural_nets import posterior_nn

ap = argparse.ArgumentParser()
ap.add_argument('--k', type=int, default=5); ap.add_argument('--n', type=int, default=100000)
ap.add_argument('--sbc', type=int, default=200); ap.add_argument('--hidden', type=int, default=100)
args = ap.parse_args()
P = OUT / 'npe_ensemble'; P.mkdir(exist_ok=True)
W = PROC / 'npe_ensemble'; W.mkdir(exist_ok=True)

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
z0 = np.load(PROC / ('npe_sims_v%d_n%d.npz' % (design['version'], args.n + args.sbc)))
theta, X = z0['theta'], z0['X']
x_obs = real_summaries(design)
keep = X[:args.n].std(0) > 1e-9
mu, sd = X[:args.n, keep].mean(0), X[:args.n, keep].std(0)
zf = lambda x: torch.tensor(np.clip((x[..., keep] - mu) / sd, -10, 10), dtype=torch.float32)
U = to_u(theta)
umu, usd = U[:args.n].mean(0), U[:args.n].std(0)
prior = torch.distributions.Independent(torch.distributions.Normal(torch.tensor(umu, dtype=torch.float32),
                                                                   torch.tensor(3 * usd, dtype=torch.float32)), 1)
posts = []
for k in range(args.k):
    f = W / f'flow_{k}.pt'
    torch.manual_seed(100 + k)
    inf = NPE(prior=prior, density_estimator=posterior_nn('nsf', hidden_features=args.hidden, num_transforms=5),
              show_progress_bars=False)
    inf.append_simulations(torch.tensor(U[:args.n], dtype=torch.float32), zf(X[:args.n]))
    if f.exists():
        est = inf._build_neural_net(torch.tensor(U[:200], dtype=torch.float32), zf(X[:200]))
        est.load_state_dict(torch.load(f)); est.eval()
        print('flow %d loaded' % k)
    else:
        t0 = time.time()
        est = inf.train(max_num_epochs=300, training_batch_size=256)
        torch.save(est.state_dict(), f)
        print('flow %d trained in %.0f s' % (k, time.time() - t0))
    posts.append(inf.build_posterior(est))


def sample(x, n):
    per = int(np.ceil(n / len(posts)))
    s = np.concatenate([p.sample((per,), x=zf(x), show_progress_bars=False).numpy() for p in posts])
    return s[:n]


ranks, cover = [], []
levels = np.linspace(0.05, 0.95, 19)
for i in range(args.n, args.n + args.sbc):
    s = sample(X[i], 500)
    ranks.append((s < U[i]).sum(0))
    lo = np.quantile(s, 0.5 - levels / 2, axis=0); hi = np.quantile(s, 0.5 + levels / 2, axis=0)
    cover.append((lo <= U[i]) & (U[i] <= hi))
ranks, cover = np.array(ranks), np.array(cover).mean(0)
res = []
for j, nme in enumerate(NAMES):
    h = np.histogram(ranks[:, j], bins=20, range=(0, 501))[0]   # plan Sec. 4: 20 bins, fail if p < 0.01
    res.append({'parameter': nme, 'chi2_p': stats.chisquare(h).pvalue, 'mean_rank': ranks[:, j].mean() / 500,
                'coverage_50': cover[np.argmin(abs(levels - .5)), j], 'coverage_90': cover[np.argmin(abs(levels - .9)), j]})
res = pd.DataFrame(res)
res.to_csv(P / 'sbc_coverage.csv', index=False, float_format='%.3f')
print(res.round(3).to_string(index=False))
np.save(W / 'sbc_ranks.npy', ranks)
print('pass (p >= 0.01, 20 bins): %d of %d | p > 0.05: %d | 90%% coverage range %.2f-%.2f'
      % ((res.chi2_p >= 0.01).sum(), len(res), (res.chi2_p > 0.05).sum(), res.coverage_90.min(), res.coverage_90.max()))
draws = from_u(sample(x_obs, 4000).astype(float))
np.savez_compressed(W / 'real_posterior_draws.npz', theta=draws, names=np.array(NAMES))
pd.DataFrame({'parameter': NAMES, 'median': np.median(draws, 0), 'q05': np.quantile(draws, .05, 0),
              'q95': np.quantile(draws, .95, 0)}).to_csv(P / 'real_data_posterior_TRIAL.csv', index=False, float_format='%.4g')

# ---- posterior predictive check on the real summaries (same rule as the prior check in src/12)
from sim.model import simulate, summary_names
rng = np.random.default_rng(15)
sims = np.array([simulate(t, design, rng) for t in draws[rng.choice(len(draws), 500, replace=False)]])
lo, hi = np.quantile(sims, [.005, .995], axis=0)
out = (x_obs < lo) | (x_obs > hi)
ppc = pd.DataFrame({'summary': summary_names(design), 'observed': x_obs, 'q005': lo, 'q50': np.median(sims, 0), 'q995': hi,
                    'outside_99': out})
ppc.to_csv(P / 'posterior_predictive_TRIAL.csv', index=False, float_format='%.4g')
print('posterior predictive: %d of %d real summaries outside the 99%% range' % (out.sum(), len(out)))
print(ppc[out].round(3).to_string(index=False))
