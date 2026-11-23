# Analysis plan (pre-registration draft)

Version 0.1, 2026-11-23. **Status: draft.** It will be frozen (git tag `plan-v1.0`) before the full
model is fitted to the real data. Changes after freezing go in the deviations log at the end, with
date and reason.

What has already been seen: descriptive statistics of the program data, the attribute-prediction
analyses (src/03–05), the detection-by-method tables (src/08), the policy cross-check (src/09), and a
pooled pilot model (src/11, D25). The design choices below were made with that knowledge. What has
**not** been run: the selection model, the hierarchical emission model, simulation-based inference
on the real data, and the decision analysis. Their outputs are the results this plan governs.

## 1. Questions and estimands

| ID | Question | Estimand (posterior median and 90% interval) |
|---|---|---|
| Q1 | After correcting for selection and detection limits, what is the distribution of pre-plugging methane emissions from documented orphaned wells? | E1: share leaking, median leak rate, P(E > 1, 10, 100 g/h), by state and pooled, for (a) plugged wells and (b) unplugged 2022-inventory wells. E2: total pre-plugging emission rate of plugged wells, and its ratio to the naive estimates (non-detects as zero; quantified wells only; the DOI method of applying the measured-well average change to unmeasured wells). |
| Q2 | How much of the variation across states reflects measurement method rather than the wells? | E3: share of between-state variance in observed detection rates explained by method and measurement regime, vs. by the emission model's state effects. |
| Q4 | For a fixed budget, how much more methane is avoided by screening before plugging, and at what cost does screening stop paying? | E4: avoided emissions (t CH4 and present-value $) of each policy in Section 5; E5: break-even screening cost c*(L, p_find); E6: expected value lost when pre-plugging measurement became optional (July 2025). |
| Q3 (secondary) | Do public attributes identify high emitters within states? | Reported from the existing analysis (src/04–05) plus the attribute coefficients of the emission model; no new primary estimand. |

## 2. Data and inclusion rules

Snapshot of 2026-11-22 (D1). Unit = one well. Rules already logged in `docs/decision_log.md`:

- Selection pool: USGS 2022 documented unplugged orphaned wells (D3); plugged = linked to a DOI
  record by API or within 100 m (D4–D6).
- Pre-plugging value classes (from `data/processed/method_group.csv`):

| Class | Rule | Use in the likelihood |
|---|---|---|
| Quantitative measurement | chamber, flow meter + GC, Hi-Flow, quantitative laser, federal guidelines | Value observed; zero = below that method's limit |
| Screening | OGI camera, handheld detector, standoff laser, FID | Detected / not detected at the screening limit; positive values used as detections only |
| Placeholder | no method recorded (D7, D20) | Missing (not measured) |
| Ambiguous | "ND" strings (WV) | Censored at the Hi-Flow limit (sensitivity: missing) |
| Field average | values repeated ≥3 times within NM (D21) | Excluded; the well counts as plugged |
| Quantized / rounded | LA Hi-Flow steps of ~1.15 g/h (D22); FL "1 g/h" (= "<1"); recurring 0.004, 0.01 (D10) | Interval-censored |

- Measurement regime r = state × period, periods split at 2023-07-07 (formula grants) and
  2025-07-11 (measurement made optional) (D19). Grant type is not in the data; where state reports
  give it, it overrides the period split.
- External data: Driscoll et al. chamber data (132 LA wells; 5 paired with DOI records, D17); USGS
  Kevin-Sunburst repeat visits (15 MT wells, D24); Kentucky Geological Survey data if shared.

## 3. Model

Notation: well i in state s, attributes x_i (type, depth, age, years since last production,
coordinates), true pre-plugging emission rate E_i (g/h).

**Emission model (hurdle, partially pooled).**
- Leaking: Z_i ~ Bernoulli(q_i), logit q_i = a_s + x_i'b, a_s ~ Normal(a_0, tau_a).
- Rate if leaking: log E_i ~ Normal(m_s + x_i'g, sigma_s), m_s ~ Normal(m_0, tau_m),
  log sigma_s ~ Normal(log sigma_0, tau_sigma). E_i = 0 if Z_i = 0.
- Tail: the primary family is log-normal; the alternative is log-Student-t (nu ~ Gamma(2, 0.1)).
  Choose by held-out posterior predictive coverage of the top 1% (Section 6), not by fit to the
  observed maximum. A hard cap is not used in the primary model (the pilot's 50 kg/h cap was a
  pilot device).

**Selection model.** Plugged by mid-2026: logit P(S_i = 1) = d_s + x_i'h + rho_s * Z_i + u_c(i),
with u_c a spatial-cluster effect (0.5 degree cells). rho_s is the dependence of selection on leaking.
Prior: rho_s ~ Normal(0, 1) in states whose documented rules score leaks or use measurements
(NM, LA, NY, KY, TX, CO, OH, KS, MT, AK, IN, FL, MS; `data/reference/state_policies.csv`), rho_s = 0
elsewhere; sensitivity with rho_s ~ Normal(0, 1) everywhere. Pennsylvania gets its own d and h
(plugger-chosen wells). rho_s is identified mainly through states that measured before selecting
and through the unplugged Louisiana chamber sample; we report how far its posterior moves from the
prior.

**Measurement model** (plugged wells, regime r, method m):
1. Measured at all: M_i ~ Bernoulli(pi_r), independent of E_i given regime (missing at random within
   regime). Sensitivity: allow dependence on Z_i in regimes where screening preceded quantification.
2. Leak found (two-stage methods: Hi-Flow, OGI-located): F_i ~ Bernoulli(p_find,m) if Z_i = 1;
   p_find ~ Beta(4, 2) (mean 0.67). Chamber and flow-meter methods enclose the wellhead: p_find = 1.
3. Single-visit reading: Y_i = E_i * eps_i, log eps_i ~ Normal(-s_v^2/2, s_v), s_v prior centered on the
   Montana between-visit CV of 0.65 (s_v ~ Normal+(0.59, 0.2)).
4. Detection: reported positive if F_i = 1 and Y_i > L_m, with log L_m ~ Normal(log of the midpoint, 0.5)
   from `data/reference/detection_limits.csv`; otherwise reported as zero / non-detect.
5. Rounding and quantization as interval censoring (Section 2).
6. Reporting error between the contractor's reading and the program record (D28): with probability
   omega a measured leak is recorded as zero, omega ~ Beta(1, 30); otherwise the recorded value is the
   reading times exp(r), r ~ Normal(0, s_r), s_r ~ Normal+(0, 0.3). Priors checked against the 202 New
   Mexico memo-record pairs (4 recorded as zero, most within 5%).

External calibration enters as additional observations of the same model: the 5 Driscoll pairs
(chamber reading with L = 0.001 g/h and p_find = 1, plus the DOI reading at the later visit), the
Driscoll unplugged wells as a sample from the LA unplugged population, and the Montana repeat
visits for s_v only (those wells were chosen because they leak).

## 4. Inference

- Neural posterior estimation (normalizing flow) trained on simulations from the priors, as in the
  landfill project. Summary statistics, per state x regime x method cell: counts measured, detected,
  share detected, 10/50/90% quantiles of log positive values, share plugged by attribute bins and
  spatial-cell plugging rates; plus the paired and external data in full.
- Validation before use on real data:
  - simulation-based calibration, 300 simulations, rank histograms for all hyperparameters and for
    E1/E2; failure = chi-square p < 0.01 on 20 bins for any target;
  - MCMC (NumPyro NUTS) on a reduced model (no selection, no spatial effect), compared with NPE on the
    same reduced model: posterior medians within 0.25 posterior SD;
  - leave-one-state-out checks for LA, CA, NM, CO: predicted vs observed detection share and
    positive-value quantiles; coverage of 90% intervals reported.
- If SBC fails, simplify in this order: drop the spatial effect, fix sigma_s = sigma_0, drop rho_s outside
  NM and LA. Each simplification is logged as a deviation.

## 5. Decision analysis

- Candidates: posterior predictive draws of unplugged 2022-inventory wells, by state.
- Budget: plug a fraction f of candidates (f = 0.1 primary, 0.3 secondary) per state.
- Policies: P0 the observed selection model; P1 random; P2 rank by posterior mean emission given
  attributes; P3 screen every candidate once with method (L, p_find, cost c per well), then plug in
  order of the reading, then by attributes.
- Value: avoided emission rate x present-value factor x social cost of methane. Primary: $1,900/t
  (EPA 2023, 2% discount), emissions declining 6.4%/yr, 25-year horizon, 2% discount (as in RFF 2026).
  Sensitivity: SC-CH4 $1,000 and $3,000; decline 1%; horizon 10 and 50 years.
- E5: c*(L, p_find) = (value of P3 - value of P2) / number screened, over L in 0.001–10,000 g/h and
  p_find in {0.25, 0.5, 0.75, 1}. E6: value of P3 with the formula-grant-era method mix minus value of
  P0, per plugged well, times the wells plugged after 2025-07-11.
- Plugging costs: state medians of the reported cost (several states report contract averages);
  used only for cost-per-tonne reporting, not for the ranking.

## 6. Sensitivity analyses (all pre-specified)

1. Detection-limit priors at half and twice the midpoint; wider (0.5 → 1.0 log-SD).
2. Leave-top-k-out (k = 1, 2, 5), including the two ~46 kg/h Colorado wells.
3. Spatial link radius 50 m and 250 m instead of 100 m (D6).
4. Exclude records plugged after 2025-07-11.
5. Tail family: log-normal vs log-Student-t.
6. Selection: rho_s = 0 everywhere; rho_s free everywhere.
7. p_find prior Beta(2, 2) and Beta(8, 2).
8. "ND" strings as missing instead of censored.
9. Handheld-screening zeros as missing instead of non-detects.

## 7. Reporting rules

- Every estimand as posterior median with 90% interval; no p-values.
- The naive-vs-corrected comparison is reported whatever its direction.
- If leave-one-state-out coverage is below 70% for a state, that state's estimates are reported as
  unreliable and the national estimate is shown with and without it.
- Code, frozen data and simulation seeds are released with the paper.

## Deviations log

(none yet; draft edited on 2026-11-23 before freezing to add the reporting-error term, D28)
