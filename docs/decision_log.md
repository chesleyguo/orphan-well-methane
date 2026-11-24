# Decision log

Append-only. Each entry: date, decision, reason, whether it was made before or after seeing results.
Entries D1–D12 were made during the 2026-11-22 feasibility exploration, i.e. **after** looking at
descriptive results. They are working choices, not pre-registered rules; review them before freezing
the analysis plan.

## 2026-11-22 — data snapshot

**D1. Frozen snapshot.** All analysis uses the raw files in `data/raw/` downloaded on 2026-11-22
(checksums in `MANIFEST.sha256`). The DOI layer is updated over time; a new pull goes to
`data/raw_<date>/` via `src/00_download.py` and is treated as a new snapshot, never mixed in silently.

**D2. Primary methane source.** Use the DOI Orphaned Wells Program ArcGIS feature service
(`DOI_BIL_Orphaned_Wells_Web_Map_Joined`, 13,108 records), not the DOI xlsx on doi.gov, because the
xlsx has only 8 columns and no methane, method or cost fields.

**D3. Pre-program selection pool.** Use the USGS 2022 DOW dataset (117,672 wells) as the pool from
which plugged wells were selected. The 2026 USGS database already drops most plugged wells
(only 1,011 DOI wells match it), so it cannot serve as the pool.

## 2026-11-22 — cleaning and linkage

**D4. API key.** Link on the first 10 digits of the API number (state + county + well). Records with
fewer than 10 digits get no API key.

**D5. Placeholder API.** `0000000000` is treated as missing (used by some LA records in both the DOI
layer and SONRIS). Found after an initial run had wrongly attached SONRIS attributes and Driscoll
pairs through this placeholder.

**D6. Spatial fallback.** Wells without an API match are linked to the 2022 inventory by nearest
coordinates within the same state, within 100 m, one-to-one, closest first. Result: 4,855 API +
1,515 spatial = 6,370 linked DOI wells (49%). The 100 m threshold is a judgment call; test 50 m and
250 m before freezing.

**D7. Missing-method zeros are not measurements.** A pre-plugging value is treated as a real
measurement only if `methane_detection_method` is filled and not one of N/A, NA, NR, NP, ND (blank).
6,080 of 8,499 records with a value have no method; 0.2% of them are non-zero (vs. 27% for Hi-Flow).
Kansas (2,568 zeros, all "N/A") is the largest case.

**D8. Method classes.** Free-text methods are grouped by keyword into: hi-flow, handheld detector,
federal-guidelines, flow+GC, chamber/flow-through, laser, OGI camera, calculated, other, none
(see `method_class()` in `src/01_build_linkage.py`). "Federal Methane Measurement Guidelines" is kept
as its own class because it names a protocol, not an instrument.

**D9. Handheld detectors excluded from emission models.** 585 wells measured with handheld detectors
(Sensit, Bascom Turner, Altair, RKI) have 0 detections; these instruments read concentration, not
flow. Excluded from prediction samples; to be modelled as screening with an unknown flow threshold.

**D10. Rounded values.** Values such as 0.01, 1.0 and 0.004 g/h recur. Kept as reported for now;
plan to treat them as interval-censored in the measurement model.

**D11. WV dates.** The WVDEP layer's `compdate`/`issuedate` are dates of the current (plugging)
permit, not spud dates, so WV wells get no age. An earlier run had used them as spud dates.

**D12. Age and date cleaning.** Spud/last-production dates outside 1850–2026 are set to missing;
`age_at_plug` outside 0–170 years is set to missing.

**D13. Reproducibility.** Gradient-boosting models use `random_state=0`; repeated-CV grouping uses
CRC32 instead of Python's salted `hash()`. Two full runs now give byte-identical logs
(`outputs/run_log.txt`). Before this fix, AUCs varied by about ±0.02 between runs.

## 2026-11-22 — detection limits

**D14. Refined instrument grouping (proposed, not yet applied to D1–D13 results).** Auditing the 76
method strings (`outputs/method_string_audit.csv`) showed misclassifications in `method_class`:
"Semtech Hi-Flow Meter" fell into flow+GC; RMLD/TDLAS standoff lasers (no detections) fell into the
quantitative laser class; several handheld detectors fell into "other". `src/08` adds `method_group`,
which separates quantitative rate methods from screening devices. Adopt it when freezing the plan; the
modeling sample and its numbers will then change slightly.

**D15. "ND" is ambiguous.** 70 TX/WV records list "ND" as the method; it may mean "not detected"
(a measured non-detect) rather than "no data". Kept separate; resolve with the states or DOI.

**D16. Detection modeled as a probability, not a floor.** Groups report values below their nominal
limits yet miss large leaks (Hi-Flow 0 vs chamber 40–48 g/h). Working priors for a detection
threshold by group are in `docs/detection_limits.md` (judgment, to revisit). Screening-device zeros
are not 0 g/h.

## 2026-11-23 — chamber pairing and Hi-Flow model

**D17. Driscoll-DOI pairing uses API, then serial number, then coordinates within 50 m.** The four
API-linked wells sit 0-36 m from their DOI record. Wells with the same lease name sit about 60-110 m
apart (for example SEYBERT: 66 m from a DOI record with the same lease name but a different API), so
a wider radius would pair sister wells. Under this rule the 22 Driscoll wells with a placeholder API
remain unmatched (nearest DOI record 66 m, median 534 m); a fifth pair came from DOI using the SONRIS
serial number as us_well_id. Result: 5 pairs (`outputs/driscoll_doi_pairs.csv`); nearest-record
distances for all 132 wells in `outputs/driscoll_nearest_doi.csv`.

**D18. Hi-Flow treated as Semtech Hi-Flow 2 throughout; no Bacharach transition-failure term.** No
DOI record names Bacharach; 350 name Semtech, and California's "Hi-Flow 2/II" is the Semtech product
name. Records that do not name a model (LA "Hiflow/Method 21", "High Flow Sampling") are assumed to be
the same instrument; a sensitivity run can widen the upper-tail error for them.

## 2026-11-23 — state and federal policy review (docs/state_policies.md)

**D19. Measurement regime depends on grant type and date.** Initial grants carried no methane
requirement; formula grants from July 2023 required screening and, on detection, before/after rate
measurement; from July 2025 measurement is optional and estimates are allowed. The DOI layer has no
grant-type field, so the measurement model uses state x period (before/after 2025-07-11) as the
regime indicator. Share of records with any value: 0.69 before, 0.49 after.

**D20. Zeros are not all measurements.** Values of 0 with no recorded method are treated as missing
(not measured). Zeros recorded with a handheld concentration detector (IN, PA, AZ) are non-detects
with the screening limit of that device. Zeros recorded with OGI are non-detects at 20-150 g/h.
About 4,100 records in KS, IL, NE, AR, IN, PA, AZ and MS have only zeros.

**D21. Field-averaged values are excluded from the likelihood.** New Mexico assigned field averages
to unmeasured wells (O'Brien field). Positive values repeated three or more times within a state
are flagged (`outputs/state_measurement_crosscheck.csv`); the NM repeats are dropped from the
measurement likelihood and kept only as evidence that the well was plugged.

**D22. Louisiana Hi-Flow values are interval-censored at ~1.15 g/h steps** (0.001 cfm resolution);
Florida values of exactly 1 g/h are treated as "<1 g/h".

**D23. Selection may depend on emissions.** NM used pre-plugging methane to choose wells; LA and NY
score observed gas leaks; PA selection is by pluggers (shallow, capped-cost wells). The selection
model allows dependence on emissions through leak indicators and models PA separately.

**D24. USGS Kevin-Sunburst (MT) wells used as an unplugged, measured comparison group and as a
snapshot-noise prior, not as paired data.** All 15 wells with emission data are in the 2022 orphan
inventory; none appears in the DOI plugged records (nearest DOI record 904 m). The wells were
chosen because they leak, so they are not a random sample of unplugged wells. Between-session
variation (median CV 0.65 for wells measured three or more times) informs the single-visit
measurement-noise prior. The data release swaps the latitude and longitude column names.

**D25. Pilot only (src/11, outputs/pilot/).** A pooled hurdle model (share not leaking + lognormal
leak rates, capped at 50 kg/h) with one censoring limit per method, and a screening simulation
(one noisy visit, CV 0.65; plug the top-ranked 10% or 30% vs a random choice; social cost of methane
$1,900/t, 25 years, 6.4%/yr decline, 2% discount). It ignores selection, state differences and the
Hi-Flow leak-finding stage. Its numbers are for discussion, not results; the planned analysis
replaces it. Two pilot findings shape the plan: (i) the value of screening scales roughly with the
probability that a leak is found (break-even $30k at 100%, $16k at 50%, $7.7k at 25%), so the
leak-finding stage must be modeled explicitly; (ii) fitted separately, states differ widely
(P > 10 g/h: LA, CA ~6%; CO 52%; OK 38%), and in single-method states the non-leaking share is not
identified (CA hits the boundary), so the full model must be hierarchical across states and use
external calibration (Driscoll chamber pairs, KGS, USGS MT).

**D26. Simulator skeleton (src/sim/model.py, v0.1).** Implements the analysis plan's generative
model with 13 global parameters; state effects are nuisance latents; attributes reduced to a gas-well
indicator and 0.5-degree cells; methods drawn from each state's empirical mix (period split not yet
used). 209 summaries computed by the same function for real and simulated data; 0.01 s per
simulation. Prior predictive check: 5 of 209 real summaries fall outside the 99% range, all in
states with very few measured wells (NY, OH, TN, MT). Marginal information check: the between-visit
noise s_v and the leak-finding probability p_find are barely informed by the program summaries
(|rank correlation| 0.07 and 0.25), as expected. Same day: the simulator now also generates the
Driscoll chamber sample of unplugged LA wells (from LA's emission model, chamber limit, no
leak-finding stage) and the Montana repeat visits (within-well SD of log readings), plus a contrast
summary (LA Hi-Flow detection share / chamber share above 1 g/h). Information rises to 0.90 for s_v
and 0.37 for p_find; 6 of 215 real summaries fall outside the 99% prior predictive range, all from
states with very few measured wells (NY, OH, NM screening, AL, MT). The 5 Driscoll-DOI pairs are not
yet simulated.

**D27. NPE trial run (src/14, outputs/npe_trial/; not part of the deterministic pipeline).** 20,000
simulations from the skeleton simulator, neural spline flow posterior (sbi 0.27, PyTorch 2.15 CPU),
SBC on 200 held-out simulations. 10 of 13 parameters pass (chi-square p > 0.05 on 10 rank bins);
log_sigma0 fails (p < 0.001, mean rank 0.43: slight upward bias), m0 (p = 0.011) and b_gas (p = 0.041)
are borderline. Posterior SD / prior SD: 0.27-0.34 for d0, m0, log_sigma0; 0.94-0.97 for tau_a and
tau_m (between-state spread barely learned). Real-data values are a pipeline check only; the
"average state" quantities are not comparable with the pilot's pooled-sample quantities.
Rerun with 50,000 simulations (default in outputs/npe_trial/ now): all 13 parameters pass SBC
(smallest p = 0.076, rho; mean ranks 0.46-0.54), so the log_sigma0 failure was a simulation-budget
issue, not a model defect. The between-state spreads tau_a, tau_m remain barely learned (posterior/prior
SD 0.89-0.98); the full model needs state-level external data or tighter priors for them.
The analysis plan's validation (MCMC comparison, leave-one-state-out) is still to come.

**D28. New Mexico contractor memos (src/13, outputs/nm_memo_audit.csv).** Downloaded the OCD well
files (2,256 PDFs filed since June 2022, text kept, 5 MB gzip) for all 217 NM program wells. 203 have
a contractor methane memo, mostly Well Done Foundation with a Ventbuster ultra-low flow meter
(111 wells name it). Findings:
- Field averaging confirmed: 11 O'Brien Fee 18/19 wells carry a field average (memo 1.67 g/h,
  DOI 1.69 g/h), as inferred in D21; 3 more memos mention averaging but report zero.
- 170 of 202 memo values agree with the DOI value within 5%. Of the rest: 4 wells are 0 in DOI but
  11-267 g/h in the memo (Clovis K Kendrick #001, Pitts Kendrick #001, Cortland Myers #002, #003); for
  several large emitters DOI is 20-55% below the memo; the others are small values or parser limits
  (128 memos are forms whose pre/post type the parser cannot tell).
Decision: the measurement model gets a reporting-error term (a small probability that a measured
leak is reported as zero, and a multiplicative error on reported values), with priors from these
202 pairs; the analysis plan will be updated before freezing.

**D29. Simulator v0.2 (src/sim/model.py).** Adds: method mix by state x period (periods split at
2023-07-07 and 2025-07-11; 70 undated wells assigned to the middle period); the reporting-error
mixture (omega, kappa, s_r); LA Hi-Flow quantization to 1.1504 g/h steps; the 4 usable Driscoll-DOI
pairs; and the 188 non-averaged New Mexico memo-record pairs. 16 parameters, 226 summaries,
0.013 s per simulation. Prior predictive: 6 of 226 real summaries outside the 99% range, all from
states with very few measured wells (NY, OH, AL, MT). Marginal information: omega 0.57, kappa 0.70,
s_r 0.34; p_find rises to 0.42.

**D30. NPE on simulator v0.2: SBC not yet passing for all parameters.** 50,000 simulations: 12 of 16
pass; a0 and m0 fail (p < 0.001) with opposite mean-rank shifts (0.42, 0.59). 100,000 simulations
(50,000 reused + 50,000 new, seed 2) and a wider flow (100 hidden units): a0 and m0 now pass, but d0
(p = 0.001), p_find (0.004), rho (0.013) and g_gas (0.046) fail, with mean ranks 0.47-0.54, i.e.
mis-calibrated width rather than bias, and the failing set changes between runs. Read as estimator
variance of a single trained flow. Next, before any real-data use: an ensemble of 5 flows (different
seeds), expected-coverage curves alongside SBC histograms, and, if needed, a learned embedding of the
summaries. Real-data values from these runs are not used.

**D31. Decision-analysis module (src/16).** Implements analysis plan Sec. 5 on unplugged 2022-inventory
wells (106,370 wells, 25 states). For each posterior draw, leak status is drawn conditional on the well
not having been plugged (selection with rho), so leak-rule states' remaining wells leak slightly less.
Each state plugs the top ceil(f x n_s) wells by policy score (f = 0.1, 0.3). P0 samples by the selection
weights without replacement (Efraimidis-Spirakis keys); P2 ranks by the posterior mean emission by
state x well type, averaged over draws (not the per-draw truth, which would be an oracle); P3 ranks
detected wells by their one-visit reading (noise s_v) and fills the rest by P2; an oracle ranking by the
true rate is reported as an upper bound. E5 grid: L = 10^-3 ... 10^4 g/h (15 points) x p_find
{0.25, 0.5, 0.75, 1}. E6 draws each candidate's method from its state's 2023-07 to 2025-07 method mix
(with leak finding and detection-limit uncertainty as in the simulator) and compares with P0, per
plugged well; multiplied by the 2,259 DOI wells with witness dates on or after 2025-07-11. Value per
g/h avoided: $179 (25 yr, 6.4%/yr decline, 2% discount, $1,900/t). Tested on prior draws only; real
runs wait for the NPE ensemble (D30) and are labelled TRIAL until validation passes.

**D32. Exact-likelihood check (src/17): reduced-model NPE is too wide; LOSO 13 of 16.** Reduced model
(no selection, no reporting error, no quantization, fixed detection limits, s_v = 0.59) fitted by NUTS
(2 chains x 1,000) and by an NPE trained on 20,000 simulations with summaries of the full simulator's
kind (108). Six synthetic data sets from the prior. Data set 4 is a NUTS failure (2,000 divergences;
one chain stuck) and is set aside. On the other five, NPE posterior means are close to NUTS (median
|difference| 0.2-0.4 NUTS SD for a0, m0, g_gas, p_find, log_sigma0) but NPE is wider: SD ratio median
1.3-2.1 for the main parameters and 3.6-7 for tau_a and tau_m. The pre-written criterion (|mean
difference| < 0.3 SD and SD ratio 0.7-1.4 on 5 of 6 data sets; the plan's own rule is medians within
0.25 SD) fails for every hyperparameter. Reading: the summaries lose information, above all about
between-state spread, and 20,000 simulations are few for 108 summaries; the error is toward
under-confidence, not bias. Consequence for the full model: SBC can pass while posteriors are wider than
the data justify. Options, to decide with the advisor: (a) learned summary embedding and more
simulations; (b) per-state summaries of the full positive-value distribution; (c) a hybrid, using NUTS
for the tractable emission-and-measurement part and simulation only for selection. Real records under
the reduced model (NUTS, 0 divergences): p_find 0.98 (0.94-0.995), log_sigma0 1.36 (sigma ~3.9),
tau_a 2.9, b_gas -1.7; not results (the reduced model ignores selection and reporting error, and
treats the 85 positive values below their method's limit as detections). Leave-one-state-out
(population-level prediction of a held-out state): 13 of 16 statistics inside the 90% interval; misses
are CA's 10% quantile (observed lower), NM's 90% quantile (observed 7.3 vs upper 6.3, log g/h; that fit
had 1,000 divergences) and CO's detection share (0.66 vs upper 0.55), consistent with CO being a
high-emission state. The intervals are wide (tau_a ~2.9), so passing says little.

**D33. NPE ensemble on simulator v0.2 passes SBC and coverage.** Five flows (seeds 100-104, 100,000
simulations, 100 hidden units), posteriors mixed equally. On 200 held-out simulations: SBC passes for
all 16 parameters by the plan's rule (20 rank bins, chi-square p >= 0.01; also all p > 0.05); expected
coverage of 90% intervals 0.86-0.94, of 50% intervals 0.44-0.55. Posterior predictive check on the real
summaries: 2 of 226 outside the 99% range (about what chance gives): Pennsylvania screening detections
(observed 0; PA screening records are placeholders, D20) and the Driscoll 10% quantile (observed
lower than predicted: the chamber finds more very small leaks than the log-normal allows). Real-data
posterior (median, 90%): a0 3.4 (2.2-4.7), i.e. almost every well leaks at some rate; log_sigma0 1.43
(1.35-1.53), i.e. a log-SD of about 4.2; p_find 0.78 (0.62-0.90); rho -0.81 (-1.55, 0.27); omega 0.06;
kappa 0.18. Still TRIAL: D32 shows that summary-based NPE is wider than the exact posterior on the
reduced model, so calibrated here may still mean inefficient.

**D34. Decision results are set by the upper tail; reported with caps.** With the D33 posterior and no
cap, the mean leak rate per candidate well is in the kg/h range, because a log-normal with log-SD ~4.2
puts mass far above any recorded value (DOI max 46 kg/h; Driscoll max 1.4 kg/h, 1.5% of wells > 1 kg/h).
The model's share > 10 g/h among candidates (25%, 17-35%) matches Driscoll's random LA sample (26%);
its share > 1 kg/h (4.4%) is about three times Driscoll's. Uncapped values are therefore not reported.
src/16 gains --cap; runs at 50 kg/h and 1 kg/h (100 draws each). Robust across caps: P0 (observed
selection), P1 (random) and P2 (rank by state and well type) avoid the same ~10% of emissions with a
10% budget; one screening visit then ranking by the reading (P3) avoids 67-77%; the oracle 84-98%.
Break-even cost per screened well c* is flat in the detection limit up to ~100 g/h and scales with
p_find: cap 1 kg/h $2.0k / 5.0k / 7.5k / 9.3k at p_find 0.25 / 0.5 / 0.75 / 1; cap 50 kg/h $15k / 38k
/ 62k / 85k. E6 (screening with the 2023-2025 method mix vs P0, per plugged well): $28k (15-39k) at
cap 1 kg/h and $227k (83-592k) at cap 50 kg/h; x 2,259 wells: $63M and $513M. Next: the plan's
log-Student-t tail and a held-out check of the top 1% (Sec. 3), and a two-component (low/high)
leak-rate mixture as an alternative; until then the cap is the stated sensitivity.

## Open decisions (to settle before freezing the analysis plan)

- Detection-limit priors: Hi-Flow and OGI updated from Riddick et al. 2023, Ravikumar et al. 2018 and
  Zimmerle et al. 2020 (now medium confidence); handheld and RMLD-type screening still low confidence.
- High-emitter threshold (currently 10 g/h).
- Treatment of the two ~46 kg/h Adams County, CO wells (keep; report leave-top-k-out).
- Spatial-block size for cross-validation (currently 0.5°).
- Requests to OK OCC (meaning of 0), IDNR (IL zeros), KGS (KY measured wells), LDENR (below).
- What LA's "Hiflow/Method 21" means in practice: Method 21 screening first and Hi-Flow only on hits, or both on every well (ask LDENR).
