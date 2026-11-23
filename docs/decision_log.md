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

## Open decisions (to settle before freezing the analysis plan)

- Detection-limit priors: Hi-Flow and OGI updated from Riddick et al. 2023, Ravikumar et al. 2018 and
  Zimmerle et al. 2020 (now medium confidence); handheld and RMLD-type screening still low confidence.
- High-emitter threshold (currently 10 g/h).
- Treatment of the two ~46 kg/h Adams County, CO wells (keep; report leave-top-k-out).
- Spatial-block size for cross-validation (currently 0.5°).
- Requests to OK OCC (meaning of 0), IDNR (IL zeros), KGS (KY measured wells), LDENR (below).
- What LA's "Hiflow/Method 21" means in practice: Method 21 screening first and Hi-Flow only on hits, or both on every well (ask LDENR).
