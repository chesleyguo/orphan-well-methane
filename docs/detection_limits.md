# Detection limits by method group

Working table for the measurement model. Machine-readable version: `data/reference/detection_limits.csv`.
Empirical columns come from `src/08_detection_empirical.py` on the 2026-11-22 snapshot. The proposed
priors are judgment calls from the sources below and should be revisited before the analysis plan is frozen.

## Summary

| Method group | Wells | Detect rate | Smallest non-zero (g/h) | Documented basis | Proposed threshold prior (g/h) | Confidence |
| --- | --- | --- | --- | --- | --- | --- |
| Chamber / flow-through cell | 186 | 20% | 0.0004 | MFM chamber ≈ 0.77 g/h; high-precision chambers ≈ 0.0001 g/h | 0.001–0.5 | medium |
| Flow meter + GC | 204 | 60% | 0.001 | DOI: limit must be ≤ 1 g/h | 0.001–0.1 | medium |
| Hi-Flow sampler | 818 | 27% | 0.0003 | Instrument floor ≈ 0.6 g/h, accurate at 1 g/h if all gas is captured; LA finds the leak with OGI first | 0.6–2 if found; finding inherits OGI (20–150) | medium |
| Laser analyzer (quantitative) | 92 | 49% | 1.0 | DOI lists as high sensitivity | 0.5–3 | medium |
| Federal guidelines (unspecified) | 259 | 32% | 0.7 | Non-detects recorded as < 1 g/h | 0.7–1.5 | medium |
| Screening: OGI camera | 54 | 17% | 1.0 | 50% detection ≈ 20 g/h at 6 m; 90% ≈ 135 g/h for experienced surveyors | 20–150 (50% point) | medium |
| Screening: handheld detector | 649 | 3% | 0.02 | ≤ 100 ppm at the source; non-detect assumed ≤ 1 g/h | 1–50 (effective) | low |
| Screening: standoff laser (RMLD) | 80 | 0% | — | ≤ 100 ppm-m at 1 m | 1–50 (effective) | low |

Excluded or undecided: calculated (13), FID (4), other/unspecified (60), "ND" (70, ambiguous), no method (6,010).

## What the table implies for the model

1. **Detection is not a fixed floor.** Several groups report values far below their nominal limit
   (Hi-Flow down to 0.0003 g/h), yet miss large leaks. The measurement model should use a detection
   probability that rises with the true rate, with a group-specific midpoint drawn from the prior above,
   rather than a hard cut-off.
2. **Hi-Flow is the weakest link and the largest quantitative group, but the instrument is not the
   problem.** In a controlled test at 1 g/h the Hi-Flow was within about 2% of the true rate, and its
   stated lower bound is 0.6 g/h; the authors name capturing all of the gas as the main limitation.
   In Louisiana the Hi-Flow only measures leak points an OGI camera has found, so detection inherits
   OGI's limit (50% detection ≈ 20 g/h at 6 m). Model it in two stages: was the leak point found, and
   if so, was it quantified. This fits the two paired wells (Hi-Flow 0, chamber 40–48 g/h) and the
   23% vs 96% detection gap reported by LSU. A fifth pair, found through the SONRIS serial number,
   adds a well that both methods put near zero.

3. **The Hi-Flow instrument is the Semtech Hi-Flow 2, not the Bacharach.** Of 975 DOI records
   that mention a high-flow sampler, 350 name Semtech (LA 251, WV 59, OK 21, AL 6, others) and none
   name Bacharach, which stopped being made before the program began. California's 238 records say
   "Hi-Flow 2" or "Hi-Flow II" without a maker; that is the Semtech product name. Louisiana's
   "Hiflow/Method 21" (291) and the "High Flow Sampling" entries in MS, MT, NE and AZ do not name a
   model. The Semtech uses a tunable diode laser, so the Bacharach sensor-transition failure
   (Howard 2015; Connolly 2019) does not apply. Its spec range of 0.0005-25 CFM works out to about
   0.6 g/h at the bottom, consistent with the 0.6-2 g/h prior. The Hi-Flow zeros at 40-48 g/h are
   therefore far above the instrument's floor, which points to leak-finding or capture rather than
   sensitivity. Table: `outputs/hiflow_model_by_state.csv`.
4. **Screening devices do not measure rates.** Handheld detectors and RMLD-type lasers give
   concentrations. Their zeros should enter the model as "below an unknown threshold", or be left out
   of the rate likelihood; they should never be read as 0 g/h. DOI's own rule assumes non-detects are
   ≤ 1 g/h, but 0 detections in 577 Sensit and Bascom Turner screenings is far below what other methods
   find, so that assumption looks too optimistic here.
5. **Rounding.** Recurring values (0.004, 0.01, 1.0, 1.69 g/h) should be treated as interval-censored.

## Sources

Read for this table:

- U.S. DOI Orphaned Wells Program Office (2023). *Assessing Methane Emissions from Orphaned Wells …
  Methane Measurement Guidelines*, July 2023 version. Detection requirements (≤ 1 g/h; 100 ppm;
  100 ppm-m at 1 m; non-detect recorded as < 1 g/h); OGI as a less sensitive screening tool.
- Haase, K. B., & Gianoutsos, N. J. (2025). Evaluating mass flow meter measurements from chambers for
  greenhouse gas emissions from orphan wells and other point sources. *AMT* 18, 4207–4226
  (abstract: leaks down to 0.77 g/h methane).
- Driscoll, F. S., et al. (2025). *JGR: Atmospheres* 130(15). Chamber with Picarro, Aeris and Gas
  Rover II; bare-ground control ≈ 0.0001 g/h; 125 of 132 wells above control. Data: Zenodo 13646591.
- Santos, O. L. A., et al. (2025). *Science of the Total Environment* 1005, 180791. LA program uses
  OGI to locate the emission point, then a SEMTECH Hi-Flow to measure it.
- Upton, G. B., Jr., Agerton, M., & Gupta, I. (2025). *Orphan and Idle Wells in Louisiana*. LSU CES.
  Hi-Flow detected methane at 23% of wells vs 96% with chambers.
- Agogino, A. M., et al. (2024). IMECE2024-145130. OGI detection accuracy declines with distance
  (GasVid leaks 5.3–2,051.6 g/h at 4.6–18.6 m).

- Riddick, S. N., Mbua, M., Riddick, J. C., Houlihan, C., Hodshire, A. L., & Zimmerle, D. J. (2023).
  Uncertainty quantification of methods used to measure methane emissions of 1 g CH4 h⁻¹. *Sensors*
  23(22), 9246. Hi-Flow +2.2% at 1 g/h (lower bound 0.6 g/h); chambers within ≈ ±6%; capture is the
  main Hi-Flow limitation; earlier work found under-estimation at 40–200 g/h.
- Ravikumar, A. P., et al. (2018). "Good versus good enough?" *ES&T* 52, 2368–2374. OGI median
  detection limit ≈ 20 g/h at 6 m; detection limit grows with distance as a power law.
- Zimmerle, D., Vaughn, T., Bell, C. S., Bennett, K., Deshmukh, P., & Thoma, E. (2020). Detection
  limits of optical gas imaging for natural gas leak detection in realistic controlled conditions.
  *ES&T*. 90% detection at ≈ 7 scfh (5.6–19.5) for experienced surveyors, ≈ 10× earlier estimates.
- Howard, T., et al. (2015). Sensor transition failure in the Bacharach Hi-Flow Sampler. Under-reporting
  by up to two orders of magnitude above ≈ 0.4 scfm.

Unit note: 1 scfm of methane ≈ 19.2 g/min, so 1 scfh ≈ 19.2 g/h (7 scfh ≈ 135 g/h; 0.4 scfm ≈ 460 g/h).
Natural gas is mostly but not entirely methane, so these conversions are upper bounds.

Still to read:

- Connolly, J. I., Robinson, R. A., & Gardiner, T. D. (2019). Assessment of the Bacharach Hi Flow
  Sampler characteristics and potential failure modes. Measurement, 145, 226-233. Noisy, non-linear
  readings in the switchover region between the two sensor modes.
- Sensors, Inc. SEMTECH HI-FLOW 2 product sheet and OOOOb readiness note (2024): tunable diode laser,
  measurable leak rate 0.0005-25 CFM, marketed as the replacement for the discontinued Bacharach unit.
- Pekney et al. (2018), Hillman State Park — the Hi-Flow adaptation used in Louisiana.
- A source for handheld concentration detectors and RMLD-type lasers expressed as a flow rate
  (the rows still marked low confidence).
