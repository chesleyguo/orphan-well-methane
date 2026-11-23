# Detection limits by method group

Working table for the measurement model. Machine-readable version: `data/reference/detection_limits.csv`.
Empirical columns come from `src/08_detection_empirical.py` on the 2026-11-22 snapshot. The proposed
priors are judgment calls from the sources below and should be revisited before the analysis plan is frozen.

## Summary

| Method group | Wells | Detect rate | Smallest non-zero (g/h) | Documented basis | Proposed threshold prior (g/h) | Confidence |
| --- | --- | --- | --- | --- | --- | --- |
| Chamber / flow-through cell | 186 | 20% | 0.0004 | MFM chamber ≈ 0.77 g/h; high-precision chambers ≈ 0.0001 g/h | 0.001–0.5 | medium |
| Flow meter + GC | 204 | 60% | 0.001 | DOI: limit must be ≤ 1 g/h | 0.001–0.1 | medium |
| Hi-Flow sampler | 818 | 27% | 0.0003 | DOI ≤ 1 g/h; LA finds the leak with OGI first | 0.5–20 (effective) | low |
| Laser analyzer (quantitative) | 92 | 49% | 1.0 | DOI lists as high sensitivity | 0.5–3 | medium |
| Federal guidelines (unspecified) | 259 | 32% | 0.7 | Non-detects recorded as < 1 g/h | 0.7–1.5 | medium |
| Screening: OGI camera | 54 | 17% | 1.0 | Less sensitive; accuracy falls with distance | 5–100 | low |
| Screening: handheld detector | 649 | 3% | 0.02 | ≤ 100 ppm at the source; non-detect assumed ≤ 1 g/h | 1–50 (effective) | low |
| Screening: standoff laser (RMLD) | 80 | 0% | — | ≤ 100 ppm-m at 1 m | 1–50 (effective) | low |

Excluded or undecided: calculated (13), FID (4), other/unspecified (60), "ND" (70, ambiguous), no method (6,010).

## What the table implies for the model

1. **Detection is not a fixed floor.** Several groups report values far below their nominal limit
   (Hi-Flow down to 0.0003 g/h), yet miss large leaks. The measurement model should use a detection
   probability that rises with the true rate, with a group-specific midpoint drawn from the prior above,
   rather than a hard cut-off.
2. **Hi-Flow is the weakest link and the largest quantitative group.** In Louisiana the Hi-Flow only
   measures leak points an OGI camera has found, so its effective limit inherits OGI's. The two paired
   wells where Hi-Flow recorded 0 while a chamber measured 40–48 g/h, and the 23% vs 96% detection gap
   reported by LSU, both point the same way.
3. **Screening devices do not measure rates.** Handheld detectors and RMLD-type lasers give
   concentrations. Their zeros should enter the model as "below an unknown threshold", or be left out
   of the rate likelihood; they should never be read as 0 g/h. DOI's own rule assumes non-detects are
   ≤ 1 g/h, but 0 detections in 577 Sensit and Bascom Turner screenings is far below what other methods
   find, so that assumption looks too optimistic here.
4. **Rounding.** Recurring values (0.004, 0.01, 1.0, 1.69 g/h) should be treated as interval-censored.

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

Still to read (would firm up the low-confidence rows):

- Zimmerle, D., et al. (2020). Detection limits of optical gas imaging for natural gas leak detection
  in realistic controlled conditions (cited in the DOI guidelines) — OGI detection-probability curves.
- The 2023 *Sensors* study comparing five methods at ~1 g/h (PMC10675349) — Hi-Flow and chamber
  uncertainty at low rates.
- Pekney et al. (2018), Hillman State Park — the Hi-Flow adaptation used in Louisiana.
