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

## Open decisions (to settle before freezing the analysis plan)

- Detection-limit prior for each method class (sources: DOI 2023 guidelines, Haase & Gianoutsos 2025,
  Driscoll 2025 instruments).
- High-emitter threshold (currently 10 g/h).
- Treatment of the two ~46 kg/h Adams County, CO wells (keep; report leave-top-k-out).
- Spatial-block size for cross-validation (currently 0.5°).
- Whether to coordinate-match the 22 Driscoll chamber wells that lack a valid API.
