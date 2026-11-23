# Data notes

Snapshot date: 2026-11-22. Data-handling choices are logged in `docs/decision_log.md`.

## Sources and what was used from each

| Source | Records | Used for |
| --- | --- | --- |
| DOI Orphaned Wells Program ArcGIS layer | 13,108 plugged wells (Oct 2022 – Jun 2026) | Pre-/post-plugging methane (g/h), detection method, plugging cost, well type, location |
| USGS Documented Unplugged Orphaned Well dataset, 2022 | 117,672 wells | Pre-program inventory; linked to DOI plugged wells |
| USGS Documented Orphaned Well Database, 2026 (subset of fields) | 130,121 wells | Reference only |
| State well databases (LA, NM, CA, WV via ArcGIS; CO, NY bulk files) | — | Spud date, depth, last production |
| Driscoll et al. (2025), Zenodo | 132 Louisiana wells | Chamber measurements; pairing with DOI records |

## Linkage

- DOI plugged wells are linked to the 2022 inventory by 10-digit API number, then by nearest location
  within the same state (≤ 100 m): 4,855 + 1,515 = 6,370 wells.
- State attributes are joined by API number. West Virginia's date fields describe the current permit,
  not the spud date, so WV wells have no age.
- The placeholder API `0000000000` is treated as missing.

## Known data issues

- 6,080 of 8,499 pre-plugging values have no detection method; nearly all are zero.
- Reported detection methods are free text and are grouped into classes (see `src/01_build_linkage.py`).
- Some values are rounded (e.g. 0.01, 1.0 g/h).
- Oklahoma's public well list has no depth or spud date; Texas attributes are not yet parsed.

## Outputs

- `data/processed/doi_wells_analysis.csv` — one row per DOI plugged well with method class, linkage
  and state attributes.
- `data/processed/usgs2022_inventory_linked.csv` — 2022 inventory with a `plugged` flag.
- `outputs/driscoll_doi_pairs.csv` — wells measured by both methods (5; column `link` = api / serial / coord).
- `outputs/driscoll_nearest_doi.csv` — nearest DOI LA record and distance for all 132 Driscoll wells.
- `data/raw/usgs_mt_kevin_sunburst/` — USGS data release doi:10.5066/P13A7RXN: repeated flow-meter measurements (56 sessions, 15 wells) and gas geochemistry (19 wells), Kevin-Sunburst Dome, Toole County, MT. Latitude/longitude column names are swapped in the geochemistry file.
- `outputs/usgs_mt_wells.csv` — per-well summary, linkage to the 2022 inventory and nearest DOI record (src/10).
- `outputs/state_measurement_crosscheck.csv` — per state: share with a value before/after July 2025, zero share, no-method share, repeated positive values (src/09).
- `data/reference/state_policies.csv` — per-state prioritization and measurement practice with sources; narrative in `docs/state_policies.md`.
- `outputs/hiflow_model_by_state.csv` — Hi-Flow maker/model named in DOI method strings, by state.
- `data/raw/states/LA_serial_coords.csv` — SONRIS coordinates for the 132 Driscoll serials (cached query).
- `outputs/detection_by_method_empirical.csv`, `outputs/method_string_audit.csv` — detection by method
  group and the raw method strings behind each group (see `docs/detection_limits.md`).
- `outputs/run_log.txt` — console output of the scripts in `src/`.
