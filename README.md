# Orphaned well methane: data linkage and exploratory analysis

Work in progress. This repository assembles public data on plugged orphaned oil and gas wells in the
United States and links them to well inventories and state well records. See `docs/data_notes.md`
for sources and known issues and `docs/decision_log.md` for data-handling choices.

## Layout

```
data/raw/            frozen 2026-11-22 snapshot (do not edit; see MANIFEST.sha256)
  doi/               DOI Orphaned Wells Program well-level layer (+ schema, + doi.gov xlsx)
  usgs/              USGS 2022 DOW dataset (selection pool) and 2026 database subset
  states/            CO and NY bulk files; LA/NM/CA/WV query caches written by src/02
  research/          Driscoll et al. (2025) Louisiana chamber data (Zenodo)
data/reference/      literature-based reference tables (detection limits)
data/processed/      outputs of src/01–08 (regenerated; *.pkl are intermediates)
outputs/             tables and run log
src/                 analysis scripts, run in numeric order from the repo root
docs/                data notes, detection limits, decision log
```

## Reproduce

```
pip install -r requirements.txt
python src/01_build_linkage.py           # DOI ↔ USGS 2022 linkage, method classes
python src/02_state_attributes.py        # state well attributes (uses cached queries)
python src/03_selection_measurement.py   # selection AUC, detection by method
python src/04_prediction_prioritization.py
python src/05_attribute_value.py
python src/06_driscoll_pairing.py
python src/07_export_tables.py
python src/08_detection_empirical.py    # detection by refined method group
```

`src/00_download.py <date>` re-downloads everything into a new, separately dated folder; it never
overwrites the frozen snapshot. The Kansas well file (49 MB) and Oklahoma's well list are not bundled;
Kansas is optional (see `src/02_state_attributes.py`).

## Sources

- DOI Orphaned Wells Program feature service:
  https://services.arcgis.com/4OV0eRKiLAYkbH2J/arcgis/rest/services/DOI_BIL_Orphaned_Wells_Web_Map_Joined/FeatureServer
- USGS 2022 DOW dataset: https://doi.org/10.5066/P91PJETI · USGS 2026 DOWDB: https://doi.org/10.5066/P13FHBYG
- Colorado ECMC wells: https://ecmc.state.co.us/documents/data/downloads/gis/WELLS_SHP.ZIP
- New York DEC wells: https://data.ny.gov/api/views/szye-wmt3
- Louisiana SONRIS, New Mexico OCD, California WellSTAR, West Virginia TAGIS ArcGIS services (URLs in `src/02`)
- Driscoll et al. (2025) data: https://doi.org/10.5281/zenodo.13646591
