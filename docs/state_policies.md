# State and federal policy: what it means for the selection and measurement models

Compiled 2026-11-23. Per-state table: `data/reference/state_policies.csv`. Data cross-check:
`src/09_policy_crosscheck.py` → `outputs/state_measurement_crosscheck.csv`.

## Federal rules changed twice, so measurement depends on grant type and date

| Period | Grant | Methane rule |
|---|---|---|
| Aug 2022 onward | Initial grants ($25M per state; $5M AR, MS) | No measurement required. DOI's 2022 FAQ says the statutory certifications do not include pre-plugging measurement. Six states (incl. CA, MS, WV) committed to before/after measurement anyway. |
| Jul 2023 – Jul 2025 | Formula grants (Phase 1, 2) | Screen every plugged well; if methane is above background, measure the rate before and after plugging (DOI July 2023 guidelines). |
| Jul 2025 onward | Revised formula and matching guidance; Phase 3 (May 2026) | Measurement no longer required unless state law requires it; estimates may be reported in place of measurements; the 2023 guidelines are optional (Oct 2025 FAQ). Existing awards keep their terms until amended. |

The DOI layer has no grant-type field, so the regime must be inferred from state and date. The
data show the change: the share of records with any methane value was 0.69 before July 2025 and
0.49 after; Texas fell from 0.58 to 0.09.

## Issues that change the analysis

1. **About 4,100 "zeros" are placeholders, not measurements.** In KS, IL, NE, AR, IN, PA, AZ and MS
   every reported value is 0. KS, IL, NE and AR record no method at all; IN, PA and AZ name only a
   handheld screening detector; MS names a Hi-Flow. Kansas used only its initial grant
   through 2025, which required no measurement. Decision D20: zeros with no recorded method are
   treated as missing, and handheld-screening zeros as non-detects with a high limit.

2. **Selection depends on the outcome in some states.** New Mexico used pre-plugging methane data to
   decide which wells to plug. New York's initial grant targeted high-priority leaking gas wells.
   Louisiana gives 32 points for a gas leak near a building. So within these states the plugged
   wells are not a random draw with respect to emissions; the selection model must allow selection
   to depend on emissions (via observed leak indicators), not only on well attributes.

3. **Pennsylvania's selection is driven by pluggers, not the state.** Its main federal program pays
   licensed pluggers per well, capped at $40,000 for wells up to 3,000 ft, first come first
   served, and the plugger screens each well before applying. Plugged PA wells will be shallow and
   cheap rather than high risk. Model PA selection separately.

4. **Some values are field averages, not measurements.** New Mexico assigned the average of
   measured wells to unmeasured wells in the same field (O'Brien: 16 of 26 measured). The value
   1.69 g/h appears 11 times, almost all on O'Brien Fee wells. These cannot be treated as
   independent measurements (D21).

5. **Louisiana Hi-Flow values are quantized.** Positive values are multiples of about 1.15 g/h
   (0.001 cfm), so each value is an interval of that width (D22). Florida's values are 0 or
   exactly 1 g/h, consistent with "<1 g/h" non-detects entered as 1.

6. **Unused measured data exist.** Kentucky Geological Survey and the Well Done Foundation measured
   80–100 Kentucky initial-grant wells (2022–24), but Kentucky's DOI records carry one value. North
   Dakota estimated pre-plugging emissions from a 2018 EPA study rather than measuring. USGS has
   published Montana (Kevin-Sunburst) orphan-well measurements. These could validate the model
   outside Louisiana.

7. **A parallel program with mandatory measurement.** The DOE/EPA Methane Emissions Reduction
   Program (MERP) plugs marginal conventional wells (not orphans) in CO, IL, KY, MI, NY, PA, TX,
   UT, VA, WV and others, and requires pre- and post-plugging measurement under DOE guidelines.
   Not in the DOI layer; a possible later data source, but a different population.

## Prioritization systems (IOGCC 2023)

Of the 31 states covered, 15 use a scoring tool. Methane or gas leaks are named or implied in AK, IN,
FL, LA, MS, KY, CO, KS, OH, NY and MT; actual measurement was used in prioritization only in KY and MT
(NM began doing so under the grant). Proximity to people and water is near-universal. AK, MI, NM, NY,
PA and WY group nearby wells to cut mobilization cost, which induces spatial clustering of plugged
wells (relevant to the spatial-block cross-validation). The July 2025 guidance adds "land use
priorities that increase domestic oil and gas production" as a ranking factor, which may shift
which wells are plugged after 2025.

## Open questions for state programs

- OK: what instrument and what a recorded 0 means; whether the three tests are averaged.
- LA: whether "Hiflow/Method 21" means Method 21 screening first and Hi-Flow only on hits.
- IL, NE, AR: whether any methane was measured, and why values are all 0.
- KGS (KY): access to the 80–100 measured wells.
- All: grant type (initial vs formula) per plugged well.

## Sources

- DOI, FAQs on initial grant guidance (Apr 2022): https://www.doi.gov/sites/doi.gov/files/faqs-on-dois-initial-grant-guidance-04.12.2022.pdf
- DOI press release, 24 states initial grants (Aug 2022): https://www.doi.gov/pressreleases/through-president-bidens-bipartisan-infrastructure-law-24-states-set-begin-plugging
- DOI/IBC formula grant webinar (Aug 2023): https://oklahoma.gov/content/dam/ok/en/iogcc/documents/webinars/2023/orphan_wells_state_formula_guidance_ppt_8.15.2023.pdf
- DOI methane measurement guidelines (Jul 2023): https://www.doi.gov/sites/default/files/orphaned-wells-methane-measurement-guidelines-july-2023-version.pdf
- DOI FY2024 report to Congress: https://www.doi.gov/sites/default/files/documents/2024-11/fy-2024-owpo-annual-congressional-reportfinal-publishing.pdf
- DOI revised state formula grant guidance (Jul 2025): https://edit.doi.gov/sites/default/files/documents/2025-07/revised-state-formula-grant-guidance508.pdf
- DOI revised formula grant FAQ (Jul 2025): https://edit.doi.gov/sites/default/files/documents/2025-07/faqowpostateformulagrants20250711.pdf
- DOI updated state grant FAQ (Oct 2025): https://www.doi.gov/sites/default/files/documents/2025-10/owpoupdated-state-grant-faqs20251022.pdf
- DOI Phase 3 formula grant guidance (May 2026): https://www.doi.gov/sites/default/files/documents/2026-05/owpophase-3-state-formula-grant-guidancefinalmay-2026.pdf
- IOGCC, State prioritization systems for orphan wells (2023): https://oklahoma.gov/content/dam/ok/en/iogcc/documents/publications/prioritization_report_7.10.23.pdf
- Stateline (Oct 2023), Oklahoma and Kansas program quotes: https://stateline.org/2023/10/12/billions-of-dollars-to-clean-up-abandoned-oil-and-gas-wells-will-only-make-a-dent/
- RRC Texas, Commissioner Wright testimony (Jul 2024): https://rrc.texas.gov/news/072524-rrc-commissioner-wright-highlights-issues-federal-orphan-well-plugging-program-testimony-congress
- Southwest Ledger, Oklahoma federal plugging: https://www.southwestledger.news/node/18099
- KGS methane project: https://scholars.uky.edu/en/projects/methane-emission-measuring-on-orphaned-oil-and-gas-wells/
- LSU CES, OSR progress report (Dec 2023): https://oes.lsu.edu/ces/publications/2023/osrprogressreport12-18-2023.pdf
- LSU CES orphan wells page: https://lsu.edu/ces/research/orphan_wells.php
- PA DEP Orphan Well Plugging Grant: https://www.pa.gov/agencies/dep/programs-and-services/oil-and-gas/legacy-wells.html
- PA DEP methane measurement: https://pa.gov/agencies/dep/programs-and-services/oil-and-gas/methane-emissions-measurement
- NM OCD well file, O'Brien Fee 19 #001 memo: https://ocdimage.emnrd.nm.gov/imaging/filestore/SantaFe/WF/20230424/30005605280000_04_24_2023_12_54_12.pdf
- NM OCD well file, Artesia Metex #038 memo: https://ocdimage.emnrd.nm.gov/Imaging/FileStore/SantaFe/WF/20230723/30015021380000_07_23_2023_01_30_43.pdf
- WVDEP initial grant abstract: https://www.highergov.com/grant/D22AP00164
- NYSDEC initial grant abstract: https://www.highergov.com/grant/D22AP00170
- NDIC initial grant narrative: https://dmr.nd.gov/dmr/sites/www/files/documents/Work%20with%20Us/Narrative%20Elements_ND.pdf
- Montana BOGC FFY23 Q1 narrative: https://bogfiles.dnrc.mt.gov//OrphanWell_IIJA_Funding\FFY23_Q1\FFY23_Q1_MT_Performance_Report_Narrative.pdf
- USGS Kevin-Sunburst data release: https://www.usgs.gov/data/data-methane-emission-measurements-and-isotope-geochemistry-samples-orphan-oil-and-gas-wells
- ADEQ press release (program partnership): https://www.azdeq.gov/press-releases/press-release-adeq-announces-new-arizona-orphaned-well-program-partnership-arizona
- ECMC MERP page: https://ecmc.colorado.gov/regulation/marginal-wells/merp
