# Cross-Check Report: GLP-1 LBM DXA Meta-Analysis Code

**Date:** 2026-04-04
**Repository:** birkanalayci-dev/glp1-lbm-dxa-meta-analysis

---

## Summary

A systematic cross-check of all code files, data sources, and R cross-validation
scripts was performed. Two critical issues were identified and fixed.

---

## Issue 1: File Content Mismatch (Upload Error)

The original git upload placed **wrong contents** into the Python script files.
The zip archive (`glp1_lbm_meta_analysis_code.zip`) contained the correct files
in a `scripts/` subdirectory, but during upload the files were flattened and
contents got swapped:

| Git repo file | Actual content (WRONG) | Correct content (from zip) |
|---|---|---|
| `01_primary_meta_analysis.py` (1,070 B) | MIT License text | Primary DL analysis script (17,159 B) |
| `02_publication_bias.py` (994 B) | Old CSV data | Publication bias script (13,921 B) |
| `03_advanced_sensitivity.py` (1,847 B) | `run_all.py` runner | Advanced sensitivity script (17,792 B) |
| `04_placebo_subgroup_advanced.py` (13,921 B) | Script 02 code | Placebo subgroup script (17,290 B) |
| `meta_analysis_data.csv` (196,371 B) | PNG image data | CSV data (994 B) |
| `download` (39 B) | requirements.txt content | Should be named `requirements.txt` |

Additionally missing from git:
- `05_tau2_estimators_rve.py` (20,109 B) — tau-squared estimators and RVE
- `run_all.py` (1,847 B) — sequential runner script

**Fix:** All files restored from the zip archive to their correct contents.

---

## Issue 2: STEP-1 Data Discrepancy (Python vs R/CSV v5)

| Source | STEP-1 MD (kg) | STEP-1 SE (kg) |
|---|---|---|
| All 5 Python scripts (zip) | -3.43 | 0.67 |
| R lean mass cross-validation | -1.79 | 0.35 |
| `meta_analysis_data_v5.csv` | -1.79 | 0.35 |
| README (DL pooled MD = -1.94) | Matches R/v5 | — |

- **Python scripts** used the original ETD directly from Wilding 2021 Suppl Table S5
  (-3.43 kg, SE 0.67)
- **R scripts and v5 CSV** use -1.79 kg (SE 0.35), described as
  "%-point ETD x pooled baseline LBM 52.1 kg"
- The **manuscript results** (DL MD = -1.94, I2 = 83.9%) match the R/v5 data,
  confirming -1.79 is the current correct value
- All other 4 studies are consistent across all sources

**Fix:** STEP-1 data updated to md=-1.79, se=0.35 in all 5 Python scripts,
with comments noting the correction.

---

## Verified Consistent Data (All Sources Agree)

| Study | MD (kg) | SE (kg) | n_int | n_ctrl |
|---|---|---|---|---|
| SURMOUNT-1 | -4.40 | 0.61 | 124 | 36 |
| STEP-1 (corrected) | -1.79 | 0.35 | 95 | 45 |
| SUSTAIN-8 | -0.78 | 0.41 | 53 | 61 |
| LEAD-2 | -1.535 | 0.473 | 30 | 32 |
| LEAD-3 | -1.508 | 0.540 | 17 | 12 |

---

## R Cross-Validation Status

Both R scripts (`R_lean_mass_cross_validation.R` and `R_fat_mass_cross_validation.R`)
are internally consistent with `meta_analysis_data_v5.csv` and the README key results.
No changes needed.

---

## Files Added/Restored

- `01_primary_meta_analysis.py` — restored from zip, STEP-1 updated
- `02_publication_bias.py` — restored from zip, STEP-1 updated
- `03_advanced_sensitivity.py` — restored from zip, STEP-1 updated
- `04_placebo_subgroup_advanced.py` — restored from zip, STEP-1 updated
- `05_tau2_estimators_rve.py` — restored from zip (was missing), STEP-1 updated
- `run_all.py` — restored from zip (was missing)
- `meta_analysis_data.csv` — restored from zip (was PNG image)
- `requirements.txt` — renamed from `download`
