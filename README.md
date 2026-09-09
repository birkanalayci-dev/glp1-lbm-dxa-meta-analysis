# DXA-measured body composition effects of GLP-1 and dual GIP/GLP-1 receptor agonists: analysis code and data

Systematic review and meta-analysis of phase 3 randomized trials with DXA-derived lean body mass and fat mass in adults with obesity or type 2 diabetes.

**Authors:** Birkan Alaycı, Öykü Zeynep Gerçek  
**Registration:** PROSPERO CRD420261323497  
**Manuscript:** under peer review at *BMC Endocrine Disorders* (submission 69e5b6d9; revised version September 2026)  
**Archive:** Zenodo concept DOI 10.5281/zenodo.19158245 (resolves to the latest version; this release = version 2.0, archived automatically from the GitHub release v2.0)

## Version 2 (September 2026): what changed

During peer review all extracted values were re-verified against the primary source tables. Three lean-mass values in version 1 were incorrect and have been corrected (see `CHANGELOG.md` and `data_raw_v1_submitted.csv` for the superseded values):

| Trial | v1 lean MD (kg) | v2 lean MD (kg) [95% CI] | Reason |
|---|---|---|---|
| STEP-1 | −1.79 | −3.43 [−4.74, −2.13] | Wilding 2021 Table S5 reports the ETD directly in kg; v1 had treated it as percentage points and rescaled by baseline lean mass |
| LEAD-2 | −1.535 | −2.816 [−3.979, −1.652] | v1 entered the within-arm change of the liraglutide 1.8 mg arm instead of the ETD vs glimepiride + metformin (CTR NN2211-1572, Table 11-51) |
| LEAD-3 | −1.508 | −0.956 [−2.697, +0.785] | same error; ETD vs glimepiride (CTR NN2211-1573, Table 11-44) |
| SUSTAIN-8 | SE 0.41 | SE 0.421 | rounding of the CI-derived SE |

Fat-mass values were correct in v1. The primary analysis (REML + Hartung–Knapp–Sidik–Jonkman) is unchanged in specification. Pooled lean-mass difference: v1 −1.96 kg [−3.64, −0.27] → v2 −2.49 kg [−4.45, −0.53]; the type 2 diabetes subgroup is no longer homogeneous (I² 5% → 73%).

The analysis is implemented in R. The Python scripts and R cross-validation scripts from the March 2026 (version 1) deposit are retained for the record but are superseded; they use the version-1 data and a DerSimonian–Laird primary model.

## Files

| File | Content |
|---|---|
| `data_raw.csv` | Trial-level data (5 primary + 2 sensitivity comparisons). `source` gives the exact source table for every value. |
| `data_raw_v1_submitted.csv` | Data as used in the originally submitted manuscript (superseded; kept for transparency). |
| `glp1_dxa_meta_v1_2.R` | Main analysis: primary HKSJ-REML models, subgroups, meta-regression, τ² estimators, profile-likelihood CIs, RVE (CR2), three-level model, leave-one-out, small-study effects (Egger, Begg, metasens LFK/Doi), expanded inclusion (k = 6, 7), descriptive lean share, forest plots, Bayesian models (brms, bayesmeta), prior sensitivity, key-results table. |
| `revision_analyses.R` | Analyses added at peer review: prediction intervals for subgroups, HKSJ leave-one-out, primary analyses excluding STEP-1, descriptive lean share. |
| `output/tables/*.csv` | All numerical outputs (key_results, tau2 comparisons, profile likelihood, RVE, three-level, LOO, small-study effects, expanded sets, Bayesian cross-validation, prior sensitivity, lean share). |
| `output/figures/*.png` | Forest plots (overall, by population, by comparator), Doi plots, bayesmeta posteriors. |
| `output/console_log.txt` | Full console log of the run that produced the deposited outputs. |
| `output/console_run_2026-09-04.txt` | Condensed record of the first corrected run (v1.1 script), including RVE and three-level results. |
| `CHANGELOG.md` | Version history. |
| `.zenodo.json`, `LICENSE` | Zenodo metadata for the archived release; licence (code MIT; data, tables and figures CC BY 4.0). |
| `01_*.py` … `05_*.py`, `run_all.py`, `requirements.txt`, `R_*_cross_validation.R`, `meta_analysis_data*.csv`, root-level `*.png`, `glp1_lbm_meta_analysis_code.zip`, `CROSS_CHECK_REPORT.md` | Version-1 material (March 2026: Python analysis with a DerSimonian–Laird primary model, R cross-validation scripts, version-1 data and figures). Superseded by the files above; retained unchanged for the record. |

## Reproducing

```r
# R >= 4.5; packages: metafor, clubSandwich, dplyr, readr, tibble, brms (Stan), posterior, bayesmeta, metasens
setwd("path/to/repo")
source("glp1_dxa_meta_v1_2.R")   # ~10-15 min (Bayesian models)
source("revision_analyses.R")
```

Outputs are written to `output/`. Bayesian summaries are subject to Monte Carlo variation (about ±0.03 kg in posterior means and ±1 percentage point in P(μ<0) between runs); `set.seed(20260508)` is used but results depend on the Stan/brms versions installed. Session information is printed at the end of the console log.

## Data sources

- SURMOUNT-1: Look et al., Diabetes Obes Metab 2025 (DXA substudy; absolute-kg estimated treatment differences in the Results text and Figure S2)
- STEP-1: Wilding et al., N Engl J Med 2021, Supplementary Table S5 (treatment-policy estimand)
- SUSTAIN-8: McCrimmon et al., Diabetologia 2020 (confirmatory on-treatment analysis)
- LEAD-2 / LEAD-3: Novo Nordisk clinical trial reports NN2211-1572 and NN2211-1573, obtained through the external researcher data access service (novonordisk-trials.com); DXA results also published by Jendle et al., Diabetes Obes Metab 2009
- S-LiTE: Lundgren et al., N Engl J Med 2021, Supplementary Table S6
- BARI-OPTIMISE: Mok et al., JAMA Surg 2023, Table 2

## License

Code: MIT License (see `LICENSE`). Extracted trial-level data, output tables and figures: CC BY 4.0. The underlying trial reports remain the property of their publishers/sponsors.
