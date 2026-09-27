# DXA-measured lean and fat mass differences with incretin-based therapy versus comparators: analysis code and data

Systematic review and meta-analysis of phase 3 randomized trials with DXA-derived lean body mass and fat mass in adults with obesity or type 2 diabetes.

**Authors:** Birkan Alaycı, Öykü Zeynep Gerçek  
**Registration:** PROSPERO CRD420261323497  
**Manuscript:** under peer review at *BMC Endocrine Disorders* (submission 69e5b6d9; second revised version, September/October 2026)  
**Archive:** Zenodo concept DOI 10.5281/zenodo.19158245 (resolves to the latest version; this release = version 2.1, archived automatically from the GitHub release v2.1; version 2.0 = 10.5281/zenodo.22528130)

## Version 2.1 (September/October 2026): what changed

Second revision for *BMC Endocrine Disorders* (Reviewer 3, round 2). No extracted value changed. The analysis was simplified and the presentation made more conservative:

- **Principal presentation is now stratified** (placebo-controlled trials in obesity, k = 2; active-comparator trials in type 2 diabetes, k = 3); the five-trial pooled estimate is reported as a secondary descriptive summary (`output/figures/forest_*_stratified.png`).
- **Withdrawn:** meta-regression, Egger and Begg tests, the Doi plot and LFK index, and the pooled lean proportion. Section 10 of the script now writes trial-level proportions only (`output/tables/lean_proportion_trial_level.csv`).
- **Leave-one-out** is the primary robustness check (`output/tables/loo_*.csv`, with a `CI_includes_0` flag); the alternative estimators, RVE, three-level and Bayesian models are described as re-using the same five estimates, not as independent confirmation.
- **Bayesian analyses** are a secondary computational sensitivity analysis; priors are stated in kilograms and prior sensitivity is reported for both outcomes (`output/tables/brms_prior_sensitivity.csv`; half-Cauchy scales 0.5/1/2 kg for lean mass, 1/2/5 kg for fat mass).
- **Risk of bias** re-assessed against the specific effect estimate synthesized (effect of assignment for four trials; effect of adhering to intervention for SUSTAIN-8). STEP-1 domain 5 changed from "some concerns" to "low". Signalling questions with source justifications: `rob2/`.
- **Searches for unpublished DXA outcomes** (ClinicalTrials.gov outcome fields, EU Clinical Trials Register, CTIS, Novo Nordisk trial-document repository; sponsor correspondence): `searches/`.
- **Independent human audit** of every value entering the meta-analysis (Dr. Nur İlayda Genç, Koç University Hospital, 27 September 2026): `audit/` (extraction table as returned and signed auditor statement; 19 values, no discrepancy).
- Script `glp1_dxa_meta_v1_3.R` replaces `glp1_dxa_meta_v1_2.R`; `refresh_forests_v1_3.R` re-renders the forest plots without re-running Stan. `revision_analyses.R` v2 drops the pooled lean-fraction block.
- `data_raw.csv`: the LEAD-3 source note now states that premature termination visits were carried forward (Jendle 2009); no value changed.

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
| `glp1_dxa_meta_v1_3.R` | Main analysis (version 2.1): primary HKSJ-REML models, stratified presentation, τ² estimators, profile-likelihood CIs, RVE (CR2), three-level model, leave-one-out with CI flag, expanded inclusion (k = 6, 7), trial-level lean proportions, forest plots (stratified Fig. 3, five-trial Figure S1), Bayesian models (brms, bayesmeta) with prior sensitivity for both outcomes, key-results table. Meta-regression, small-study tests and the pooled lean proportion were removed at version 2.1. |
| `glp1_dxa_meta_v1_2.R` | Version 2.0 script (first revision), superseded by `glp1_dxa_meta_v1_3.R`; retained for the record. |
| `refresh_forests_v1_3.R` | Re-renders all forest plots from `data_raw.csv` in seconds (identical to Section 11 of the main script; no Stan). |
| `rob2/` | RoB 2 re-assessment: estimand per trial and answers to all signalling questions with the supporting source passage. |
| `searches/` | Registry and sponsor-repository searches for unpublished DXA outcomes (21 September 2026), the raw CTIS export, and the record-level screening decisions for the update searches (PubMed, CENTRAL, Scopus, 4 September 2026; Web of Science, 6 September 2026). |
| `audit/` | Independent human audit of data extraction: the table as returned by the auditor, the signed statement and a README. |
| `revision_analyses.R` | Analyses added at the first revision: prediction intervals for subgroups, HKSJ leave-one-out, primary analyses excluding STEP-1 (Table S8). Version 2 (2.1): pooled lean-fraction block removed. |
| `output/tables/*.csv` | All numerical outputs (key_results, tau2 comparisons, profile likelihood, RVE, three-level, LOO, Bayesian two-implementation comparison, prior sensitivity for both outcomes, trial-level lean proportions). |
| `output/figures/*.png` | Forest plots (stratified, five-trial overall, by population, by comparator) and bayesmeta posteriors. |
| `output/console_log.txt` | Full console log of the run that produced the deposited outputs. |
| `output/console_run_2026-09-04.txt` | Condensed record of the first corrected run (v1.1 script), including RVE and three-level results. |
| `CHANGELOG.md` | Version history. |
| `.zenodo.json`, `LICENSE` | Zenodo metadata for the archived release; licence (code MIT; data, tables and figures CC BY 4.0). |
| `01_*.py` … `05_*.py`, `run_all.py`, `requirements.txt`, `R_*_cross_validation.R`, `meta_analysis_data*.csv`, root-level `*.png`, `glp1_lbm_meta_analysis_code.zip`, `CROSS_CHECK_REPORT.md` | Version-1 material (March 2026: Python analysis with a DerSimonian–Laird primary model, R cross-validation scripts, version-1 data and figures). Superseded by the files above; retained unchanged for the record. |

## Reproducing

```r
# R >= 4.5; packages: metafor, clubSandwich, dplyr, readr, tibble, brms (Stan), posterior, bayesmeta
setwd("path/to/repo")
source("glp1_dxa_meta_v1_3.R")   # 15-25 min (eight brms fits)
source("revision_analyses.R")
source("refresh_forests_v1_3.R") # optional: forest plots only, a few seconds
```

Outputs are written to `output/`. Bayesian summaries are subject to Monte Carlo variation (about ±0.03 kg in posterior means and ±1 percentage point in P(μ<0) between runs); `set.seed(20260508)` is used but results depend on the Stan/brms versions installed. Session information is printed at the end of the console log.

## Data sources

- SURMOUNT-1: Look et al., Diabetes Obes Metab 2025 (DXA substudy; absolute-kg estimated treatment differences in the Results text and Figure S2)
- STEP-1: Wilding et al., N Engl J Med 2021, Supplementary Table S5 (treatment-policy estimand)
- SUSTAIN-8: McCrimmon et al., Diabetologia 2020 (confirmatory on-treatment analysis)
- LEAD-2 / LEAD-3: redacted Novo Nordisk clinical trial reports NN2211-1572 and NN2211-1573, publicly downloadable from novonordisk-trials.com (the sponsor directed us to them on 10 March 2026 in reply to a request for the lean-mass variance estimates); DXA results also published in summary form by Jendle et al., Diabetes Obes Metab 2009
- S-LiTE: Lundgren et al., N Engl J Med 2021, Supplementary Table S6
- BARI-OPTIMISE: Mok et al., JAMA Surg 2023, Table 2

## License

Code: MIT License (see `LICENSE`). Extracted trial-level data, output tables and figures: CC BY 4.0. The underlying trial reports remain the property of their publishers/sponsors.
