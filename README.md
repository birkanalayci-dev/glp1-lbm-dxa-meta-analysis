# GLP-1 and Dual GIP/GLP-1 Receptor Agonists — Body Composition Meta-Analysis

## Reproducible Code Repository

**Manuscript:** "Composition of Weight Loss Associated with Incretin-Based Therapy: A Systematic Review and Meta-Analysis of DXA-Measured Lean and Fat Mass in Randomized Controlled Trials"

**Journal:** *Diabetes, Obesity and Metabolism* — Manuscript ID: DOM-26-1453-OP

**Authors:** Birkan Alayci, Öykü Zeynep Gerçek

**PROSPERO:** CRD420261323497

---

## Overview

This repository contains the analytical code, input data, R cross-validation scripts, and output figures for a systematic review and meta-analysis examining the effect of incretin-based therapies on DXA-measured lean body mass and fat mass. Five Phase 3/3b RCTs (k = 5, n = 505 DXA-evaluated participants) were included.

Given the small number of studies, we employed **14 complementary statistical models** to ensure conclusions were not model-dependent:

| Model | Purpose | Script |
|-------|---------|--------|
| DerSimonian-Laird (DL) | Primary random-effects | `01` |
| Subgroup analyses (placebo/active/GLP-1 only) | Comparator-driven heterogeneity | `01` |
| Leave-one-out sensitivity | Individual study influence | `01` |
| Doi plot + LFK index | Small-study effects (k < 10) | `02` |
| Peters' weighted regression | Publication bias | `02` |
| Begg's rank correlation | Publication bias | `02` |
| Egger's test (supplementary) | Publication bias (discounted at k < 10) | `02` |
| HKSJ correction | Conservative CI at small k | `03` |
| Prediction interval (t-based, df=k-2) | Between-study heterogeneity range | `03` |
| Bayesian random-effects | Posterior probabilities, prior sensitivity | `03`, `04` |
| REML τ² estimator | Recommended by Cochrane Handbook 6.4 | `03` |
| Paule-Mandel τ² estimator | Iterative τ² estimation | `03` |
| Profile likelihood CI | Exact CI for μ and τ² | `03` |
| Robust variance estimation (CR2) | LEAD-2/3 within-publication dependency | `03` |
| Three-level model | Within-study correlation sensitivity | `03` |

## Repository Contents

```
├── README.md
├── LICENSE
├── meta_analysis_data_v5.csv           # Study-level data (5 primary + 2 sensitivity)
├── 01_primary_meta_analysis.py         # DL, subgroups, leave-one-out, forest plots
├── 02_publication_bias.py              # Doi/LFK, Peters, Begg, Egger, funnel plot
├── 03_advanced_sensitivity.py          # HKSJ, PI, Bayesian, REML, PM, RVE, 3-level
├── 04_placebo_subgroup_advanced.py     # Placebo subgroup: HKSJ, PI, Bayesian
├── R_lean_mass_cross_validation.R      # R metafor cross-validation (lean mass)
├── R_fat_mass_cross_validation.R       # R metafor cross-validation (fat mass)
├── forest_all.png                      # Figure 1: Overall forest plot (k=5)
├── forest_placebo.png                  # Figure 2A: Placebo subgroup (k=2)
├── forest_glp1.png                     # Figure S1: GLP-1 RA only (k=4)
├── funnel_plot.png                     # Figure S2: Funnel plot
├── doi_plot.png                        # Figure 4 / S3: Doi plot with LFK index
├── peters_test.png                     # Figure S4: Peters' regression
├── loo_sensitivity.png                 # Figure 3 / S5: Leave-one-out
├── forest_model_comparison.png         # Figure S6: Multi-model comparison (overall)
├── bayesian_posteriors.png             # Figure S8: Bayesian posteriors (overall)
├── forest_placebo_advanced.png         # Figure S7: Multi-model comparison (placebo)
├── bayesian_placebo.png                # Figure S9: Bayesian posteriors (placebo)
└── forest_tau2_comparison.png          # Figure S10: τ² estimator comparison
```

## Data

`meta_analysis_data_v5.csv` contains the following fields for each study comparison:

| Field | Description |
|-------|-------------|
| `study_id` | Study identifier |
| `drug` | Intervention agent |
| `dose_mg` | Dose(s) evaluated |
| `comparator` | Placebo or active comparator |
| `duration_wks` | Treatment duration (weeks) |
| `n_int` / `n_ctrl` | DXA subpopulation sample sizes |
| `lean_md_kg` | Mean difference in lean mass (kg) |
| `lean_se_kg` | Standard error of lean mass MD |
| `fat_md_kg` | Mean difference in fat mass (kg) |
| `fat_se_kg` | Standard error of fat mass MD |
| `analysis` | Primary (k=5) or sensitivity (k=6–7) |
| `source` | Data source (publication or clinical trial report) |

### Data Sources

- **SURMOUNT-1:** Look et al. (DOM, 2025) — ETD
- **STEP-1:** Wilding et al. (NEJM, 2021) — Table S5, treatment policy estimand
- **SUSTAIN-8:** McCrimmon et al. (Diabetologia, 2020) — Table 2 ANCOVA
- **LEAD-2:** Novo Nordisk CTR NN2211-1572
- **LEAD-3:** Novo Nordisk CTR NN2211-1573

## Key Results

### Lean Mass (Primary Outcome)

| Analysis | MD (kg) | 95% CI/CrI | Significance |
|----------|---------|------------|-------------|
| Overall (DL, k=5) | −1.94 | [−2.96; −0.93] | p = 0.0002 |
| Placebo (DL, k=2) | −3.05 | [−5.60; −0.49] | p = 0.020 |
| Active (DL, k=3) | −1.20 | [−1.73; −0.68] | p < 0.0001 |
| Overall HKSJ | −1.94 | [−3.62; −0.27] | p = 0.032 |
| 95% Prediction interval | −1.94 | [−5.66; +1.77] | Crosses zero |
| Overall Bayesian HC(0,2) | −1.94 | [−3.51; −0.43] | P(Δ<0) = 98.9% |
| Profile likelihood | −1.95 | [−2.97; −0.92] | Excludes zero |
| RVE (CR2, m=4) | −1.96 | [−4.01; +0.09] | p = 0.058 |

### Fat Mass (Secondary Outcome)

| Analysis | MD (kg) | 95% CI | Significance |
|----------|---------|--------|-------------|
| Overall (DL, k=5) | −5.27 | [−8.81; −1.72] | p = 0.004 |
| Placebo (DL, k=2) | −9.65 | [−14.86; −4.45] | p = 0.0003 |
| Active (DL, k=3) | −2.45 | [−4.42; −0.47] | p = 0.015 |

### Lean-to-Fat Ratio

| Subgroup | Lean % | Fat % |
|----------|--------|-------|
| Overall (k=5) | 27.0% | 73.0% |
| Placebo (k=2) | 24.0% | 76.0% |
| Active (k=3) | 32.9% | 67.1% |

## R Cross-Validation

Python results were independently cross-validated using R 4.5 with `metafor` and `clubSandwich`:

```r
# Lean mass
Rscript R_lean_mass_cross_validation.R

# Fat mass
Rscript R_fat_mass_cross_validation.R
```

All models (DL, REML, PM, HKSJ, PI, RVE, three-level) were verified to match within numerical precision (Δ < 0.01).

## Requirements

```
Python >= 3.8
numpy >= 1.20
scipy >= 1.7
matplotlib >= 3.4
```

## Usage

Run individual scripts:

```bash
python 01_primary_meta_analysis.py
python 02_publication_bias.py
python 03_advanced_sensitivity.py
python 04_placebo_subgroup_advanced.py
```

## Statistical Notes

- **Primary analysis:** DerSimonian-Laird random-effects
- **Prediction intervals:** t-distribution (df = k−2), consistent with R metafor `predict()`
- **HKSJ at k=2:** Non-informative (df=1, t₀.₉₇₅ = 12.706) — reflects frequentist limitation, not evidence against effect
- **Bayesian priors:** μ ~ N(0, 10²); τ ~ Half-Cauchy(0, scale); sensitivity across scale = {0.5, 1.0, 2.0}
- **RVE:** CR2 bias-reduced linearization; LEAD-2 + LEAD-3 treated as single cluster (same publication)

## Citation

> Alayci B, Gerçek ÖZ. Composition of weight loss associated with incretin-based therapy: a systematic review and meta-analysis of DXA-measured lean and fat mass in randomized controlled trials. *Diabetes, Obesity and Metabolism*. 2026. Manuscript ID: DOM-26-1453-OP.

## Permanent Archive

**Zenodo DOI:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19158245.svg)](https://doi.org/10.5281/zenodo.19158245)

## License

MIT License — see [LICENSE](LICENSE)
