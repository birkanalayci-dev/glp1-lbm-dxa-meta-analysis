# Changelog

## Version 2 — September 2026 (revision for BMC Endocrine Disorders)

### Data (`data_raw.csv`)
- STEP-1 lean MD: −1.79 (SE 0.35) → **−3.43 kg (SE 0.666)**. Source: Wilding 2021 NEJM Supplementary Table S5, treatment-policy estimand, change in total lean body mass in kg (semaglutide −5.26, placebo −1.83, ETD −3.43 [−4.74, −2.13]). The kg ETD had been mistaken for a percentage-point value and multiplied by baseline lean mass (52.1 kg).
- LEAD-2 lean MD: −1.535 (SE 0.473) → **−2.816 kg (SE 0.594)**. Source: CTR NN2211-1572 Table 11-51 (ANCOVA, LOCF, week 26), ETD liraglutide 1.8 mg + metformin vs glimepiride + metformin (−2,816 g [−3,979; −1,652], Dunnett-adjusted CI). The within-arm LS-mean change had been entered instead of the ETD.
- LEAD-3 lean MD: −1.508 (SE 0.540) → **−0.956 kg (SE 0.888)**. Source: CTR NN2211-1573 Table 11-44 (ANCOVA, week 52), ETD liraglutide 1.8 mg vs glimepiride (−955.8 g [−2,696.9; 785.3]). Same error.
- SUSTAIN-8 lean SE: 0.41 → 0.4209 (derived from the published CI −1.61 to 0.04).
- All `source` fields now name the exact source table. Fat-mass values unchanged. The version-1 file is retained as `data_raw_v1_submitted.csv`.

### Code
- `glp1_dxa_meta_v1_2.R` replaces `glp1_dxa_meta_v1_1.R` (file corrupted on disk during the revision; v1.2 reconstructs the same pipeline): parametric simulation for the lean-fraction interval removed (ratio reported descriptively only); LFK index and Doi plots via `metasens`; RVE, three-level, profile-likelihood, small-study tests and prediction intervals written to CSV; console log fixed; data-consistency assertions added.
- `revision_analyses.R` added (subgroup prediction intervals, HKSJ leave-one-out, k = 4 sensitivity excluding STEP-1).

### Results affected
- Lean overall (HKSJ-REML): −1.96 [−3.64, −0.27] → −2.49 [−4.45, −0.53]; τ² 1.52 → 2.11; PI −6.97 to +2.00.
- Lean T2D subgroup: −1.21 [−2.33, −0.09], I² 5% → −1.53 [−4.38, +1.32], I² 73%.
- Lean obesity subgroup: −3.05 [−19.6, +13.5] → −3.95 [−10.10, +2.19], I² 13%.
- Expanded sets: k = 6 −2.55 [−4.05, −1.05]; k = 7 −2.64 [−3.88, −1.40].
- RVE lean −2.49 [−4.82, −0.16]; three-level −2.49 [−3.88, −1.09].
- Bayesian lean: brms −2.44 [−4.17, −0.59]; bayesmeta −2.47 [−4.13, −0.80].
- Lean share of pooled tissue-mass difference 27% → 32%.
- LFK (metasens): lean −0.03, fat −2.28 (version 1 reported +0.19 / +0.34 from a different implementation).
- Fat-mass analyses unchanged.

### Post-release note (7 September 2026, main branch only)
- `data_raw.csv` and `README.md`: the source note for SURMOUNT-1 now points to the Results text and Figure S2 of Look 2025 (absolute-kg treatment differences) instead of "Fig. 1" (which shows percentage changes). Values unchanged. The Zenodo v2.0 archive carries the earlier wording.

## Version 1 — March 2026
- Initial deposit (Python analysis with DerSimonian–Laird primary model and R cross-validation scripts; data as extracted in March 2026).
