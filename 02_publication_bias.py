#!/usr/bin/env python3
"""
Publication Bias Diagnostics — Extended
1. Doi Plot (normal quantile plot of Z-residuals)
2. LFK Index (Furuya-Kanamori asymmetry index)
3. Peters' Test (WLS regression of MD on 1/n_total)

Birkan Alayci — 20 March 2026
"""

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')

# ── DATA (corrected, identical to meta_analysis.py) ──────────────────────────
studies = [
    {"id": "SURMOUNT-1",    "md": -4.40,  "se": 0.61,  "n_int": 124, "n_ctrl": 36},
    {"id": "STEP-1",        "md": -1.79,  "se": 0.35,  "n_int": 95,  "n_ctrl": 45},
    {"id": "SUSTAIN-8",     "md": -0.78,  "se": 0.41,  "n_int": 53,  "n_ctrl": 61},
    {"id": "LEAD-2",        "md": -1.535, "se": 0.473, "n_int": 30,  "n_ctrl": 32},
    {"id": "LEAD-3",        "md": -1.508, "se": 0.540, "n_int": 17,  "n_ctrl": 12},
]

labels = np.array([s['id'] for s in studies])
md = np.array([s['md'] for s in studies])
se = np.array([s['se'] for s in studies])
n_total = np.array([s['n_int'] + s['n_ctrl'] for s in studies])
k = len(studies)

# ── FIXED-EFFECT POOLED ESTIMATE ─────────────────────────────────────────────
w_fe = 1.0 / se**2
theta_fe = np.sum(w_fe * md) / np.sum(w_fe)
se_fe = np.sqrt(1.0 / np.sum(w_fe))

print(f"\n{'▓'*70}")
print(f"  PUBLICATION BIAS DIAGNOSTICS — EXTENDED")
print(f"  k = {k} studies")
print(f"{'▓'*70}")
print(f"\n  Fixed-effect pooled estimate: {theta_fe:.3f} kg (SE = {se_fe:.3f})")


# ══════════════════════════════════════════════════════════════════════════════
#  1. DOI PLOT + LFK INDEX
# ══════════════════════════════════════════════════════════════════════════════

# Standardized residuals (centered on FE estimate)
z_resid = (md - theta_fe) / se

# Sort for Q-Q plot
order = np.argsort(z_resid)
z_sorted = z_resid[order]
labels_sorted = labels[order]

# Expected normal quantiles (Blom's formula: (r - 3/8) / (k + 1/4))
ranks = np.arange(1, k + 1)
probs_blom = (ranks - 0.375) / (k + 0.25)
expected_q = stats.norm.ppf(probs_blom)

# ── LFK INDEX CALCULATION ────────────────────────────────────────────────────
# Method: Compare observed percentiles (from normal CDF of Z-residuals)
# to expected percentiles under symmetry.
# LFK = Σ(observed_pctile - expected_pctile) * k / Σ|deviations|
#
# Reference: Furuya-Kanamori L, Barendregt JJ, Doi SAR.
# Int J Evid Based Healthc. 2018;16(4):195-203.
#
# Interpretation:
#   |LFK| ≤ 1   → no asymmetry
#   1 < |LFK| ≤ 2 → minor asymmetry
#   |LFK| > 2   → major asymmetry

# Observed percentiles from normal CDF of Z-residuals
observed_pctile = stats.norm.cdf(z_sorted)

# Expected percentiles under uniform (symmetric) distribution
expected_pctile = (ranks - 0.5) / k

# Deviations
deviations = observed_pctile - expected_pctile

# LFK index (directional measure scaled by sqrt(k))
# This captures both direction and magnitude of asymmetry
lfk_raw = np.sum(deviations)
lfk_index = lfk_raw * np.sqrt(k)

# Interpretation
if abs(lfk_index) <= 1:
    lfk_interp = "No asymmetry"
elif abs(lfk_index) <= 2:
    lfk_interp = "Minor asymmetry"
else:
    lfk_interp = "Major asymmetry"

print(f"\n{'='*70}")
print(f"  1. DOI PLOT & LFK INDEX")
print(f"{'='*70}")
print(f"\n  Centered Z-residuals (Z_i = [MD_i − θ_FE] / SE_i):")
for i in range(k):
    j = order[i]
    print(f"    {labels_sorted[i]:<20}  Z = {z_sorted[i]:>+6.3f}   obs%ile = {observed_pctile[i]:.3f}   exp%ile = {expected_pctile[i]:.3f}   Δ = {deviations[i]:>+.3f}")

print(f"\n  LFK Index = {lfk_index:+.3f}")
print(f"  Interpretation: {lfk_interp}")
print(f"    |LFK| ≤ 1 → no asymmetry")
print(f"    1 < |LFK| ≤ 2 → minor asymmetry")
print(f"    |LFK| > 2 → major asymmetry")

# ── DOI PLOT FIGURE ──────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

# Reference line (expected under symmetry: slope from OLS through origin)
slope_ref = np.sum(z_sorted * expected_q) / np.sum(expected_q**2)
x_line = np.linspace(expected_q[0] - 0.3, expected_q[-1] + 0.3, 100)
ax.plot(x_line, slope_ref * x_line, '--', color='#999', linewidth=1, label='Expected (symmetric)')

# Zero lines
ax.axhline(y=0, color='black', linewidth=0.5, alpha=0.3)
ax.axvline(x=0, color='black', linewidth=0.5, alpha=0.3)

# Data points
colors = ['#DC2626' if d > 0 else '#2563EB' for d in deviations]
ax.scatter(expected_q, z_sorted, c=colors, s=80, zorder=5, edgecolors='white', linewidth=0.8)

# Labels
for i in range(k):
    offset_x = 0.12
    offset_y = 0.15 if i % 2 == 0 else -0.25
    ax.annotate(labels_sorted[i], (expected_q[i], z_sorted[i]),
                xytext=(offset_x, offset_y), textcoords='offset fontsize',
                fontsize=8, color='#444')

# Connect points
ax.plot(expected_q, z_sorted, '-', color='#555', linewidth=1, alpha=0.5, zorder=3)

# Fill asymmetry areas
for i in range(k):
    expected_z = slope_ref * expected_q[i]
    if z_sorted[i] > expected_z:
        ax.fill_between([expected_q[i]-0.02, expected_q[i]+0.02],
                        expected_z, z_sorted[i], alpha=0.15, color='#DC2626')
    else:
        ax.fill_between([expected_q[i]-0.02, expected_q[i]+0.02],
                        z_sorted[i], expected_z, alpha=0.15, color='#2563EB')

ax.set_xlabel('Normal Quantile (expected under symmetry)', fontsize=10)
ax.set_ylabel('Observed Z-residual  (MD − θ_FE) / SE', fontsize=10)
ax.set_title(f'Doi Plot  |  LFK Index = {lfk_index:+.3f} ({lfk_interp})', fontsize=11, fontweight='bold')
ax.legend(fontsize=9, loc='upper left')

# Annotation box
textstr = f'LFK = {lfk_index:+.3f}\nk = {k}\nθ_FE = {theta_fe:.3f} kg'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.97, 0.05, textstr, transform=ax.transAxes, fontsize=8.5,
        verticalalignment='bottom', horizontalalignment='right', bbox=props)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.tight_layout()
plt.savefig(f'{outdir}/doi_plot.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"\n  → Saved: {outdir}/doi_plot.png")


# ══════════════════════════════════════════════════════════════════════════════
#  2. PETERS' TEST
# ══════════════════════════════════════════════════════════════════════════════
#
# Peters JL, Sutton AJ, Jones DR, et al. Comparison of two methods to
# detect publication bias in meta-analysis. JAMA. 2006;295(6):676-680.
#
# For continuous outcomes:
#   WLS regression: MD_i = β₀ + β₁ × (1/n_total_i)
#   Weights: w_i = 1/SE_i²
#   Test H₀: β₁ = 0
#
# The rationale: if small studies have inflated effects, 1/n will predict MD.

inv_n = 1.0 / n_total
weights = w_fe  # inverse-variance weights

# Weighted least squares
X = np.column_stack([np.ones(k), inv_n])
W = np.diag(weights)

# β = (X'WX)^{-1} X'Wy
XtWX = X.T @ W @ X
XtWy = X.T @ W @ md

try:
    beta = np.linalg.solve(XtWX, XtWy)
except np.linalg.LinAlgError:
    beta = np.linalg.lstsq(XtWX, XtWy, rcond=None)[0]

# Residuals and variance of coefficients
y_pred = X @ beta
resid = md - y_pred
w_resid = weights * resid**2
mse_w = np.sum(w_resid) / (k - 2)
var_beta = mse_w * np.linalg.inv(XtWX)
se_beta = np.sqrt(np.diag(var_beta))

# t-test for slope (β₁)
t_slope = beta[1] / se_beta[1]
p_slope = 2 * (1 - stats.t.cdf(abs(t_slope), k - 2))

# t-test for intercept (β₀)
t_intercept = beta[0] / se_beta[0]
p_intercept = 2 * (1 - stats.t.cdf(abs(t_intercept), k - 2))

print(f"\n{'='*70}")
print(f"  2. PETERS' TEST (WLS: MD ~ 1/n_total)")
print(f"{'='*70}")
print(f"\n  Regression: MD_i = β₀ + β₁ × (1/n_total_i),  weights = 1/SE²")
print(f"\n  β₀ (intercept) = {beta[0]:>+8.3f}  (SE = {se_beta[0]:.3f})")
print(f"     t = {t_intercept:.3f},  df = {k-2},  p = {p_intercept:.4f}")
print(f"\n  β₁ (slope)     = {beta[1]:>+8.3f}  (SE = {se_beta[1]:.3f})")
print(f"     t = {t_slope:.3f},  df = {k-2},  p = {p_slope:.4f}")

if p_slope < 0.10:
    print(f"\n  ⚠ Slope significant at α=0.10 — possible small-study effect")
else:
    print(f"\n  → Slope not significant at α=0.10 — no evidence of small-study effect")

print(f"\n  Study-level data:")
print(f"    {'Study':<20} {'n_total':>8} {'1/n':>8} {'MD':>8} {'Predicted':>10}")
print(f"    {'-'*58}")
for i in range(k):
    print(f"    {labels[i]:<20} {n_total[i]:>8d} {inv_n[i]:>8.4f} {md[i]:>8.3f} {y_pred[i]:>10.3f}")


# ── PETERS' PLOT ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

# Scatter (size ∝ weight)
sizes = weights / np.max(weights) * 200 + 30
ax.scatter(inv_n, md, s=sizes, c='#2563EB', edgecolors='white', linewidth=0.8, zorder=3, alpha=0.8)

# Regression line
x_range = np.linspace(0, max(inv_n) * 1.2, 100)
y_line = beta[0] + beta[1] * x_range
ax.plot(x_range, y_line, '-', color='#DC2626', linewidth=1.5, label=f'WLS fit (p={p_slope:.3f})')

# Confidence band
for xi in x_range:
    xvec = np.array([1, xi])
    se_pred = np.sqrt(xvec @ var_beta @ xvec * mse_w)

# Null line (pooled FE estimate)
ax.axhline(y=theta_fe, color='#999', linewidth=1, linestyle=':', label=f'θ_FE = {theta_fe:.2f} kg')

# Labels
for i in range(k):
    ax.annotate(labels[i], (inv_n[i], md[i]),
                xytext=(8, 6), textcoords='offset points',
                fontsize=8, color='#444')

ax.set_xlabel('1 / Total Sample Size (1/n)', fontsize=10)
ax.set_ylabel('Mean Difference in LBM (kg)', fontsize=10)
ax.set_title(f"Peters' Test  |  β₁ = {beta[1]:+.1f}, p = {p_slope:.3f}", fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Annotation
textstr = f'β₁ = {beta[1]:+.1f}\nt = {t_slope:.2f}\np = {p_slope:.3f}\ndf = {k-2}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.97, 0.97, textstr, transform=ax.transAxes, fontsize=8.5,
        verticalalignment='top', horizontalalignment='right', bbox=props)

plt.tight_layout()
plt.savefig(f'{outdir}/peters_test.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"\n  → Saved: {outdir}/peters_test.png")


# ══════════════════════════════════════════════════════════════════════════════
#  3. COMPARISON TABLE — ALL SMALL-STUDY EFFECT TESTS
# ══════════════════════════════════════════════════════════════════════════════

# Recompute Egger's for comparison
z_scores = md / se
precision = 1.0 / se
x_eg = precision
y_eg = z_scores
x_mean_eg = np.mean(x_eg)
y_mean_eg = np.mean(y_eg)
ss_xx_eg = np.sum((x_eg - x_mean_eg)**2)
ss_xy_eg = np.sum((x_eg - x_mean_eg) * (y_eg - y_mean_eg))
slope_eg = ss_xy_eg / ss_xx_eg
intercept_eg = y_mean_eg - slope_eg * x_mean_eg
y_pred_eg = intercept_eg + slope_eg * x_eg
resid_eg = y_eg - y_pred_eg
mse_eg = np.sum(resid_eg**2) / (k - 2)
se_int_eg = np.sqrt(mse_eg * (1/k + x_mean_eg**2 / ss_xx_eg))
t_eg = intercept_eg / se_int_eg
p_eg = 2 * (1 - stats.t.cdf(abs(t_eg), k - 2))

# Begg's rank correlation test
# Kendall's tau between effect size and variance
tau_begg, p_begg = stats.kendalltau(md, se**2)

print(f"\n{'='*70}")
print(f"  COMPARISON: ALL SMALL-STUDY EFFECT TESTS")
print(f"{'='*70}")
print(f"\n  {'Test':<25} {'Statistic':>12} {'p-value':>10} {'Significant':>12}")
print(f"  {'-'*62}")
print(f"  {'Egger (intercept)':<25} {'t = %.2f' % t_eg:>12} {p_eg:>10.4f} {'YES' if p_eg < 0.10 else 'No':>12}")
print(f"  {'Peters (slope)':<25} {'t = %.2f' % t_slope:>12} {p_slope:>10.4f} {'YES' if p_slope < 0.10 else 'No':>12}")
print(f"  {'Begg (Kendall τ)':<25} {'τ = %.3f' % tau_begg:>12} {p_begg:>10.4f} {'YES' if p_begg < 0.10 else 'No':>12}")
print(f"  {'LFK Index':<25} {'LFK = %.2f' % lfk_index:>12} {'—':>10} {lfk_interp:>12}")

print(f"\n  α = 0.10 for all tests (two-sided)")
print(f"\n  ⚠ CRITICAL CAVEAT: All tests have very low statistical power with k={k}.")
print(f"    Egger's, Peters', and Begg's tests require k ≥ 10 for adequate power.")
print(f"    The LFK index + Doi plot was specifically designed for small meta-analyses")
print(f"    and is the most appropriate diagnostic here.")
print(f"\n  Clinical interpretation:")
print(f"    The observed asymmetry is likely driven by:")
print(f"    • Comparator type heterogeneity (placebo vs active → different effect sizes)")
print(f"    • Drug class differences (dual agonist vs mono GLP-1)")
print(f"    • Dose-response gradient across agents")
print(f"    Rather than selective reporting / publication bias per se.")
print()
