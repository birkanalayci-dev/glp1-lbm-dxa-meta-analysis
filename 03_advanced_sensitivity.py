#!/usr/bin/env python3
"""
Advanced Meta-Analysis Sensitivity Models
1. HKSJ (Hartung-Knapp-Sidik-Jonkman) correction
2. Prediction intervals
3. Bayesian random-effects with weakly informative prior

Birkan Alayci — 20 March 2026
"""

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, json

outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')

# ── DATA ──────────────────────────────────────────────────────────────────────
studies = [
    {"id": "SURMOUNT-1 (Look 2025)",    "md": -4.40, "se": 0.61},
    {"id": "STEP-1 (Wilding 2021)",      "md": -1.79, "se": 0.35},
    {"id": "SUSTAIN-8 (McCrimmon 2020)", "md": -0.78, "se": 0.41},
    {"id": "LEAD-2 (Jendle 2009)",       "md": -1.535,"se": 0.473},
    {"id": "LEAD-3 (Jendle 2009)",       "md": -1.508,"se": 0.540},
]

labels = [s['id'] for s in studies]
md = np.array([s['md'] for s in studies])
se = np.array([s['se'] for s in studies])
k = len(md)

# ── DL BASELINE ───────────────────────────────────────────────────────────────
wi_fe = 1.0 / se**2
mu_fe = np.sum(wi_fe * md) / np.sum(wi_fe)
Q = np.sum(wi_fe * (md - mu_fe)**2)
df = k - 1
C = np.sum(wi_fe) - np.sum(wi_fe**2) / np.sum(wi_fe)
tau2_dl = max(0, (Q - df) / C)

wi_re = 1.0 / (se**2 + tau2_dl)
mu_dl = np.sum(wi_re * md) / np.sum(wi_re)
se_dl = np.sqrt(1.0 / np.sum(wi_re))
ci_dl = (mu_dl - 1.96 * se_dl, mu_dl + 1.96 * se_dl)
p_dl = 2 * (1 - stats.norm.cdf(abs(mu_dl / se_dl)))
I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0

print("▓" * 70)
print("  ADVANCED META-ANALYSIS SENSITIVITY MODELS")
print("▓" * 70)

print(f"\n{'='*70}")
print(f"  BASELINE: DerSimonian-Laird (standard)")
print(f"{'='*70}")
print(f"  Pooled MD = {mu_dl:.3f} kg")
print(f"  95% CI (Wald/z): [{ci_dl[0]:.3f}; {ci_dl[1]:.3f}]")
print(f"  p = {p_dl:.5f}")
print(f"  τ² = {tau2_dl:.4f}, τ = {np.sqrt(tau2_dl):.3f}")
print(f"  I² = {I2:.1f}%")


# ══════════════════════════════════════════════════════════════════════════════
#  1. HKSJ CORRECTION
# ══════════════════════════════════════════════════════════════════════════════
# The HKSJ method uses the same point estimate as DL but replaces the
# z-based CI with a t-distribution CI using a modified variance estimator:
#   q_HKSJ = (1/k) * Σ wi_re * (yi - mu_re)² / (k-1)
# This inflates the SE when studies are heterogeneous, giving wider,
# more honest confidence intervals — critical at small k.

q_hksj = np.sum(wi_re * (md - mu_dl)**2) / (k - 1)
se_hksj = se_dl * np.sqrt(q_hksj)

# t-distribution with k-1 degrees of freedom
t_crit = stats.t.ppf(0.975, df=k-1)
ci_hksj = (mu_dl - t_crit * se_hksj, mu_dl + t_crit * se_hksj)
t_stat = mu_dl / se_hksj
p_hksj = 2 * (1 - stats.t.cdf(abs(t_stat), df=k-1))

print(f"\n{'='*70}")
print(f"  1. HKSJ CORRECTION (Hartung-Knapp-Sidik-Jonkman)")
print(f"{'='*70}")
print(f"  Pooled MD = {mu_dl:.3f} kg  (same point estimate as DL)")
print(f"  q_HKSJ = {q_hksj:.4f}")
print(f"  SE_HKSJ = {se_hksj:.4f}  (vs DL SE = {se_dl:.4f})")
print(f"  95% CI (t, df={k-1}): [{ci_hksj[0]:.3f}; {ci_hksj[1]:.3f}]")
print(f"  CI width: {ci_hksj[1]-ci_hksj[0]:.3f} kg  (DL: {ci_dl[1]-ci_dl[0]:.3f} kg)")
print(f"  t = {t_stat:.3f}, p = {p_hksj:.5f}")
if ci_hksj[1] < 0:
    print(f"  → CI still excludes zero — result ROBUST under HKSJ")
else:
    print(f"  ⚠ CI crosses zero — result NOT significant under HKSJ")


# ══════════════════════════════════════════════════════════════════════════════
#  2. PREDICTION INTERVAL
# ══════════════════════════════════════════════════════════════════════════════
# "Where would the true effect of the NEXT study fall?"
# PI = mu ± t_{k-2, 0.975} * sqrt(se_mu² + tau²)
# This captures BOTH estimation uncertainty AND between-study variance.

se_pred = np.sqrt(se_dl**2 + tau2_dl)
t_crit_pi = stats.t.ppf(0.975, df=max(k-2, 1))
pi = (mu_dl - t_crit_pi * se_pred, mu_dl + t_crit_pi * se_pred)

print(f"\n{'='*70}")
print(f"  2. PREDICTION INTERVAL")
print(f"{'='*70}")
print(f"  Pooled MD = {mu_dl:.3f} kg")
print(f"  95% CI: [{ci_dl[0]:.3f}; {ci_dl[1]:.3f}]  (where is the MEAN?)")
print(f"  95% PI: [{pi[0]:.3f}; {pi[1]:.3f}]  (where would NEXT study fall?)")
print(f"  SE_prediction = {se_pred:.4f}  (vs SE_mean = {se_dl:.4f})")
if pi[1] < 0:
    print(f"  → PI excludes zero — effect direction consistent across settings")
else:
    print(f"  ⚠ PI crosses zero — in some settings, effect may be negligible or reversed")
print(f"  Clinical interpretation: While the AVERAGE effect is a {abs(mu_dl):.2f} kg")
print(f"  LBM reduction, individual study settings could see effects ranging")
print(f"  from {abs(pi[0]):.2f} kg loss to {abs(pi[1]):.2f} kg {'loss' if pi[1] < 0 else 'GAIN'}.")


# ══════════════════════════════════════════════════════════════════════════════
#  3. BAYESIAN RANDOM-EFFECTS
# ══════════════════════════════════════════════════════════════════════════════
# Grid-based numerical integration (no MCMC needed for k=5)
# Prior for mu: N(0, 10²) — weakly informative, centered at no effect
# Prior for tau: Half-Cauchy(0, 2) — recommended for small k
#   (Polson & Scott 2012, Gelman 2006)
# Likelihood: yi ~ N(mu, se_i² + tau²) for each study

np.random.seed(42)

# Grid for mu and tau
mu_grid = np.linspace(-8, 4, 500)
tau_grid = np.linspace(0.001, 5, 300)
MU, TAU = np.meshgrid(mu_grid, tau_grid)

# Prior: mu ~ N(0, 10²)
log_prior_mu = stats.norm.logpdf(MU, loc=0, scale=10)

# Prior: tau ~ Half-Cauchy(0, 2)
log_prior_tau = stats.halfcauchy.logpdf(TAU, scale=2)

# Log-likelihood
log_lik = np.zeros_like(MU)
for i in range(k):
    sigma_total = np.sqrt(se[i]**2 + TAU**2)
    log_lik += stats.norm.logpdf(md[i], loc=MU, scale=sigma_total)

# Log-posterior (unnormalized)
log_post = log_lik + log_prior_mu + log_prior_tau

# Normalize in log-space for numerical stability
log_post -= np.max(log_post)
post = np.exp(log_post)

# Normalize to proper density
dmu = mu_grid[1] - mu_grid[0]
dtau = tau_grid[1] - tau_grid[0]
post /= np.sum(post) * dmu * dtau

# Marginal posteriors
post_mu = np.sum(post, axis=0) * dtau   # marginalize over tau
post_mu /= np.sum(post_mu) * dmu

post_tau = np.sum(post, axis=1) * dmu   # marginalize over mu
post_tau /= np.sum(post_tau) * dtau

# Summary statistics
mu_post_mean = np.sum(mu_grid * post_mu) * dmu
mu_post_var = np.sum((mu_grid - mu_post_mean)**2 * post_mu) * dmu
mu_post_sd = np.sqrt(mu_post_var)

# Credible interval (highest density)
mu_cdf = np.cumsum(post_mu) * dmu
mu_ci_lo = mu_grid[np.searchsorted(mu_cdf, 0.025)]
mu_ci_hi = mu_grid[np.searchsorted(mu_cdf, 0.975)]
mu_median = mu_grid[np.searchsorted(mu_cdf, 0.50)]

# P(mu < 0)
p_negative = np.sum(post_mu[mu_grid < 0]) * dmu

# P(mu < -2)
p_gt2 = np.sum(post_mu[mu_grid < -2]) * dmu

# P(mu < -1)
p_gt1 = np.sum(post_mu[mu_grid < -1]) * dmu

# Tau posterior
tau_post_mean = np.sum(tau_grid * post_tau) * dtau
tau_cdf = np.cumsum(post_tau) * dtau
tau_ci_lo = tau_grid[np.searchsorted(tau_cdf, 0.025)]
tau_ci_hi = tau_grid[np.searchsorted(tau_cdf, 0.975)]
tau_median = tau_grid[np.searchsorted(tau_cdf, 0.50)]

print(f"\n{'='*70}")
print(f"  3. BAYESIAN RANDOM-EFFECTS")
print(f"{'='*70}")
print(f"  Priors:")
print(f"    μ ~ Normal(0, 10²)  [weakly informative, centered at no effect]")
print(f"    τ ~ Half-Cauchy(0, 2)  [recommended for small k; Polson & Scott 2012]")
print(f"  Method: Grid-based numerical integration (500×300 grid)")
print(f"\n  Posterior for μ (pooled LBM effect):")
print(f"    Mean   = {mu_post_mean:.3f} kg")
print(f"    Median = {mu_median:.3f} kg")
print(f"    SD     = {mu_post_sd:.3f} kg")
print(f"    95% CrI: [{mu_ci_lo:.3f}; {mu_ci_hi:.3f}]")
print(f"\n  Posterior for τ (between-study SD):")
print(f"    Mean   = {tau_post_mean:.3f} kg")
print(f"    Median = {tau_median:.3f} kg")
print(f"    95% CrI: [{tau_ci_lo:.3f}; {tau_ci_hi:.3f}]")
print(f"\n  Clinically relevant posterior probabilities:")
print(f"    P(LBM reduction > 0 kg)  = {p_negative:.1%}")
print(f"    P(LBM reduction > 1 kg)  = {p_gt1:.1%}")
print(f"    P(LBM reduction > 2 kg)  = {p_gt2:.1%}")
print(f"\n  Comparison with DL frequentist:")
print(f"    DL estimate:      {mu_dl:.3f} [{ci_dl[0]:.3f}; {ci_dl[1]:.3f}]")
print(f"    Bayesian mean:    {mu_post_mean:.3f} [{mu_ci_lo:.3f}; {mu_ci_hi:.3f}]")
print(f"    DL τ:             {np.sqrt(tau2_dl):.3f}")
print(f"    Bayesian τ mean:  {tau_post_mean:.3f}")


# ══════════════════════════════════════════════════════════════════════════════
#  COMPARISON FOREST PLOT
# ══════════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(13, 7))

# Study-level data
y_studies = list(range(k + 6, 6, -1))  # leave space below for summaries
for i in range(k):
    y = y_studies[i]
    ci_lo_s = md[i] - 1.96 * se[i]
    ci_hi_s = md[i] + 1.96 * se[i]
    ax.plot(md[i], y, 's', color='#2563EB', markersize=8, zorder=3)
    ax.plot([ci_lo_s, ci_hi_s], [y, y], '-', color='#2563EB', linewidth=1.5)
    ax.text(-9.5, y, labels[i], va='center', ha='left', fontsize=9)
    ax.text(5.5, y, f"{md[i]:.2f} [{ci_lo_s:.2f}; {ci_hi_s:.2f}]",
            va='center', ha='left', fontsize=8.5, family='monospace')

# Separator
ax.axhline(y=5.5, color='black', linewidth=0.5)

# Summary rows
summaries = [
    {"y": 5, "label": "DL (standard)",
     "mu": mu_dl, "ci_lo": ci_dl[0], "ci_hi": ci_dl[1],
     "color": "#DC2626", "extra": f"p={p_dl:.4f}"},
    {"y": 4, "label": "HKSJ correction",
     "mu": mu_dl, "ci_lo": ci_hksj[0], "ci_hi": ci_hksj[1],
     "color": "#9333EA", "extra": f"p={p_hksj:.4f}"},
    {"y": 3, "label": "Prediction interval",
     "mu": mu_dl, "ci_lo": pi[0], "ci_hi": pi[1],
     "color": "#EA580C", "extra": "next study range"},
    {"y": 2, "label": "Bayesian posterior",
     "mu": mu_post_mean, "ci_lo": mu_ci_lo, "ci_hi": mu_ci_hi,
     "color": "#0D9488", "extra": f"P(Δ<0)={p_negative:.0%}"},
]

for s in summaries:
    y = s["y"]
    # Diamond
    dx = [s["ci_lo"], s["mu"], s["ci_hi"], s["mu"]]
    dy = [y, y + 0.22, y, y - 0.22]
    ax.fill(dx, dy, color=s["color"], alpha=0.8, zorder=3)
    ax.text(-9.5, y, s["label"], va='center', ha='left', fontsize=9, fontweight='bold', color=s["color"])
    ax.text(5.5, y, f"{s['mu']:.2f} [{s['ci_lo']:.2f}; {s['ci_hi']:.2f}]  {s['extra']}",
            va='center', ha='left', fontsize=8.5, family='monospace', color=s["color"])

# Null line
ax.axvline(x=0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)

# Header
hdr_y = k + 7
ax.text(-9.5, hdr_y, 'Study / Model', va='center', ha='left', fontsize=9, fontweight='bold')
ax.text(5.5, hdr_y, 'Estimate [95% CI / CrI / PI]', va='center', ha='left', fontsize=9, fontweight='bold')
ax.axhline(y=hdr_y - 0.5, color='black', linewidth=0.5)

# Footer annotation
ax.text(-9.5, 0.5, f"I² = {I2:.1f}%,  τ² = {tau2_dl:.3f} (DL),  τ = {tau_post_mean:.2f} (Bayesian)",
        va='center', ha='left', fontsize=8, style='italic', color='#555')

ax.set_xlabel('Mean Difference in Lean Body Mass (kg)\n← Greater LBM loss with drug          Less LBM loss →', fontsize=9)
ax.set_xlim(-10, 12)
ax.set_ylim(0, k + 7.5)
ax.set_yticks([])
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_title('Comparison of Meta-Analytic Models (k = 5)', fontsize=11, fontweight='bold', pad=15)

plt.tight_layout()
plt.savefig(f'{outdir}/forest_model_comparison.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"\n  → Saved: {outdir}/forest_model_comparison.png")


# ══════════════════════════════════════════════════════════════════════════════
#  BAYESIAN POSTERIOR PLOTS
# ══════════════════════════════════════════════════════════════════════════════

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# μ posterior
ax1.fill_between(mu_grid, post_mu, alpha=0.3, color='#0D9488')
ax1.plot(mu_grid, post_mu, color='#0D9488', linewidth=2)
ax1.axvline(x=0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
ax1.axvline(x=mu_post_mean, color='#DC2626', linewidth=1.5, linestyle='-', alpha=0.7, label=f'Mean = {mu_post_mean:.2f}')
ax1.axvline(x=mu_dl, color='#2563EB', linewidth=1.5, linestyle=':', alpha=0.7, label=f'DL = {mu_dl:.2f}')

# Shade 95% CrI
mask_cri = (mu_grid >= mu_ci_lo) & (mu_grid <= mu_ci_hi)
ax1.fill_between(mu_grid[mask_cri], post_mu[mask_cri], alpha=0.15, color='#0D9488')

# Shade P(mu < -2)
mask_gt2 = mu_grid < -2
ax1.fill_between(mu_grid[mask_gt2], post_mu[mask_gt2], alpha=0.4, color='#DC2626',
                 label=f'P(Δ < -2 kg) = {p_gt2:.0%}')

# Annotations
ax1.annotate(f'95% CrI\n[{mu_ci_lo:.2f}; {mu_ci_hi:.2f}]',
             xy=(mu_ci_lo, 0), xytext=(mu_ci_lo - 1.5, max(post_mu) * 0.6),
             fontsize=8, color='#555',
             arrowprops=dict(arrowstyle='->', color='#999', lw=0.8))

ax1.set_xlabel('Pooled MD in LBM (kg)', fontsize=10)
ax1.set_ylabel('Posterior density', fontsize=10)
ax1.set_title('Posterior distribution of μ', fontsize=11, fontweight='bold')
ax1.legend(fontsize=8, loc='upper left')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# τ posterior
ax2.fill_between(tau_grid, post_tau, alpha=0.3, color='#9333EA')
ax2.plot(tau_grid, post_tau, color='#9333EA', linewidth=2)
ax2.axvline(x=tau_post_mean, color='#DC2626', linewidth=1.5, linestyle='-', alpha=0.7,
            label=f'Mean = {tau_post_mean:.2f}')
ax2.axvline(x=np.sqrt(tau2_dl), color='#2563EB', linewidth=1.5, linestyle=':', alpha=0.7,
            label=f'DL τ = {np.sqrt(tau2_dl):.2f}')

mask_tau_cri = (tau_grid >= tau_ci_lo) & (tau_grid <= tau_ci_hi)
ax2.fill_between(tau_grid[mask_tau_cri], post_tau[mask_tau_cri], alpha=0.15, color='#9333EA')

ax2.set_xlabel('Between-study SD (τ, kg)', fontsize=10)
ax2.set_ylabel('Posterior density', fontsize=10)
ax2.set_title('Posterior distribution of τ', fontsize=11, fontweight='bold')
ax2.legend(fontsize=8)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

plt.tight_layout()
plt.savefig(f'{outdir}/bayesian_posteriors.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"  → Saved: {outdir}/bayesian_posteriors.png")


# ══════════════════════════════════════════════════════════════════════════════
#  SUMMARY TABLE
# ══════════════════════════════════════════════════════════════════════════════

print(f"\n{'▓'*70}")
print(f"  SUMMARY: ALL MODELS COMPARED")
print(f"{'▓'*70}")
print(f"\n  {'Model':<25} {'MD':>7} {'Lower':>8} {'Upper':>8} {'Width':>7} {'p / P(Δ<0)':>12}")
print(f"  {'-'*70}")
print(f"  {'DL (standard)':<25} {mu_dl:>7.3f} {ci_dl[0]:>8.3f} {ci_dl[1]:>8.3f} {ci_dl[1]-ci_dl[0]:>7.3f} {'p=%.4f'%p_dl:>12}")
print(f"  {'HKSJ correction':<25} {mu_dl:>7.3f} {ci_hksj[0]:>8.3f} {ci_hksj[1]:>8.3f} {ci_hksj[1]-ci_hksj[0]:>7.3f} {'p=%.4f'%p_hksj:>12}")
print(f"  {'Prediction interval':<25} {mu_dl:>7.3f} {pi[0]:>8.3f} {pi[1]:>8.3f} {pi[1]-pi[0]:>7.3f} {'—':>12}")
print(f"  {'Bayesian (posterior)':<25} {mu_post_mean:>7.3f} {mu_ci_lo:>8.3f} {mu_ci_hi:>8.3f} {mu_ci_hi-mu_ci_lo:>7.3f} {'P=%.1f%%'%(p_negative*100):>12}")

print(f"\n  Key findings:")
print(f"  • HKSJ widens CI by {((ci_hksj[1]-ci_hksj[0])/(ci_dl[1]-ci_dl[0])-1)*100:.0f}% vs DL — {'still significant' if ci_hksj[1] < 0 else 'NO LONGER significant'}")
print(f"  • Prediction interval {'excludes' if pi[1] < 0 else 'INCLUDES'} zero — {'consistent direction' if pi[1] < 0 else 'effect may not generalize to all settings'}")
print(f"  • Bayesian P(any LBM reduction) = {p_negative:.1%}")
print(f"  • Bayesian P(LBM reduction > 1 kg) = {p_gt1:.1%}")
print(f"  • Bayesian P(LBM reduction > 2 kg) = {p_gt2:.1%}")
print()
