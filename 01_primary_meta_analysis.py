#!/usr/bin/env python3
"""
Meta-Analysis: GLP-1/Dual Agonist Effects on Lean Body Mass (kg)
DerSimonian-Laird Random-Effects Model
Birkan Alayci — 20 March 2026

CORRECTED DATA (verified from primary sources, matched to manuscript v8):
- STEP-1: %-point ETD x baseline LBM from Wilding 2021 Table S5 = -1.79 kg (SE 0.35)
- SURMOUNT-1: Direct ETD from Look 2025 Figure S2 = -4.4 kg
- LEAD-2 comp corrected to Active (Glimepiride+Metformin)
"""

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import os

# ── DATA ──────────────────────────────────────────────────────────────────────
studies = [
    {"id": "SURMOUNT-1 (Look 2025)",     "drug": "Tirzepatide pooled",  "md": -4.40,  "se": 0.61,  "n_int": 124, "n_ctrl": 36,  "comp": "Placebo", "pop": "Obesity",  "agent": "Dual"},
    {"id": "STEP-1 (Wilding 2021)",       "drug": "Semaglutide 2.4mg",   "md": -1.79,  "se": 0.35,  "n_int": 95,  "n_ctrl": 45,  "comp": "Placebo", "pop": "Obesity",  "agent": "GLP-1"},
    {"id": "SUSTAIN-8 (McCrimmon 2020)",  "drug": "Semaglutide 1.0mg",   "md": -0.78,  "se": 0.41,  "n_int": 53,  "n_ctrl": 61,  "comp": "Active",  "pop": "T2DM",     "agent": "GLP-1"},
    {"id": "LEAD-2 (Jendle 2009)",        "drug": "Liraglutide 1.8mg",   "md": -1.535, "se": 0.473, "n_int": 30,  "n_ctrl": 32,  "comp": "Active",  "pop": "T2DM",     "agent": "GLP-1"},
    {"id": "LEAD-3 (Jendle 2009)",        "drug": "Liraglutide 1.8mg",   "md": -1.508, "se": 0.540, "n_int": 17,  "n_ctrl": 12,  "comp": "Active",  "pop": "T2DM",     "agent": "GLP-1"},
]

# ── DerSimonian-Laird RANDOM EFFECTS ─────────────────────────────────────────
def dl_meta(md_arr, se_arr, study_labels=None):
    """DerSimonian-Laird random-effects meta-analysis."""
    md = np.array(md_arr, dtype=float)
    se = np.array(se_arr, dtype=float)
    k = len(md)
    wi = 1.0 / se**2  # inverse-variance weights (fixed)

    # Fixed-effect estimate
    mu_fe = np.sum(wi * md) / np.sum(wi)

    # Cochran's Q
    Q = np.sum(wi * (md - mu_fe)**2)
    df = k - 1
    p_Q = 1.0 - stats.chi2.cdf(Q, df) if df > 0 else 1.0

    # Tau² (DL estimator)
    C = np.sum(wi) - np.sum(wi**2) / np.sum(wi)
    tau2 = max(0, (Q - df) / C)
    tau = np.sqrt(tau2)

    # I²
    I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0

    # Random-effects weights
    wi_re = 1.0 / (se**2 + tau2)
    mu_re = np.sum(wi_re * md) / np.sum(wi_re)
    se_re = np.sqrt(1.0 / np.sum(wi_re))
    ci_lo = mu_re - 1.96 * se_re
    ci_hi = mu_re + 1.96 * se_re
    z = mu_re / se_re
    p_z = 2 * (1 - stats.norm.cdf(abs(z)))

    # Per-study weights (%)
    wt_pct = wi_re / np.sum(wi_re) * 100

    return {
        "mu": mu_re, "se": se_re, "ci_lo": ci_lo, "ci_hi": ci_hi,
        "z": z, "p": p_z,
        "tau2": tau2, "tau": tau, "I2": I2,
        "Q": Q, "Q_df": df, "Q_p": p_Q,
        "k": k, "weights": wt_pct,
        "md": md, "se_arr": se, "labels": study_labels
    }

def print_result(res, title=""):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")
    print(f"  k = {res['k']} studies")
    print(f"  Pooled MD  = {res['mu']:.3f} kg  [{res['ci_lo']:.3f}; {res['ci_hi']:.3f}]")
    print(f"  z = {res['z']:.3f},  p = {res['p']:.4f}")
    print(f"  tau² = {res['tau2']:.4f},  tau = {res['tau']:.3f}")
    print(f"  Q = {res['Q']:.2f}  (df={res['Q_df']}, p={res['Q_p']:.4f})")
    print(f"  I² = {res['I2']:.1f}%")
    if res['labels'] is not None:
        print(f"\n  {'Study':<35} {'MD':>7} {'SE':>7} {'Weight%':>8}")
        print(f"  {'-'*60}")
        for i, lab in enumerate(res['labels']):
            print(f"  {lab:<35} {res['md'][i]:>7.3f} {res['se_arr'][i]:>7.3f} {res['weights'][i]:>7.1f}%")
    print()


# ── EGGER'S TEST ──────────────────────────────────────────────────────────────
def egger_test(md_arr, se_arr):
    """Egger's regression test for funnel plot asymmetry."""
    md = np.array(md_arr, dtype=float)
    se = np.array(se_arr, dtype=float)
    # Standard normal deviate (z_i = md_i / se_i) regressed on 1/se_i
    # Equivalent: regress md_i/se_i on 1/se_i, test intercept
    y = md / se  # z-scores
    x = 1.0 / se  # precision
    k = len(md)
    # OLS
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    ss_xx = np.sum((x - x_mean)**2)
    ss_xy = np.sum((x - x_mean) * (y - y_mean))
    slope = ss_xy / ss_xx
    intercept = y_mean - slope * x_mean
    y_pred = intercept + slope * x
    resid = y - y_pred
    mse = np.sum(resid**2) / (k - 2)
    se_intercept = np.sqrt(mse * (1/k + x_mean**2 / ss_xx))
    t_val = intercept / se_intercept
    p_val = 2 * (1 - stats.t.cdf(abs(t_val), k - 2))
    return {"intercept": intercept, "se": se_intercept, "t": t_val, "p": p_val, "df": k-2}


# ── LEAVE-ONE-OUT ─────────────────────────────────────────────────────────────
def leave_one_out(md_arr, se_arr, labels):
    md = np.array(md_arr, dtype=float)
    se = np.array(se_arr, dtype=float)
    results = []
    for i in range(len(md)):
        idx = [j for j in range(len(md)) if j != i]
        r = dl_meta(md[idx], se[idx])
        results.append({
            "omitted": labels[i],
            "mu": r["mu"], "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"],
            "I2": r["I2"], "tau2": r["tau2"]
        })
    return results


# ── FOREST PLOT ───────────────────────────────────────────────────────────────
def forest_plot(res, title, filename, subgroup_label=None):
    k = res['k']
    labels = res['labels']
    md = res['md']
    se = res['se_arr']
    ci_lo_s = md - 1.96 * se
    ci_hi_s = md + 1.96 * se
    wt = res['weights']

    fig, ax = plt.subplots(figsize=(12, max(4, 1.2 * k + 2.5)))

    y_positions = list(range(k, 0, -1))
    y_summary = 0

    # Study-level
    for i in range(k):
        y = y_positions[i]
        marker_size = max(4, min(20, wt[i] * 0.8))
        ax.plot(md[i], y, 's', color='#2563EB', markersize=marker_size, zorder=3)
        ax.plot([ci_lo_s[i], ci_hi_s[i]], [y, y], '-', color='#2563EB', linewidth=1.5, zorder=2)
        ax.text(-10.5, y, labels[i], va='center', ha='left', fontsize=9)
        ax.text(5.0, y, f"{md[i]:.2f} [{ci_lo_s[i]:.2f}; {ci_hi_s[i]:.2f}]", va='center', ha='left', fontsize=8.5, family='monospace')
        ax.text(9.8, y, f"{wt[i]:.1f}%", va='center', ha='right', fontsize=8.5, family='monospace')

    # Summary diamond
    diamond_x = [res['ci_lo'], res['mu'], res['ci_hi'], res['mu']]
    diamond_y = [y_summary, y_summary + 0.25, y_summary, y_summary - 0.25]
    ax.fill(diamond_x, diamond_y, color='#DC2626', alpha=0.8, zorder=3)
    ax.text(-10.5, y_summary, f"RE Model (k={k})", va='center', ha='left', fontsize=9, fontweight='bold')
    ax.text(5.0, y_summary, f"{res['mu']:.2f} [{res['ci_lo']:.2f}; {res['ci_hi']:.2f}]", va='center', ha='left', fontsize=8.5, family='monospace', fontweight='bold')

    # Null line
    ax.axvline(x=0, color='black', linewidth=0.8, linestyle='--', alpha=0.5)

    # Heterogeneity annotation
    het_text = f"I² = {res['I2']:.1f}%,  τ² = {res['tau2']:.4f},  Q = {res['Q']:.2f} (p = {res['Q_p']:.3f})"
    ax.text(-10.5, y_summary - 1.0, het_text, va='center', ha='left', fontsize=8, style='italic', color='#555')

    # Axis labels
    ax.set_xlabel('Mean Difference in Lean Body Mass (kg)\n← Favours drug (more LBM loss)          Favours comparator →', fontsize=9)

    # Headers
    header_y = k + 1
    ax.text(-10.5, header_y, 'Study', va='center', ha='left', fontsize=9, fontweight='bold')
    ax.text(5.0, header_y, 'MD [95% CI]', va='center', ha='left', fontsize=9, fontweight='bold')
    ax.text(9.8, header_y, 'Weight', va='center', ha='right', fontsize=9, fontweight='bold')

    # Separator line
    ax.axhline(y=k + 0.5, color='black', linewidth=0.5)
    ax.axhline(y=0.5, color='black', linewidth=0.5)

    ax.set_xlim(-11, 10)
    ax.set_ylim(y_summary - 1.5, k + 1.5)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title(title, fontsize=11, fontweight='bold', pad=15)

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  → Saved: {filename}")


# ── FUNNEL PLOT ───────────────────────────────────────────────────────────────
def funnel_plot(res, filename):
    md = res['md']
    se = res['se_arr']
    mu = res['mu']

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(md, se, s=60, c='#2563EB', edgecolors='white', zorder=3)

    # Pseudo-CI triangle
    se_range = np.linspace(0, max(se) * 1.3, 100)
    ax.plot(mu - 1.96 * se_range, se_range, '--', color='gray', linewidth=0.8)
    ax.plot(mu + 1.96 * se_range, se_range, '--', color='gray', linewidth=0.8)
    ax.axvline(x=mu, color='#DC2626', linewidth=1, linestyle='-', alpha=0.7)

    ax.invert_yaxis()
    ax.set_xlabel('Mean Difference (kg)', fontsize=10)
    ax.set_ylabel('Standard Error', fontsize=10)
    ax.set_title('Funnel Plot', fontsize=11, fontweight='bold')

    # Label points
    for i, lab in enumerate(res['labels']):
        short = lab.split('(')[0].strip()
        ax.annotate(short, (md[i], se[i]), textcoords="offset points",
                    xytext=(8, -5), fontsize=7, color='#555')

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  → Saved: {filename}")


# ── LEAVE-ONE-OUT PLOT ────────────────────────────────────────────────────────
def loo_plot(loo_results, overall_mu, overall_ci, filename):
    k = len(loo_results)
    fig, ax = plt.subplots(figsize=(10, max(4, k * 0.8 + 2)))

    for i, r in enumerate(loo_results):
        y = k - i
        ax.plot(r['mu'], y, 'D', color='#2563EB', markersize=6, zorder=3)
        ax.plot([r['ci_lo'], r['ci_hi']], [y, y], '-', color='#2563EB', linewidth=1.5)
        short = r['omitted'].split('(')[0].strip()
        ax.text(ax.get_xlim()[0] if ax.get_xlim()[0] != 0 else -7, y,
                f"Omitting {short}", va='center', ha='left', fontsize=9)

    # Overall band
    ax.axvspan(overall_ci[0], overall_ci[1], alpha=0.1, color='#DC2626')
    ax.axvline(x=overall_mu, color='#DC2626', linewidth=1, linestyle='--')
    ax.axvline(x=0, color='black', linewidth=0.5, linestyle=':')

    ax.set_xlabel('Pooled MD (kg) when study omitted', fontsize=10)
    ax.set_yticks(range(1, k+1))
    ax.set_yticklabels([r['omitted'].split('(')[0].strip() for r in reversed(loo_results)], fontsize=9)
    ax.set_title('Leave-One-Out Sensitivity Analysis', fontsize=11, fontweight='bold')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  → Saved: {filename}")


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')
os.makedirs(outdir, exist_ok=True)

labels_all = [s['id'] for s in studies]
md_all = [s['md'] for s in studies]
se_all = [s['se'] for s in studies]

# 1. OVERALL (all 5 studies)
print("\n" + "▓"*70)
print("  GLP-1 / DUAL AGONIST — LEAN BODY MASS META-ANALYSIS")
print("  Corrected Data — 20 March 2026")
print("▓"*70)

res_all = dl_meta(md_all, se_all, labels_all)
print_result(res_all, "OVERALL ANALYSIS (k=5, all studies)")
forest_plot(res_all, "All Studies — LBM Mean Difference (kg)", f"{outdir}/forest_all.png")

# 2. PLACEBO-CONTROLLED ONLY
idx_pbo = [i for i, s in enumerate(studies) if s['comp'] == 'Placebo']
if len(idx_pbo) >= 2:
    labels_pbo = [studies[i]['id'] for i in idx_pbo]
    md_pbo = [studies[i]['md'] for i in idx_pbo]
    se_pbo = [studies[i]['se'] for i in idx_pbo]
    res_pbo = dl_meta(md_pbo, se_pbo, labels_pbo)
    print_result(res_pbo, "SUBGROUP: Placebo-controlled only")
    forest_plot(res_pbo, "Placebo-Controlled — LBM Mean Difference (kg)", f"{outdir}/forest_placebo.png")

# 3. ACTIVE COMPARATOR ONLY
idx_act = [i for i, s in enumerate(studies) if s['comp'] == 'Active']
if len(idx_act) >= 2:
    labels_act = [studies[i]['id'] for i in idx_act]
    md_act = [studies[i]['md'] for i in idx_act]
    se_act = [studies[i]['se'] for i in idx_act]
    res_act = dl_meta(md_act, se_act, labels_act)
    print_result(res_act, "SUBGROUP: Active comparator only")
    forest_plot(res_act, "Active Comparator — LBM Mean Difference (kg)", f"{outdir}/forest_active.png")

# 4. BY AGENT TYPE
idx_glp1 = [i for i, s in enumerate(studies) if s['agent'] == 'GLP-1']
idx_dual = [i for i, s in enumerate(studies) if s['agent'] == 'Dual']

if len(idx_glp1) >= 2:
    labels_g = [studies[i]['id'] for i in idx_glp1]
    md_g = [studies[i]['md'] for i in idx_glp1]
    se_g = [studies[i]['se'] for i in idx_glp1]
    res_glp1 = dl_meta(md_g, se_g, labels_g)
    print_result(res_glp1, "SUBGROUP: GLP-1 RA only (mono)")
    forest_plot(res_glp1, "GLP-1 RA Only — LBM Mean Difference (kg)", f"{outdir}/forest_glp1.png")

print(f"\n  Note: Dual agonist subgroup has k={len(idx_dual)} (SURMOUNT-1 only) — no pooling possible.\n")

# 5. BY POPULATION
idx_ob = [i for i, s in enumerate(studies) if s['pop'] == 'Obesity']
idx_t2 = [i for i, s in enumerate(studies) if s['pop'] == 'T2DM']

if len(idx_ob) >= 2:
    labels_ob = [studies[i]['id'] for i in idx_ob]
    md_ob = [studies[i]['md'] for i in idx_ob]
    se_ob = [studies[i]['se'] for i in idx_ob]
    res_ob = dl_meta(md_ob, se_ob, labels_ob)
    print_result(res_ob, "SUBGROUP: Obesity population")

if len(idx_t2) >= 2:
    labels_t2 = [studies[i]['id'] for i in idx_t2]
    md_t2 = [studies[i]['md'] for i in idx_t2]
    se_t2 = [studies[i]['se'] for i in idx_t2]
    res_t2 = dl_meta(md_t2, se_t2, labels_t2)
    print_result(res_t2, "SUBGROUP: T2DM population")

# 6. LEAVE-ONE-OUT
print(f"\n{'='*70}")
print(f"  LEAVE-ONE-OUT SENSITIVITY ANALYSIS")
print(f"{'='*70}")
loo = leave_one_out(md_all, se_all, labels_all)
print(f"\n  {'Omitted':<35} {'Pooled MD':>9} {'95% CI':>20} {'I²':>7}")
print(f"  {'-'*75}")
for r in loo:
    short = r['omitted']
    print(f"  {short:<35} {r['mu']:>9.3f} [{r['ci_lo']:.3f}; {r['ci_hi']:.3f}] {r['I2']:>6.1f}%")

loo_plot(loo, res_all['mu'], (res_all['ci_lo'], res_all['ci_hi']), f"{outdir}/loo_sensitivity.png")

# 7. FUNNEL PLOT + EGGER'S
funnel_plot(res_all, f"{outdir}/funnel_plot.png")

egger = egger_test(md_all, se_all)
print(f"\n{'='*70}")
print(f"  EGGER'S TEST FOR FUNNEL PLOT ASYMMETRY")
print(f"{'='*70}")
print(f"  Intercept = {egger['intercept']:.3f} (SE = {egger['se']:.3f})")
print(f"  t = {egger['t']:.3f},  df = {egger['df']},  p = {egger['p']:.4f}")
if egger['p'] < 0.10:
    print(f"  ⚠ Significant asymmetry at α=0.10 — possible publication bias")
else:
    print(f"  → No significant asymmetry detected (α=0.10)")
print(f"\n  NOTE: Egger's test has low power with k={res_all['k']} studies.")
print(f"        Interpret with extreme caution.\n")

# 8. SUMMARY TABLE
print(f"\n{'▓'*70}")
print(f"  SUMMARY OF ALL ANALYSES")
print(f"{'▓'*70}")
print(f"\n  {'Analysis':<30} {'k':>3} {'MD':>7} {'95% CI':>20} {'I²':>7} {'p':>8}")
print(f"  {'-'*80}")

all_results = [
    ("Overall", res_all),
    ("Placebo-controlled", res_pbo if len(idx_pbo) >= 2 else None),
    ("Active comparator", res_act if len(idx_act) >= 2 else None),
    ("GLP-1 RA only", res_glp1 if len(idx_glp1) >= 2 else None),
    ("Obesity population", res_ob if len(idx_ob) >= 2 else None),
    ("T2DM population", res_t2 if len(idx_t2) >= 2 else None),
]

for name, r in all_results:
    if r:
        print(f"  {name:<30} {r['k']:>3} {r['mu']:>7.3f} [{r['ci_lo']:.3f}; {r['ci_hi']:.3f}] {r['I2']:>6.1f}% {r['p']:>8.4f}")

print(f"\n  All analyses: DerSimonian-Laird random-effects model")
print(f"  Effect measure: Mean Difference (kg), 95% CI")
print(f"  Primary estimand: Treatment policy (ANCOVA)")
print()
