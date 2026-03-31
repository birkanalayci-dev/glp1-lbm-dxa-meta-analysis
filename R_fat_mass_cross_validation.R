#!/usr/bin/env Rscript
# ============================================================
# FAT MASS META-ANALYSIS — R CROSS-VALIDATION
# Alayci & Gercek, 2026 — DOM-26-1453-OP
# Run in R 4.x with: metafor, clubSandwich
# All 'Manuscript' values correspond to v8
# ============================================================

if (!require("metafor")) install.packages("metafor", repos="https://cran.r-project.org")
if (!require("clubSandwich")) install.packages("clubSandwich", repos="https://cran.r-project.org")

library(metafor)
library(clubSandwich)

# ============================================================
# DATA
# ============================================================
dat <- data.frame(
  study   = c("SURMOUNT-1", "STEP-1", "SUSTAIN-8", "LEAD-2", "LEAD-3"),
  md      = c(-12.30, -6.99, -0.79, -3.532, -3.367),
  ci_lo   = c(-15.1, -9.79, -2.10, -5.181, -5.744),
  ci_hi   = c(-9.6, -4.19, 0.51, -1.884, -0.989),
  n_int   = c(124, 95, 53, 30, 17),
  n_ctrl  = c(36, 45, 61, 32, 12),
  comparator = c("placebo", "placebo", "active", "active", "active")
)
dat$se <- (dat$ci_hi - dat$ci_lo) / 3.92
dat$vi <- dat$se^2

cat("============================================================\n")
cat("FAT MASS DATA\n")
cat("============================================================\n")
print(dat[, c("study", "md", "se", "ci_lo", "ci_hi", "comparator")])

# ============================================================
# 1. DerSimonian-Laird
# ============================================================
cat("\n\n--- 1. DerSimonian-Laird ---\n")
res_dl <- rma(yi = md, vi = vi, data = dat, method = "DL")
print(summary(res_dl))

# ============================================================
# 2. REML
# ============================================================
cat("\n\n--- 2. REML ---\n")
res_reml <- rma(yi = md, vi = vi, data = dat, method = "REML")
print(summary(res_reml))

# ============================================================
# 3. Paule-Mandel
# ============================================================
cat("\n\n--- 3. Paule-Mandel ---\n")
res_pm <- rma(yi = md, vi = vi, data = dat, method = "PM")
print(summary(res_pm))

# ============================================================
# 4. HKSJ (DL)
# ============================================================
cat("\n\n--- 4. HKSJ with DL tau2 ---\n")
res_hksj_dl <- rma(yi = md, vi = vi, data = dat, method = "DL", test = "knha")
print(summary(res_hksj_dl))

# ============================================================
# 5. HKSJ (REML)
# ============================================================
cat("\n\n--- 5. HKSJ with REML tau2 ---\n")
res_hksj_reml <- rma(yi = md, vi = vi, data = dat, method = "REML", test = "knha")
print(summary(res_hksj_reml))

# ============================================================
# 6. Prediction Interval (t-based, df=k-2=3)
# ============================================================
cat("\n\n--- 6. Prediction Interval ---\n")
pi <- predict(res_dl, level = 0.95)
cat(sprintf("  Manuscript: PI = [-18.90; 8.36]\n"))
cat(sprintf("  R:          PI = [%.2f; %.2f]\n", pi$pi.lb, pi$pi.ub))

# ============================================================
# 7. Profile Likelihood (mu CI)
# ============================================================
cat("\n\n--- 7. Profile Likelihood (mu) ---\n")
res_ml <- rma(yi = md, vi = vi, data = dat, method = "ML")
cat(sprintf("  MD = %.2f [%.2f; %.2f] (ML Wald)\n", coef(res_ml), res_ml$ci.lb, res_ml$ci.ub))

# ============================================================
# 8. Profile Likelihood (tau2 CI)
# ============================================================
cat("\n\n--- 8. Profile Likelihood (tau2) ---\n")
pl_tau2 <- confint(res_reml)
cat(sprintf("  Manuscript: tau2 PL CI = [5.64; 161.63]\n"))
cat(sprintf("  R:          tau2 = %.3f [%.3f; %.3f]\n", res_reml$tau2, pl_tau2$random[1,2], pl_tau2$random[1,3]))
cat(sprintf("  I2   = %.1f%% [%.1f%%; %.1f%%]\n", res_reml$I2, pl_tau2$random[3,2], pl_tau2$random[3,3]))

# ============================================================
# 9. RVE (CR2 + Satterthwaite)
# ============================================================
cat("\n\n--- 9. Robust Variance Estimation (CR2) ---\n")
rve <- coef_test(res_dl, vcov = "CR2", cluster = dat$study, test = "Satterthwaite")
cat(sprintf("  Manuscript: p = 0.055, df = 4.0\n"))
cat(sprintf("  R:          MD = %.2f, SE = %.3f\n", rve$beta, rve$SE))
cat(sprintf("  t = %.4f, df = %.1f, p = %.6f\n", rve$tstat, rve$df, rve$p_Satt))
ci_rve_lo <- rve$beta - qt(0.975, rve$df) * rve$SE
ci_rve_hi <- rve$beta + qt(0.975, rve$df) * rve$SE
cat(sprintf("  95%% CI = [%.2f; %.2f]\n", ci_rve_lo, ci_rve_hi))

# ============================================================
# 10. Three-level model (rho = 0)
# ============================================================
cat("\n\n--- 10. Three-level model (rho=0) ---\n")
dat$cluster_id <- 1:nrow(dat)
res_3lvl <- rma.mv(yi = md, V = vi, random = ~ 1 | cluster_id, data = dat, method = "REML")
print(summary(res_3lvl))

# ============================================================
# 11. Leave-One-Out
# ============================================================
cat("\n\n--- 11. Leave-One-Out ---\n")
loo <- leave1out(res_dl)
print(loo)

# ============================================================
# SUBGROUP: Placebo-controlled
# ============================================================
cat("\n\n--- SUBGROUP: Placebo-controlled (k=2) ---\n")
dat_p <- dat[dat$comparator == "placebo", ]
res_p <- rma(yi = md, vi = vi, data = dat_p, method = "DL")
print(summary(res_p))

# ============================================================
# SUBGROUP: Active-comparator
# ============================================================
cat("\n\n--- SUBGROUP: Active-comparator (k=3) ---\n")
dat_a <- dat[dat$comparator == "active", ]
res_a <- rma(yi = md, vi = vi, data = dat_a, method = "DL")
print(summary(res_a))

# ============================================================
# PUBLICATION BIAS
# ============================================================
cat("\n\n--- Publication Bias ---\n")

cat("Begg's rank correlation:\n")
rt <- ranktest(res_dl)
print(rt)

cat("\nEgger's regression:\n")
reg <- regtest(res_dl, model = "lm")
print(reg)

cat("\nPeters' test:\n")
dat$n_total <- dat$n_int + dat$n_ctrl
dat$inv_n <- 1 / dat$n_total
reg_peters <- lm(md/se ~ inv_n, data = dat, weights = 1/dat$vi)
peters_summary <- summary(reg_peters)
cat(sprintf("  Peters' z = %.2f, p = %.4f\n",
    coef(peters_summary)[2,3], coef(peters_summary)[2,4]))

# ============================================================
# SUMMARY TABLE
# ============================================================
cat("\n\n============================================================\n")
cat("SUMMARY — FAT MASS CROSS-VALIDATION (Manuscript v8)\n")
cat("============================================================\n")
cat(sprintf("  DL:           MD = %.2f [%.2f; %.2f], p = %.6f, tau2 = %.3f\n",
    coef(res_dl), res_dl$ci.lb, res_dl$ci.ub, res_dl$pval, res_dl$tau2))
cat(sprintf("  REML:         MD = %.2f [%.2f; %.2f], p = %.6f, tau2 = %.3f\n",
    coef(res_reml), res_reml$ci.lb, res_reml$ci.ub, res_reml$pval, res_reml$tau2))
cat(sprintf("  PM:           MD = %.2f [%.2f; %.2f], p = %.6f, tau2 = %.3f\n",
    coef(res_pm), res_pm$ci.lb, res_pm$ci.ub, res_pm$pval, res_pm$tau2))
cat(sprintf("  HKSJ(DL):     MD = %.2f [%.2f; %.2f], p = %.6f\n",
    coef(res_hksj_dl), res_hksj_dl$ci.lb, res_hksj_dl$ci.ub, res_hksj_dl$pval))
cat(sprintf("  HKSJ(REML):   MD = %.2f [%.2f; %.2f], p = %.6f\n",
    coef(res_hksj_reml), res_hksj_reml$ci.lb, res_hksj_reml$ci.ub, res_hksj_reml$pval))
cat(sprintf("  PI:           [%.2f; %.2f]\n", pi$pi.lb, pi$pi.ub))
cat(sprintf("  RVE(CR2):     MD = %.2f [%.2f; %.2f], p = %.6f, df = %.1f\n",
    rve$beta, ci_rve_lo, ci_rve_hi, rve$p_Satt, rve$df))
cat(sprintf("  Placebo(k=2): MD = %.2f [%.2f; %.2f], p = %.6f, I2 = %.1f%%\n",
    coef(res_p), res_p$ci.lb, res_p$ci.ub, res_p$pval, res_p$I2))
cat(sprintf("  Active(k=3):  MD = %.2f [%.2f; %.2f], p = %.6f, I2 = %.1f%%\n",
    coef(res_a), res_a$ci.lb, res_a$ci.ub, res_a$pval, res_a$I2))

# ============================================================
# LEAN/FAT RATIO (using manuscript v8 lean mass pooled MDs)
# ============================================================
cat("\n  LEAN/FAT RATIO (from pooled MDs):\n")
lean_all <- -1.94; fat_all <- coef(res_dl)
cat(sprintf("    Overall: Lean %.1f%% / Fat %.1f%%\n",
    abs(lean_all)/abs(lean_all+fat_all)*100,
    abs(fat_all)/abs(lean_all+fat_all)*100))
lean_p <- -3.05; fat_p <- coef(res_p)
cat(sprintf("    Placebo: Lean %.1f%% / Fat %.1f%%\n",
    abs(lean_p)/abs(lean_p+fat_p)*100,
    abs(fat_p)/abs(lean_p+fat_p)*100))
lean_a <- -1.20; fat_a <- coef(res_a)
cat(sprintf("    Active:  Lean %.1f%% / Fat %.1f%%\n",
    abs(lean_a)/abs(lean_a+fat_a)*100,
    abs(fat_a)/abs(lean_a+fat_a)*100))

cat("\nDone. Cross-validation complete.\n")
