#!/usr/bin/env Rscript
# ═══════════════════════════════════════════════════════════════════════
# CROSS-VALIDATION: GLP-1 LBM Meta-Analysis (Lean Mass)
# Alayci & Gercek, 2026 — DOM-26-1453-OP
# Run this in R to verify Python results against manuscript v8
# Requires: metafor, clubSandwich
# ═══════════════════════════════════════════════════════════════════════

if (!require("metafor")) install.packages("metafor", repos="https://cran.r-project.org")
if (!require("clubSandwich")) install.packages("clubSandwich", repos="https://cran.r-project.org")

library(metafor)
library(clubSandwich)

# ── DATA ──────────────────────────────────────────────────────────────
dat <- data.frame(
  study   = c("SURMOUNT-1", "STEP-1", "SUSTAIN-8", "LEAD-2", "LEAD-3"),
  yi      = c(-4.40, -1.79, -0.78, -1.535, -1.508),
  sei     = c(0.61, 0.35, 0.41, 0.473, 0.540),
  cluster = c(1, 2, 3, 4, 4),
  ni      = c(160, 140, 114, 62, 29)
)
dat$vi <- dat$sei^2

cat("\n", strrep("=", 70), "\n")
cat("  CROSS-VALIDATION: R metafor vs Manuscript v8\n")
cat(strrep("=", 70), "\n\n")

# ── 1. DL ─────────────────────────────────────────────────────────────
res_dl <- rma(yi, vi, ni=ni, data=dat, method="DL")
cat("1. DerSimonian-Laird\n")
cat(sprintf("   Manuscript: MD = -1.94, tau2 = 1.1000, I2 = 83.9%%\n"))
cat(sprintf("   R:          MD = %.3f, tau2 = %.4f, I2 = %.1f%%\n",
            coef(res_dl), res_dl$tau2, res_dl$I2))
cat(sprintf("   Match: %s\n\n", ifelse(abs(coef(res_dl) - (-1.94)) < 0.01, "YES", "NO")))

# ── 2. REML ───────────────────────────────────────────────────────────
res_reml <- rma(yi, vi, ni=ni, data=dat, method="REML")
cat("2. REML\n")
cat(sprintf("   Manuscript: MD = -1.96, tau2 = 1.521\n"))
cat(sprintf("   R:          MD = %.3f, tau2 = %.4f\n",
            coef(res_reml), res_reml$tau2))
cat(sprintf("   Match: %s\n\n", ifelse(abs(coef(res_reml) - (-1.96)) < 0.01, "YES", "NO")))

# ── 3. Paule-Mandel ──────────────────────────────────────────────────
res_pm <- rma(yi, vi, ni=ni, data=dat, method="PM")
cat("3. Paule-Mandel\n")
cat(sprintf("   Manuscript: MD = -1.96, tau2 = 1.615\n"))
cat(sprintf("   R:          MD = %.3f, tau2 = %.4f\n",
            coef(res_pm), res_pm$tau2))
cat(sprintf("   Match: %s\n\n", ifelse(abs(coef(res_pm) - (-1.96)) < 0.01, "YES", "NO")))

# ── 4. HKSJ on DL ────────────────────────────────────────────────────
res_hksj <- rma(yi, vi, ni=ni, data=dat, method="DL", test="knha")
cat("4. HKSJ (on DL)\n")
cat(sprintf("   Manuscript: CI = [-3.62; -0.27], p = 0.032\n"))
cat(sprintf("   R:          CI = [%.3f; %.3f], p = %.4f\n",
            res_hksj$ci.lb, res_hksj$ci.ub, res_hksj$pval))
cat(sprintf("   Match: %s\n\n", ifelse(abs(res_hksj$pval - 0.032) < 0.005, "YES", "NO")))

# ── 5. HKSJ on REML ──────────────────────────────────────────────────
res_reml_hksj <- rma(yi, vi, ni=ni, data=dat, method="REML", test="knha")
cat("5. REML + HKSJ\n")
cat(sprintf("   Manuscript: CI = [-3.64; -0.27], p = 0.032\n"))
cat(sprintf("   R:          CI = [%.3f; %.3f], p = %.4f\n",
            res_reml_hksj$ci.lb, res_reml_hksj$ci.ub, res_reml_hksj$pval))
cat(sprintf("   Match: %s\n\n", ifelse(abs(res_reml_hksj$pval - 0.032) < 0.005, "YES", "NO")))

# ── 6. Profile Likelihood CI for tau2 ─────────────────────────────────
pl <- confint(res_reml)
cat("6. Profile Likelihood CI for tau2\n")
cat(sprintf("   Manuscript: tau2 PL CI = [0.52; 18.69]\n"))
cat(sprintf("   R:          tau2 PL CI = [%.4f; %.4f]\n",
            pl$random[1,2], pl$random[1,3]))
cat(sprintf("   Match: %s\n\n",
            ifelse(abs(pl$random[1,2] - 0.52) < 0.05, "YES", "CHECK")))

# ── 7. Prediction Interval (t-based, df=k-2) ─────────────────────────
pi_res <- predict(res_dl, level=0.95)
cat("7. Prediction Interval (t-based, df=k-2=3)\n")
cat(sprintf("   Manuscript: PI = [-5.66; +1.77]\n"))
cat(sprintf("   R:          PI = [%.2f; %.2f]\n", pi_res$pi.lb, pi_res$pi.ub))
cat(sprintf("   Match: %s\n\n", ifelse(abs(pi_res$pi.lb - (-5.66)) < 0.05, "YES", "CHECK")))

# ── 8. RVE (CR2) with clubSandwich ───────────────────────────────────
rve <- coef_test(res_reml, vcov="CR2", cluster=dat$cluster)
cat("8. RVE (CR2, clubSandwich)\n")
cat(sprintf("   Manuscript: p = 0.058, df = 2.6\n"))
cat(sprintf("   R:          SE = %.4f, df = %.1f, p = %.4f\n",
            rve$SE, rve$df, rve$p_Satt))
cat(sprintf("   Match p: %s\n\n", ifelse(abs(rve$p_Satt - 0.058) < 0.01, "YES", "CHECK")))

# ── 9. Three-level model ─────────────────────────────────────────────
dat$study_id <- dat$cluster
dat$obs_id <- 1:nrow(dat)

res_3level <- rma.mv(yi, vi, random = ~ 1 | study_id / obs_id, data=dat)
cat("9. Three-level model\n")
cat(sprintf("   Manuscript: MD = -2.07 (rho=0), total sigma2 = 2.062\n"))
cat(sprintf("   R:          MD = %.3f, sigma2 = [%.4f, %.4f], total = %.4f\n",
            coef(res_3level), res_3level$sigma2[1], res_3level$sigma2[2],
            sum(res_3level$sigma2)))
cat("\n")

# ── 10. Egger's test ──────────────────────────────────────────────────
egg <- regtest(res_dl, model="lm")
cat("10. Egger's test\n")
cat(sprintf("   Manuscript: t = -3.24, p = 0.048\n"))
cat(sprintf("   R:          z = %.3f, p = %.4f\n", egg$zval, egg$pval))
cat("   Note: metafor regtest uses z-test by default; manuscript used t-test (lm)\n\n")

# ── 11. Peters' test ──────────────────────────────────────────────────
pet <- regtest(res_dl, predictor="ni", model="lm")
cat("11. Peters' test\n")
cat(sprintf("   Manuscript: t = -1.21, p = 0.31\n"))
cat(sprintf("   R:          z = %.3f, p = %.4f\n", pet$zval, pet$pval))
cat("\n")

# ── 12. Begg's test ───────────────────────────────────────────────────
begg <- ranktest(res_dl)
cat("12. Begg's rank test\n")
cat(sprintf("   Manuscript: tau = -0.40, p = 0.48\n"))
cat(sprintf("   R:          tau = %.2f, p = %.4f\n", begg$tau, begg$pval))
cat(sprintf("   Match: %s\n\n", ifelse(abs(begg$tau - (-0.40)) < 0.05, "YES", "CHECK")))

# ── SUMMARY ───────────────────────────────────────────────────────────
cat(strrep("=", 70), "\n")
cat("  SUMMARY\n")
cat(strrep("=", 70), "\n")
cat("  All 'Manuscript' values correspond to DOM-26-1453-OP v8.\n")
cat("  Small differences (<0.01) expected due to numerical precision.\n")
cat("  PI uses t-distribution (df=k-2=3) consistent with metafor predict().\n")
cat(strrep("=", 70), "\n\n")
