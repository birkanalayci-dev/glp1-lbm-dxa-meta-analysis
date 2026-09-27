# =============================================================================
# revision_analyses.R (v2, second revision, 22 September 2026)
# Companion to glp1_dxa_meta_v1_3.R. Run from the project folder (C:/R DXA).
# v2 change: the pooled lean-fraction rows (R5) and the optional simulation (R5b)
# were removed, because the pooled proportion was withdrawn at the second
# revision (Reviewer 3, point 6); per-trial fractions are produced by Section 10
# of glp1_dxa_meta_v1_3.R (output/tables/lean_proportion_trial_level.csv).
#
# Produces, into output/revision/:
#   R1  t2d_prediction_interval.csv      (R3 major #3)
#   R2  loo_hksj_lean.csv / loo_hksj_fat.csv   (R3 major #4; supplement LOO
#                                          tables currently show DL values)
#   R3  (consistency assertions only; the STEP-1 conversion no longer exists,
#        see 02_extraction_verification.md)
#   R4  primary_excluding_step1.csv      (R3 major #4: k = 4 sensitivity)
#
# Nothing here changes the primary analysis. If the updated search adds a
# trial, add it to data_raw.csv with analysis = "primary" and re-run BOTH
# scripts.
# =============================================================================

suppressPackageStartupMessages({
  library(metafor)
  library(dplyr)
})

dir.create("output/revision", showWarnings = FALSE, recursive = TRUE)

dat <- read.csv("data_raw.csv", stringsAsFactors = FALSE) %>%
  mutate(
    lean_vi = lean_se_kg^2,
    fat_vi  = fat_se_kg^2,
    population = case_when(
      study_id %in% c("SURMOUNT-1", "STEP-1", "S-LiTE", "BARI-OPTIMISE") ~ "obesity",
      TRUE ~ "T2D"
    )
  )
dat_primary <- dat %>% filter(analysis == "primary")
stopifnot(nrow(dat_primary) >= 4)

fit_hksj <- function(d, yi, vi) {
  rma(yi = d[[yi]], vi = d[[vi]], slab = d$study_id,
      method = "REML", test = "knha")
}

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# -----------------------------------------------------------------------------
# R1. Prediction intervals for the T2D subgroup (and overall, for reference)
#     metafor::predict() uses t with k-1 df (k minus the number of model
#     coefficients, here 1) for the PI under test = "knha"; for k = 5 the
#     multiplier is t(0.975, 4) = 2.776. (Comment corrected 22 Sept 2026; the
#     earlier comment said k-2, the Methods text has always said k-1.)
#     With k = 3 the PI has 1 df and will be very wide; that is the point.
# -----------------------------------------------------------------------------
cat("\n==== R1: prediction intervals ====\n")
pi_row <- function(m, label) {
  p <- predict(m, level = 95)
  data.frame(analysis = label, k = m$k,
             MD = p$pred, CI_lb = p$ci.lb, CI_ub = p$ci.ub,
             PI_lb = p$pi.lb, PI_ub = p$pi.ub,
             tau2 = m$tau2, I2 = m$I2)
}
d_t2d <- dat_primary %>% filter(population == "T2D")
m_lean_all <- fit_hksj(dat_primary, "lean_md_kg", "lean_vi")
m_fat_all  <- fit_hksj(dat_primary, "fat_md_kg",  "fat_vi")
m_lean_t2d <- fit_hksj(d_t2d, "lean_md_kg", "lean_vi")
m_fat_t2d  <- fit_hksj(d_t2d, "fat_md_kg",  "fat_vi")

r1 <- bind_rows(
  pi_row(m_lean_all, "Lean, overall (HKSJ-REML)"),
  pi_row(m_fat_all,  "Fat, overall (HKSJ-REML)"),
  pi_row(m_lean_t2d, "Lean, T2D subgroup (HKSJ-REML)"),
  pi_row(m_fat_t2d,  "Fat, T2D subgroup (HKSJ-REML)")
)
print(r1, digits = 3)
write.csv(r1, "output/revision/t2d_prediction_interval.csv", row.names = FALSE)

# Also: tau2 CI for the T2D subgroup, to support the sentence that I2 at k=3
# is itself imprecise (R3 major #3).
cat("\nT2D lean: confint() for tau2 / I2\n")
print(confint(m_lean_t2d))

# -----------------------------------------------------------------------------
# R2. Leave-one-out under the PRIMARY model (HKSJ-REML), both outcomes.
#     Supplement Figure S1 tables currently report DL values without saying so.
#     Replace them with these (label the model) or add these as a second table.
# -----------------------------------------------------------------------------
cat("\n==== R2: leave-one-out, HKSJ-REML ====\n")
loo_tab <- function(m, outcome) {
  l <- leave1out(m)
  data.frame(outcome = outcome, omitted = m$slab,
             MD = l$estimate, CI_lb = l$ci.lb, CI_ub = l$ci.ub,
             p_HKSJ = l$pval, tau2 = l$tau2, I2 = l$I2)
}
loo_lean <- loo_tab(m_lean_all, "Lean")
loo_fat  <- loo_tab(m_fat_all,  "Fat")
print(loo_lean, digits = 3); print(loo_fat, digits = 3)
write.csv(loo_lean, "output/revision/loo_hksj_lean.csv", row.names = FALSE)
write.csv(loo_fat,  "output/revision/loo_hksj_fat.csv",  row.names = FALSE)

# -----------------------------------------------------------------------------
# R3. STEP-1 conversion: NO LONGER APPLICABLE (3 Sept 2026)
#     Re-verification showed Wilding 2021 Table S5 reports the LBM change in kg
#     directly (ETD -3.43 [-4.74; -2.13]); the submitted value (-1.79) came from
#     treating that kg ETD as percentage points. data_raw.csv now carries the
#     kg ETD, so there is no conversion and no baseline uncertainty to propagate.
#     LEAD-2 / LEAD-3 lean values were also corrected (ETD instead of within-arm
#     change). See SNAPP/02_extraction_verification.md. The block below only
#     asserts that the corrected values are in the CSV.
# -----------------------------------------------------------------------------
cat("\n==== R3: data consistency assertions ====\n")
chk <- dat_primary %>% select(study_id, lean_md_kg, lean_se_kg, fat_md_kg, fat_se_kg)
print(chk)
stopifnot(abs(chk$lean_md_kg[chk$study_id == "STEP-1"] - (-3.43)) < 1e-6,
          abs(chk$lean_md_kg[chk$study_id == "LEAD-2"] - (-2.816)) < 1e-6,
          abs(chk$lean_md_kg[chk$study_id == "LEAD-3"] - (-0.956)) < 1e-6)
cat("Corrected lean values present.\n")

# -----------------------------------------------------------------------------
# R4. Sensitivity excluding STEP-1 (k = 4), both outcomes, HKSJ-REML.
#     Same numbers as the LOO row, but reported as its own named analysis
#     because the reviewer asked for it explicitly. Obesity subgroup becomes
#     k = 1 (SURMOUNT-1 only): report descriptively, do not pool.
# -----------------------------------------------------------------------------
cat("\n==== R4: primary excluding STEP-1 ====\n")
d_no_step1 <- dat_primary %>% filter(study_id != "STEP-1")
m_lean_ns <- fit_hksj(d_no_step1, "lean_md_kg", "lean_vi")
m_fat_ns  <- fit_hksj(d_no_step1, "fat_md_kg",  "fat_vi")
r4 <- bind_rows(
  pi_row(m_lean_all, "Lean, k=5 (as submitted)"),
  pi_row(m_lean_ns,  "Lean, k=4 (STEP-1 excluded)"),
  pi_row(m_fat_all,  "Fat, k=5 (as submitted)"),
  pi_row(m_fat_ns,   "Fat, k=4 (STEP-1 excluded)")
)
r4$p_HKSJ <- c(m_lean_all$pval, m_lean_ns$pval, m_fat_all$pval, m_fat_ns$pval)
print(r4, digits = 3)
write.csv(r4, "output/revision/primary_excluding_step1.csv", row.names = FALSE)
# DL for the same, so the response letter can say "point estimate unchanged;
# HKSJ CI widens at k=4, DL CI still excludes null" if that is what it shows:
cat("\nDL, STEP-1 excluded (lean / fat):\n")
print(rma(yi = d_no_step1$lean_md_kg, vi = d_no_step1$lean_vi, method = "DL"))
print(rma(yi = d_no_step1$fat_md_kg,  vi = d_no_step1$fat_vi,  method = "DL"))

cat("\nDone. Files in output/revision/\n")
sessionInfo()
