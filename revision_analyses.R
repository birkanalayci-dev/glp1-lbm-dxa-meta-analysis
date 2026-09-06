# =============================================================================
# revision_analyses.R  —  BMC Endocrine Disorders revision (Sept 2026)
# Companion to glp1_dxa_meta_v1_1.R. Run from the project folder (C:/R DXA).
#
# Produces, into output/revision/:
#   R1  t2d_prediction_interval.csv      (R3 major #3)
#   R2  loo_hksj_lean.csv / loo_hksj_fat.csv   (R3 major #4; supplement LOO
#                                          tables currently show DL values)
#   R3  (consistency assertions only — the STEP-1 conversion no longer exists,
#        see 02_extraction_verification.md)
#   R4  primary_excluding_step1.csv      (R3 major #4: k = 4 sensitivity)
#   R5  lean_fraction_descriptive.csv    (R3 major #9: descriptive ratio only)
#   R5b lean_fraction_simulation.csv     (OPTIONAL; off by default)
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
#     metafor::predict() uses t with k-2 df for the PI under test = "knha".
#     With k = 3 the PI has 1 df and will be very wide — that is the point.
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
# R3. STEP-1 conversion — NO LONGER APPLICABLE (3 Sept 2026)
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
#     k = 1 (SURMOUNT-1 only) — report descriptively, do not pool.
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

# -----------------------------------------------------------------------------
# R5. Lean fraction — DESCRIPTIVE ONLY (R3 major #9).
#     Ratio of pooled point estimates; no CI. State in Methods that it is a
#     ratio of trial-level pooled estimates, not a patient-level quantity, and
#     that lean and fat estimates are correlated within trials.
#     Per-trial fractions are given too so the reader sees the spread.
# -----------------------------------------------------------------------------
cat("\n==== R5: lean fraction, descriptive ====\n")
frac <- function(l, f) abs(l) / (abs(l) + abs(f)) * 100
r5 <- bind_rows(
  data.frame(level = "Pooled, overall k=5",   lean_pct = frac(m_lean_all$b[1], m_fat_all$b[1])),
  data.frame(level = "Pooled, T2D k=3",       lean_pct = frac(m_lean_t2d$b[1], m_fat_t2d$b[1])),
  dat_primary %>% transmute(level = paste0("Trial: ", study_id),
                            lean_pct = frac(lean_md_kg, fat_md_kg))
) %>% mutate(fat_per_kg_lean = (100 - lean_pct) / lean_pct)
print(r5, digits = 3)
write.csv(r5, "output/revision/lean_fraction_descriptive.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# R5b. OPTIONAL — only if you decide to KEEP an interval for the lean fraction.
#      This is what the submitted "bootstrap" actually was: a parametric
#      Monte Carlo simulation from the fitted random-effects models. If kept,
#      Methods must say exactly that, and it must be labelled post hoc.
#      Adds a within-trial correlation rho between lean and fat errors
#      (sensitivity over rho), which the submitted version ignored.
# -----------------------------------------------------------------------------
RUN_R5B <- FALSE
if (RUN_R5B) {
  set.seed(20260903)
  B <- 10000
  k <- nrow(dat_primary)
  mu_l <- m_lean_all$b[1]; t2_l <- m_lean_all$tau2
  mu_f <- m_fat_all$b[1];  t2_f <- m_fat_all$tau2
  vi_l <- dat_primary$lean_vi; vi_f <- dat_primary$fat_vi
  sim_one <- function(rho) {
    fr <- numeric(B)
    for (i in seq_len(B)) {
      z1 <- rnorm(k); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(k)  # correlated within-trial errors
      u1 <- rnorm(k); u2 <- rho * u1 + sqrt(1 - rho^2) * rnorm(k)  # correlated random effects
      yl <- mu_l + sqrt(t2_l) * u1 + sqrt(vi_l) * z1
      yf <- mu_f + sqrt(t2_f) * u2 + sqrt(vi_f) * z2
      ml <- rma(yi = yl, vi = vi_l, method = "REML")$b[1]
      mf <- rma(yi = yf, vi = vi_f, method = "REML")$b[1]
      fr[i] <- frac(ml, mf)
    }
    c(rho = rho, mean = mean(fr), lb = unname(quantile(fr, .025)), ub = unname(quantile(fr, .975)))
  }
  r5b <- as.data.frame(do.call(rbind, lapply(c(0, 0.5, 0.8), sim_one)))
  print(r5b, digits = 3)
  write.csv(r5b, "output/revision/lean_fraction_simulation.csv", row.names = FALSE)
}

cat("\nDone. Files in output/revision/\n")
sessionInfo()
