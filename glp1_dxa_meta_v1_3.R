# ===============================================================
# GLP-1 / dual GIP-GLP-1 agonists and DXA body composition
# Systematic review and meta-analysis: analysis script v1.3
# Authors: Birkan Alaycı, Öykü Zeynep Gerçek
# PROSPERO CRD420261323497
#
# v1.3 (2026-09-22), second revision for BMC Endocrine Disorders
#   (Reviewer 3, round 2: 'simpler and more conservative'):
#   - meta-regression removed (k = 5: moderator effects not estimable)
#   - Egger/Begg tests, LFK index and Doi plots removed (not
#     interpretable at k = 5; reporting bias stated as not assessable)
#   - pooled 'lean share' removed; trial-level proportions kept as a
#     descriptive CSV only
#   - Bayesian analyses reduced to a secondary computational
#     sensitivity analysis: two implementations, no 'corroboration'
#     claim; prior sensitivity now reported for BOTH outcomes with
#     scales chosen in the units of each outcome
#   - leave-one-out table flags whether each HKSJ interval includes 0
#   - forest plots: axis limits cover all intervals; heterogeneity
#     annotation no longer clipped
#   - comparator/population-specific results are the principal
#     presentation in the manuscript; the k = 5 pooled estimate is a
#     secondary descriptive summary (script order unchanged)
#
# v1.2 (2026-09-04), first revision:
#   - reconstructed after v1.1 file corruption; same structure/outputs
#   - data_raw.csv corrected (STEP-1, LEAD-2, LEAD-3 lean; SUSTAIN-8 SE)
#   - parametric simulation for the lean fraction removed (ratio is
#     descriptive only); LFK index via metasens::lfkindex
#   - RVE, three-level, profile-likelihood, small-study tests and
#     prediction intervals now written to CSV
#   - console log fixed (sink closed at the end, not via on.exit)
# Run from the project folder (C:/R DXA). All outputs -> ./output/
# ===============================================================

required_pkgs <- c("metafor", "clubSandwich", "dplyr", "readr", "tibble",
                   "brms", "posterior", "bayesmeta")
missing <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)

suppressPackageStartupMessages({
  library(metafor)
  library(clubSandwich)
  library(dplyr)
  library(readr)
  library(tibble)
  library(brms)
  library(posterior)
  library(bayesmeta)
})

n_cores <- if (requireNamespace("rstudioapi", quietly = TRUE)) 4 else 1

dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)
dir.create("output/tables", showWarnings = FALSE)
dir.create("output/bayes", showWarnings = FALSE)

set.seed(20260508)

con <- file("output/console_log.txt", open = "wt")
sink(con, split = TRUE)

cat("================================================================\n")
cat("GLP-1/Dual Agonist DXA Meta-Analysis, v1.3 (second revision)\n")
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat("brms cores:", n_cores, "\n")
cat("================================================================\n\n")

# ---- 1. DATA LOAD ----

dat_raw <- read_csv("data_raw.csv", show_col_types = FALSE)

dat <- dat_raw %>%
  mutate(
    population = case_when(
      study_id %in% c("SURMOUNT-1", "STEP-1", "S-LiTE", "BARI-OPTIMISE") ~ "obesity",
      TRUE ~ "T2D"),
    drug_class = ifelse(drug == "Tirzepatide", "dual_GIP_GLP1", "mono_GLP1"),
    comparator_type = ifelse(comparator == "Placebo", "placebo", "active"),
    lean_vi = lean_se_kg^2,
    fat_vi  = fat_se_kg^2
  )

dat_primary <- dat %>% filter(analysis == "primary")
dat_k6      <- dat %>% filter(analysis %in% c("primary", "sensitivity_A"))
dat_k7      <- dat

# Guard against the extraction errors found in Sept 2026
stopifnot(abs(dat_primary$lean_md_kg[dat_primary$study_id == "STEP-1"] - (-3.43)) < 1e-6,
          abs(dat_primary$lean_md_kg[dat_primary$study_id == "LEAD-2"] - (-2.816)) < 1e-6,
          abs(dat_primary$lean_md_kg[dat_primary$study_id == "LEAD-3"] - (-0.956)) < 1e-6)

cat("=== DATA ===\n")
cat("Primary studies (k=", nrow(dat_primary), "):\n", sep = "")
print(dat_primary %>% select(study_id, drug, comparator, duration_wks,
                             n_int, n_ctrl, lean_md_kg, lean_se_kg,
                             fat_md_kg, fat_se_kg, population, drug_class))
cat("\n")

fit_hksj <- function(d, yi, vi) {
  rma(yi = d[[yi]], vi = d[[vi]], slab = d$study_id,
      method = "REML", test = "knha")
}
pi_of <- function(m) {
  p <- predict(m, level = 95)
  c(PI_lb = unname(p$pi.lb), PI_ub = unname(p$pi.ub))
}

# ---- 2. PRIMARY (HKSJ-REML primary; DL sensitivity) ----

cat("\n================================================================\n")
cat("SECTION 2: PRIMARY ANALYSES (Frequentist)\n")
cat("================================================================\n")

m_lean_hksj <- fit_hksj(dat_primary, "lean_md_kg", "lean_vi")
m_lean_dl   <- rma(yi = lean_md_kg, vi = lean_vi, data = dat_primary,
                   slab = study_id, method = "DL")
cat("\n--- LEAN MASS: Overall (k=5) ---\n")
cat("[PRIMARY] HKSJ-REML:\n"); print(m_lean_hksj)
cat("\n[Sensitivity] DerSimonian-Laird:\n"); print(m_lean_dl)
pi_lean <- pi_of(m_lean_hksj)
cat(sprintf("\n95%% PI (lean, t with k-1 df as implemented in metafor): [%.3f, %.3f] kg\n",
            pi_lean["PI_lb"], pi_lean["PI_ub"]))

m_fat_hksj <- fit_hksj(dat_primary, "fat_md_kg", "fat_vi")
m_fat_dl   <- rma(yi = fat_md_kg, vi = fat_vi, data = dat_primary,
                  slab = study_id, method = "DL")
cat("\n--- FAT MASS: Overall (k=5) ---\n")
cat("[PRIMARY] HKSJ-REML:\n"); print(m_fat_hksj)
cat("\n[Sensitivity] DerSimonian-Laird:\n"); print(m_fat_dl)
pi_fat <- pi_of(m_fat_hksj)
cat(sprintf("\n95%% PI (fat): [%.3f, %.3f] kg\n", pi_fat["PI_lb"], pi_fat["PI_ub"]))

# ---- 3. SUBGROUP ----

cat("\n================================================================\n")
cat("SECTION 3: SUBGROUP ANALYSES\n")
cat("================================================================\n")

run_sub <- function(dat_sub, outcome, label) {
  yi_col <- if (outcome == "lean") "lean_md_kg" else "fat_md_kg"
  vi_col <- if (outcome == "lean") "lean_vi" else "fat_vi"
  cat(sprintf("\n--- %s | %s (k=%d) ---\n", outcome, label, nrow(dat_sub)))
  if (nrow(dat_sub) < 2) { cat("  k<2, descriptive only.\n"); return(NULL) }
  m <- fit_hksj(dat_sub, yi_col, vi_col)
  print(m)
  if (m$k >= 3) {
    p <- pi_of(m)
    cat(sprintf("  95%% PI: [%.3f, %.3f]\n", p["PI_lb"], p["PI_ub"]))
    cat("  tau2 / I2 confidence intervals (profile likelihood):\n")
    print(tryCatch(confint(rma(yi = dat_sub[[yi_col]], vi = dat_sub[[vi_col]], method = "REML")),
                   error = function(e) "confint failed"))
  }
  invisible(m)
}

cat("\n>>> POPULATION (pre-specified primary subgroup)\n")
m_lean_obesity <- run_sub(dat_primary %>% filter(population == "obesity"), "lean", "obesity")
m_lean_t2d     <- run_sub(dat_primary %>% filter(population == "T2D"),     "lean", "T2D")
m_fat_obesity  <- run_sub(dat_primary %>% filter(population == "obesity"), "fat",  "obesity")
m_fat_t2d      <- run_sub(dat_primary %>% filter(population == "T2D"),     "fat",  "T2D")

cat("\n>>> COMPARATOR TYPE (fully confounded with population in this set)\n")
m_lean_pbo <- run_sub(dat_primary %>% filter(comparator_type == "placebo"), "lean", "placebo-controlled")
m_lean_act <- run_sub(dat_primary %>% filter(comparator_type == "active"),  "lean", "active-comparator")
m_fat_pbo  <- run_sub(dat_primary %>% filter(comparator_type == "placebo"), "fat",  "placebo-controlled")
m_fat_act  <- run_sub(dat_primary %>% filter(comparator_type == "active"),  "fat",  "active-comparator")

cat("\n>>> DRUG CLASS\n")
m_lean_mono <- run_sub(dat_primary %>% filter(drug_class == "mono_GLP1"), "lean", "mono GLP-1")
cat("\nDual GIP/GLP-1: only SURMOUNT-1 (k=1), descriptive:\n")
print(dat_primary %>% filter(drug_class == "dual_GIP_GLP1") %>%
        select(study_id, lean_md_kg, lean_se_kg, fat_md_kg, fat_se_kg))

# ---- 4. TAU^2 ESTIMATORS + PROFILE LIKELIHOOD ----

cat("\n================================================================\n")
cat("SECTION 4: TAU^2 ESTIMATORS\n")
cat("================================================================\n")

tau_methods <- c("DL", "REML", "PM", "ML", "EB", "SJ")
tau_compare <- function(dat_sub, outcome) {
  yi_col <- if (outcome == "lean") "lean_md_kg" else "fat_md_kg"
  vi_col <- if (outcome == "lean") "lean_vi" else "fat_vi"
  do.call(rbind, lapply(tau_methods, function(meth) {
    m <- rma(yi = dat_sub[[yi_col]], vi = dat_sub[[vi_col]], method = meth)
    data.frame(method = meth, tau2 = round(m$tau2, 4), tau = round(sqrt(m$tau2), 4),
               mu = round(as.numeric(m$b), 3), ci_lb = round(m$ci.lb, 3),
               ci_ub = round(m$ci.ub, 3), p = signif(m$pval, 2))
  }))
}
tau_table_lean <- tau_compare(dat_primary, "lean")
cat("\n[Lean] (Wald CIs)\n"); print(tau_table_lean)
write_csv(tau_table_lean, "output/tables/tau2_comparison_lean.csv")
tau_table_fat <- tau_compare(dat_primary, "fat")
cat("\n[Fat] (Wald CIs)\n"); print(tau_table_fat)
write_csv(tau_table_fat, "output/tables/tau2_comparison_fat.csv")

pl_lean <- confint(rma(yi = lean_md_kg, vi = lean_vi, data = dat_primary, method = "REML"))
pl_fat  <- confint(rma(yi = fat_md_kg,  vi = fat_vi,  data = dat_primary, method = "REML"))
cat("\n[Profile likelihood: Lean]\n"); print(pl_lean)
cat("\n[Profile likelihood: Fat]\n");  print(pl_fat)
pl_tab <- rbind(
  data.frame(outcome = "Lean", quantity = rownames(pl_lean$random), pl_lean$random, row.names = NULL),
  data.frame(outcome = "Fat",  quantity = rownames(pl_fat$random),  pl_fat$random,  row.names = NULL))
write_csv(pl_tab, "output/tables/profile_likelihood_tau2.csv")

# ---- 5. RVE (CR2, LEAD-2/LEAD-3 as one cluster) ----

cat("\n================================================================\n")
cat("SECTION 5: RVE WITH CR2 (LEAD-2/3 cluster)\n")
cat("================================================================\n")

dat_rve <- dat_primary %>%
  mutate(cluster = ifelse(study_id %in% c("LEAD-2", "LEAD-3"), "Jendle2009", study_id))
rve_row <- function(fit, outcome) {
  r <- robust(fit, cluster = dat_rve$cluster, clubSandwich = TRUE)
  print(r)
  data.frame(outcome = outcome, MD = as.numeric(r$beta), SE = as.numeric(r$se),
             df_Satterthwaite = if (!is.null(r$dfs)) as.numeric(r$dfs) else NA,
             p = as.numeric(r$pval), CI_lb = as.numeric(r$ci.lb), CI_ub = as.numeric(r$ci.ub),
             clusters = 4)
}
cat("\n[Lean] RVE CR2:\n")
rve_lean <- rve_row(rma(yi = lean_md_kg, vi = lean_vi, data = dat_rve, method = "REML"), "Lean")
cat("\n[Fat] RVE CR2:\n")
rve_fat  <- rve_row(rma(yi = fat_md_kg,  vi = fat_vi,  data = dat_rve, method = "REML"), "Fat")
write_csv(rbind(rve_lean, rve_fat), "output/tables/rve_cr2.csv")

# ---- 6. THREE-LEVEL MODEL ----

cat("\n================================================================\n")
cat("SECTION 6: THREE-LEVEL MODEL\n")
cat("================================================================\n")

dat_3l <- dat_primary %>%
  mutate(cluster_id = ifelse(study_id %in% c("LEAD-2", "LEAD-3"), "Jendle2009", study_id),
         obs_id = seq_len(n()))
three_level <- function(yi, vi, outcome) {
  fit <- tryCatch(rma.mv(yi = dat_3l[[yi]], V = dat_3l[[vi]],
                         random = ~ 1 | cluster_id/obs_id, data = dat_3l, method = "REML"),
                  error = function(e) { cat("  three-level model failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)
  print(fit)
  data.frame(outcome = outcome, MD = as.numeric(fit$b), SE = as.numeric(fit$se),
             p = as.numeric(fit$pval), CI_lb = as.numeric(fit$ci.lb), CI_ub = as.numeric(fit$ci.ub),
             sigma2_cluster = fit$sigma2[1], sigma2_within = fit$sigma2[2])
}
cat("\n[Lean]\n"); tl_lean <- three_level("lean_md_kg", "lean_vi", "Lean")
cat("\n[Fat]\n");  tl_fat  <- three_level("fat_md_kg",  "fat_vi",  "Fat")
tl_tab <- rbind(tl_lean, tl_fat)
if (!is.null(tl_tab)) write_csv(tl_tab, "output/tables/three_level.csv")

# ---- 7. LEAVE-ONE-OUT (HKSJ-REML) ----

cat("\n================================================================\n")
cat("SECTION 7: LEAVE-ONE-OUT (HKSJ-REML)\n")
cat("================================================================\n")
loo_tab <- function(m, outcome) {
  l <- leave1out(m)
  data.frame(outcome = outcome, omitted = m$slab, MD = l$estimate, SE = l$se,
             CI_lb = l$ci.lb, CI_ub = l$ci.ub, p_HKSJ = l$pval, tau2 = l$tau2, I2 = l$I2,
             CI_includes_0 = (l$ci.lb < 0 & l$ci.ub > 0))
}
loo_lean <- loo_tab(m_lean_hksj, "Lean"); cat("\n[Lean]\n"); print(loo_lean)
loo_fat  <- loo_tab(m_fat_hksj,  "Fat");  cat("\n[Fat]\n");  print(loo_fat)
write_csv(loo_lean, "output/tables/loo_lean.csv")
write_csv(loo_fat,  "output/tables/loo_fat.csv")
cat(sprintf("\nLean: HKSJ interval includes 0 in %d of %d leave-one-out analyses (%s).\n",
            sum(loo_lean$CI_includes_0), nrow(loo_lean),
            paste(loo_lean$omitted[loo_lean$CI_includes_0], collapse = ", ")))
cat(sprintf("Fat:  HKSJ interval includes 0 in %d of %d leave-one-out analyses.\n",
            sum(loo_fat$CI_includes_0), nrow(loo_fat)))
cat("This instability is reported in the manuscript as the main robustness finding.\n")

# ---- 8. REPORTING BIAS ----
cat("\n================================================================\n")
cat("SECTION 8: REPORTING BIAS\n")
cat("================================================================\n")
cat("With k = 5, funnel-plot-based tests and asymmetry indices have no useful\n")
cat("diagnostic performance and are not computed. Reporting bias is reported as\n")
cat("'not formally assessable' (v1.2 computed Egger, Begg and LFK; removed in v1.3).\n")

# ---- 9. EXPANDED INCLUSION ----

cat("\n================================================================\n")
cat("SECTION 9: EXPANDED INCLUSION\n")
cat("================================================================\n")
run_exp <- function(dat_sub, label) {
  cat(sprintf("\n>>> %s (k=%d)\n", label, nrow(dat_sub)))
  m_l <- fit_hksj(dat_sub, "lean_md_kg", "lean_vi")
  m_f <- fit_hksj(dat_sub, "fat_md_kg",  "fat_vi")
  cat("\n[Lean]\n"); print(m_l); cat("\n[Fat]\n"); print(m_f)
  list(lean = m_l, fat = m_f)
}
m_k6_set <- run_exp(dat_k6, "k=6 (+ S-LiTE)")
m_k7_set <- run_exp(dat_k7, "k=7 (+ S-LiTE + BARI-OPTIMISE)")

# ---- 10. TRIAL-LEVEL LEAN PROPORTION (descriptive only; no pooled value) ----

cat("\n================================================================\n")
cat("SECTION 10: TRIAL-LEVEL LEAN PROPORTION OF THE TISSUE-MASS DIFFERENCE\n")
cat("================================================================\n")
frac <- function(l, f) abs(l) / (abs(l) + abs(f)) * 100
share_tab <- dat_primary %>%
  transmute(study_id, lean_md_kg, fat_md_kg,
            lean_pct_of_tissue_difference = round(frac(lean_md_kg, fat_md_kg), 1))
print(share_tab)
write_csv(share_tab, "output/tables/lean_proportion_trial_level.csv")
cat("Note: trial-level ratios of the two treatment differences; descriptive only.\n")
cat("A pooled proportion is NOT computed (v1.2 reported 32%; removed in v1.3 because the\n")
cat("two meta-analyses weight trials differently and the within-trial covariance is unknown).\n")

# ---- 11. FOREST PLOTS ----

cat("\n================================================================\n")
cat("SECTION 11: FOREST PLOTS\n")
cat("================================================================\n")
# Axis limits and tick positions are fixed per outcome so that the ticks are
# round numbers. Polygons whose interval exceeds the axis limits are truncated
# by metafor at the limits; the numerical interval is always given in the
# right-hand annotation column.
fmt2  <- function(x) sprintf("%.2f", unname(x))
fmt1p <- function(x) sprintf("%.1f%%", unname(x))
het_label <- function(m, prefix = "") {
  t2 <- fmt2(m$tau2); i2 <- fmt1p(m$I2)
  if (m$k >= 3) {
    p <- pi_of(m)
    pi_txt <- sprintf("; 95%% PI [%s, %s]", fmt2(p["PI_lb"]), fmt2(p["PI_ub"]))
    bquote(paste(.(prefix), tau^2 == .(t2), "; ", I^2 == .(i2), .(pi_txt)))
  } else {
    bquote(paste(.(prefix), tau^2 == .(t2), "; ", I^2 == .(i2), "; PI not estimable (k = 2)"))
  }
}
axis_lean       <- list(alim = c(-8, 3),   at = seq(-8, 2, by = 2))
axis_lean_strat <- list(alim = c(-11, 3),  at = seq(-10, 2, by = 2))  # wide enough for the k = 2 stratum interval
axis_fat        <- list(alim = c(-20, 10), at = seq(-20, 10, by = 5))  # k = 2 fat stratum interval is truncated
# Plot limits: label column = 60% of the axis width to the left, annotation column = 32% to the right,
# so that study and stratum labels never overlap the axis area (metafor draws polygons up to alim).
xlim_of <- function(ax) c(ax$alim[1] - 0.60 * diff(ax$alim), ax$alim[2] + 0.32 * diff(ax$alim))

forest_one <- function(m, outfile, title_text, xlab_text, ax, width = 2400, height = 1500) {
  png(outfile, width = width, height = height, res = 220)
  par(mar = c(6, 4, 3, 2))
  forest(m, header = c("Study", paste0(xlab_text, " [95% CI]")), xlab = xlab_text,
         refline = 0, mlab = "Pooled (REML, HKSJ)", digits = c(2, 0), alim = ax$alim, at = ax$at, xlim = xlim_of(ax))
  title(title_text)
  mtext(het_label(m), side = 1, line = 4.5, adj = 0.5, cex = 0.8)
  dev.off(); cat("  Saved:", outfile, "\n")
}
forest_one(m_lean_hksj, "output/figures/forest_lean_overall.png",
           "Lean body mass, incretin therapy vs comparator (k = 5)", "Mean difference (kg)", axis_lean)
forest_one(m_fat_hksj, "output/figures/forest_fat_overall.png",
           "Fat mass, incretin therapy vs comparator (k = 5)", "Mean difference (kg)", axis_fat)

forest_two <- function(m_a, m_b, lab_a, lab_b, outfile, xlab_text, ax) {
  png(outfile, width = 2400, height = 2000, res = 220)
  par(mfrow = c(2, 1), mar = c(5, 4, 3, 2))
  forest(m_a, header = c("Study", paste0(xlab_text, " [95% CI]")), xlab = xlab_text,
         refline = 0, mlab = "Pooled (REML, HKSJ)", digits = c(2, 0), alim = ax$alim, at = ax$at, xlim = xlim_of(ax)); title(lab_a)
  forest(m_b, header = c("Study", paste0(xlab_text, " [95% CI]")), xlab = xlab_text,
         refline = 0, mlab = "Pooled (REML, HKSJ)", digits = c(2, 0), alim = ax$alim, at = ax$at, xlim = xlim_of(ax)); title(lab_b)
  dev.off(); cat("  Saved:", outfile, "\n")
}
if (!is.null(m_lean_t2d) && !is.null(m_lean_obesity))
  forest_two(m_lean_t2d, m_lean_obesity, "A. Type 2 diabetes (k = 3)", "B. Obesity (k = 2)",
             "output/figures/forest_lean_by_population.png", "Mean difference in lean mass (kg)", axis_lean)
if (!is.null(m_fat_t2d) && !is.null(m_fat_obesity))
  forest_two(m_fat_t2d, m_fat_obesity, "A. Type 2 diabetes (k = 3)", "B. Obesity (k = 2)",
             "output/figures/forest_fat_by_population.png", "Mean difference in fat mass (kg)", axis_fat)

# Stratum summary row: metafor polygon when the interval fits the axis; otherwise a
# diamond at the point estimate with the interval drawn as a line and arrows at the
# axis limits (the k = 2 fat-mass stratum has a HKSJ interval of about -43 to +24 kg).
strat_row <- function(m, row, label, alim) {
  est <- as.numeric(m$b); lb <- m$ci.lb; ub <- m$ci.ub
  if (lb >= alim[1] && ub <= alim[2]) {
    addpoly(m, row = row, mlab = label, digits = 2)
    return(invisible(NULL))
  }
  usr <- par("usr"); aw <- 0.014 * (usr[2] - usr[1]); ah <- 0.22
  segments(max(lb, alim[1]), row, min(ub, alim[2]), row, lwd = 1.3)
  if (lb < alim[1]) polygon(c(alim[1], alim[1] + aw, alim[1] + aw), c(row, row + ah, row - ah), col = "black")
  else segments(lb, row - 0.15, lb, row + 0.15)
  if (ub > alim[2]) polygon(c(alim[2], alim[2] - aw, alim[2] - aw), c(row, row + ah, row - ah), col = "black")
  else segments(ub, row - 0.15, ub, row + 0.15)
  pw <- 0.015 * diff(alim); ph <- 0.30
  polygon(c(est - pw, est, est + pw, est), c(row, row + ph, row, row - ph), col = "black")
  text(usr[1], row, label, pos = 4, cex = 1)
  text(usr[2], row, sprintf("%.2f [%.2f, %.2f]", est, lb, ub), pos = 2, cex = 1)
}
# Fig. 3 of the manuscript (v1.3): stratified forest plot with the five-trial
# summary shown below the strata as a secondary descriptive summary.
# Row layout (ylim upper limit L = n_ob + n_t2 + 9): header text at L - 1,
# header line at L - 2, obesity label at L - 3, obesity trials below it, obesity
# polygon, T2D label, T2D trials, T2D polygon at 1.5, separator at 0, five-trial
# summary at -1.
forest_strat <- function(outcome, outfile, xlab_text, ax) {
  yi <- if (outcome == "lean") "lean_md_kg" else "fat_md_kg"
  vi <- if (outcome == "lean") "lean_vi" else "fat_vi"
  d <- dat_primary %>% arrange(desc(population == "obesity"), study_id)
  m_all <- fit_hksj(d, yi, vi)
  m_ob  <- fit_hksj(d %>% filter(population == "obesity"), yi, vi)
  m_t2d <- fit_hksj(d %>% filter(population == "T2D"), yi, vi)
  n_ob <- m_ob$k; n_t2 <- m_t2d$k
  rows <- c((n_ob + n_t2 + 5):(n_t2 + 6), (n_t2 + 2):3)   # obesity block on top, T2D below
  png(outfile, width = 2400, height = 1700, res = 220)
  par(mar = c(6, 4, 3, 2))
  forest(m_all, rows = rows, ylim = c(-2, n_ob + n_t2 + 9), alim = ax$alim, at = ax$at, xlim = xlim_of(ax),
         header = c("Study", paste0(xlab_text, " [95% CI]")), xlab = xlab_text,
         refline = 0, mlab = "All five trials (secondary summary)", digits = c(2, 0))
  strat_row(m_ob,  n_t2 + 4.5, sprintf("Obesity, placebo-controlled (k = %d)", n_ob), ax$alim)
  strat_row(m_t2d, 1.5, sprintf("Type 2 diabetes, active comparator (k = %d)", n_t2), ax$alim)
  text(par("usr")[1], n_ob + n_t2 + 6, "Obesity, placebo-controlled", pos = 4, font = 2, cex = 0.9)
  text(par("usr")[1], n_t2 + 3, "Type 2 diabetes, active comparator", pos = 4, font = 2, cex = 0.9)
  mtext(het_label(m_all, "Five-trial summary: "), side = 1, line = 4.5, adj = 0.5, cex = 0.8)
  dev.off(); cat("  Saved:", outfile, "\n")
}
forest_strat("lean", "output/figures/forest_lean_stratified.png", "Mean difference in lean body mass (kg)", axis_lean_strat)
forest_strat("fat",  "output/figures/forest_fat_stratified.png",  "Mean difference in fat mass (kg)", axis_fat)

if (!is.null(m_lean_act) && !is.null(m_lean_pbo))
  forest_two(m_lean_act, m_lean_pbo, "A. Active comparator (k = 3)", "B. Placebo (k = 2)",
             "output/figures/forest_lean_by_comparator.png", "Mean difference in lean mass (kg)", axis_lean)

# ---- 12. BAYESIAN: brms ----

cat("\n================================================================\n")
cat("SECTION 12: BAYESIAN RANDOM-EFFECTS, SECONDARY COMPUTATIONAL SENSITIVITY (brms / Stan)\n")
cat("================================================================\n")

# Priors are stated in the units of the outcome (kg). mu ~ Normal(0, 10^2): a
# weakly informative prior placing 95% mass within +/- 20 kg. tau ~ half-Cauchy(0, 1):
# prior median 1 kg for the between-trial SD; heavy tail allows tau of several kg.
# Because fat-mass heterogeneity is much larger (REML tau about 4 kg), prior
# sensitivity uses scales 0.5, 1, 2 for lean and 1, 2, 5 for fat (Section 12b).
priors <- c(prior(normal(0, 10), class = "Intercept"),
            prior(cauchy(0, 1),  class = "sd"))
fit_brms <- function(formula, seed_add = 0) {
  brm(formula, data = dat_primary, prior = priors,
      iter = 4000, warmup = 1000, chains = 4, cores = n_cores,
      seed = 20260508 + seed_add, control = list(adapt_delta = 0.99), refresh = 0)
}
summ_brms <- function(fit, label, thresholds) {
  print(summary(fit))
  mu <- as_draws_matrix(fit)[, "b_Intercept"]
  cat(sprintf("\nbrms %s posterior:\n  Mean: %.3f kg\n  95%% CrI: [%.3f, %.3f]\n",
              label, mean(mu), quantile(mu, .025), quantile(mu, .975)))
  for (th in thresholds) cat(sprintf("  P(mu < %g) = %.4f\n", th, mean(mu < th)))
  mu
}
mu_post_l <- NA; mu_post_f <- NA; bayes_lean <- NULL; bayes_fat <- NULL
cat("\n[brms Lean: fitting...]\n")
bayes_lean <- tryCatch(fit_brms(lean_md_kg | se(lean_se_kg) ~ 1 + (1 | study_id)),
                       error = function(e) { cat("brms lean failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(bayes_lean)) {
  saveRDS(bayes_lean, "output/bayes/brms_bayes_lean.rds")
  mu_post_l <- summ_brms(bayes_lean, "Lean", c(0, -1, -2))
}
cat("\n[brms Fat: fitting...]\n")
bayes_fat <- tryCatch(fit_brms(fat_md_kg | se(fat_se_kg) ~ 1 + (1 | study_id), seed_add = 1),
                      error = function(e) { cat("brms fat failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(bayes_fat)) {
  saveRDS(bayes_fat, "output/bayes/brms_bayes_fat.rds")
  mu_post_f <- summ_brms(bayes_fat, "Fat", c(0, -3, -5))
}

cat("\n[12b. brms prior sensitivity, both outcomes]\n")
prior_sens <- function(outcome, scales, seed_add) {
  yi <- if (outcome == "Lean") "lean_md_kg" else "fat_md_kg"
  se <- if (outcome == "Lean") "lean_se_kg" else "fat_se_kg"
  f  <- as.formula(sprintf("%s | se(%s) ~ 1 + (1 | study_id)", yi, se))
  out <- data.frame()
  for (sc in scales) {
    cat(sprintf("  %s, half-Cauchy scale = %.1f ...\n", outcome, sc))
    pr_sens <- c(prior(normal(0, 10), class = "Intercept"),
                 prior_string(sprintf("cauchy(0, %s)", sc), class = "sd"))
    fit_s <- tryCatch(brm(f, data = dat_primary, prior = pr_sens, iter = 4000, warmup = 1000,
                          chains = 4, cores = n_cores, seed = 20260508 + seed_add,
                          control = list(adapt_delta = 0.99), refresh = 0),
                      error = function(e) NULL)
    if (!is.null(fit_s)) {
      mu_s  <- as_draws_matrix(fit_s)[, "b_Intercept"]
      tau_s <- as_draws_matrix(fit_s)[, "sd_study_id__Intercept"]
      out <- rbind(out, data.frame(
        outcome = outcome, tau_prior_scale = sc,
        mu_mean = round(mean(mu_s), 3),
        crI_lb = round(unname(quantile(mu_s, .025)), 3), crI_ub = round(unname(quantile(mu_s, .975)), 3),
        p_mu_lt_0 = round(mean(mu_s < 0), 4),
        tau_post_median = round(median(tau_s), 3)))
    }
  }
  out
}
prior_sens_results <- rbind(prior_sens("Lean", c(0.5, 1.0, 2.0), 0),
                            prior_sens("Fat",  c(1.0, 2.0, 5.0), 1))
if (nrow(prior_sens_results) > 0) {
  cat("\nPrior sensitivity table (both outcomes):\n"); print(prior_sens_results)
  write_csv(prior_sens_results, "output/tables/brms_prior_sensitivity.csv")
}

# ---- 13. BAYESIAN: bayesmeta (second implementation, same model) ----

cat("\n================================================================\n")
cat("SECTION 13: BAYESIAN, SECOND IMPLEMENTATION (bayesmeta; same data, priors and model)\n")
cat("================================================================\n")
fit_bm <- function(y, s, label) {
  bm <- tryCatch(bayesmeta(y = y, sigma = s, labels = dat_primary$study_id,
                           mu.prior.mean = 0, mu.prior.sd = 10,
                           tau.prior = function(t) dhalfcauchy(t, scale = 1)),
                 error = function(e) { cat("bayesmeta failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(bm)) return(NULL)
  print(bm)
  cat(sprintf("\nbayesmeta %s:\n  Posterior mean (mu): %.3f kg\n  95%% CrI: [%.3f, %.3f]\n  Posterior tau mean: %.3f\n  P(mu < 0) = %.4f\n",
              label, bm$summary["mean", "mu"], bm$summary["95% lower", "mu"],
              bm$summary["95% upper", "mu"], bm$summary["mean", "tau"], bm$pposterior(mu = 0)))
  saveRDS(bm, sprintf("output/bayes/bayesmeta_%s.rds", tolower(label)))
  png(sprintf("output/figures/bayesmeta_%s.png", tolower(label)), width = 2400, height = 1800, res = 220)
  tryCatch(plot(bm, which = 3, main = paste("bayesmeta:", label)),
           error = function(e) tryCatch(plot(bm), error = function(e2) plot.new()))
  dev.off()
  bm
}
cat("\n[bayesmeta Lean]\n"); bm_lean <- fit_bm(dat_primary$lean_md_kg, dat_primary$lean_se_kg, "Lean")
cat("\n[bayesmeta Fat]\n");  bm_fat  <- fit_bm(dat_primary$fat_md_kg,  dat_primary$fat_se_kg,  "Fat")

if (!is.null(bayes_lean) && !is.null(bm_lean) && !is.null(bayes_fat) && !is.null(bm_fat)) {
  cv_table <- data.frame(
    Outcome = c("Lean", "Fat"),
    brms_mean = round(c(mean(mu_post_l), mean(mu_post_f)), 3),
    brms_CrI = c(sprintf("[%.3f, %.3f]", quantile(mu_post_l, .025), quantile(mu_post_l, .975)),
                 sprintf("[%.3f, %.3f]", quantile(mu_post_f, .025), quantile(mu_post_f, .975))),
    bayesmeta_mean = round(c(bm_lean$summary["mean", "mu"], bm_fat$summary["mean", "mu"]), 3),
    bayesmeta_CrI = c(sprintf("[%.3f, %.3f]", bm_lean$summary["95% lower", "mu"], bm_lean$summary["95% upper", "mu"]),
                      sprintf("[%.3f, %.3f]", bm_fat$summary["95% lower", "mu"],  bm_fat$summary["95% upper", "mu"])))
  cat("\n[Two implementations of the same model; agreement is a computational check,\n")
  cat(" not independent evidence]\n"); print(cv_table)
  write_csv(cv_table, "output/tables/bayesian_two_implementations.csv")
}

# ---- 14. KEY RESULTS SUMMARY ----

cat("\n================================================================\n")
cat("SECTION 14: KEY RESULTS SUMMARY\n")
cat("================================================================\n")
extract <- function(m, outcome, group) {
  if (is.null(m)) return(NULL)
  pi <- if (m$k >= 3) pi_of(m) else c(PI_lb = NA, PI_ub = NA)
  data.frame(Outcome = outcome, Subgroup = group, k = m$k,
             MD_kg = round(as.numeric(m$b), 3), CI_lb = round(m$ci.lb, 3), CI_ub = round(m$ci.ub, 3),
             p = round(m$pval, 4), I2 = round(m$I2, 1), tau2 = round(m$tau2, 3),
             PI_lb = round(unname(pi["PI_lb"]), 3), PI_ub = round(unname(pi["PI_ub"]), 3))
}
bayes_row <- function(outcome, label, mu_vec, tau2 = NA) {
  data.frame(Outcome = outcome, Subgroup = label, k = 5, MD_kg = round(mean(mu_vec), 3),
             CI_lb = round(unname(quantile(mu_vec, .025)), 3), CI_ub = round(unname(quantile(mu_vec, .975)), 3),
             p = round(mean(mu_vec < 0), 4), I2 = NA, tau2 = tau2, PI_lb = NA, PI_ub = NA)
}
key_rows <- list(
  extract(m_lean_hksj,    "Lean", "Overall (HKSJ-REML primary)"),
  extract(m_lean_dl,      "Lean", "Overall (DL sensitivity)"),
  extract(m_lean_obesity, "Lean", "Obesity"),
  extract(m_lean_t2d,     "Lean", "T2D"),
  extract(m_lean_mono,    "Lean", "GLP-1 mono-agonist trials"),
  extract(m_fat_hksj,     "Fat",  "Overall (HKSJ-REML primary)"),
  extract(m_fat_dl,       "Fat",  "Overall (DL sensitivity)"),
  extract(m_fat_obesity,  "Fat",  "Obesity"),
  extract(m_fat_t2d,      "Fat",  "T2D"),
  extract(m_k6_set$lean,  "Lean", "Expanded k=6 (+S-LiTE)"),
  extract(m_k7_set$lean,  "Lean", "Expanded k=7 (+S-LiTE+BARI)"),
  extract(m_k6_set$fat,   "Fat",  "Expanded k=6 (+S-LiTE)"),
  extract(m_k7_set$fat,   "Fat",  "Expanded k=7 (+S-LiTE+BARI)")
)
if (!is.null(bayes_lean)) key_rows[[length(key_rows) + 1]] <- bayes_row("Lean", "Bayesian brms HC(0,1)", mu_post_l)
if (!is.null(bayes_fat))  key_rows[[length(key_rows) + 1]] <- bayes_row("Fat",  "Bayesian brms HC(0,1)", mu_post_f)
if (!is.null(bm_lean)) key_rows[[length(key_rows) + 1]] <- data.frame(
  Outcome = "Lean", Subgroup = "Bayesian bayesmeta HC(0,1)", k = 5,
  MD_kg = round(bm_lean$summary["mean", "mu"], 3), CI_lb = round(bm_lean$summary["95% lower", "mu"], 3),
  CI_ub = round(bm_lean$summary["95% upper", "mu"], 3), p = round(bm_lean$pposterior(mu = 0), 4),
  I2 = NA, tau2 = round(bm_lean$summary["mean", "tau"]^2, 3), PI_lb = NA, PI_ub = NA)
if (!is.null(bm_fat)) key_rows[[length(key_rows) + 1]] <- data.frame(
  Outcome = "Fat", Subgroup = "Bayesian bayesmeta HC(0,1)", k = 5,
  MD_kg = round(bm_fat$summary["mean", "mu"], 3), CI_lb = round(bm_fat$summary["95% lower", "mu"], 3),
  CI_ub = round(bm_fat$summary["95% upper", "mu"], 3), p = round(bm_fat$pposterior(mu = 0), 4),
  I2 = NA, tau2 = round(bm_fat$summary["mean", "tau"]^2, 3), PI_lb = NA, PI_ub = NA)
key_results <- do.call(rbind, key_rows)
rownames(key_results) <- NULL
cat("\n"); print(key_results)
write_csv(key_results, "output/tables/key_results.csv")

# ---- 15. SESSION INFO ----

cat("\n================================================================\n")
cat("SECTION 15: SESSION INFO\n")
cat("================================================================\n")
print(sessionInfo())
cat("\nALL ANALYSES COMPLETE. Output: ./output/\n")

sink(); close(con)
