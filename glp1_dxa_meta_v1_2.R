# ===============================================================
# GLP-1 / dual GIP-GLP-1 agonists and DXA body composition
# Systematic review and meta-analysis — analysis script v1.2
# Authors: Birkan Alaycı, Öykü Zeynep Gerçek
# PROSPERO CRD420261323497
#
# v1.2 (2026-09-04), revision for BMC Endocrine Disorders:
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
                   "brms", "posterior", "bayesmeta", "metasens")
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
cat("GLP-1/Dual Agonist DXA Meta-Analysis — v1.2 (revision)\n")
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
cat("SECTION 2 — PRIMARY ANALYSES (Frequentist)\n")
cat("================================================================\n")

m_lean_hksj <- fit_hksj(dat_primary, "lean_md_kg", "lean_vi")
m_lean_dl   <- rma(yi = lean_md_kg, vi = lean_vi, data = dat_primary,
                   slab = study_id, method = "DL")
cat("\n--- LEAN MASS — Overall (k=5) ---\n")
cat("[PRIMARY] HKSJ-REML:\n"); print(m_lean_hksj)
cat("\n[Sensitivity] DerSimonian-Laird:\n"); print(m_lean_dl)
pi_lean <- pi_of(m_lean_hksj)
cat(sprintf("\n95%% PI (lean, t with k-1 df as implemented in metafor): [%.3f, %.3f] kg\n",
            pi_lean["PI_lb"], pi_lean["PI_ub"]))

m_fat_hksj <- fit_hksj(dat_primary, "fat_md_kg", "fat_vi")
m_fat_dl   <- rma(yi = fat_md_kg, vi = fat_vi, data = dat_primary,
                  slab = study_id, method = "DL")
cat("\n--- FAT MASS — Overall (k=5) ---\n")
cat("[PRIMARY] HKSJ-REML:\n"); print(m_fat_hksj)
cat("\n[Sensitivity] DerSimonian-Laird:\n"); print(m_fat_dl)
pi_fat <- pi_of(m_fat_hksj)
cat(sprintf("\n95%% PI (fat): [%.3f, %.3f] kg\n", pi_fat["PI_lb"], pi_fat["PI_ub"]))

# ---- 3. SUBGROUP ----

cat("\n================================================================\n")
cat("SECTION 3 — SUBGROUP ANALYSES\n")
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

# ---- 4. META-REGRESSION (exploratory, k=5) ----

cat("\n================================================================\n")
cat("SECTION 4 — META-REGRESSION (Exploratory, k=5)\n")
cat("================================================================\n")
cat("\n[Lean] ~ duration_wks\n")
print(rma(yi = lean_md_kg, vi = lean_vi, mods = ~ duration_wks, data = dat_primary, method = "REML", test = "knha"))
cat("\n[Lean] ~ drug_class\n")
print(rma(yi = lean_md_kg, vi = lean_vi, mods = ~ drug_class, data = dat_primary, method = "REML", test = "knha"))
cat("\n[Lean] ~ population\n")
print(rma(yi = lean_md_kg, vi = lean_vi, mods = ~ population, data = dat_primary, method = "REML", test = "knha"))
cat("\n[Fat] ~ duration_wks\n")
print(rma(yi = fat_md_kg, vi = fat_vi, mods = ~ duration_wks, data = dat_primary, method = "REML", test = "knha"))

# ---- 5. TAU^2 ESTIMATORS + PROFILE LIKELIHOOD ----

cat("\n================================================================\n")
cat("SECTION 5 — TAU^2 ESTIMATORS\n")
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
cat("\n[Profile likelihood — Lean]\n"); print(pl_lean)
cat("\n[Profile likelihood — Fat]\n");  print(pl_fat)
pl_tab <- rbind(
  data.frame(outcome = "Lean", quantity = rownames(pl_lean$random), pl_lean$random, row.names = NULL),
  data.frame(outcome = "Fat",  quantity = rownames(pl_fat$random),  pl_fat$random,  row.names = NULL))
write_csv(pl_tab, "output/tables/profile_likelihood_tau2.csv")

# ---- 6. RVE (CR2, LEAD-2/LEAD-3 as one cluster) ----

cat("\n================================================================\n")
cat("SECTION 6 — RVE WITH CR2 (LEAD-2/3 cluster)\n")
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

# ---- 7. THREE-LEVEL MODEL ----

cat("\n================================================================\n")
cat("SECTION 7 — THREE-LEVEL MODEL\n")
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

# ---- 8. LEAVE-ONE-OUT (HKSJ-REML) ----

cat("\n================================================================\n")
cat("SECTION 8 — LEAVE-ONE-OUT (HKSJ-REML)\n")
cat("================================================================\n")
loo_tab <- function(m, outcome) {
  l <- leave1out(m)
  data.frame(outcome = outcome, omitted = m$slab, MD = l$estimate, SE = l$se,
             CI_lb = l$ci.lb, CI_ub = l$ci.ub, p_HKSJ = l$pval, tau2 = l$tau2, I2 = l$I2)
}
loo_lean <- loo_tab(m_lean_hksj, "Lean"); cat("\n[Lean]\n"); print(loo_lean)
loo_fat  <- loo_tab(m_fat_hksj,  "Fat");  cat("\n[Fat]\n");  print(loo_fat)
write_csv(loo_lean, "output/tables/loo_lean.csv")
write_csv(loo_fat,  "output/tables/loo_fat.csv")

# ---- 9. SMALL-STUDY EFFECTS (k=5: report as not assessable) ----

cat("\n================================================================\n")
cat("SECTION 9 — SMALL-STUDY EFFECTS (k=5, underpowered; descriptive only)\n")
cat("================================================================\n")

ssb <- function(m_dl, yi, vi, outcome) {
  eg <- regtest(m_dl, model = "lm"); bg <- ranktest(m_dl)
  cat(sprintf("\n[%s] Egger:\n", outcome)); print(eg)
  cat(sprintf("\n[%s] Begg:\n", outcome));  print(bg)
  lfk <- NA; lfk_txt <- "metasens not available"
  if (requireNamespace("metasens", quietly = TRUE)) {
    res <- tryCatch(metasens::lfkindex(TE = dat_primary[[yi]], seTE = sqrt(dat_primary[[vi]])),
                    error = function(e) NULL)
    if (!is.null(res)) { lfk <- res$lfkindex; lfk_txt <- res$interpretation; print(res) }
    png(sprintf("output/figures/doi_%s.png", tolower(outcome)), width = 1600, height = 1400, res = 200)
    tryCatch(metasens::doiplot(TE = dat_primary[[yi]], seTE = sqrt(dat_primary[[vi]]),
                               main = paste("Doi plot —", outcome)),
             error = function(e) plot.new())
    dev.off()
  }
  data.frame(outcome = outcome, egger_t = as.numeric(eg$zval), egger_p = eg$pval,
             begg_tau = as.numeric(bg$tau), begg_p = bg$pval,
             LFK = as.numeric(lfk), LFK_interpretation = lfk_txt)
}
ssb_tab <- rbind(ssb(m_lean_dl, "lean_md_kg", "lean_vi", "Lean"),
                 ssb(m_fat_dl,  "fat_md_kg",  "fat_vi",  "Fat"))
print(ssb_tab)
write_csv(ssb_tab, "output/tables/small_study_effects.csv")

# ---- 10. EXPANDED INCLUSION ----

cat("\n================================================================\n")
cat("SECTION 10 — EXPANDED INCLUSION\n")
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

# ---- 11. LEAN SHARE OF POOLED LOSS (descriptive only) ----

cat("\n================================================================\n")
cat("SECTION 11 — LEAN SHARE OF POOLED LOSS (descriptive; no interval)\n")
cat("================================================================\n")
frac <- function(l, f) abs(l) / (abs(l) + abs(f)) * 100
share_tab <- bind_rows(
  data.frame(level = "Pooled overall (k=5)", lean_pct = frac(m_lean_hksj$b[1], m_fat_hksj$b[1])),
  data.frame(level = "Pooled T2D (k=3)",     lean_pct = frac(m_lean_t2d$b[1],  m_fat_t2d$b[1])),
  data.frame(level = "Pooled obesity (k=2)", lean_pct = frac(m_lean_obesity$b[1], m_fat_obesity$b[1])),
  dat_primary %>% transmute(level = paste0("Trial: ", study_id), lean_pct = frac(lean_md_kg, fat_md_kg))
) %>% mutate(fat_per_kg_lean = (100 - lean_pct) / lean_pct)
print(share_tab)
write_csv(share_tab, "output/tables/lean_share_descriptive.csv")
cat("Note: ratios of trial-level pooled point estimates; not patient-level; no interval reported.\n")

# ---- 12. FOREST PLOTS ----

cat("\n================================================================\n")
cat("SECTION 12 — FOREST PLOTS\n")
cat("================================================================\n")
forest_one <- function(m, outfile, title_text, xlab_text, width = 2400, height = 1400) {
  png(outfile, width = width, height = height, res = 220)
  forest(m, header = c("Study", paste0(xlab_text, " [95% CI]")), xlab = xlab_text,
         refline = 0, mlab = "Pooled (REML, HKSJ)", digits = 2)
  title(title_text)
  mtext(sprintf("tau^2 = %.2f; I^2 = %.1f%%; 95%% PI [%.2f, %.2f]",
                m$tau2, m$I2, pi_of(m)["PI_lb"], pi_of(m)["PI_ub"]),
        side = 1, line = 3.2, adj = 0, cex = 0.8)
  dev.off(); cat("  Saved:", outfile, "\n")
}
forest_one(m_lean_hksj, "output/figures/forest_lean_overall.png",
           "Lean body mass, incretin therapy vs comparator (k = 5)", "Mean difference (kg)")
forest_one(m_fat_hksj, "output/figures/forest_fat_overall.png",
           "Fat mass, incretin therapy vs comparator (k = 5)", "Mean difference (kg)")

forest_two <- function(m_a, m_b, lab_a, lab_b, outfile, xlab_text) {
  png(outfile, width = 2400, height = 1800, res = 220)
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
  forest(m_a, header = c("Study", paste0(xlab_text, " [95% CI]")), xlab = xlab_text,
         refline = 0, mlab = "Pooled (REML, HKSJ)", digits = 2); title(lab_a)
  forest(m_b, header = c("Study", paste0(xlab_text, " [95% CI]")), xlab = xlab_text,
         refline = 0, mlab = "Pooled (REML, HKSJ)", digits = 2); title(lab_b)
  dev.off(); cat("  Saved:", outfile, "\n")
}
if (!is.null(m_lean_t2d) && !is.null(m_lean_obesity))
  forest_two(m_lean_t2d, m_lean_obesity, "A. Type 2 diabetes (k = 3)", "B. Obesity (k = 2)",
             "output/figures/forest_lean_by_population.png", "Mean difference in lean mass (kg)")
if (!is.null(m_fat_t2d) && !is.null(m_fat_obesity))
  forest_two(m_fat_t2d, m_fat_obesity, "A. Type 2 diabetes (k = 3)", "B. Obesity (k = 2)",
             "output/figures/forest_fat_by_population.png", "Mean difference in fat mass (kg)")
if (!is.null(m_lean_act) && !is.null(m_lean_pbo))
  forest_two(m_lean_act, m_lean_pbo, "A. Active comparator (k = 3)", "B. Placebo (k = 2)",
             "output/figures/forest_lean_by_comparator.png", "Mean difference in lean mass (kg)")

# ---- 13. BAYESIAN — brms ----

cat("\n================================================================\n")
cat("SECTION 13 — BAYESIAN RANDOM-EFFECTS (brms / Stan)\n")
cat("================================================================\n")

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
cat("\n[brms Lean — fitting...]\n")
bayes_lean <- tryCatch(fit_brms(lean_md_kg | se(lean_se_kg) ~ 1 + (1 | study_id)),
                       error = function(e) { cat("brms lean failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(bayes_lean)) {
  saveRDS(bayes_lean, "output/bayes/brms_bayes_lean.rds")
  mu_post_l <- summ_brms(bayes_lean, "Lean", c(0, -1, -2))
}
cat("\n[brms Fat — fitting...]\n")
bayes_fat <- tryCatch(fit_brms(fat_md_kg | se(fat_se_kg) ~ 1 + (1 | study_id), seed_add = 1),
                      error = function(e) { cat("brms fat failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(bayes_fat)) {
  saveRDS(bayes_fat, "output/bayes/brms_bayes_fat.rds")
  mu_post_f <- summ_brms(bayes_fat, "Fat", c(0, -3, -5))
}

cat("\n[brms Prior sensitivity — Lean]\n")
prior_sens_results <- data.frame()
for (sc in c(0.5, 1.0, 2.0)) {
  cat(sprintf("\n  scale = %.1f ...\n", sc))
  pr_sens <- c(prior(normal(0, 10), class = "Intercept"),
               prior_string(sprintf("cauchy(0, %s)", sc), class = "sd"))
  fit_s <- tryCatch(brm(lean_md_kg | se(lean_se_kg) ~ 1 + (1 | study_id), data = dat_primary,
                        prior = pr_sens, iter = 4000, warmup = 1000, chains = 4, cores = n_cores,
                        seed = 20260508, control = list(adapt_delta = 0.99), refresh = 0),
                    error = function(e) NULL)
  if (!is.null(fit_s)) {
    mu_s <- as_draws_matrix(fit_s)[, "b_Intercept"]
    prior_sens_results <- rbind(prior_sens_results, data.frame(
      tau_prior_scale = sc, mu_mean = round(mean(mu_s), 3),
      crI_lb = round(unname(quantile(mu_s, .025)), 3), crI_ub = round(unname(quantile(mu_s, .975)), 3),
      p_mu_lt_0 = round(mean(mu_s < 0), 4)))
  }
}
if (nrow(prior_sens_results) > 0) {
  cat("\nPrior sensitivity table (Lean):\n"); print(prior_sens_results)
  write_csv(prior_sens_results, "output/tables/brms_prior_sensitivity_lean.csv")
}

# ---- 14. BAYESIAN — bayesmeta (cross-validation) ----

cat("\n================================================================\n")
cat("SECTION 14 — BAYESIAN CROSS-VALIDATION (bayesmeta)\n")
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
  tryCatch(plot(bm, which = 3, main = paste("bayesmeta —", label)),
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
  cv_table$abs_diff_pct <- round(abs(cv_table$brms_mean - cv_table$bayesmeta_mean) / abs(cv_table$brms_mean) * 100, 2)
  cat("\n[Cross-validation: brms vs bayesmeta]\n"); print(cv_table)
  write_csv(cv_table, "output/tables/bayesian_cross_validation.csv")
}

# ---- 15. KEY RESULTS SUMMARY ----

cat("\n================================================================\n")
cat("SECTION 15 — KEY RESULTS SUMMARY\n")
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
  extract(m_lean_mono,    "Lean", "Mono GLP-1 RA"),
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

# ---- 16. SESSION INFO ----

cat("\n================================================================\n")
cat("SECTION 16 — SESSION INFO\n")
cat("================================================================\n")
print(sessionInfo())
cat("\nALL ANALYSES COMPLETE. Output: ./output/\n")

sink(); close(con)
