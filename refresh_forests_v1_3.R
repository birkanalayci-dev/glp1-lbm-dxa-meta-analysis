# refresh_forests_v1_3.R
# Re-renders the forest plots of glp1_dxa_meta_v1_3.R (Section 11) without
# re-running the Bayesian models. Run from C:/R DXA:  source("refresh_forests_v1_3.R")
# The code of Section 11 is identical to the main script.
suppressPackageStartupMessages({ library(metafor); library(dplyr); library(readr) })
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
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
fit_hksj <- function(d, yi, vi) {
  rma(yi = d[[yi]], vi = d[[vi]], slab = d$study_id, method = "REML", test = "knha")
}
pi_of <- function(m) {
  p <- predict(m, level = 95)
  c(PI_lb = unname(p$pi.lb), PI_ub = unname(p$pi.ub))
}
m_lean_hksj <- fit_hksj(dat_primary, "lean_md_kg", "lean_vi")
m_fat_hksj  <- fit_hksj(dat_primary, "fat_md_kg",  "fat_vi")
sub_fit <- function(d, yi, vi) if (nrow(d) >= 2) fit_hksj(d, yi, vi) else NULL
m_lean_t2d     <- sub_fit(dat_primary %>% filter(population == "T2D"),     "lean_md_kg", "lean_vi")
m_lean_obesity <- sub_fit(dat_primary %>% filter(population == "obesity"), "lean_md_kg", "lean_vi")
m_fat_t2d      <- sub_fit(dat_primary %>% filter(population == "T2D"),     "fat_md_kg",  "fat_vi")
m_fat_obesity  <- sub_fit(dat_primary %>% filter(population == "obesity"), "fat_md_kg",  "fat_vi")
m_lean_act     <- sub_fit(dat_primary %>% filter(comparator_type == "active"),  "lean_md_kg", "lean_vi")
m_lean_pbo     <- sub_fit(dat_primary %>% filter(comparator_type == "placebo"), "lean_md_kg", "lean_vi")

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

