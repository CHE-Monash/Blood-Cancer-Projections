# ================
# Blood Cancer Projections (lymphoma): Out-of-sample incidence holdout validation
# ================
# Primary out-of-sample validation of the incidence projections
# (QL major comment #1; myeloma-template approach, Luo et al. 2024).
#
# Each of the 10 APC models is refitted on a TRUNCATED window that ends
# at `holdout_split` (2011), then projected forward over the held-out
# years 2012-2021 and compared against OBSERVED registry incidence
# (case counts + ASR). This is genuine out-of-sample validation against
# hard registry ground truth, distinct from the AIHW comparison (which
# is a forecast-vs-forecast "comparison with prior forecasts", not
# validation).
#
# Split (see _notes.md "Incidence holdout validation"):
#   Aggregate tier (NHL, HL): fit 1982-2011 (30 yr) -> hold out 2012-2021.
#   Subtype tier (DLBCL, FL, MCL): fit 2003-2011 (9 yr) -> hold out 2012-2021.
# Both tiers share the 2012-2021 holdout so aggregate and subtype are
# validated on the same out-of-sample years. Projection uses the same
# damping engine, seed (20260507) and Monte Carlo machinery as the
# production model; only the fit window and the projection target years
# (which draw HISTORICAL population, not ABS projections) differ.
#
# Library + driver pattern:
#   source("code/holdout.R")   -> defines functions; no side effects
#   run_holdout()              -> runs validation, writes table + figure
#   Rscript code/holdout.R     -> auto-runs via the bottom guard
# ================

set.seed(20260507)
stopifnot(
  requireNamespace("MASS"), requireNamespace("Epi"),
  requireNamespace("dplyr"), requireNamespace("readr"),
  requireNamespace("purrr"), requireNamespace("tidyr"),
  requireNamespace("ggplot2")
)
suppressPackageStartupMessages({
  library(tidyverse); library(Epi); library(MASS, exclude = "select")
})

# Shared constants + the APC library (apc_basis, apc_components,
# sample_beta_draws, prep_*_data, build_std_pop, calc_hist_asr, save_fig).
# Sourcing apc_model.R defines functions only (guard is sys.nframe()==0).
source("code/_setup.R")
source("code/apc_model.R")

# ----------------------
# Holdout configuration
# ----------------------
holdout_split <- 2011           # last year in every truncated fit
holdout_years <- (holdout_split + 1):hist_end   # 2012-2021, compared vs observed

# ----------------------
# Historical population matrix for arbitrary years (rows = years, cols = 9 age groups)
# ----------------------
# Mirrors build_pop_mat() but pulls OBSERVED ERP from pop_hist for the
# held-out years, rather than ABS projections from pop_proj.

build_pop_mat_years <- function(sex_name, years, pop_hist) {
  wide <- pop_hist |>
    dplyr::filter(sex == sex_name, year %in% years,
                  age_group %in% names(age_mid)) |>
    dplyr::select(year, age_group, population) |>
    tidyr::pivot_wider(names_from = age_group, values_from = population) |>
    dplyr::arrange(year)
  if (nrow(wide) != length(years) || any(is.na(wide[, names(age_mid)]))) {
    stop("pop_hist missing entries for sex=", sex_name,
         " across holdout years x 9 age groups")
  }
  m <- as.matrix(wide[, names(age_mid)])
  rownames(m) <- as.character(wide$year)
  m
}

# ----------------------
# Damped projection over the holdout years (per draw)
# ----------------------
# Identical construction to apc_model.R::project_one_draw(), but the
# projection anchor `fit_end` and the target `proj_years` are explicit
# arguments (rather than the hist_end/proj_start/proj_end globals) so the
# production path is left untouched.

project_one_draw_holdout <- function(beta, basis, pop_mat, fit_end, proj_years) {
  comp   <- apc_components(beta, basis)
  age_rr <- comp$Age
  per_rr <- comp$Per
  coh_rr <- comp$Coh

  last_per_RR <- per_rr[length(per_rr)]
  last_coh_C  <- as.numeric(names(coh_rr))[length(coh_rr)]
  last_coh_RR <- coh_rr[length(coh_rr)]

  coh_interp <- approxfun(as.numeric(names(coh_rr)), as.numeric(coh_rr),
                          method = "linear", rule = 2)

  P_grid <- proj_years
  A_grid <- age_mid

  t_per       <- P_grid - fit_end
  per_rr_proj <- 1 + (last_per_RR - 1) * damping ^ t_per

  C_mat <- outer(P_grid, A_grid, FUN = "-")
  t_coh <- pmax(C_mat - last_coh_C, 0)
  coh_rr_within <- matrix(coh_interp(C_mat), nrow = length(P_grid),
                          ncol = length(A_grid))
  coh_rr_beyond <- 1 + (last_coh_RR - 1) * damping ^ t_coh
  coh_rr_proj   <- ifelse(C_mat <= last_coh_C, coh_rr_within, coh_rr_beyond)

  rate_mat <- outer(per_rr_proj, age_rr, "*") * coh_rr_proj
  cases    <- rate_mat / 1e5 * pop_mat
  dimnames(cases) <- list(as.character(P_grid), names(A_grid))
  cases
}

# ----------------------
# Fit one model on the truncated window and project the holdout years
# ----------------------
# Returns per-(year, draw) matrices of projected cases and ASR, plus the
# observed series over both the fit window and the holdout window.

run_holdout_model <- function(tier, name, sex_name, fit_start,
                              inc_subtype, inc_agg, pop_hist, std_pop, B) {
  df_full <- if (tier == "agg") {
    prep_agg_data(name, sex_name, inc_agg, pop_hist)
  } else {
    prep_subtype_data(name, sex_name, inc_subtype, pop_hist)
  }
  df <- df_full |> dplyr::filter(P >= fit_start, P <= holdout_split)

  knots <- if (tier == "agg") agg_knots else subtype_knots
  fit   <- apc.fit(data = df[, c("A", "P", "D", "Y")],
                   model = "ns", dist = "poisson",
                   npar = knots, print.AOV = FALSE)
  basis <- apc_basis(fit, df)
  beta_draws <- sample_beta_draws(fit, B)

  pop_mat <- build_pop_mat_years(sex_name, holdout_years, pop_hist)

  # per-draw projected cases (years x age x B)
  cases_arr <- array(NA_real_,
    dim = c(length(holdout_years), length(age_mid), B),
    dimnames = list(as.character(holdout_years), names(age_mid),
                    paste0("d", seq_len(B))))
  for (b in seq_len(B)) {
    cases_arr[, , b] <- project_one_draw_holdout(
      beta_draws[b, ], basis, pop_mat, holdout_split, holdout_years)
  }

  # cases per (year, draw); ASR per (year, draw) on the 2001 ASP
  cases_yr_b <- apply(cases_arr, c(1, 3), sum)
  pop_arr    <- array(rep(pop_mat, B), dim = dim(cases_arr))
  rate_arr   <- cases_arr / pop_arr * 1e5
  std_vec    <- std_pop[dimnames(cases_arr)[[2]]]
  asr_yr_b   <- apply(rate_arr, c(1, 3),
                      function(r) sum(r * std_vec) / sum(std_pop))

  cs <- summarise_draws(cases_yr_b, 1)
  as_ <- summarise_draws(asr_yr_b, 1)

  # observed over the whole available series (fit + holdout), for the figure
  obs_inc <- if (tier == "agg") {
    lbl <- switch(name, nhl = "Non-Hodgkin lymphoma", hodgkin = "Hodgkin lymphoma")
    inc_agg |> dplyr::filter(cancer_group == lbl, sex == sex_name) |>
      dplyr::select(year, age_group, count)
  } else {
    inc_subtype |> dplyr::filter(subtype == name, sex == sex_name) |>
      dplyr::select(year, age_group, count)
  }
  # cap at hist_end: inc_agg carries AIHW projection rows (2022+) that have
  # no matching historical population and would otherwise render as ASR = 0.
  obs <- calc_hist_asr(obs_inc, sex_name, pop_hist, std_pop) |>
    dplyr::filter(year >= fit_start, year <= hist_end)

  proj_tbl <- tibble::tibble(
    year       = holdout_years,
    cases_mid  = cs$mid,  cases_p025 = cs$p025, cases_p975 = cs$p975,
    asr_mid    = as_$mid, asr_p025   = as_$p025, asr_p975  = as_$p975
  ) |>
    dplyr::left_join(obs |> dplyr::select(year, cases_obs, asr_obs), by = "year") |>
    dplyr::mutate(subtype = name, sex = sex_name, tier = tier,
                  fit_start = fit_start)

  list(proj = proj_tbl, obs = obs |> dplyr::mutate(subtype = name, sex = sex_name),
       dev_df = round(fit$Model$deviance / fit$Model$df.residual, 2))
}

# ----------------------
# Summaries: per-model error metrics
# ----------------------

summarise_holdout <- function(by_year) {
  by_year |>
    dplyr::group_by(subtype, sex, tier) |>
    dplyr::summarise(
      n_years        = dplyr::n(),
      cases_mape     = mean(abs(cases_mid - cases_obs) / cases_obs) * 100,
      asr_mape       = mean(abs(asr_mid   - asr_obs)   / asr_obs)   * 100,
      cases_pe_2021  = (cases_mid[year == hist_end] - cases_obs[year == hist_end]) /
                        cases_obs[year == hist_end] * 100,
      asr_pe_2021    = (asr_mid[year == hist_end]   - asr_obs[year == hist_end]) /
                        asr_obs[year == hist_end]   * 100,
      coverage_cases = mean(cases_obs >= cases_p025 & cases_obs <= cases_p975) * 100,
      cases_obs_2021  = cases_obs[year == hist_end],
      cases_mid_2021  = cases_mid[year == hist_end],
      cases_p025_2021 = cases_p025[year == hist_end],
      cases_p975_2021 = cases_p975[year == hist_end],
      .groups = "drop"
    )
}

# ----------------------
# Figure: observed vs projected ASR over the holdout window
# ----------------------
# Observed ASR points across the whole available series; the held-out
# projection (2012-2021) drawn as a line + 95% ribbon; dotted rule at the
# 2011/2012 split. Faceted subtype (rows) x sex (cols).

build_holdout_figure <- function(by_year, obs_all) {
  lvl <- c("NHL", "HL", "DLBCL", "FL", "MCL")
  obs_all  <- obs_all  |> dplyr::mutate(label = factor(label_map[subtype], levels = lvl))
  by_year  <- by_year  |> dplyr::mutate(label = factor(label_map[subtype], levels = lvl))

  ggplot() +
    geom_vline(xintercept = holdout_split + 0.5, linetype = "dotted",
               colour = "grey55", linewidth = 0.4) +
    geom_point(data = obs_all,
               aes(x = year, y = asr_obs), colour = "grey30", size = 0.7) +
    geom_ribbon(data = by_year,
                aes(x = year, ymin = asr_p025, ymax = asr_p975, fill = label),
                alpha = 0.25) +
    geom_line(data = by_year,
              aes(x = year, y = asr_mid, colour = label), linewidth = 0.7) +
    facet_grid(label ~ sex, scales = "free_y",
               labeller = labeller(sex = sex_labels)) +
    scale_colour_manual(values = line_colours, guide = "none") +
    scale_fill_manual(values = line_colours, guide = "none") +
    labs(x = "Year", y = "Age-standardised rate (per 100,000)",
         subtitle = paste0("Grey points: observed. Coloured line + band: ",
                           "out-of-sample projection (fit to ", holdout_split,
                           "), 95% credible interval.")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 10, colour = "grey30"))
}

# ----------------------
# Driver
# ----------------------

run_holdout <- function(B = 1000, save_dir = "output",
                        seed = 20260507, verbose = TRUE) {
  set.seed(seed)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  inc_subtype <- readr::read_csv("data/incidence_subtype.csv", show_col_types = FALSE)
  inc_agg     <- readr::read_csv("data/incidence_agg.csv",     show_col_types = FALSE)
  pop_hist    <- readr::read_csv("data/pop_hist.csv",          show_col_types = FALSE)
  std_pop     <- build_std_pop(pop_hist)

  spec <- c(
    setNames(rep("agg", length(agg_groups)), agg_groups),        # nhl, hodgkin
    setNames(rep("subtype", length(nhl_subtypes)), nhl_subtypes) # dlbcl, follicular, mantle_cell
  )

  results <- list()
  for (name in names(spec)) {
    tier      <- spec[[name]]
    fit_start <- if (tier == "agg") agg_start else subtype_start
    for (sx in sexes) {
      key <- paste(name, sx, sep = "_")
      if (verbose) cat(sprintf("[holdout %s] %s (fit %d-%d)\n",
                               tier, key, fit_start, holdout_split))
      results[[key]] <- run_holdout_model(
        tier, name, sx, fit_start,
        inc_subtype, inc_agg, pop_hist, std_pop, B)
    }
  }

  by_year <- purrr::map_dfr(results, "proj")
  obs_all <- purrr::map_dfr(results, "obs") |>
    dplyr::mutate(label = label_map[subtype])
  summary_tbl <- summarise_holdout(by_year)

  readr::write_csv(by_year,     file.path(save_dir, "holdout_by_year.csv"))
  readr::write_csv(summary_tbl, file.path(save_dir, "table_holdout_validation.csv"))

  fig <- build_holdout_figure(by_year, obs_all)
  save_fig(fig, file.path(save_dir, "figure_holdout_validation"),
           width = 9, height = 11)

  if (verbose) {
    cat("\n=== Holdout validation summary (fit ≤ 2011, project 2012-2021) ===\n")
    print(as.data.frame(summary_tbl |>
      dplyr::select(subtype, sex, cases_mape, asr_mape,
                    cases_pe_2021, coverage_cases) |>
      dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 1)))))
    cat("\nWrote table_holdout_validation.csv, holdout_by_year.csv,",
        "figure_holdout_validation.{pdf,png}\n")
  }

  invisible(list(by_year = by_year, summary = summary_tbl, results = results))
}

# ----------------------
# Auto-run guard
# ----------------------
if (!interactive() && sys.nframe() == 0) {
  invisible(run_holdout())
}
