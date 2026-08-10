# ================
# Blood Cancer Projections (lymphoma): damping-rate sensitivity (aggregate tier)
# ================
# Varies the period/cohort damping rate (the `damping` constant in
# _setup.R; base 0.92 = 8%/yr decay of the excess RR toward 1.0) and
# reports, for the aggregate tier (NHL, HL):
#   (i)  out-of-sample holdout error (refit to 2011, predict 2012-2021)
#        -- how well each rate calibrates against observed incidence;
#   (ii) the 2045 projection (cases, ASR) -- how much the headline moves.
#
# Motivation: the holdout shows the 1982-start aggregate NHL over-projects,
# and knot selection does not fix it (see knot_selection.R / _notes.md).
# This isolates the damping rate as a candidate lever. Result (see
# _notes.md "Damping sensitivity"): the 2045 projection is ROBUST to the
# damping rate (NHL 2045 moves <5% across 0-12%/yr), the over-projection
# is intrinsic to the long-series APC trend (present even undamped), and
# the base 8%/yr is retained (literature convention + far-horizon
# stability). Reported as a robustness/sensitivity analysis, not a tuning
# knob.
#
# Diagnostic driver, run separately:
#   Rscript code/damping_sensitivity.R  -> output/table_s_damping_sensitivity.csv
# ================

set.seed(20260507)
suppressPackageStartupMessages({
  library(tidyverse); library(Epi); library(MASS, exclude = "select")
})
source("code/_setup.R")
source("code/apc_model.R")
source("code/holdout.R")   # build_pop_mat_years, project_one_draw_holdout, holdout_split/years

# Evaluate `expr` with the global `damping` temporarily set to `rate`.
with_damping <- function(rate, expr) {
  old <- get("damping", envir = .GlobalEnv)
  assign("damping", rate, envir = .GlobalEnv)
  on.exit(assign("damping", old, envir = .GlobalEnv))
  force(expr)
}

run_damping_sensitivity <- function(rates_pct = c(0, 4, 8, 12, 16, 20),
                                     B = 200, save_dir = "output",
                                     seed = 20260507, verbose = TRUE) {
  inc_agg  <- readr::read_csv("data/incidence_agg.csv", show_col_types = FALSE)
  pop_hist <- readr::read_csv("data/pop_hist.csv",      show_col_types = FALSE)
  pop_proj <- readr::read_csv("data/pop_proj.csv",      show_col_types = FALSE)
  std_pop  <- build_std_pop(pop_hist)
  pmbs <- list(males = build_pop_mat("males", pop_proj),
               females = build_pop_mat("females", pop_proj))

  # observed held-out counts per (group, sex)
  obs_hold <- function(grp, sx) {
    lbl <- switch(grp, nhl = "Non-Hodgkin lymphoma", hodgkin = "Hodgkin lymphoma")
    d <- inc_agg |> dplyr::filter(cancer_group == lbl, sex == sx) |>
      dplyr::select(year, age_group, count)
    calc_hist_asr(d, sx, pop_hist, std_pop) |>
      dplyr::filter(year %in% holdout_years) |> dplyr::arrange(year)
  }

  rows <- list()
  for (pct in rates_pct) {
    rate <- 1 - pct / 100
    for (grp in agg_groups) for (sx in sexes) {
      df_full <- prep_agg_data(grp, sx, inc_agg, pop_hist)

      # (i) holdout: fit <= 2011, project the held-out decade at this rate
      df_tr <- df_full |> dplyr::filter(P <= holdout_split)
      ft <- suppressWarnings(suppressMessages(apc.fit(
        data = df_tr[, c("A", "P", "D", "Y")], model = "ns", dist = "poisson",
        npar = agg_knots, print.AOV = FALSE)))
      bs <- apc_basis(ft, df_tr)
      pm <- build_pop_mat_years(sx, holdout_years, pop_hist)
      proj_h <- with_damping(rate,
        rowSums(project_one_draw_holdout(coef(ft$Model), bs, pm, holdout_split, holdout_years)))
      obs <- obs_hold(grp, sx)
      mape   <- mean(abs(proj_h - obs$cases_obs) / obs$cases_obs) * 100
      pe2021 <- (proj_h[length(proj_h)] - obs$cases_obs[obs$year == hist_end]) /
                 obs$cases_obs[obs$year == hist_end] * 100

      # (ii) production: full-window fit projected to 2045 at this rate
      set.seed(seed)
      ff <- suppressWarnings(suppressMessages(apc.fit(
        data = df_full[, c("A", "P", "D", "Y")], model = "ns", dist = "poisson",
        npar = agg_knots, print.AOV = FALSE)))
      bf <- apc_basis(ff, df_full)
      cd <- with_damping(rate,
        project_all_draws(sample_beta_draws(ff, B), bf, pmbs[[sx]]))
      entry <- list(list(subtype = grp, sex = sx, tier = "agg", cases_draws = cd))
      pa <- with_damping(rate, build_proj_annual(entry, std_pop, pmbs)) |>
        dplyr::filter(P == 2045)

      rows[[length(rows) + 1]] <- tibble::tibble(
        damping_pct = pct, subtype = grp, sex = sx,
        holdout_mape = mape, holdout_pe_2021 = pe2021,
        cases_2045 = pa$cases_mid, asr_2045 = pa$asr_mid)
    }
    if (verbose) cat(sprintf("  damping %2d%%/yr done\n", pct))
  }
  tab <- dplyr::bind_rows(rows) |> dplyr::arrange(subtype, sex, damping_pct)
  readr::write_csv(tab, file.path(save_dir, "table_s_damping_sensitivity.csv"))

  if (verbose) {
    base_pct <- round((1 - get("damping", envir = .GlobalEnv)) * 100)
    cat(sprintf("\n=== Damping sensitivity (aggregate tier); base = %d%%/yr ===\n", base_pct))
    show <- tab |>
      dplyr::mutate(M = paste(toupper(dplyr::recode(subtype, hodgkin = "hl")),
                              ifelse(sex == "males", "M", "F"))) |>
      dplyr::select(M, damping_pct, holdout_pe_2021, holdout_mape, cases_2045) |>
      dplyr::mutate(dplyr::across(c(holdout_pe_2021, holdout_mape), ~round(., 1)),
                    cases_2045 = round(cases_2045))
    print(as.data.frame(show), row.names = FALSE)
    cat("\nWrote table_s_damping_sensitivity.csv\n")
  }
  invisible(tab)
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_damping_sensitivity())
}
