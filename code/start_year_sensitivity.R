# ================
# Blood Cancer Projections (lymphoma): aggregate-tier start-year sensitivity
# ================
# Varies the aggregate-tier fit-start year (the `agg_start` constant in
# _setup.R; base 1990) and reports, for NHL and HL:
#   (i)  out-of-sample holdout error (refit to 2011, predict 2012-2021)
#        -- how well each start year calibrates against observed incidence;
#   (ii) the 2045 projection (cases, growth) -- the headline sensitivity.
#
# Motivation and result (see _notes.md "Start-year selection"): fitting the
# aggregate tier from the full observed series (1982) over-projects NHL,
# because it extrapolates the transient pre-1990 NHL surge (HIV/AIDS-
# associated NHL + pre-WHO reclassification). Out-of-sample error falls
# sharply and then PLATEAUS for any start year >= ~1990, so 1990 is chosen
# as the earliest well-validated year (longest series, 32 yr). Subtype-tier
# data begin in 2003, so only the aggregate tier is affected.
#
# Diagnostic driver, run separately:
#   Rscript code/start_year_sensitivity.R -> output/table_s_start_year_sensitivity.csv
# ================

set.seed(20260507)
suppressPackageStartupMessages({
  library(tidyverse); library(Epi); library(MASS, exclude = "select")
})
source("code/_setup.R")
source("code/apc_model.R")
source("code/holdout.R")   # build_pop_mat_years, project_one_draw_holdout, holdout_split/years

# Evaluate `expr` with the global `agg_start` temporarily set to `yr`
# (agg_start governs the P >= agg_start filter in prep_agg_data()).
with_agg_start <- function(yr, expr) {
  old <- get("agg_start", envir = .GlobalEnv)
  assign("agg_start", yr, envir = .GlobalEnv)
  on.exit(assign("agg_start", old, envir = .GlobalEnv))
  force(expr)
}

run_start_year_sensitivity <- function(starts = c(1982, 1985, 1990, 1995, 2000, 2003),
                                       B = 200, save_dir = "output",
                                       seed = 20260507, verbose = TRUE) {
  inc_agg  <- readr::read_csv("data/incidence_agg.csv", show_col_types = FALSE)
  pop_hist <- readr::read_csv("data/pop_hist.csv",      show_col_types = FALSE)
  pop_proj <- readr::read_csv("data/pop_proj.csv",      show_col_types = FALSE)
  std_pop  <- build_std_pop(pop_hist)
  pmbs <- list(males = build_pop_mat("males", pop_proj),
               females = build_pop_mat("females", pop_proj))

  obs_of <- function(grp, sx) {
    lbl <- switch(grp, nhl = "Non-Hodgkin lymphoma", hodgkin = "Hodgkin lymphoma")
    d <- inc_agg |> dplyr::filter(cancer_group == lbl, sex == sx) |>
      dplyr::select(year, age_group, count)
    calc_hist_asr(d, sx, pop_hist, std_pop)
  }

  rows <- list()
  for (yr in starts) {
    for (grp in agg_groups) for (sx in sexes) {
      row <- with_agg_start(yr, {
        # prep_agg_data now filters P >= agg_start (= yr)
        df_full <- prep_agg_data(grp, sx, inc_agg, pop_hist)

        # (i) holdout: fit yr..2011, predict 2012-2021
        df_tr <- df_full |> dplyr::filter(P <= holdout_split)
        ft <- suppressWarnings(suppressMessages(apc.fit(
          data = df_tr[, c("A", "P", "D", "Y")], model = "ns", dist = "poisson",
          npar = agg_knots, print.AOV = FALSE)))
        bs <- apc_basis(ft, df_tr)
        pm <- build_pop_mat_years(sx, holdout_years, pop_hist)
        ph <- rowSums(project_one_draw_holdout(coef(ft$Model), bs, pm, holdout_split, holdout_years))
        obs <- obs_of(grp, sx); oh <- obs |> dplyr::filter(year %in% holdout_years) |> dplyr::arrange(year)
        mape   <- mean(abs(ph - oh$cases_obs) / oh$cases_obs) * 100
        pe2021 <- (ph[length(ph)] - oh$cases_obs[oh$year == hist_end]) / oh$cases_obs[oh$year == hist_end] * 100

        # (ii) production: fit yr..2021, project to 2045
        set.seed(seed)
        ff <- suppressWarnings(suppressMessages(apc.fit(
          data = df_full[, c("A", "P", "D", "Y")], model = "ns", dist = "poisson",
          npar = agg_knots, print.AOV = FALSE)))
        bf <- apc_basis(ff, df_full)
        cd <- project_all_draws(sample_beta_draws(ff, B), bf, pmbs[[sx]])
        pa <- build_proj_annual(list(list(subtype = grp, sex = sx, tier = "agg", cases_draws = cd)),
                                std_pop, pmbs) |> dplyr::filter(P == 2045)
        c2021 <- obs$cases_obs[obs$year == hist_end]

        tibble::tibble(
          start_year = yr, n_years = 2011 - yr + 1, subtype = grp, sex = sx,
          holdout_mape = mape, holdout_pe_2021 = pe2021,
          cases_2045 = pa$cases_mid, growth_pct = (pa$cases_mid - c2021) / c2021 * 100,
          dev_df = ff$Model$deviance / ff$Model$df.residual)
      })
      rows[[length(rows) + 1]] <- row
    }
    if (verbose) cat(sprintf("  start %d done\n", yr))
  }
  tab <- dplyr::bind_rows(rows) |> dplyr::arrange(subtype, sex, start_year)
  readr::write_csv(tab, file.path(save_dir, "table_s_start_year_sensitivity.csv"))

  if (verbose) {
    cat(sprintf("\n=== Start-year sensitivity (aggregate tier); base = %d ===\n",
                get("agg_start", envir = .GlobalEnv)))
    show <- tab |>
      dplyr::mutate(M = paste(toupper(dplyr::recode(subtype, hodgkin = "hl")),
                              ifelse(sex == "males", "M", "F"))) |>
      dplyr::select(M, start_year, n_years, holdout_pe_2021, holdout_mape, cases_2045, growth_pct) |>
      dplyr::mutate(dplyr::across(c(holdout_pe_2021, holdout_mape, growth_pct), ~round(., 1)),
                    cases_2045 = round(cases_2045))
    print(as.data.frame(show), row.names = FALSE)
    cat("\nWrote table_s_start_year_sensitivity.csv\n")
  }
  invisible(tab)
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_start_year_sensitivity())
}
