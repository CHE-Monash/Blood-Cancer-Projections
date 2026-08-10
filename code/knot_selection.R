# ================
# Blood Cancer Projections (lymphoma): APC knot-selection grid
# ================
# Refits every model over a grid of natural-cubic-spline knot counts for
# age, period and cohort, and scores each configuration by (i) OUT-OF-SAMPLE
# holdout error — refit to 2011, predict 2012-2021 vs observed counts, the
# same criterion as code/holdout.R — and (ii) in-sample fit on the full
# window (deviance/df, AIC). Selection follows the myeloma template
# (Luo 2024): the parsimonious configuration with strong out-of-sample
# validation (QL major comment #2). Motivates the knot counts hard-coded as
# agg_knots / subtype_knots in _setup.R and documented in _notes.md.
#
# This is a diagnostic/selection script, run separately from the main
# apc -> prev -> supplement chain:
#   Rscript code/knot_selection.R   -> writes output/table_s_knot_selection.csv
# Runtime ~2-3 min (grid x 10 models x {truncated + full} fits).
# ================

set.seed(20260507)
suppressPackageStartupMessages({
  library(tidyverse); library(Epi); library(MASS, exclude = "select")
})
source("code/_setup.R")
source("code/apc_model.R")   # apc_basis, apc_components, prep_*_data, build_std_pop, calc_hist_asr
source("code/holdout.R")     # build_pop_mat_years, project_one_draw_holdout, holdout_split/years

# Knot grid. Minimum 4: apc.fit()/apc_basis need >=2 columns per A/C block,
# which a 3-knot ns basis does not always provide.
knot_grid_vals <- 4:7

quiet_apc_fit <- function(df, cfg) tryCatch(
  suppressWarnings(suppressMessages(
    apc.fit(data = df[, c("A", "P", "D", "Y")], model = "ns", dist = "poisson",
            npar = c(A = cfg[1], P = cfg[2], C = cfg[3]), print.AOV = FALSE))),
  error = function(e) NULL)

# Score one (model, config): out-of-sample holdout cases-MAPE + in-sample fit.
score_knot_config <- function(tier, name, sex_name, fit_start, cfg,
                              inc_subtype, inc_agg, pop_hist, std_pop) {
  df_full <- if (tier == "agg") prep_agg_data(name, sex_name, inc_agg, pop_hist)
             else               prep_subtype_data(name, sex_name, inc_subtype, pop_hist)
  df_tr <- df_full |> dplyr::filter(P >= fit_start, P <= holdout_split)

  ft <- quiet_apc_fit(df_tr, cfg); if (is.null(ft)) return(NULL)
  bs <- tryCatch(apc_basis(ft, df_tr), error = function(e) NULL); if (is.null(bs)) return(NULL)
  pm <- build_pop_mat_years(sex_name, holdout_years, pop_hist)
  cases <- tryCatch(project_one_draw_holdout(coef(ft$Model), bs, pm, holdout_split, holdout_years),
                    error = function(e) NULL)
  if (is.null(cases) || any(!is.finite(cases))) return(NULL)
  proj <- rowSums(cases)

  obs_inc <- if (tier == "agg") {
    lbl <- switch(name, nhl = "Non-Hodgkin lymphoma", hodgkin = "Hodgkin lymphoma")
    inc_agg |> dplyr::filter(cancer_group == lbl, sex == sex_name) |>
      dplyr::select(year, age_group, count)
  } else {
    inc_subtype |> dplyr::filter(subtype == name, sex == sex_name) |>
      dplyr::select(year, age_group, count)
  }
  obs <- calc_hist_asr(obs_inc, sex_name, pop_hist, std_pop) |>
    dplyr::filter(year %in% holdout_years) |> dplyr::arrange(year)

  ff <- quiet_apc_fit(df_full, cfg); if (is.null(ff)) return(NULL)
  m <- ff$Model
  tibble::tibble(
    tier = tier, subtype = name, sex = sex_name,
    n_a = cfg[1], n_p = cfg[2], n_c = cfg[3], tot_knots = sum(cfg),
    cases_mape = mean(abs(proj - obs$cases_obs) / obs$cases_obs) * 100,
    pe_2021    = (proj[length(proj)] - obs$cases_obs[obs$year == hist_end]) /
                  obs$cases_obs[obs$year == hist_end] * 100,
    n_coef = length(coef(m)),
    dev_df_full = m$deviance / m$df.residual, aic_full = AIC(m)
  )
}

run_knot_selection <- function(save_dir = "output", verbose = TRUE) {
  inc_subtype <- readr::read_csv("data/incidence_subtype.csv", show_col_types = FALSE)
  inc_agg     <- readr::read_csv("data/incidence_agg.csv",     show_col_types = FALSE)
  pop_hist    <- readr::read_csv("data/pop_hist.csv",          show_col_types = FALSE)
  std_pop     <- build_std_pop(pop_hist)

  tasks <- c(
    unlist(lapply(agg_groups,   function(g) lapply(sexes, function(sx)
      list(tier = "agg", name = g, sex = sx, fit_start = agg_start))), recursive = FALSE),
    unlist(lapply(nhl_subtypes, function(s) lapply(sexes, function(sx)
      list(tier = "subtype", name = s, sex = sx, fit_start = subtype_start))), recursive = FALSE)
  )
  grid <- as.matrix(expand.grid(n_a = knot_grid_vals, n_p = knot_grid_vals, n_c = knot_grid_vals))

  rows <- list()
  for (ci in seq_len(nrow(grid))) {
    cfg <- grid[ci, ]
    for (tk in tasks) {
      r <- score_knot_config(tk$tier, tk$name, tk$sex, tk$fit_start, cfg,
                             inc_subtype, inc_agg, pop_hist, std_pop)
      if (!is.null(r)) rows[[length(rows) + 1]] <- r
    }
    if (verbose && ci %% 8 == 0) cat(sprintf("  scored %d / %d configs\n", ci, nrow(grid)))
  }
  detail <- dplyr::bind_rows(rows)

  # Per-tier mean over that tier's models (tier-mean out-of-sample error).
  n_expected <- c(agg = length(agg_groups) * length(sexes),
                  subtype = length(nhl_subtypes) * length(sexes))
  tier_summary <- detail |>
    dplyr::group_by(tier, n_a, n_p, n_c, tot_knots) |>
    dplyr::summarise(
      n_models      = dplyr::n(),
      mean_mape     = mean(cases_mape),
      max_mape      = max(cases_mape),
      mean_abs_pe21 = mean(abs(pe_2021)),
      mean_dev_df   = mean(dev_df_full),
      mean_aic      = mean(aic_full),
      .groups = "drop"
    ) |>
    dplyr::filter(n_models == n_expected[tier]) |>
    dplyr::arrange(tier, mean_mape)

  readr::write_csv(detail,       file.path(save_dir, "knot_selection_detail.csv"))
  readr::write_csv(tier_summary, file.path(save_dir, "table_s_knot_selection.csv"))

  if (verbose) {
    for (tt in c("agg", "subtype")) {
      cat(sprintf("\n=== Tier %s — top 6 configs by out-of-sample cases-MAPE ===\n", tt))
      st <- tier_summary |> dplyr::filter(tier == tt)
      print(as.data.frame(head(st, 6) |>
        dplyr::select(n_a, n_p, n_c, mean_mape, max_mape, mean_abs_pe21, mean_dev_df) |>
        dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 2)))), row.names = FALSE)
      sel <- if (tt == "agg") agg_knots else subtype_knots
      srow <- st |> dplyr::filter(n_a == sel["A"], n_p == sel["P"], n_c == sel["C"])
      cat(sprintf("  selected (%d,%d,%d): MAPE %.2f, rank %d/%d\n",
                  sel["A"], sel["P"], sel["C"], srow$mean_mape,
                  which(st$n_a == sel["A"] & st$n_p == sel["P"] & st$n_c == sel["C"]), nrow(st)))
    }
    cat("\nWrote table_s_knot_selection.csv (tier summary) + knot_selection_detail.csv (per model)\n")
  }
  invisible(list(detail = detail, tier_summary = tier_summary))
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_knot_selection())
}
