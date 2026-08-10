# ================
# Blood Cancer Projections (lymphoma): COVID (2020-21) exclusion sensitivity
# ================
# Refits every model EXCLUDING the 2020-2021 pandemic years and compares the
# 2045 projection against the base (all data, fit to 2021). Purpose (QL): test
# whether the pandemic years distort the projections; if not, keep all data
# and report exclusion as a sensitivity.
#
# Base:      fit <start>-2021, project 2022-2045 (the production model; read
#            from output/incidence_projections.csv).
# Excluded:  fit <start>-2019, project 2020-2045 (2020-21 treated as projected),
#            computed here with the same damping engine, seed and knots.
# Population for projected years is observed ERP where available (2020-21, from
# pop_hist) and ABS projections thereafter (pop_proj).
#
# Context (see _notes.md "COVID sensitivity"): Australian lymphoma incidence
# shows no 2020 dip (2020 is often the high point), so there is no pandemic
# artifact to remove; the test is whether the recent years bias the fit.
#
# Diagnostic driver, run separately:
#   Rscript code/covid_sensitivity.R -> output/table_s_covid_sensitivity.csv
# ================

set.seed(20260507)
suppressPackageStartupMessages({
  library(tidyverse); library(Epi); library(MASS, exclude = "select")
})
source("code/_setup.R")
source("code/apc_model.R")
source("code/holdout.R")   # build_pop_mat_years, project_one_draw_holdout

# Population matrix for arbitrary future years: observed ERP (pop_hist) for
# years <= hist_end, ABS projections (pop_proj) thereafter.
pop_mat_stitch <- function(sex_name, years, pop_hist, pop_proj) {
  h <- pop_hist |> dplyr::filter(sex == sex_name, year %in% years[years <= hist_end])
  p <- pop_proj |> dplyr::filter(sex == sex_name, year %in% years[years >  hist_end])
  w <- dplyr::bind_rows(h, p) |>
    dplyr::filter(age_group %in% names(age_mid)) |>
    dplyr::select(year, age_group, population) |>
    tidyr::pivot_wider(names_from = age_group, values_from = population) |>
    dplyr::arrange(year)
  if (nrow(w) != length(years) || any(is.na(w[, names(age_mid)]))) {
    stop("population missing for sex=", sex_name, " across requested years")
  }
  m <- as.matrix(w[, names(age_mid)]); rownames(m) <- as.character(w$year); m
}

# MC 2045 projection for a model fitted to `fit_end`, projecting to 2045.
covid_project_2045 <- function(tier, name, sex_name, fit_end, B,
                               inc_subtype, inc_agg, pop_hist, pop_proj, std_pop) {
  start <- if (tier == "agg") agg_start else subtype_start
  knots <- if (tier == "agg") agg_knots else subtype_knots
  df_full <- if (tier == "agg") prep_agg_data(name, sex_name, inc_agg, pop_hist)
             else               prep_subtype_data(name, sex_name, inc_subtype, pop_hist)
  df <- df_full |> dplyr::filter(P <= fit_end)

  ft <- suppressWarnings(suppressMessages(apc.fit(
    data = df[, c("A", "P", "D", "Y")], model = "ns", dist = "poisson",
    npar = knots, print.AOV = FALSE)))
  bs <- apc_basis(ft, df)
  proj_years <- (fit_end + 1):proj_end
  pm <- pop_mat_stitch(sex_name, proj_years, pop_hist, pop_proj)
  bd <- sample_beta_draws(ft, B)

  cases <- numeric(B); asr <- numeric(B)
  for (b in seq_len(B)) {
    cm <- project_one_draw_holdout(bd[b, ], bs, pm, fit_end, proj_years)
    r  <- cm["2045", ] / pm["2045", ] * 1e5
    cases[b] <- sum(cm["2045", ]); asr[b] <- sum(r * std_pop[names(r)]) / sum(std_pop)
  }
  list(cases = cases, asr = asr)
}

run_covid_sensitivity <- function(B = 1000, save_dir = "output",
                                  seed = 20260507, verbose = TRUE) {
  inc_subtype <- readr::read_csv("data/incidence_subtype.csv", show_col_types = FALSE)
  inc_agg     <- readr::read_csv("data/incidence_agg.csv",     show_col_types = FALSE)
  pop_hist    <- readr::read_csv("data/pop_hist.csv",          show_col_types = FALSE)
  pop_proj    <- readr::read_csv("data/pop_proj.csv",          show_col_types = FALSE)
  std_pop     <- build_std_pop(pop_hist)

  # Base 2045 (median + 95% CrI) from the production projections.
  base <- readr::read_csv(file.path(save_dir, "incidence_projections.csv"),
                          show_col_types = FALSE) |>
    dplyr::filter(P == 2045) |>
    dplyr::select(subtype, sex, tier,
                  base_mid = cases_mid, base_p025 = cases_p025, base_p975 = cases_p975)

  models <- c(lapply(agg_groups,   function(g) list(tier = "agg",     name = g)),
              lapply(nhl_subtypes, function(s) list(tier = "subtype", name = s)))
  rows <- list()
  for (m in models) for (sx in sexes) {
    set.seed(seed)
    ex <- covid_project_2045(m$tier, m$name, sx, hist_end - 2, B,
                             inc_subtype, inc_agg, pop_hist, pop_proj, std_pop)
    rows[[length(rows) + 1]] <- tibble::tibble(
      subtype = m$name, sex = sx, tier = m$tier,
      excl_mid = median(ex$cases),
      excl_p025 = unname(quantile(ex$cases, 0.025)),
      excl_p975 = unname(quantile(ex$cases, 0.975)))
    if (verbose) cat(sprintf("  %s %s done\n", m$name, sx))
  }
  tab <- dplyr::bind_rows(rows) |>
    dplyr::left_join(base, by = c("subtype", "sex", "tier")) |>
    dplyr::mutate(
      diff_pct       = (excl_mid - base_mid) / base_mid * 100,
      excl_within_ci = excl_mid >= base_p025 & excl_mid <= base_p975) |>
    dplyr::select(subtype, sex, tier, base_mid, base_p025, base_p975,
                  excl_mid, diff_pct, excl_within_ci) |>
    dplyr::arrange(tier, subtype, sex)
  readr::write_csv(tab, file.path(save_dir, "table_s_covid_sensitivity.csv"))

  if (verbose) {
    cat("\n=== COVID (2020-21) exclusion: 2045 cases, base vs excluded ===\n")
    show <- tab |>
      dplyr::mutate(M = paste(toupper(dplyr::recode(subtype, hodgkin = "hl", nhl = "nhl")),
                              ifelse(sex == "males", "M", "F")),
                    base = sprintf("%d (%d-%d)", round(base_mid), round(base_p025), round(base_p975))) |>
      dplyr::select(M, base, excl_mid, diff_pct, excl_within_ci) |>
      dplyr::mutate(excl_mid = round(excl_mid), diff_pct = round(diff_pct, 1))
    print(as.data.frame(show), row.names = FALSE)
    cat(sprintf("\nAll excluded-model 2045 medians within the base 95%% CrI: %s\n",
                all(tab$excl_within_ci)))
    cat("Wrote table_s_covid_sensitivity.csv\n")
  }
  invisible(tab)
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_covid_sensitivity())
}
