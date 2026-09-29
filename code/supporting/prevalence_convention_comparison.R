# ================
# Blood Cancer Projections (lymphoma): comparison of prevalence counting
# conventions, and diagnoses at constant age-specific rates
# ================
# Supports one sentence of the Discussion, which compares this study's
# aggregate NHL prevalence with that of Luo et al. (Lancet Reg Health West Pac
# 2026;73:101946). It answers two questions from this study's OWN inputs:
#
# 1. How much of a difference in level is due to the counting convention?
#      end_of_year  (this study, and the AIHW): people alive at the end of
#                   year Y who were diagnosed in years Y-K+1 to Y. The cohort
#                   diagnosed in year Y is weighted by S(0) = sqrt(S(1)).
#      during_year  (Luo et al. 2026, supplement section 3.3): cohorts
#                   diagnosed 1 to K years before year Y, each weighted by the
#                   mean of survival at k-1 and k years, with S(0) = 1. The
#                   cohort diagnosed in year Y itself is not counted.
#    Both are applied to the same diagnoses and the same survival, so any
#    difference between them is the convention and nothing else.
#
# 2. How fast would diagnoses grow from population change alone?
#    Age- and sex-specific rates are held at their 2017-2021 mean and applied
#    to the ABS Series B population.
#
# READ-ONLY with respect to the analysis. Sourcing apc_model.R and
# prev_model.R defines functions only (their drivers sit behind
# sys.nframe() == 0 guards). Nothing is fitted or sampled. Diagnoses for
# 2022 onwards are the medians of the stored draws in output/apc_results.rds;
# survival is the deterministic point curve of build_projected_surv_pt(), so
# the end_of_year column is a point calculation and agrees with the published
# draw-level medians only to within about 0.1%.
#
# Writes output/table_prevalence_convention.csv and
#        output/table_constant_rate_diagnoses.csv
#
#   Rscript code/supporting/prevalence_convention_comparison.R
# ================

grDevices::pdf(NULL)
suppressPackageStartupMessages({
  source("code/apc_model.R")
  source("code/prev_model.R")
})

.nhl_cases <- function(sx, inputs, apc_results) {
  hist <- inputs$inc_agg |>
    dplyr::filter(cancer_group == "Non-Hodgkin lymphoma", sex == sx, year <= hist_end) |>
    dplyr::select(year, age_group, count)
  med <- apply(apc_results[[paste0("nhl_", sx)]]$cases_draws, c(1, 2), median)
  proj <- as.data.frame(as.table(med), stringsAsFactors = FALSE)
  names(proj) <- c("year", "age_group", "count")
  proj$year <- as.integer(proj$year)
  dplyr::bind_rows(hist, proj[proj$year >= proj_start, ])
}

.prevalence <- function(sx, Y, K, convention, cases, surv_obs, improvement) {
  surv <- function(dx) build_projected_surv_pt("nhl", sx, dx, improve_cap,
                                                surv_obs, improvement) |>
    dplyr::select(age_group, t = years_since_dx, s = surv_mid)
  total <- 0
  if (convention == "end_of_year") {
    for (t in 0:(K - 1)) {
      sp <- surv(Y - t)
      st <- if (t == 0) transform(sp[sp$t == 1, ], s = sqrt(s)) else sp[sp$t == t, ]
      m  <- merge(cases[cases$year == Y - t, ], st, by = "age_group")
      stopifnot(nrow(m) == length(age_groups))
      total <- total + sum(m$count * m$s)
    }
  } else {
    for (k in 1:K) {
      sp <- surv(Y - k)
      a  <- if (k == 1) transform(sp[sp$t == 1, ], s = 1) else sp[sp$t == k - 1, ]
      m  <- merge(merge(cases[cases$year == Y - k, ], a, by = "age_group"),
                  sp[sp$t == k, ], by = "age_group", suffixes = c(".a", ".b"))
      stopifnot(nrow(m) == length(age_groups))
      total <- total + sum(m$count * (m$s.a + m$s.b) / 2)
    }
  }
  total
}

run_prevalence_convention_comparison <- function(save_dir = "output",
                                                 years = c(2021, 2025, 2045),
                                                 durations = c(1, 5)) {
  inputs      <- load_prev_inputs("data")
  surv_obs    <- build_surv_obs(inputs$surv_raw)
  improvement <- build_improvement_final(surv_obs)
  apc_results <- readRDS(file.path(save_dir, "apc_results.rds"))
  pop_hist    <- readr::read_csv("data/pop_hist.csv", show_col_types = FALSE)
  pop_proj    <- readr::read_csv("data/pop_proj.csv", show_col_types = FALSE)

  cases <- lapply(stats::setNames(sexes, sexes), .nhl_cases,
                  inputs = inputs, apc_results = apc_results)
  grid <- expand.grid(duration = durations, year = years, stringsAsFactors = FALSE)
  conv <- purrr::pmap_dfr(grid, function(duration, year) {
    v <- sapply(c("end_of_year", "during_year"), function(cv)
      sum(sapply(sexes, function(sx)
        .prevalence(sx, year, duration, cv, cases[[sx]], surv_obs, improvement))))
    tibble::tibble(subtype = "nhl", sex = "persons", duration = duration, year = year,
                   end_of_year = round(v[["end_of_year"]]),
                   during_year = round(v[["during_year"]]),
                   difference_pct = round((v[["during_year"]] / v[["end_of_year"]] - 1) * 100, 1))
  }) |> dplyr::arrange(duration, year)

  rate <- purrr::map_dfr(sexes, function(sx) {
    d <- inputs$inc_agg |>
      dplyr::filter(cancer_group == "Non-Hodgkin lymphoma", sex == sx, year %in% 2017:2021) |>
      dplyr::group_by(age_group) |> dplyr::summarise(d = sum(count), .groups = "drop")
    y <- pop_hist |> dplyr::filter(sex == sx, year %in% 2017:2021) |>
      dplyr::group_by(age_group) |> dplyr::summarise(y = sum(population), .groups = "drop")
    dplyr::inner_join(d, y, by = "age_group") |> dplyr::mutate(sex = sx, rate = d / y)
  })
  stopifnot(nrow(rate) == 2 * length(age_groups))
  const <- pop_proj |>
    dplyr::filter(year %in% c(2024, 2025, 2045, 2049, 2050), sex %in% sexes) |>
    dplyr::inner_join(rate[, c("sex", "age_group", "rate")], by = c("sex", "age_group")) |>
    dplyr::group_by(year) |>
    dplyr::summarise(diagnoses_constant_rate = round(sum(rate * population)), .groups = "drop") |>
    dplyr::mutate(subtype = "nhl", sex = "persons",
                  rates = "mean observed age- and sex-specific rate, 2017-2021",
                  population = "ABS Series B") |>
    dplyr::select(subtype, sex, year, diagnoses_constant_rate, rates, population)

  readr::write_csv(conv,  file.path(save_dir, "table_prevalence_convention.csv"))
  readr::write_csv(const, file.path(save_dir, "table_constant_rate_diagnoses.csv"))
  print(as.data.frame(conv)); print(as.data.frame(const[, 1:4]))
  g <- function(a, b) (const$diagnoses_constant_rate[const$year == b] /
                       const$diagnoses_constant_rate[const$year == a] - 1) * 100
  cat(sprintf("Constant-rate growth in NHL diagnoses: 2025-2045 %+.1f%%; 2025-2050 %+.1f%%; 2024-2049 %+.1f%%\n",
              g(2025, 2045), g(2025, 2050), g(2024, 2049)))
  invisible(list(convention = conv, constant_rate = const))
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_prevalence_convention_comparison())
  invisible(grDevices::dev.off())
}
