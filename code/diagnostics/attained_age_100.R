# ================
# Blood Cancer Projections (lymphoma): diagnostic of attained age above 100
# ================
# Question. Luo et al. (Lancet Reg Health West Pac 2026;73:101946) age
# prevalent survivors through single years of age and drop anyone whose
# attained age passes 100. This study has no attained-age dimension: survival
# is by age at diagnosis in ten-year bands, and beyond year 5 the years-3-to-5
# conditional survival is held constant to year 40. How many of the survivors
# this study counts would Luo's cut-off remove, in 2021, 2031 and 2045?
#
# Method. A weighted analogue of compute_prevalence(). Each term
#   cases[dy, g, b] x S[dy, g, t, b]
# is multiplied by w_g(dy, t), the share of band g's diagnoses in year dy whose
# single-year age at diagnosis x satisfies x + t <= 100 (t = 0 is the
# diagnosis-year cohort). Everything else is the published model, draw for draw.
#
# The single-year age distribution of diagnoses within a band is not observed,
# so it is ASSUMED, in three versions:
#   population (main)  - proportional to the ABS single-year population of that
#                        sex and year within the band. ERP (3101.0 Table 59) to
#                        2021, Series B (3222.0 Table B9) from 2022, parsed as
#                        import_data.R parses them; "100 and over" is age 100.
#                        Population is a PROXY for the age distribution of
#                        diagnoses within each band. Whether it over- or
#                        understates the truncation in the 85+ band is unknown,
#                        since NHL rates within that band are not observed. The
#                        lower-bound assumption gives the least truncation; the
#                        population and midpoint assumptions give similar
#                        results at 40 years.
#   lower              - every diagnosis at the band's lower age (85 for 85+,
#                        75 for 75-84): the least truncation possible.
#   midpoint           - every diagnosis at age_mid from _setup.R (90 for 85+).
#
# Reading the result. This is a sensitivity analysis under the stated
# within-band age allocations, not an estimate of the over-count: the cut-off
# does not make mortality rise between diagnosis and 100, and it also removes
# people who do survive past 100.
#
# READ-ONLY with respect to the analysis. Sourcing apc_model.R and
# prev_model.R defines functions only (their drivers sit behind
# sys.nframe() == 0 guards). The saved fit output/apc_results.rds is used; no
# model is refitted. The draws are those of run_prev_model() (B = 1000, seed
# 20260507, same call sequence). build_cases_array() and build_surv_array()
# consume no random numbers, so the conservative and optimistic arrays are not
# needed to keep the stream aligned and only the base scenario is built. Two
# checks run before anything is reported: with all weights 1 the weighted
# function must equal compute_prevalence() on the same arrays, and the
# untruncated medians are compared with output/prevalence_projections.csv and
# output/table_3_prevalence_combined.csv. The back-estimate is written to a
# temporary file, not to data/. No existing file is modified.
#
# Needs data/raw/ABS/ (git-ignored; see README.md and import_data.R).
#
# Writes output/diagnostic_attained_age_100.csv
#
#   Rscript code/diagnostics/attained_age_100.R        (from the repository root)
# ================

grDevices::pdf(NULL)
suppressPackageStartupMessages({
  source("code/apc_model.R")
  source("code/prev_model.R")
})

AGE_CUT <- 100
band_lower <- c("0–14" = 0, "15–24" = 15, "25–34" = 25, "35–44" = 35, "45–54" = 45,
                "55–64" = 55, "65–74" = 65, "75–84" = 75, "85+" = 85)
stopifnot(identical(names(band_lower), age_groups))

# ----------------------
# Single-year population, 1982-2045, by sex (ages 0-99 and 100 = "100 and over")
# ----------------------
.read_abs_single_year <- function(path, years) {
  nm  <- suppressMessages(names(readxl::read_excel(path, sheet = "Data1", skip = 0)))
  raw <- suppressMessages(readxl::read_excel(path, sheet = "Data1", skip = 9,
                                             col_names = FALSE))
  names(raw) <- nm
  raw |>
    dplyr::rename(date = 1) |>
    dplyr::mutate(year = suppressWarnings(as.integer(format(
      as.Date(as.numeric(date), origin = "1899-12-30"), "%Y")))) |>
    dplyr::filter(year %in% years) |>
    dplyr::select(-date) |>
    tidyr::pivot_longer(-year, names_to = "description", values_to = "population") |>
    dplyr::mutate(
      population = as.numeric(population),
      sex = dplyr::case_when(
        stringr::str_detect(description, "Male")   ~ "males",
        stringr::str_detect(description, "Female") ~ "females",
        TRUE ~ NA_character_),
      age = stringr::str_extract(description, "(?<=;  )[0-9]+(?= ;)"),
      age = dplyr::if_else(stringr::str_detect(description, "100 and over"), "100", age),
      age = as.integer(age)) |>
    dplyr::filter(!is.na(sex), !is.na(age)) |>
    dplyr::select(year, sex, age, population)
}

load_single_year_pop <- function() {
  pop <- dplyr::bind_rows(
    .read_abs_single_year("data/raw/ABS/3101059.xlsx", 1982:hist_end),
    .read_abs_single_year("data/raw/ABS/3222_Table_B9.xlsx", proj_start:proj_end)
  ) |> dplyr::mutate(age_group = as.character(cut(
    age, breaks = c(band_lower - 1, Inf), labels = age_groups)))

  # The single years must add up to the banded denominators the model uses.
  banded <- dplyr::bind_rows(
    readr::read_csv("data/pop_hist.csv", show_col_types = FALSE),
    readr::read_csv("data/pop_proj.csv", show_col_types = FALSE)) |>
    dplyr::filter(year %in% 1982:proj_end)
  chk <- pop |>
    dplyr::group_by(year, sex, age_group) |>
    dplyr::summarise(p1 = sum(population), .groups = "drop") |>
    dplyr::full_join(banded, by = c("year", "sex", "age_group"))
  stopifnot(nrow(chk) == 2 * 9 * length(1982:proj_end),
            !anyNA(chk$p1), !anyNA(chk$population),
            max(abs(chk$p1 - chk$population)) < 1e-6)
  pop
}

# w[dx_year, age_group, t]: share of the band's diagnoses with x + t <= AGE_CUT
build_weights <- function(assumption = c("population", "lower", "midpoint", "ones"),
                          pop = NULL, sx = NULL) {
  assumption <- match.arg(assumption)
  years <- 1982:proj_end
  w <- array(1, dim = c(length(years), length(age_groups), 40),
             dimnames = list(year = as.character(years), age = age_groups,
                             t = as.character(0:39)))
  if (assumption == "ones") return(w)
  if (assumption %in% c("lower", "midpoint")) {
    x <- if (assumption == "lower") band_lower else age_mid[age_groups]
    for (g in age_groups) for (t in 0:39) w[, g, t + 1] <- as.numeric(x[[g]] + t <= AGE_CUT)
    return(w)
  }
  for (yr in years) for (g in age_groups) {
    d <- pop[pop$year == yr & pop$sex == sx & pop$age_group == g, ]
    share <- d$population / sum(d$population)
    for (t in 0:39) w[as.character(yr), g, t + 1] <- sum(share[d$age + t <= AGE_CUT])
  }
  stopifnot(all(w >= 0), all(w <= 1 + 1e-12))
  w
}

# ----------------------
# Weighted analogue of compute_prevalence(): returns the year x draw matrix
# ----------------------
compute_prevalence_w <- function(cases_array, surv_array, calc_years, duration, weights) {
  B <- dim(cases_array)[3]
  stopifnot(dim(surv_array)[4] == B)
  min_dx_year <- max(min(as.numeric(dimnames(cases_array)[[1]])),
                     min(as.numeric(dimnames(surv_array)[[1]])))
  out <- matrix(0, nrow = length(calc_years), ncol = B,
                dimnames = list(as.character(calc_years), NULL))
  for (ci in seq_along(calc_years)) {
    for (t in 0:(duration - 1)) {
      dy <- calc_years[ci] - t
      if (dy < min_dx_year) next
      dy_str <- as.character(dy)
      cm <- cases_array[dy_str, , ]
      sm <- if (t == 0) sqrt(surv_array[dy_str, , "1", ])
            else surv_array[dy_str, , as.character(t), ]
      # weights[dy, , t] has one value per age band and recycles down the rows
      out[ci, ] <- out[ci, ] + colSums(cm * sm * weights[dy_str, , t + 1])
    }
  }
  out
}

# ----------------------
# Driver
# ----------------------
run_attained_age_100 <- function(B = 1000, seed = 20260507, save_dir = "output",
                                 report_years = c(2021, 2031, 2045),
                                 durations = c(5, 10, 40)) {
  t_start <- Sys.time()
  # Same call sequence as run_prev_model()
  set.seed(seed)
  apc_results <- readRDS(file.path(save_dir, "apc_results.rds"))
  inputs      <- load_prev_inputs()
  inc_back    <- run_back_estimation(inputs$inc_subtype, inputs$inc_agg,
                                     save_path = tempfile(fileext = ".csv"))
  surv_obs       <- build_surv_obs(inputs$surv_raw)
  improvement_pt <- build_improvement_final(surv_obs)
  surv_draws_tbl        <- sample_surv_draws(surv_obs, B = B, seed = seed)
  improvement_draws_tbl <- sample_improvement_draws(surv_draws_tbl, improvement_pt, B = B)

  pop <- load_single_year_pop()
  assumptions <- c("population", "lower", "midpoint")
  dkey <- function(st, sx, a, dur) paste(st, sx, a, dur, sep = "|")
  store <- list()
  ones_maxdiff <- 0

  for (st in prev_series) for (sx in sexes) {
    cat(sprintf("  %s / %s\n", st, sx))
    cases_array <- build_cases_array(st, sx, B,
                                     inc_subtype = inputs$inc_subtype,
                                     inc_agg     = inputs$inc_agg,
                                     inc_back    = inc_back,
                                     apc_results = apc_results)
    surv_array <- build_surv_array(st, sx, dx_years = 1982:proj_end,
                                   cap_year = improve_cap,
                                   surv_draws_tbl, improvement_draws_tbl, B = B)
    W <- list(ones       = build_weights("ones"),
              population = build_weights("population", pop, sx),
              lower      = build_weights("lower"),
              midpoint   = build_weights("midpoint"))
    for (dur in durations) {
      ref  <- compute_prevalence(cases_array, surv_array,
                                 calc_years = report_years, duration = dur)$draws
      ones <- compute_prevalence_w(cases_array, surv_array, report_years, dur, W$ones)
      # Check 1: all weights 1 must reproduce compute_prevalence()
      rel <- max(abs(ones / ref - 1))
      stopifnot(rel < 1e-12)
      ones_maxdiff <- max(ones_maxdiff, rel)
      store[[dkey(st, sx, "untruncated", dur)]] <- ref
      for (a in assumptions) {
        store[[dkey(st, sx, a, dur)]] <-
          compute_prevalence_w(cases_array, surv_array, report_years, dur, W[[a]])
      }
    }
    rm(cases_array, surv_array); gc(verbose = FALSE)
  }

  # Persons and the headline total (HL + aggregate NHL) at the draw level, as in
  # run_prev_model() step 4b. The subtypes are never summed into the total.
  versions <- c("untruncated", assumptions)
  for (a in versions) for (dur in durations) {
    for (st in prev_series) {
      store[[dkey(st, "persons", a, dur)]] <-
        store[[dkey(st, "males", a, dur)]] + store[[dkey(st, "females", a, dur)]]
    }
    for (sx in c(sexes, "persons")) {
      store[[dkey("total", sx, a, dur)]] <-
        store[[dkey("hodgkin", sx, a, dur)]] + store[[dkey("nhl", sx, a, dur)]]
    }
  }

  # ---- Check 2: untruncated medians against the published outputs ----
  pub <- readr::read_csv(file.path(save_dir, "prevalence_projections.csv"),
                         show_col_types = FALSE)
  worst <- 0
  for (st in c(prev_series, "total")) for (sx in c(sexes, "persons"))
    for (dur in durations) for (yr in report_years) {
      mine   <- median(store[[dkey(st, sx, "untruncated", dur)]][as.character(yr), ])
      theirs <- pub$prev_mid[pub$subtype == st & pub$sex == sx &
                             pub$duration == dur & pub$year == yr]
      stopifnot(length(theirs) == 1)
      worst <- max(worst, abs(mine / theirs - 1))
    }
  t3 <- readr::read_csv(file.path(save_dir, "table_3_prevalence_combined.csv"),
                        show_col_types = FALSE)
  t3_ok <- TRUE
  for (st in c(prev_series, "total")) for (dur in durations) {
    m   <- store[[dkey(st, "persons", "untruncated", dur)]]
    row <- t3[t3$subtype == st & t3$duration == dur, ]
    stopifnot(nrow(row) == 1)
    t3_ok <- t3_ok && round(median(m["2021", ])) == row$prev_2021 &&
                      round(median(m["2045", ])) == row$prev_2045
  }

  # ---- Output table ----
  q <- function(x, p) unname(stats::quantile(x, p))
  rows <- list()
  for (st in c("total", "nhl", "hodgkin", nhl_subtypes))
    for (sx in c("females", "males", "persons")) for (dur in durations)
      for (yr in report_years) {
        y  <- as.character(yr)
        m0 <- store[[dkey(st, sx, "untruncated", dur)]][y, ]
        r  <- tibble::tibble(subtype = st, sex = sx, duration = dur, year = yr,
                             untruncated = median(m0),
                             untruncated_p025 = q(m0, .025), untruncated_p975 = q(m0, .975))
        for (a in assumptions) {
          m <- store[[dkey(st, sx, a, dur)]][y, ]
          r[[paste0(a, "")]]            <- median(m)
          r[[paste0(a, "_p025")]]       <- q(m, .025)
          r[[paste0(a, "_p975")]]       <- q(m, .975)
          # change formed within each draw, then the median across draws
          r[[paste0(a, "_pct_change")]] <- median((m / m0 - 1) * 100)
        }
        rows[[length(rows) + 1]] <- r
      }
  res <- dplyr::bind_rows(rows)
  out <- res |>
    dplyr::mutate(dplyr::across(dplyr::ends_with("_pct_change"), ~ round(.x, 2)),
                  dplyr::across(-c(subtype, sex, duration, year,
                                   dplyr::ends_with("_pct_change")), ~ round(.x)))
  out_path <- file.path(save_dir, "diagnostic_attained_age_100.csv")
  readr::write_csv(out, out_path)

  # ---- 2021 40-year NHL and HL against the AIHW (as Table S9) ----
  aihw <- inputs$prev_agg |>
    dplyr::filter(cancer_group %in% c("Non-Hodgkin lymphoma", "Hodgkin lymphoma"),
                  year == 2021, duration == 40, sex %in% sexes) |>
    dplyr::mutate(subtype = ifelse(cancer_group == "Hodgkin lymphoma", "hodgkin", "nhl")) |>
    dplyr::group_by(subtype, sex) |>
    dplyr::summarise(aihw_40yr = sum(prevalence, na.rm = TRUE), .groups = "drop")
  val <- res |>
    dplyr::filter(duration == 40, year == 2021) |>
    dplyr::inner_join(aihw, by = c("subtype", "sex")) |>
    dplyr::transmute(subtype, sex, aihw_40yr,
                     untruncated_pct = untruncated / aihw_40yr * 100,
                     population_pct  = population  / aihw_40yr * 100,
                     lower_pct       = lower       / aihw_40yr * 100,
                     midpoint_pct    = midpoint    / aihw_40yr * 100)
  stopifnot(nrow(val) == 4)

  # ---- Console report ----
  op <- options(width = 200); on.exit(options(op))
  show <- function(d, digits = 1) {
    d <- as.data.frame(d); num <- vapply(d, is.numeric, logical(1))
    d[num] <- lapply(d[num], round, digits); print(d, row.names = FALSE)
  }
  cat(sprintf("\nCheck 1 passed: with all weights 1 the weighted function equals compute_prevalence() (max relative difference %.1e)\n", ones_maxdiff))
  cat(sprintf("Check 2: untruncated medians against prevalence_projections.csv, max relative difference %.1e (%s); table_3 2021 and 2045 values %s\n",
              worst, if (worst < 1e-9) "exact match" else "NOT an exact match: RNG alignment",
              if (t3_ok) "match" else "DO NOT match"))
  cat("\nPopulation is a proxy for the age distribution of diagnoses within each band.\n")
  cat("Whether it over- or understates the truncation in the 85+ band is unknown, since\n")
  cat("NHL rates within that band are not observed. The lower-bound assumption gives\n")
  cat("the least truncation; the population and midpoint assumptions give similar\n")
  cat("results at 40 years.\n")
  cat("\n=== % change from untruncated (median of within-draw change), total and aggregate NHL, persons ===\n")
  show(res |>
         dplyr::filter(subtype %in% c("total", "nhl"), sex == "persons") |>
         dplyr::select(subtype, duration, year, untruncated, population,
                       population_pct_change, lower_pct_change, midpoint_pct_change) |>
         dplyr::mutate(untruncated = round(untruncated), population = round(population)), 2)
  cat("\n=== 2021 40-year prevalence as % of the AIHW ===\n")
  show(val)
  big <- res$population_pct_change[res$subtype == "total" & res$sex == "persons" &
                                   res$duration == 40 & res$year == 2045]
  cat(sprintf("\n2045 40-year total, population assumption: %.2f%% %s\n", big,
              if (abs(big) > 3) "** EXCEEDS 3% **" else "(within 3%)"))
  cat(sprintf("Wrote %s (%d rows). Runtime %.0f s.\n", out_path, nrow(out),
              as.numeric(Sys.time() - t_start, units = "secs")))

  invisible(list(results = res, validation = val))
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_attained_age_100())
  invisible(grDevices::dev.off())
}
