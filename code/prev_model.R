# ================
# Blood Cancer Projections (lymphoma): Prevalence Projection (MC-only)
# ================
# Projects 2-, 3-, 5-, 10-, and 40-year prevalence to 2045 with 95% credible
# intervals from parameter-level Monte Carlo simulation, in a two-tier
# structure that is consistent across incidence and prevalence:
#   Headline tier      - Hodgkin lymphoma (HL) and aggregate NHL. These two
#                        are disjoint and together are ALL lymphoma, so the
#                        headline TOTAL lymphoma = HL + NHL (built at the draw
#                        level, step 4b: total_b = HL_b + NHL_b, then median +
#                        percentiles; HL and NHL are independent APC fits).
#   Decomposition tier - DLBCL, FL, MCL, a SUBSET of aggregate NHL (they cover
#                        ~53-73% of it, depending on year/duration). They
#                        decompose NHL and are NEVER summed into the total
#                        (that would double-count). In Figures 2/3 they are
#                        drawn dashed, sitting within the solid NHL line.
# Aggregate NHL is a headline result. Its validation against AIHW NHL
# prevalence (QL #3) is retained but reclassified as validation OF a headline
# result, not the reason the series exists.
#
# Method: modified counting method
#   Prevalence(Y) = sum_{t=0..T-1} Incidence(Y-t) x S(t)
#   with S(0) := sqrt(S(1)) (constant-hazard half-year approximation)
#
# Library + driver pattern:
#   source("code/prev_model.R")    -> defines functions; no side effects
#   run_prev_model()               -> runs full pipeline
#   Rscript code/prev_model.R      -> auto-runs via the bottom guard
# ================

set.seed(20260507)
stopifnot(
  requireNamespace("dplyr"),
  requireNamespace("tidyr"),
  requireNamespace("readr"),
  requireNamespace("purrr"),
  requireNamespace("ggplot2"),
  requireNamespace("cowplot"),
  requireNamespace("grid"),
  requireNamespace("gtable")
)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

# Shared constants, colours, grob makers, save_fig. See code/_setup.R.
source("code/_setup.R")

# ----------------------
# Data loaders
# ----------------------

load_prev_inputs <- function(data_dir = "data") {
  list(
    surv_raw     = read_csv(file.path(data_dir, "survival.csv"),       show_col_types = FALSE),
    inc_subtype  = read_csv(file.path(data_dir, "incidence_subtype.csv"), show_col_types = FALSE),
    inc_agg      = read_csv(file.path(data_dir, "incidence_agg.csv"),  show_col_types = FALSE),
    prev_subtype = read_csv(file.path(data_dir, "prevalence_subtype.csv"), show_col_types = FALSE),
    prev_agg     = read_csv(file.path(data_dir, "prevalence_agg.csv"), show_col_types = FALSE)
  )
}

# ----------------------
# Survival cleaning (with explicit donor_age_group)
# ----------------------
# Replaces float-equality donor lookup in sample_surv_draws() with an
# explicit donor identity. Filled cells (n.p. in AIHW) get a
# `donor_age_group` column pointing to the cell whose values they share.

.fill_nearest <- function(df) {
  if (!any(is.na(df$survival_pct))) {
    return(df |> mutate(was_filled = FALSE,
                        donor_age_group = age_group))
  }
  df <- df |> mutate(age_pos = match(age_group, age_groups))
  available <- df |> filter(!is.na(survival_pct)) |>
    mutate(was_filled = FALSE, donor_age_group = age_group)
  missing   <- df |> filter(is.na(survival_pct)) |>
    mutate(was_filled = TRUE, donor_age_group = NA_character_)
  if (nrow(available) == 0) {
    return(df |> mutate(was_filled = TRUE,
                        donor_age_group = NA_character_) |>
             select(-age_pos))
  }
  for (i in seq_len(nrow(missing))) {
    distances <- abs(missing$age_pos[i] - available$age_pos)
    nearest   <- available |> slice(which.min(distances))
    missing$survival_pct[i]    <- nearest$survival_pct
    missing$ci_lower[i]        <- nearest$ci_lower
    missing$ci_upper[i]        <- nearest$ci_upper
    missing$donor_age_group[i] <- nearest$age_group
  }
  bind_rows(available, missing) |> arrange(age_pos) |> select(-age_pos)
}

build_surv_obs <- function(surv_raw) {
  surv_obs <- surv_raw |>
    filter(survival_type == "Observed",
           sex != "persons",
           period != "2007–2021") |>
    mutate(
      survival_pct = as.numeric(survival_pct),
      ci_lower     = as.numeric(ci_lower),
      ci_upper     = as.numeric(ci_upper)
    ) |>
    select(period, sex, age_group, years_since_dx, subtype,
           survival_pct, ci_lower, ci_upper) |>
    group_by(subtype, sex, period, years_since_dx) |>
    group_modify(~ .fill_nearest(.x)) |>
    ungroup()
  stopifnot(sum(is.na(surv_obs$survival_pct)) == 0)
  surv_obs
}

# ----------------------
# Improvement rates (point estimate)
# ----------------------
# Used by the source classification (held fixed across draws), and
# directly by the point-estimate survival builder below.

build_improvement_final <- function(surv_obs) {
  filled_lookup <- surv_obs |>
    filter(period == "2017–2021", years_since_dx == 5) |>
    select(subtype, sex, age_group, was_filled) |>
    distinct()

  surv_wide <- surv_obs |>
    filter(!was_filled) |>
    mutate(period_num = period_map[period]) |>
    select(period_num, sex, age_group, years_since_dx, subtype, survival_pct) |>
    pivot_wider(names_from = period_num,
                values_from = survival_pct,
                names_prefix = "s") |>
    mutate(annual_rate = (s3 - s1) / 10)

  mono_check <- surv_wide |>
    filter(years_since_dx == 5) |>
    mutate(is_improving = (s3 > s1)) |>
    select(sex, age_group, subtype, is_improving)

  age_rates <- surv_wide |>
    select(sex, age_group, years_since_dx, subtype, annual_rate) |>
    left_join(mono_check, by = c("sex", "age_group", "subtype"))

  surv_fallback <- age_rates |>
    filter(is_improving) |>
    group_by(sex, subtype, years_since_dx) |>
    summarise(fallback_rate = mean(annual_rate), .groups = "drop")

  full_grid <- expand_grid(
    sex = sexes, age_group = age_groups,
    subtype = prev_series, years_since_dx = 1:5
  )

  full_grid |>
    left_join(age_rates,     by = c("sex", "age_group", "subtype", "years_since_dx")) |>
    left_join(surv_fallback, by = c("sex", "subtype", "years_since_dx")) |>
    left_join(filled_lookup, by = c("sex", "age_group", "subtype")) |>
    mutate(
      fallback_rate = replace_na(fallback_rate, 0),
      source = case_when(
        was_filled                       ~ "fallback (filled cell)",
        is.na(is_improving)              ~ "fallback (no data)",
        is_improving                     ~ "age-specific",
        TRUE                             ~ "fallback (non-improving)"
      ),
      annual_rate = if_else(source == "age-specific", annual_rate, fallback_rate)
    ) |>
    select(sex, age_group, years_since_dx, subtype, annual_rate, source)
}

# ----------------------
# Point-estimate survival projection (used only by Figure S2)
# ----------------------
# Mirrors the original prev_model.R::build_projected_surv() exactly.
# Kept here so supplement.R can build Figure S2 without rebuilding the
# MC array (which would be wasteful for one curve).

assign_period <- function(dx_year) {
  dplyr::case_when(dx_year <= 2011 ~ "2007–2011",
                   dx_year <= 2016 ~ "2012–2016",
                   TRUE            ~ "2017–2021")
}

build_projected_surv_pt <- function(subtype_name, sex_name, dx_year,
                                    cap_year, surv_obs, improvement_final) {
  dx_period <- if (dx_year <= 2021) assign_period(dx_year) else "2017–2021"
  midpoint  <- as.numeric(period_midpoints[dx_period])
  yrs_from_mid <- if (dx_year <= 2021) dx_year - midpoint
                  else min(dx_year - midpoint, cap_year - midpoint)

  base <- surv_obs |>
    filter(subtype == subtype_name, sex == sex_name, period == dx_period) |>
    select(age_group, years_since_dx,
           surv_mid = survival_pct, surv_lo = ci_lower, surv_hi = ci_upper)

  if (yrs_from_mid != 0) {
    rates <- improvement_final |>
      filter(subtype == subtype_name, sex == sex_name) |>
      select(age_group, years_since_dx, annual_rate)
    base <- base |>
      left_join(rates, by = c("age_group", "years_since_dx")) |>
      mutate(
        surv_mid = pmin(pmax(surv_mid + annual_rate * yrs_from_mid, 0), 1),
        surv_lo  = pmin(pmax(surv_lo  + annual_rate * yrs_from_mid, 0), 1),
        surv_hi  = pmin(pmax(surv_hi  + annual_rate * yrs_from_mid, 0), 1)
      ) |> select(-annual_rate)
  }

  anchor <- base |>
    filter(years_since_dx %in% c(3, 5)) |>
    pivot_wider(names_from = years_since_dx,
                values_from = c(surv_mid, surv_lo, surv_hi)) |>
    mutate(
      cond_mid = if_else(surv_mid_3 > 0, pmin(sqrt(surv_mid_5 / surv_mid_3), 1), 0),
      cond_lo  = if_else(surv_lo_3  > 0, pmin(sqrt(surv_lo_5  / surv_lo_3 ), 1), 0),
      cond_hi  = if_else(surv_hi_3  > 0, pmin(sqrt(surv_hi_5  / surv_hi_3 ), 1), 0)
    )

  extrap <- anchor |>
    crossing(years_since_dx = 6:40) |>
    mutate(
      surv_mid = surv_mid_5 * cond_mid^(years_since_dx - 5),
      surv_lo  = surv_lo_5  * cond_lo^(years_since_dx - 5),
      surv_hi  = surv_hi_5  * cond_hi^(years_since_dx - 5)
    ) |>
    mutate(across(starts_with("surv_"), ~ pmax(.x, 0))) |>
    select(age_group, years_since_dx, surv_mid, surv_lo, surv_hi)

  bind_rows(base |> filter(years_since_dx <= 5), extrap) |>
    arrange(age_group, years_since_dx)
}

# ----------------------
# MC sampling: per-cell survival draws (logit-normal)
# ----------------------
# Filled cells share their donor's draw vector via reference equality
# (same list element). Donor identity comes from the `donor_age_group`
# column added by .fill_nearest() above.

sample_surv_draws <- function(surv_obs, B = 1000, seed = 20260507) {
  set.seed(seed)
  obs <- surv_obs |> mutate(.row = row_number())

  # Resolve donor row index: by name (donor_age_group) within slice.
  donor_row <- integer(nrow(obs))
  for (i in seq_len(nrow(obs))) {
    if (!obs$was_filled[i]) {
      donor_row[i] <- obs$.row[i]
    } else {
      cand <- which(
        obs$subtype        == obs$subtype[i]        &
        obs$sex            == obs$sex[i]            &
        obs$period         == obs$period[i]         &
        obs$years_since_dx == obs$years_since_dx[i] &
        obs$age_group      == obs$donor_age_group[i]
      )
      if (length(cand) != 1L) {
        stop("Filled cell ", i, " has ", length(cand),
             " donor candidates (expected 1)")
      }
      donor_row[i] <- cand
    }
  }

  unique_donors <- sort(unique(donor_row))
  draws_by_row  <- vector("list", nrow(obs))
  for (k in unique_donors) {
    p  <- pmin(pmax(obs$survival_pct[k], SURV_LO), SURV_HI)
    lo <- pmin(pmax(obs$ci_lower[k],     SURV_LO), SURV_HI)
    hi <- pmin(pmax(obs$ci_upper[k],     SURV_LO), SURV_HI)
    mu <- qlogis(p)
    sigma <- (qlogis(hi) - qlogis(lo)) / (2 * 1.96)
    if (!is.finite(sigma) || sigma < 0) sigma <- 0
    eta <- if (sigma == 0) rep(mu, B) else rnorm(B, mean = mu, sd = sigma)
    draws_by_row[[k]] <- plogis(eta)
  }
  for (i in seq_len(nrow(obs))) {
    if (obs$was_filled[i]) draws_by_row[[i]] <- draws_by_row[[donor_row[i]]]
  }

  obs |> mutate(surv_draws = draws_by_row) |> select(-.row)
}

# ----------------------
# Improvement rate draws
# ----------------------
# Per-draw rate = (S3_b - S1_b) / 10 for "age-specific" cells; for
# fallback cells, per-draw rate is the mean across improving cells in
# the (subtype, sex, ysd) group. Source classification held fixed at
# point estimate (simpler-alternative path the plan permits).

sample_improvement_draws <- function(surv_draws_tbl, improvement_final, B = 1000) {
  s1 <- surv_draws_tbl |>
    filter(period == "2007–2011") |>
    select(subtype, sex, age_group, years_since_dx, s1_draws = surv_draws)
  s3 <- surv_draws_tbl |>
    filter(period == "2017–2021") |>
    select(subtype, sex, age_group, years_since_dx, s3_draws = surv_draws)
  paired <- s1 |> inner_join(s3, by = c("subtype", "sex", "age_group", "years_since_dx"))
  paired$rate_age <- mapply(function(a, b) (b - a) / 10, paired$s1_draws, paired$s3_draws,
                            SIMPLIFY = FALSE)

  is_improving_pt <- improvement_final |>
    filter(years_since_dx == 5) |>
    mutate(is_improving = source == "age-specific") |>
    select(subtype, sex, age_group, is_improving)

  paired_improving <- paired |>
    inner_join(is_improving_pt, by = c("subtype", "sex", "age_group")) |>
    filter(is_improving)

  fallback_lookup <- paired_improving |>
    group_by(subtype, sex, years_since_dx) |>
    summarise(
      fallback_draws = list({
        m <- do.call(rbind, rate_age)
        if (is.null(m)) rep(0, B) else colMeans(m)
      }),
      .groups = "drop"
    )

  out <- improvement_final |>
    left_join(paired |> select(subtype, sex, age_group, years_since_dx, rate_age),
              by = c("subtype", "sex", "age_group", "years_since_dx")) |>
    left_join(fallback_lookup, by = c("subtype", "sex", "years_since_dx"))

  rate_draws <- vector("list", nrow(out))
  for (i in seq_len(nrow(out))) {
    if (out$source[i] == "age-specific" && !is.null(out$rate_age[[i]])) {
      rate_draws[[i]] <- out$rate_age[[i]]
    } else {
      fb <- out$fallback_draws[[i]]
      rate_draws[[i]] <- if (is.null(fb)) rep(0, B) else fb
    }
  }
  out |> mutate(rate_draws = rate_draws) |>
    select(subtype, sex, age_group, years_since_dx, source, rate_draws)
}

# ----------------------
# Build per-(subtype, sex) 4D survival array
# ----------------------
# dim: (length(dx_years), 9, 40, B); axes (dx_year, age, ysd, draw).
# Per-draw flooring at 0 and clipping at 1; long-tail extrapolation via
# cond^(ysd-5) with cond capped at 1 per draw.

build_surv_array <- function(subtype_name, sex_name,
                             dx_years, cap_year = improve_cap,
                             surv_draws_tbl, improvement_draws_tbl,
                             B = NULL) {
  base_tbl <- surv_draws_tbl |>
    filter(subtype == subtype_name, sex == sex_name) |>
    arrange(match(period, period_levels),
            match(age_group, age_groups),
            years_since_dx)

  if (is.null(B)) B <- length(base_tbl$surv_draws[[1]])

  S_base <- array(NA_real_, dim = c(3, 9, 5, B),
                  dimnames = list(period = period_levels,
                                  age    = age_groups,
                                  ysd    = as.character(1:5),
                                  draw   = NULL))
  for (i in seq_len(nrow(base_tbl))) {
    p <- match(base_tbl$period[i], period_levels)
    a <- match(base_tbl$age_group[i], age_groups)
    y <- base_tbl$years_since_dx[i]
    S_base[p, a, y, ] <- base_tbl$surv_draws[[i]]
  }
  stopifnot(!any(is.na(S_base)))

  rate_tbl <- improvement_draws_tbl |>
    filter(subtype == subtype_name, sex == sex_name) |>
    arrange(match(age_group, age_groups), years_since_dx)
  rates <- array(NA_real_, dim = c(9, 5, B),
                 dimnames = list(age = age_groups,
                                 ysd = as.character(1:5),
                                 draw = NULL))
  for (i in seq_len(nrow(rate_tbl))) {
    a <- match(rate_tbl$age_group[i], age_groups)
    y <- rate_tbl$years_since_dx[i]
    rates[a, y, ] <- rate_tbl$rate_draws[[i]]
  }
  stopifnot(!any(is.na(rates)))

  n_dx <- length(dx_years)
  out <- array(NA_real_, dim = c(n_dx, 9, 40, B),
               dimnames = list(dx_year = as.character(dx_years),
                               age     = age_groups,
                               ysd     = as.character(1:40),
                               draw    = NULL))

  for (k in seq_len(n_dx)) {
    dx <- dx_years[k]
    if (dx <= 2021) {
      dx_period <- if (dx <= 2011) "2007–2011"
                   else if (dx <= 2016) "2012–2016"
                   else "2017–2021"
    } else {
      dx_period <- "2017–2021"
    }
    midpoint <- period_midpoints[dx_period]
    yrs_from_mid <- if (dx <= 2021) {
      dx - as.numeric(midpoint)
    } else {
      min(dx - as.numeric(midpoint), cap_year - as.numeric(midpoint))
    }
    period_idx <- match(dx_period, period_levels)

    S_adj <- S_base[period_idx, , , ] + rates * yrs_from_mid
    S_adj[S_adj < 0] <- 0
    S_adj[S_adj > 1] <- 1

    out[k, , 1:5, ] <- S_adj

    S_3 <- S_adj[, 3, , drop = FALSE]
    S_5 <- S_adj[, 5, , drop = FALSE]
    S_3_floored <- pmax(S_3, SURV_FLOOR)
    cond <- sqrt(S_5 / S_3_floored)
    cond[cond > 1] <- 1
    cond[!is.finite(cond)] <- 0

    S_5m  <- array(S_5,  dim = c(9, B))
    condm <- array(cond, dim = c(9, B))
    for (yy in 6:40) {
      vals <- S_5m * condm^(yy - 5)
      vals[vals < 0] <- 0
      out[k, , yy, ] <- vals
    }
  }
  out
}

# ----------------------
# Cases array: historical broadcast + MC projection
# ----------------------

build_cases_array <- function(subtype_name, sex_name, B,
                              inc_subtype, inc_agg, inc_back,
                              apc_results) {
  years <- 1982:2045
  n_age <- length(age_groups)
  n_yr  <- length(years)

  if (subtype_name %in% nhl_subtypes) {
    early <- inc_back |>
      filter(year <= 2002, subtype == subtype_name, sex == sex_name) |>
      select(year, age_group, count)
    late <- inc_subtype |>
      filter(year >= 2003 & year <= 2021,
             subtype == subtype_name, sex == sex_name) |>
      select(year, age_group, count)
    hist_cases <- bind_rows(early, late)
  } else if (subtype_name == "hodgkin") {
    hist_cases <- inc_agg |>
      filter(cancer_group == "Hodgkin lymphoma", year <= 2021,
             sex == sex_name) |>
      select(year, age_group, count)
  } else if (subtype_name == "nhl") {
    # Aggregate NHL: observed 1982-2021 directly (no back-estimation needed,
    # unlike the NHL subtypes); projection from the aggregate-tier NHL APC fit.
    hist_cases <- inc_agg |>
      filter(cancer_group == "Non-Hodgkin lymphoma", year <= 2021,
             sex == sex_name) |>
      select(year, age_group, count)
  } else {
    stop("unknown subtype: ", subtype_name)
  }

  key <- paste0(subtype_name, "_", sex_name)
  proj_arr <- apc_results[[key]]$cases_draws
  if (is.null(proj_arr)) {
    stop("apc_results has no cases_draws for key: ", key)
  }
  if (dim(proj_arr)[3] != B) {
    stop("B mismatch: cases_draws has ", dim(proj_arr)[3], " draws, expected ", B)
  }

  cases_array <- array(0, dim = c(n_yr, n_age, B),
                       dimnames = list(year = as.character(years),
                                       age = age_groups,
                                       draw = NULL))
  for (i in seq_len(nrow(hist_cases))) {
    yr <- as.character(hist_cases$year[i])
    ag <- hist_cases$age_group[i]
    cases_array[yr, ag, ] <- hist_cases$count[i]
  }
  proj_years <- as.character(2022:2045)
  cases_array[proj_years, , ] <- proj_arr[proj_years, age_groups, ]

  cases_array
}

# ----------------------
# Per-draw prevalence calculation
# ----------------------
# Reported prev_mid is the median (matches a deterministic point
# estimate under the log/logit links). prev_mean is also returned.
# Periodic in-loop gc() keeps transient working set under ~60 MB at
# B = 1000 (slice copies in R accumulate otherwise).

compute_prevalence <- function(cases_array, surv_array,
                                calc_years, duration) {
  B <- dim(cases_array)[3]
  if (dim(surv_array)[4] != B) {
    stop("B mismatch: cases has ", B, " draws, surv has ", dim(surv_array)[4])
  }
  n_age <- dim(cases_array)[2]
  cases_years <- as.numeric(dimnames(cases_array)[[1]])
  surv_years  <- as.numeric(dimnames(surv_array)[[1]])
  min_dx_year <- max(min(cases_years), min(surv_years))

  prev_draws <- matrix(0, nrow = length(calc_years), ncol = B)
  for (ci in seq_along(calc_years)) {
    cy <- calc_years[ci]
    for (t in 0:(duration - 1)) {
      dy <- cy - t
      if (dy < min_dx_year) next
      dy_str <- as.character(dy)
      cm <- cases_array[dy_str, , ]
      sm <- if (t == 0) sqrt(surv_array[dy_str, , "1", ])
            else surv_array[dy_str, , as.character(t), ]
      prev_draws[ci, ] <- prev_draws[ci, ] + colSums(cm * sm)
    }
    if (ci %% 5 == 0) gc(verbose = FALSE)
  }

  s <- summarise_draws(prev_draws, 1)
  rownames(prev_draws) <- as.character(calc_years)
  # Return BOTH the per-year summary and the (year x draw) matrix. The
  # driver keeps the draws so the headline total (HL + NHL) and the
  # persons (both-sex) rows can be formed at the DRAW level: total_b =
  # HL_b + NHL_b per draw, then median + 2.5/97.5 percentiles. HL and NHL
  # come from independent APC fits, so draw-wise summation propagates
  # their variance correctly (they are paired only by draw index).
  list(
    summary = tibble::tibble(
      year      = calc_years,
      prev_mid  = s$mid,
      prev_p025 = s$p025,
      prev_p975 = s$p975,
      prev_mean = s$mean
    ),
    draws = prev_draws
  )
}

# ----------------------
# NHL subtype back-estimation (1982-2002)
# ----------------------
# Sources: 1982-2002 from aggregate NHL via 2003-2007 subtype
# proportions; HL used directly. Idempotent - re-running overwrites
# data/incidence_subtype_back.csv with deterministic content.

run_back_estimation <- function(inc_subtype, inc_agg,
                                save_path = "data/incidence_subtype_back.csv",
                                prop_years = 2003:2007,
                                nhl_subtypes_arg = nhl_subtypes) {
  subtype_avg <- inc_subtype |>
    filter(subtype %in% nhl_subtypes_arg, year %in% prop_years) |>
    group_by(sex, age_group, subtype) |>
    summarise(avg_count = mean(count, na.rm = TRUE), .groups = "drop")

  nhl_total_avg <- inc_agg |>
    filter(cancer_group == "Non-Hodgkin lymphoma", year %in% prop_years) |>
    group_by(sex, age_group) |>
    summarise(avg_nhl = mean(count, na.rm = TRUE), .groups = "drop")

  subtype_props <- subtype_avg |>
    left_join(nhl_total_avg, by = c("sex", "age_group")) |>
    mutate(raw_prop = avg_count / avg_nhl)

  nhl_agg_hist <- inc_agg |>
    filter(cancer_group == "Non-Hodgkin lymphoma", year <= 2002)

  nhl_back <- subtype_props |>
    select(sex, age_group, subtype, raw_prop) |>
    left_join(nhl_agg_hist,
              by = c("sex", "age_group"),
              relationship = "many-to-many") |>
    mutate(count = round(count * raw_prop)) |>
    select(year, sex, age_group, subtype, count)

  hodgkin_back <- inc_agg |>
    filter(cancer_group == "Hodgkin lymphoma", year <= 2002) |>
    mutate(subtype = "hodgkin") |>
    select(year, sex, age_group, subtype, count)

  inc_back <- bind_rows(nhl_back, hodgkin_back) |>
    arrange(year, sex, age_group, subtype)

  write_csv(inc_back, save_path)
  invisible(inc_back)
}

# ----------------------
# Figure 2 / 3: prevalence projections (with MC ribbons)
# ----------------------
# Uses shared make_line_grob / make_point_grob / make_text_grob /
# make_spacer / line_colours / label_map / sex_labels from _setup.R.

# The legend is the shared, box-free lymphoma legend defined in _setup.R
# (make_lymphoma_legend / attach_lymphoma_legend), used by Figures 1, 2, 3 and
# SI S3c-S3e so they all render identically.

build_prev_fig <- function(prev_df, duration_value,
                           prev_subtype_aihw = NULL,
                           prev_agg_aihw     = NULL) {
  # Headline tier (HL, aggregate NHL; solid) plus the NHL decomposition
  # (DLBCL/FL/MCL; dashed) so the subtypes read as sitting WITHIN NHL.
  draw_series <- c("hodgkin", "nhl", nhl_subtypes)
  lvl <- c("HL", "NHL", "DLBCL", "FL", "MCL")
  plot_data <- prev_df |>
    filter(duration == duration_value, year >= 2012,
           subtype %in% draw_series, sex %in% sexes) |>
    mutate(
      label    = factor(label_map[subtype], levels = lvl),
      tier_lty = if_else(subtype %in% c("hodgkin", "nhl"), "headline", "decomposition"),
      prev_mid = prev_mid / 1000,
      prev_lo  = prev_p025 / 1000,
      prev_hi  = prev_p975 / 1000
    )

  aihw_plot <- NULL
  include_aihw <- duration_value %in% c(5, 10) && !is.null(prev_agg_aihw)
  if (include_aihw) {
    aihw_agg <- prev_agg_aihw |>
      filter(duration == duration_value,
             cancer_group %in% c("Hodgkin lymphoma", "Non-Hodgkin lymphoma")) |>
      mutate(subtype = if_else(cancer_group == "Hodgkin lymphoma", "hodgkin", "nhl")) |>
      group_by(year, subtype, sex) |>
      summarise(prev_aihw = sum(prevalence, na.rm = TRUE), .groups = "drop")
    aihw_sub <- if (!is.null(prev_subtype_aihw)) {
      prev_subtype_aihw |>
        filter(duration == duration_value, subtype %in% nhl_subtypes) |>
        group_by(year, subtype, sex) |>
        summarise(prev_aihw = sum(prevalence, na.rm = TRUE), .groups = "drop")
    } else NULL
    aihw_plot <- bind_rows(aihw_agg, aihw_sub) |>
      mutate(label = factor(label_map[subtype], levels = lvl),
             prev_aihw = prev_aihw / 1000)
  }

  p <- ggplot(plot_data, aes(x = year, colour = label, fill = label)) +
    geom_ribbon(
      data = ~ filter(.x, year > hist_end),
      aes(ymin = prev_lo, ymax = prev_hi),
      alpha = 0.15, colour = NA
    ) +
    geom_line(data = ~ filter(.x, tier_lty == "decomposition"),
              aes(y = prev_mid), linewidth = unname(lymphoma_lwd["decomposition"])) +
    geom_line(data = ~ filter(.x, tier_lty == "headline"),
              aes(y = prev_mid), linewidth = unname(lymphoma_lwd["headline"])) +
    geom_vline(xintercept = hist_end + 0.5, linetype = "dotted",
               colour = "grey50", linewidth = 0.4) +
    facet_wrap(~ sex, ncol = 2, labeller = labeller(sex = sex_labels)) +
    scale_colour_manual(values = line_colours) +
    scale_fill_manual(values = line_colours) +
    labs(x = "Year", y = "Prevalence (thousands)") +
    theme_bw(base_size = 13) +
    theme(legend.position  = "none",
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(colour = NA, fill = NA))

  if (include_aihw) {
    p <- p + geom_point(
      data = aihw_plot,
      aes(x = year, y = prev_aihw, colour = label),
      size = 1.2, shape = 16, show.legend = FALSE, inherit.aes = FALSE
    )
  }

  attach_lymphoma_legend(p)
}

# ----------------------
# Driver
# ----------------------

run_prev_model <- function(B = 1000,
                           save_dir = "output",
                           seed = 20260507,
                           verbose = TRUE) {
  set.seed(seed)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  apc_results_path <- file.path(save_dir, "apc_results.rds")
  if (!file.exists(apc_results_path)) {
    stop(apc_results_path, " missing - run apc_model.R first")
  }
  apc_results <- readRDS(apc_results_path)

  inputs <- load_prev_inputs()

  # 1. Back-estimation (writes data/incidence_subtype_back.csv)
  inc_back <- run_back_estimation(inputs$inc_subtype, inputs$inc_agg)

  # 2. Survival cleaning + improvement rates (point estimate)
  surv_obs       <- build_surv_obs(inputs$surv_raw)
  improvement_pt <- build_improvement_final(surv_obs)

  # 3. Survival + improvement draws (sampled once, shared across chunks)
  surv_draws_tbl        <- sample_surv_draws(surv_obs, B = B, seed = seed)
  improvement_draws_tbl <- sample_improvement_draws(surv_draws_tbl,
                                                    improvement_pt, B = B)

  # 4. Loop over (subtype, sex, scenario, duration)
  scenarios  <- c(base = improve_cap, conservative = 2021, optimistic = 2036)
  durations  <- c(2, 3, 5, 10, 40)
  calc_years <- 2012:proj_end

  results <- list()
  # Per-draw prevalence matrices (year x B), keyed "series|sex|scenario|dur".
  # Retained so the headline TOTAL (HL + NHL) and the persons (both-sex)
  # rows are built at the draw level (see step 4b). ~41 MB at B = 1000.
  draw_store <- list()
  dkey <- function(st, sx, sc, dur) paste(st, sx, sc, dur, sep = "|")

  total_t0 <- Sys.time()
  for (st in prev_series) for (sx in sexes) {
    chunk_t0 <- Sys.time()
    if (verbose) cat(sprintf("\n=== %s / %s ===\n", st, sx))

    cases_array <- build_cases_array(st, sx, B,
                                     inc_subtype = inputs$inc_subtype,
                                     inc_agg     = inputs$inc_agg,
                                     inc_back    = inc_back,
                                     apc_results = apc_results)

    for (sc_name in names(scenarios)) {
      cap <- scenarios[[sc_name]]
      surv_array <- build_surv_array(st, sx,
                                     dx_years = 1982:proj_end,
                                     cap_year = cap,
                                     surv_draws_tbl,
                                     improvement_draws_tbl,
                                     B = B)

      for (dur in durations) {
        cp <- compute_prevalence(cases_array, surv_array,
                                 calc_years = calc_years, duration = dur)
        prev <- cp$summary
        prev$subtype  <- st
        prev$sex      <- sx
        prev$duration <- dur
        prev$scenario <- sc_name
        results[[length(results) + 1]] <- prev
        draw_store[[dkey(st, sx, sc_name, dur)]] <- cp$draws
      }
      rm(surv_array); gc(verbose = FALSE)
    }
    rm(cases_array); gc(verbose = FALSE)
    if (verbose) cat(sprintf("  done in %.1f s\n",
                              as.numeric(Sys.time() - chunk_t0, units = "secs")))
  }
  total_elapsed <- as.numeric(Sys.time() - total_t0, units = "secs")
  if (verbose) cat(sprintf("\nTotal: %.1f s\n", total_elapsed))

  # ------------------------------------------------------------------
  # 4b. DRAW-LEVEL persons (both-sex) rows and the headline TOTAL lymphoma
  #     (HL + aggregate NHL). HL and aggregate NHL are disjoint and together
  #     are all lymphoma; DLBCL/FL/MCL are a SUBSET of NHL and are NEVER
  #     summed into the total (that would double-count). Each combination is
  #     formed by adding the relevant per-draw matrices element-wise (paired
  #     by draw index), then taking the median and 2.5/97.5 percentiles.
  # ------------------------------------------------------------------
  summarise_draw_mat <- function(m, st_label, sex_label, sc_name, dur) {
    s <- summarise_draws(m, 1)
    tibble::tibble(
      year = as.numeric(rownames(m)),
      prev_mid = s$mid, prev_p025 = s$p025, prev_p975 = s$p975, prev_mean = s$mean,
      duration = dur, subtype = st_label, sex = sex_label, scenario = sc_name
    )
  }
  extra <- list()
  for (sc_name in names(scenarios)) for (dur in durations) {
    # persons (both sexes) for every series, incl. aggregate NHL
    persons_mat <- list()
    for (st in prev_series) {
      pm <- draw_store[[dkey(st, "males",   sc_name, dur)]] +
            draw_store[[dkey(st, "females", sc_name, dur)]]
      persons_mat[[st]] <- pm
      extra[[length(extra) + 1]] <- summarise_draw_mat(pm, st, "persons", sc_name, dur)
    }
    # Total lymphoma = HL + NHL, per draw, by sex and persons
    for (sx in sexes) {
      tot_sx <- draw_store[[dkey("hodgkin", sx, sc_name, dur)]] +
                draw_store[[dkey("nhl",     sx, sc_name, dur)]]
      extra[[length(extra) + 1]] <- summarise_draw_mat(tot_sx, "total", sx, sc_name, dur)
    }
    tot_persons <- persons_mat[["hodgkin"]] + persons_mat[["nhl"]]
    extra[[length(extra) + 1]] <- summarise_draw_mat(tot_persons, "total", "persons", sc_name, dur)
  }

  prev_full <- bind_rows(results, extra)

  # 5. Mirror analytical pipeline outputs (drop _mc suffix)
  prev_base <- prev_full |>
    filter(scenario == "base") |>
    select(year, prev_mid, prev_p025, prev_p975, prev_mean,
           duration, subtype, sex)

  prev_sensitivity <- prev_full |>
    select(scenario, year, subtype, sex, duration,
           prev_mid, prev_p025, prev_p975, prev_mean)

  baseline_2021 <- prev_base |>
    filter(year == 2021) |>
    select(subtype, sex, duration, base_mid = prev_mid)
  prev_summary <- prev_base |>
    filter(year %in% c(2021, 2030, 2040, 2045)) |>
    left_join(baseline_2021, by = c("subtype", "sex", "duration")) |>
    mutate(pct_change = round((prev_mid - base_mid) / base_mid * 100, 1)) |>
    select(subtype, sex, duration, year,
           prev_mid, prev_p025, prev_p975, pct_change) |>
    arrange(duration, subtype, sex, year)

  # 10-year validation file derived from main results (was a standalone
  # compute pass before durations was extended to include 10).
  prev_10yr_validation <- prev_full |>
    filter(scenario == "base", duration == 10, year <= 2021) |>
    select(year, prev_mid, prev_p025, prev_p975, prev_mean,
           duration, subtype, sex)

  # ------------------------------------------------------------------
  # table_3: combined-duration prevalence (persons, base case), with the
  # headline TOTAL lymphoma (HL + aggregate NHL). Built directly from the
  # per-draw persons matrices so the % change carries a draw-level CrI too.
  # Subtype rows (DLBCL/FL/MCL) are flagged subset_of_nhl = TRUE and are a
  # decomposition of the NHL row, never summed into the total.
  # ------------------------------------------------------------------
  t3_series <- c("hodgkin", "nhl", "dlbcl", "follicular", "mantle_cell", "total")
  q <- function(x, p) unname(quantile(x, p))
  t3_rows <- list()
  for (dur in durations) {
    pb <- list()
    for (st in prev_series) {
      pb[[st]] <- draw_store[[dkey(st, "males", "base", dur)]] +
                  draw_store[[dkey(st, "females", "base", dur)]]
    }
    pb[["total"]] <- pb[["hodgkin"]] + pb[["nhl"]]
    for (st in t3_series) {
      m   <- pb[[st]]
      d21 <- m[match("2021", rownames(m)), ]
      d45 <- m[match("2045", rownames(m)), ]
      pc  <- (d45 / d21 - 1) * 100
      t3_rows[[length(t3_rows) + 1]] <- tibble::tibble(
        subtype = st, subset_of_nhl = st %in% nhl_subtypes, duration = dur,
        prev_2021 = median(d21), prev_2021_p025 = q(d21, .025), prev_2021_p975 = q(d21, .975),
        prev_2045 = median(d45), prev_2045_p025 = q(d45, .025), prev_2045_p975 = q(d45, .975),
        pct_change = median(pc), pct_change_p025 = q(pc, .025), pct_change_p975 = q(pc, .975)
      )
    }
  }
  table_3 <- bind_rows(t3_rows) |>
    mutate(across(starts_with("prev_"), ~ round(.x)),
           across(starts_with("pct_"),  ~ round(.x, 1))) |>
    arrange(duration, match(subtype, t3_series))

  write_csv(table_3,             file.path(save_dir, "table_3_prevalence_combined.csv"))
  write_csv(prev_base,            file.path(save_dir, "prevalence_projections.csv"))
  write_csv(prev_sensitivity,     file.path(save_dir, "prevalence_sensitivity.csv"))
  write_csv(prev_summary,         file.path(save_dir, "table_2_prevalence_summary.csv"))
  write_csv(prev_10yr_validation, file.path(save_dir, "prevalence_validation_10yr.csv"))

  # 6. Render Figures 2 and 3
  fig2 <- build_prev_fig(prev_base, duration_value = 5,
                         prev_subtype_aihw = inputs$prev_subtype,
                         prev_agg_aihw     = inputs$prev_agg)
  fig3 <- build_prev_fig(prev_base, duration_value = 40)

  save_fig(fig2, file.path(save_dir, "figure_2_prevalence_5yr"),
           width = 12, height = 6.5)
  save_fig(fig3, file.path(save_dir, "figure_3_prevalence_40yr"),
           width = 12, height = 6.5)

  if (verbose) cat(sprintf("Wrote 5 CSVs + Figs 2 & 3 to %s\n", save_dir))

  res <- list(base = prev_base, sensitivity = prev_sensitivity,
              summary = prev_summary, combined = table_3,
              validation_10yr = prev_10yr_validation,
              surv_draws_tbl = surv_draws_tbl,
              improvement_draws_tbl = improvement_draws_tbl,
              surv_obs = surv_obs,
              improvement_final = improvement_pt)
  attr(res, "elapsed_s") <- total_elapsed
  invisible(res)
}

# ----------------------
# Auto-run guard
# ----------------------

if (!interactive() && sys.nframe() == 0) {
  invisible(run_prev_model())
}
