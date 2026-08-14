# ================
# Blood Cancer Projections (lymphoma): Age-Period-Cohort Incidence Modelling (MC-only)
# ================
# Two-tier projections to 2045 with parameter-level Monte Carlo CIs:
#   Tier 1 (aggregate): NHL and HL fitted on 1982-2021 (4 models)
#   Tier 2 (subtype):   DLBCL, FL, MCL fitted on 2003-2021 (6 models)
#
# Method: Epi::apc.fit() with natural cubic splines (Sasieni/Carstensen).
# Uncertainty: B = 1000 parameter draws from MVN(coef, vcov) per fit;
# year-level totals computed by summing within draws and taking
# percentiles (avoids the perfect-correlation assumption of the previous
# bound-propagation method).
#
# Library + driver pattern:
#   source("code/apc_model.R")    -> defines functions; no side effects
#   run_apc_model()               -> runs full pipeline
#   Rscript code/apc_model.R      -> auto-runs via the bottom guard
# ================

set.seed(20260507)
stopifnot(
  requireNamespace("MASS"),
  requireNamespace("Epi"),
  requireNamespace("dplyr"),
  requireNamespace("readr"),
  requireNamespace("purrr"),
  requireNamespace("tidyr"),
  requireNamespace("ggplot2"),
  requireNamespace("cowplot"),
  requireNamespace("grid"),
  requireNamespace("gtable")
)

suppressPackageStartupMessages({
  library(tidyverse)
  library(Epi)
  library(MASS, exclude = "select")
  library(cowplot)
  library(grid)
})

# Shared constants, colours, grob makers, save_fig. See code/_setup.R.
source("code/_setup.R")

# ----------------------
# Data prep helpers
# ----------------------

prep_subtype_data <- function(subtype_name, sex_name, inc_subtype, pop_hist) {
  inc_subtype |>
    dplyr::filter(subtype == subtype_name, sex == sex_name) |>
    dplyr::left_join(
      pop_hist |> dplyr::filter(sex == sex_name),
      by = c("year", "sex", "age_group")
    ) |>
    dplyr::mutate(
      A = age_mid[age_group], P = year, D = count, Y = population / 1e5
    ) |>
    dplyr::filter(!is.na(A), !is.na(Y), Y > 0, P >= subtype_start) |>
    dplyr::select(A, P, D, Y, age_group)
}

prep_agg_data <- function(group_name, sex_name, inc_agg, pop_hist) {
  cancer_group_label <- switch(group_name,
                               nhl     = "Non-Hodgkin lymphoma",
                               hodgkin = "Hodgkin lymphoma")
  inc_agg |>
    dplyr::filter(cancer_group == cancer_group_label, sex == sex_name) |>
    dplyr::left_join(
      pop_hist |> dplyr::filter(sex == sex_name),
      by = c("year", "sex", "age_group")
    ) |>
    dplyr::mutate(
      A = age_mid[age_group], P = year, D = count, Y = population / 1e5
    ) |>
    dplyr::filter(!is.na(A), !is.na(Y), Y > 0, P >= agg_start) |>
    dplyr::select(A, P, D, Y, age_group)
}

prepare_fit_data <- function(key, tier, inc_subtype, inc_agg, pop_hist) {
  parts        <- strsplit(key, "_")[[1]]
  sex_name     <- parts[length(parts)]
  subtype_name <- paste(parts[-length(parts)], collapse = "_")
  if (tier == "subtype") {
    prep_subtype_data(subtype_name, sex_name, inc_subtype, pop_hist)
  } else {
    prep_agg_data(subtype_name, sex_name, inc_agg, pop_hist)
  }
}

# ----------------------
# 2001 ASP weights (per 100,000)
# ----------------------
# Sums across sexes; expresses as proportion x 1e5. Includes the 0-4 age
# band (rate = 0 for lymphoma) so the denominator matches the AIHW
# convention.

build_std_pop <- function(pop_hist) {
  s <- pop_hist |>
    dplyr::filter(year == 2001) |>
    dplyr::group_by(age_group) |>
    dplyr::summarise(population = sum(population), .groups = "drop") |>
    dplyr::mutate(weight = population / sum(population) * 1e5)
  setNames(s$weight, s$age_group)
}

# ----------------------
# Population matrix (used by per-draw projection)
# ----------------------
# Built once per sex inside the driver, passed into project_one_draw().

build_pop_mat <- function(sex_name, pop_proj_df) {
  P_grid <- proj_start:proj_end
  wide <- pop_proj_df |>
    dplyr::filter(sex == sex_name,
                  year %in% P_grid,
                  age_group %in% names(age_mid)) |>
    dplyr::select(year, age_group, population) |>
    tidyr::pivot_wider(names_from = age_group, values_from = population) |>
    dplyr::arrange(year)
  if (nrow(wide) != length(P_grid)) {
    stop("pop_proj_df missing entries for sex=", sex_name,
         " across proj_start..proj_end")
  }
  pop_mat <- as.matrix(wide[, names(age_mid)])
  rownames(pop_mat) <- as.character(wide$year)
  if (any(is.na(pop_mat))) {
    stop("pop_proj_df missing entries for sex=", sex_name,
         " across proj_start..proj_end x 9 age groups")
  }
  pop_mat
}

# ----------------------
# APC basis + per-draw component reconstruction
# ----------------------
# `apc_basis(fit, df)` precomputes the design-matrix slices needed to
# reconstruct (Age, Per, Coh) effects from any beta vector. Variable
# parameter count (11 or 12) is handled by regex column matching.
#
# `apc_components(beta, basis)` returns Age/Per/Coh as named vectors at
# unique data points. Bitwise-equal to fit$Age/Per/Coh when beta = coef.

apc_basis <- function(fit, df) {
  if (!inherits(fit, "apc")) stop("`fit` must be an apc.fit object")
  if (!all(c("A", "P") %in% names(df))) stop("`df` must have columns A and P")
  mm <- model.matrix(fit$Model)
  if (nrow(mm) != nrow(df)) {
    stop("model.matrix has ", nrow(mm), " rows but df has ", nrow(df))
  }

  cn       <- colnames(mm)
  ma_cols  <- grep("^MA([0-9]+)?$",  cn, value = TRUE)
  mpr_cols <- grep("^MPr[0-9]+$",    cn, value = TRUE)
  mcr_cols <- grep("^MCr([0-9]+)?$", cn, value = TRUE)
  if (length(ma_cols) < 2 || length(mpr_cols) < 1 || length(mcr_cols) < 2) {
    stop("design matrix is missing expected blocks - is parm = 'ACP', model = 'ns'?")
  }
  if (length(c(ma_cols, mpr_cols, mcr_cols)) != ncol(mm)) {
    stop("column groups don't cover the full design matrix")
  }

  A_pt  <- sort(unique(df$A));     A_pos <- match(A_pt, df$A)
  P_pt  <- sort(unique(df$P));     P_pos <- match(P_pt, df$P)
  C     <- df$P - df$A
  C_pt  <- sort(unique(C));        C_pos <- match(C_pt, C)

  structure(
    list(
      ma_cols = ma_cols, mpr_cols = mpr_cols, mcr_cols = mcr_cols,
      A_pt = A_pt, P_pt = P_pt, C_pt = C_pt,
      ctr_A = mm[A_pos, ma_cols,  drop = FALSE],
      ctr_P = mm[P_pos, mpr_cols, drop = FALSE],
      ctr_C = mm[C_pos, mcr_cols, drop = FALSE],
      n_params = length(coef(fit$Model))
    ),
    class = "apc_basis"
  )
}

apc_components <- function(beta, basis) {
  if (!inherits(basis, "apc_basis")) {
    stop("`basis` must be an apc_basis object")
  }
  if (length(beta) != basis$n_params) {
    stop("beta length ", length(beta), " mismatch (expected ", basis$n_params, ")")
  }
  list(
    Age = setNames(as.numeric(exp(basis$ctr_A %*% beta[basis$ma_cols])),
                   basis$A_pt),
    Per = setNames(as.numeric(exp(basis$ctr_P %*% beta[basis$mpr_cols])),
                   basis$P_pt),
    Coh = setNames(as.numeric(exp(basis$ctr_C %*% beta[basis$mcr_cols])),
                   basis$C_pt)
  )
}

# ----------------------
# Damped projection (per draw)
# ----------------------
# Returns a (24 x 9) matrix of projected case counts for one beta draw.
# `pop_mat` is the precomputed (24 x 9) population matrix for this sex.

project_one_draw <- function(beta, basis, pop_mat) {
  comp   <- apc_components(beta, basis)
  age_rr <- comp$Age
  per_rr <- comp$Per
  coh_rr <- comp$Coh

  last_per_P  <- as.numeric(names(per_rr))[length(per_rr)]
  last_per_RR <- per_rr[length(per_rr)]
  last_coh_C  <- as.numeric(names(coh_rr))[length(coh_rr)]
  last_coh_RR <- coh_rr[length(coh_rr)]

  coh_interp <- approxfun(as.numeric(names(coh_rr)), as.numeric(coh_rr),
                          method = "linear", rule = 2)

  P_grid <- proj_start:proj_end
  A_grid <- age_mid

  t_per <- P_grid - hist_end
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
# Beta sampling and batch projection
# ----------------------

sample_beta_draws <- function(fit, B = 1000) {
  mu    <- coef(fit$Model)
  Sigma <- vcov(fit$Model)
  draws <- MASS::mvrnorm(n = B, mu = mu, Sigma = Sigma)
  if (B == 1) draws <- matrix(draws, nrow = 1)
  colnames(draws) <- names(mu)
  draws
}

project_all_draws <- function(beta_draws, basis, pop_mat) {
  B   <- nrow(beta_draws)
  P_grid <- proj_start:proj_end
  cases_array <- array(
    NA_real_,
    dim = c(length(P_grid), length(age_mid), B),
    dimnames = list(year      = as.character(P_grid),
                    age_group = names(age_mid),
                    draw      = paste0("d", seq_len(B)))
  )
  for (b in seq_len(B)) {
    cases_array[, , b] <- project_one_draw(beta_draws[b, ], basis, pop_mat)
  }
  cases_array
}

# Cell-level (year x age_group) summary across draws. Reported
# `cases_mid` is the median (matches a deterministic point estimate
# under the log link); `cases_mean` is retained for bias-aware analyses.
build_cases_summary <- function(cases_array) {
  P_grid <- as.numeric(dimnames(cases_array)[[1]])
  A_grid <- dimnames(cases_array)[[2]]
  s <- summarise_draws(cases_array, c(1, 2))
  expand.grid(P = P_grid, age_group = A_grid, KEEP.OUT.ATTRS = FALSE,
              stringsAsFactors = FALSE) |>
    dplyr::mutate(
      cases_mid  = as.vector(s$mid),
      cases_p025 = as.vector(s$p025),
      cases_p975 = as.vector(s$p975),
      cases_mean = as.vector(s$mean)
    ) |>
    dplyr::as_tibble()
}

# ----------------------
# Year-level aggregation (sum WITHIN draws first, then percentile)
# ----------------------
# Replaces the analytical compile_proj_annual()'s incorrect bound-sum
# (which assumed perfect correlation across age groups). Returns one
# row per (subtype, sex, year).

build_proj_annual <- function(apc_results, std_pop, pop_mat_by_sex) {
  purrr::map_dfr(apc_results, function(entry) {
    sex_name    <- entry$sex
    cases_draws <- entry$cases_draws  # (24 x 9 x B)
    pop_mat     <- pop_mat_by_sex[[sex_name]]

    # Cases per (year, draw)
    cases_yr_b <- apply(cases_draws, c(1, 3), sum)  # (24 x B)

    # Rates per (year, age, draw)
    pop_arr     <- array(rep(pop_mat, dim(cases_draws)[3]), dim = dim(cases_draws))
    rates_draws <- cases_draws / pop_arr * 1e5

    # ASR per (year, draw): sum_age rate * std_weight / sum(std_weights)
    std_vec   <- std_pop[dimnames(cases_draws)[[2]]]
    sum_std   <- sum(std_pop)
    asr_yr_b  <- apply(rates_draws, c(1, 3),
                       function(r) sum(r * std_vec) / sum_std)

    cs <- summarise_draws(cases_yr_b, 1)
    as_ <- summarise_draws(asr_yr_b,  1)
    tibble::tibble(
      P          = as.numeric(dimnames(cases_draws)[[1]]),
      cases_mid  = cs$mid,
      cases_p025 = cs$p025,
      cases_p975 = cs$p975,
      asr_mid    = as_$mid,
      asr_p025   = as_$p025,
      asr_p975   = as_$p975,
      subtype    = entry$subtype,
      sex        = entry$sex,
      tier       = entry$tier
    )
  })
}

# ----------------------
# Draw-level combined incidence total (HL + aggregate NHL)
# ----------------------
# HL and aggregate NHL are disjoint and together are all lymphoma, so the
# headline total = HL + NHL. Formed at the DRAW level (sum the age x draw
# case arrays, then compute cases and ASR per draw and take percentiles),
# by sex and persons. The NHL subtypes are a subset of NHL and are NOT
# summed in. Persons ASR uses both-sex population (not the sum of by-sex
# ASRs). Returns rows with subtype = "total".

build_incidence_total <- function(apc_results, std_pop, pop_mat_by_sex) {
  ref     <- apc_results[["nhl_males"]]$cases_draws
  years   <- as.numeric(dimnames(ref)[[1]])
  std_vec <- std_pop[dimnames(ref)[[2]]]
  sum_std <- sum(std_pop)
  cd <- function(grp, sx) apc_results[[paste(grp, sx, sep = "_")]]$cases_draws

  annual <- function(cases_age_b, pop_mat) {  # (year x age x B), pop (year x age)
    cases_yr_b <- apply(cases_age_b, c(1, 3), sum)
    pop_arr    <- array(rep(pop_mat, dim(cases_age_b)[3]), dim = dim(cases_age_b))
    rates      <- cases_age_b / pop_arr * 1e5
    asr_yr_b   <- apply(rates, c(1, 3), function(r) sum(r * std_vec) / sum_std)
    cs <- summarise_draws(cases_yr_b, 1); as_ <- summarise_draws(asr_yr_b, 1)
    tibble::tibble(P = years,
                   cases_mid = cs$mid, cases_p025 = cs$p025, cases_p975 = cs$p975,
                   asr_mid = as_$mid, asr_p025 = as_$p025, asr_p975 = as_$p975)
  }

  rows <- lapply(sexes, function(sx)
    annual(cd("hodgkin", sx) + cd("nhl", sx), pop_mat_by_sex[[sx]]) |>
      dplyr::mutate(subtype = "total", sex = sx, tier = "total"))
  tot_p <- cd("hodgkin", "males") + cd("hodgkin", "females") +
           cd("nhl", "males")     + cd("nhl", "females")
  pop_p <- pop_mat_by_sex[["males"]] + pop_mat_by_sex[["females"]]
  rows[[length(rows) + 1]] <- annual(tot_p, pop_p) |>
    dplyr::mutate(subtype = "total", sex = "persons", tier = "total")
  dplyr::bind_rows(rows)
}

# Observed (deterministic) baseline total (HL + NHL), by sex and persons.
build_hist_total <- function(inc_agg, pop_hist, std_pop) {
  lymph <- c("Non-Hodgkin lymphoma", "Hodgkin lymphoma")
  by_sex <- purrr::map_dfr(sexes, function(sx) {
    inc_tot <- inc_agg |>
      dplyr::filter(cancer_group %in% lymph, sex == sx) |>
      dplyr::group_by(year, age_group) |>
      dplyr::summarise(count = sum(count), .groups = "drop")
    calc_hist_asr(inc_tot, sx, pop_hist, std_pop) |>
      dplyr::mutate(subtype = "total", sex = sx, tier = "total")
  })
  pop_p <- pop_hist |>
    dplyr::group_by(year, age_group) |>
    dplyr::summarise(population = sum(population), .groups = "drop")
  persons <- inc_agg |>
    dplyr::filter(cancer_group %in% lymph) |>
    dplyr::group_by(year, age_group) |>
    dplyr::summarise(count = sum(count), .groups = "drop") |>
    dplyr::left_join(pop_p, by = c("year", "age_group")) |>
    dplyr::mutate(rate = count / population * 1e5, std_wt = std_pop[age_group]) |>
    dplyr::group_by(year) |>
    dplyr::summarise(cases_obs = sum(count),
                     asr_obs = sum(rate * std_wt, na.rm = TRUE) / sum(std_pop),
                     .groups = "drop") |>
    dplyr::mutate(subtype = "total", sex = "persons", tier = "total")
  dplyr::bind_rows(by_sex, persons)
}

# ----------------------
# Historical observed ASR (deterministic)
# ----------------------

calc_hist_asr <- function(inc_df, sex_name, pop_hist, std_pop) {
  inc_df |>
    dplyr::left_join(
      pop_hist |> dplyr::filter(sex == sex_name) |>
        dplyr::select(year, age_group, population),
      by = c("year", "age_group")
    ) |>
    dplyr::mutate(
      rate   = count / population * 1e5,
      std_wt = std_pop[age_group]
    ) |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      cases_obs = sum(count),
      asr_obs   = sum(rate * std_wt, na.rm = TRUE) / sum(std_pop),
      .groups   = "drop"
    )
}

compile_hist_annual <- function(apc_results) {
  purrr::map_dfr(apc_results, function(x) {
    x$hist |> dplyr::mutate(subtype = x$subtype, sex = x$sex, tier = x$tier)
  })
}

# ----------------------
# Fitted historical ASR (model-smoothed; deterministic point estimate)
# ----------------------
# Uses point-estimate effects so Figure 1 has a smooth fitted line over
# the historical period. (CIs on this line are not displayed.)

calc_fitted_asr <- function(fit, sex_name, label, fit_start, std_pop) {
  age_df <- as.data.frame(fit$Age) |>
    dplyr::rename(A = Age, age_rr = Rate)
  per_df <- as.data.frame(fit$Per) |>
    dplyr::rename(P = Per, per_rr = `P-RR`)
  coh_df <- as.data.frame(fit$Coh) |>
    dplyr::rename(C = Coh, coh_rr = `C-RR`)

  per_interp <- approxfun(per_df$P, per_df$per_rr, method = "linear", rule = 2)
  coh_interp <- approxfun(coh_df$C, coh_df$coh_rr, method = "linear", rule = 2)

  hist_grid <- expand.grid(A = age_mid, P = fit_start:hist_end) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      age_group = names(age_mid)[match(A, age_mid)],
      C = P - A
    ) |>
    dplyr::left_join(age_df |> dplyr::select(A, age_rr), by = "A") |>
    dplyr::mutate(
      per_rr      = per_interp(P),
      coh_rr      = coh_interp(C),
      rate_fitted = age_rr * per_rr * coh_rr
    )

  hist_grid |>
    dplyr::group_by(P) |>
    dplyr::summarise(
      asr_fitted = sum(rate_fitted * std_pop[age_group], na.rm = TRUE) /
                   sum(std_pop),
      .groups = "drop"
    ) |>
    dplyr::mutate(subtype = label, sex = sex_name)
}

# ----------------------
# APC effects from beta draws (replaces Wald CIs from Epi::apc.fit)
# ----------------------
# Long-format tibble of effect, value, rr, rr_p025, rr_p975 per fit.
# Rebuilds (Age, Per, Coh) for every draw and takes percentiles. Single
# CI method end-to-end (consistent with everything else MC).

build_apc_effects <- function(apc_results, prepared_data) {
  purrr::map_dfr(apc_results, function(entry) {
    key   <- paste0(entry$subtype, "_", entry$sex)
    df    <- prepared_data[[key]]
    basis <- apc_basis(entry$fit, df)

    # Vectorise across all draws: a single matrix multiply per effect
    # block replaces the previous B-call apc_components() loop.
    bd  <- entry$beta_draws
    age_mat <- t(exp(basis$ctr_A %*% t(bd[, basis$ma_cols,  drop = FALSE])))
    per_mat <- t(exp(basis$ctr_P %*% t(bd[, basis$mpr_cols, drop = FALSE])))
    coh_mat <- t(exp(basis$ctr_C %*% t(bd[, basis$mcr_cols, drop = FALSE])))
    colnames(age_mat) <- as.character(basis$A_pt)
    colnames(per_mat) <- as.character(basis$P_pt)
    colnames(coh_mat) <- as.character(basis$C_pt)

    .effect_tbl <- function(M, effect_name) {
      s <- summarise_draws(M, 2)
      tibble::tibble(
        effect  = effect_name,
        value   = as.numeric(colnames(M)),
        rr      = s$mid,
        rr_p025 = s$p025,
        rr_p975 = s$p975
      )
    }
    dplyr::bind_rows(
      .effect_tbl(age_mat, "age"),
      .effect_tbl(per_mat, "period"),
      .effect_tbl(coh_mat, "cohort")
    ) |>
      dplyr::mutate(subtype = entry$subtype,
                    sex     = entry$sex,
                    tier    = entry$tier)
  })
}

# ----------------------
# Fit statistics (deterministic; uses the underlying glm)
# ----------------------

build_fit_stats <- function(apc_results) {
  purrr::map_dfr(apc_results, function(x) {
    mdl <- x$fit$Model
    tibble::tibble(
      subtype      = x$subtype,
      sex          = x$sex,
      deviance     = if (!is.null(mdl)) mdl$deviance     else NA_real_,
      df_resid     = if (!is.null(mdl)) mdl$df.residual  else NA_integer_,
      dev_df_ratio = if (!is.null(mdl)) round(mdl$deviance / mdl$df.residual, 2)
                     else NA_real_
    )
  })
}

# ----------------------
# Figure 1: incidence ASR, two sex panels, one line per lymphoma
# ----------------------
# Matches the convention of Figures 2 & 3: two panels (Females, Males), one
# line per lymphoma coloured by type, with the headline tier (HL, aggregate
# NHL) solid and the decomposition tier (DLBCL/FL/MCL, a subset of NHL)
# dashed. Faint observed points, AIHW projection triangles, projection 95%
# CrI ribbon, dotted line at the last observed year. The combined total
# (subtype "total") is a table/abstract quantity, not plotted here.

build_figure_1 <- function(proj_annual, fitted_annual, hist_annual, aihw_proj_df) {
  lvl         <- c("HL", "NHL", "DLBCL", "FL", "MCL")
  draw_series <- c("hodgkin", "nhl", nhl_subtypes)
  .lab  <- function(st) factor(label_map[st], levels = lvl)
  .tier <- function(st) dplyr::if_else(st %in% c("hodgkin", "nhl"),
                                       "headline", "decomposition")

  model_line <- dplyr::bind_rows(
    fitted_annual |>
      dplyr::transmute(year = P, asr = asr_fitted, subtype, sex,
                       asr_p025 = NA_real_, asr_p975 = NA_real_),
    proj_annual |>
      dplyr::transmute(year = P, asr = asr_mid, subtype, sex, asr_p025, asr_p975)
  ) |>
    dplyr::filter(subtype %in% draw_series, sex %in% sexes) |>
    dplyr::mutate(label = .lab(subtype), tier = .tier(subtype))

  obs <- hist_annual |>
    dplyr::filter(subtype %in% draw_series, asr_obs > 0, sex %in% sexes,
                  year >= dplyr::if_else(subtype %in% c("hodgkin", "nhl"),
                                         agg_start, subtype_start)) |>
    dplyr::mutate(label = .lab(subtype))

  aihw_proj <- aihw_proj_df |>
    dplyr::mutate(
      sex     = tolower(sex),
      subtype = dplyr::case_when(
        cancer_group == "Non-Hodgkin lymphoma" ~ "nhl",
        cancer_group == "Hodgkin lymphoma"     ~ "hodgkin")) |>
    dplyr::filter(!is.na(subtype)) |>
    dplyr::mutate(label = .lab(subtype))

  # One legend keyed on lymphoma showing both colour AND line style (solid for
  # headline HL/NHL, dashed for the subtypes); colour + linetype scales share
  # the `label` mapping so ggplot merges them. Observed points, AIHW triangles
  # and the ribbon are kept out of the legend (show.legend = FALSE). No in-image
  # caption — the solid/dashed key, AIHW and CrI are explained in the .docx caption.
  lty_vals <- c(HL = "solid", NHL = "solid", DLBCL = "22", FL = "22", MCL = "22")
  ggplot(model_line, aes(year, colour = label)) +
    geom_ribbon(data = ~ dplyr::filter(.x, year > hist_end),
                aes(ymin = asr_p025, ymax = asr_p975, fill = label),
                alpha = 0.15, colour = NA, show.legend = FALSE) +
    geom_point(data = obs, aes(year, asr_obs), size = 0.6, alpha = 0.45,
               show.legend = FALSE) +
    geom_line(aes(y = asr, linetype = label), linewidth = 0.8) +
    geom_point(data = aihw_proj, aes(year, asr_2001), shape = 2,
               size = 1.3, stroke = 0.5, show.legend = FALSE) +
    geom_vline(xintercept = hist_end + 0.5, linetype = "dotted",
               colour = "grey50", linewidth = 0.4) +
    facet_wrap(~ sex, ncol = 2, labeller = labeller(sex = sex_labels)) +
    scale_colour_manual(values = line_colours, name = NULL) +
    scale_fill_manual(values = line_colours, guide = "none") +
    scale_linetype_manual(values = lty_vals, name = NULL) +
    labs(x = "Year", y = "Age-standardised rate (per 100,000)") +
    theme_bw(base_size = 13) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank(),
          plot.background  = element_rect(colour = NA, fill = NA))
}

# ----------------------
# Driver
# ----------------------

run_apc_model <- function(B = 1000,
                          save_dir = "output",
                          seed = 20260507,
                          verbose = TRUE) {
  set.seed(seed)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. Load data
  inc_subtype <- readr::read_csv("data/incidence_subtype.csv", show_col_types = FALSE)
  inc_agg     <- readr::read_csv("data/incidence_agg.csv",     show_col_types = FALSE)
  pop_hist    <- readr::read_csv("data/pop_hist.csv",          show_col_types = FALSE)
  pop_proj    <- readr::read_csv("data/pop_proj.csv",          show_col_types = FALSE)
  aihw_proj   <- readr::read_csv("data/incidence_agg_proj.csv", show_col_types = FALSE)

  # 2. ASP weights and population matrices
  std_pop <- build_std_pop(pop_hist)
  pop_mat_by_sex <- list(
    males   = build_pop_mat("males",   pop_proj),
    females = build_pop_mat("females", pop_proj)
  )

  # 3. Fit + sample beta draws + project per draw
  results <- list()
  prepared_data <- list()  # cached per-key df for rebuilding bases later
  total_t0 <- Sys.time()

  for (grp in agg_groups) for (sx in sexes) {
    key <- paste(grp, sx, sep = "_")
    if (verbose) cat(sprintf("[Tier 1 agg] %s\n", key))
    df <- prep_agg_data(grp, sx, inc_agg, pop_hist)
    prepared_data[[key]] <- df

    fit <- apc.fit(data = df[, c("A", "P", "D", "Y")],
                   model = "ns", dist = "poisson",
                   npar = agg_knots, print.AOV = FALSE)

    basis        <- apc_basis(fit, df)
    beta_draws   <- sample_beta_draws(fit, B)
    cases_draws  <- project_all_draws(beta_draws, basis, pop_mat_by_sex[[sx]])
    cases_summary <- build_cases_summary(cases_draws)

    cancer_group_label <- switch(grp,
                                 nhl     = "Non-Hodgkin lymphoma",
                                 hodgkin = "Hodgkin lymphoma")
    inc_for_hist <- inc_agg |>
      dplyr::filter(cancer_group == cancer_group_label, sex == sx) |>
      dplyr::select(year, age_group, count)
    hist <- calc_hist_asr(inc_for_hist, sx, pop_hist, std_pop)

    results[[key]] <- list(
      subtype       = grp, sex = sx, tier = "agg",
      fit           = fit, hist = hist,
      beta_draws    = beta_draws,
      cases_draws   = cases_draws,
      cases_summary = cases_summary
    )
  }

  for (s in nhl_subtypes) for (sx in sexes) {
    key <- paste(s, sx, sep = "_")
    if (verbose) cat(sprintf("[Tier 2 subtype] %s\n", key))
    df <- prep_subtype_data(s, sx, inc_subtype, pop_hist)
    prepared_data[[key]] <- df

    fit <- apc.fit(data = df[, c("A", "P", "D", "Y")],
                   model = "ns", dist = "poisson",
                   npar = subtype_knots, print.AOV = FALSE)

    basis        <- apc_basis(fit, df)
    beta_draws   <- sample_beta_draws(fit, B)
    cases_draws  <- project_all_draws(beta_draws, basis, pop_mat_by_sex[[sx]])
    cases_summary <- build_cases_summary(cases_draws)

    inc_for_hist <- inc_subtype |>
      dplyr::filter(subtype == s, sex == sx) |>
      dplyr::select(year, age_group, count)
    hist <- calc_hist_asr(inc_for_hist, sx, pop_hist, std_pop)

    results[[key]] <- list(
      subtype       = s, sex = sx, tier = "subtype",
      fit           = fit, hist = hist,
      beta_draws    = beta_draws,
      cases_draws   = cases_draws,
      cases_summary = cases_summary
    )
  }

  total_elapsed <- as.numeric(Sys.time() - total_t0, units = "secs")
  if (verbose) cat(sprintf("Fits + MC sampling complete: %.1f s\n", total_elapsed))

  # 4. Annual aggregations (+ draw-level headline total: HL + aggregate NHL)
  proj_annual   <- dplyr::bind_rows(
    build_proj_annual(results, std_pop, pop_mat_by_sex),
    build_incidence_total(results, std_pop, pop_mat_by_sex))
  hist_annual   <- dplyr::bind_rows(
    compile_hist_annual(results),
    build_hist_total(inc_agg, pop_hist, std_pop))
  fitted_annual <- purrr::map_dfr(results, function(x) {
    fit_start <- if (x$tier == "agg") agg_start else subtype_start
    calc_fitted_asr(x$fit, x$sex, x$subtype, fit_start, std_pop)
  })

  # 5. Effect CIs (from beta draws) and fit statistics
  apc_effects   <- build_apc_effects(results, prepared_data)
  apc_fit_stats <- build_fit_stats(results)

  # 6. Summary table for Table 1
  summary_table <- proj_annual |>
    dplyr::filter(P %in% c(2030, 2040, 2045)) |>
    dplyr::left_join(
      hist_annual |>
        dplyr::filter(year == hist_end) |>
        dplyr::select(subtype, sex, cases_2021 = cases_obs, asr_2021 = asr_obs),
      by = c("subtype", "sex")
    ) |>
    dplyr::mutate(
      pct_change_cases = (cases_mid - cases_2021) / cases_2021 * 100,
      pct_change_asr   = (asr_mid   - asr_2021)   / asr_2021   * 100
    )

  # 7. Write CSVs
  readr::write_csv(proj_annual,   file.path(save_dir, "incidence_projections.csv"))
  readr::write_csv(hist_annual,   file.path(save_dir, "incidence_historical_asr.csv"))
  readr::write_csv(summary_table, file.path(save_dir, "table_1_incidence_summary.csv"))
  readr::write_csv(apc_effects,   file.path(save_dir, "apc_effects.csv"))
  readr::write_csv(apc_fit_stats, file.path(save_dir, "table_s1_apc_fit_stats.csv"))

  # 8. Render Figure 1
  fig1 <- build_figure_1(proj_annual, fitted_annual, hist_annual, aihw_proj)
  save_fig(fig1, file.path(save_dir, "figure_1_combined_projections"),
           width = 11, height = 6.5)

  # 9. Save the augmented apc_results object for downstream consumers
  saveRDS(results, file.path(save_dir, "apc_results.rds"), compress = "xz")

  if (verbose) {
    sz <- file.info(file.path(save_dir, "apc_results.rds"))$size / 1024^2
    cat(sprintf("Wrote %s (%.1f MB) and 5 CSVs + Figure 1\n",
                file.path(save_dir, "apc_results.rds"), sz))
    cat(sprintf("Total run time: %.1f s\n",
                as.numeric(Sys.time() - total_t0, units = "secs")))
  }

  attr(results, "elapsed_s") <- as.numeric(Sys.time() - total_t0, units = "secs")
  invisible(results)
}

# ----------------------
# Auto-run guard
# ----------------------
# `Rscript code/apc_model.R` from the repo root runs the pipeline.
# Sourcing the file in interactive R is silent.

if (!interactive() && sys.nframe() == 0) {
  invisible(run_apc_model())
}
