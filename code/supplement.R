# ================
# Blood Cancer Projections (lymphoma): Supplementary Information Outputs (MC-only)
# ================
# Produces the supplementary tables and figures that appear in the
# Supporting Information. Reads inputs from disk - no in-memory state
# required from prior scripts (see disk-input pattern in §4.3 of the
# transition plan).
#
# Tables (CSVs), flat SI numbering (round 3):
#   S1   Disease definitions            - this file
#   S2   APC fit statistics             - written by apc_model.R (no code here)
#   S3   Back-estimation validation     - this file
#   S4   Incidence vs AIHW projections  - this file (forecast comparison)
#   S5   Subtype/aggregate consistency  - this file
#   S6   Suppressed survival cells      - this file
#   S7   Survival improvement rates     - this file
#   S8   5- and 10-yr prevalence vs AIHW - this file
#   S9   40-yr prevalence vs AIHW       - this file
#   S11  Prevalence by sex, 2/3/10-yr   - this file (from table_2_change_cri.csv)
#   S12  Sensitivity scenarios          - this file
#
# Figures:
#   S1/S2/S3  APC effect plots (age, period, cohort) - this file
#   S4        Survival curves used in prevalence     - this file
#   S5        Worked-example survival curve          - this file
#   S6        Population pyramid 2021 vs 2045        - this file
#   S7/S8/S9  Additional duration prevalence trajectories (2/3/10-yr) - this file
#   S10/S11   Survival by diagnosis year (5-yr, 10-yr) - this file
#   S12/S13   Sensitivity prevalence trajectories (5-yr, 40-yr) - this file
#
# Library + driver pattern:
#   source("code/supplement.R")    -> defines functions; no side effects
#   run_supplement()               -> reads inputs, writes outputs
#   Rscript code/supplement.R      -> auto-runs via the bottom guard
# ================

set.seed(20260507)
stopifnot(
  requireNamespace("dplyr"),
  requireNamespace("tidyr"),
  requireNamespace("readr"),
  requireNamespace("purrr"),
  requireNamespace("ggplot2"),
  requireNamespace("cowplot"),
  requireNamespace("scales")
)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(cowplot)
})

# Shared constants, colours, grob makers, save_fig. See code/_setup.R.
source("code/_setup.R")

# ----------------------
# Table builders
# ----------------------

# Table S4 (was S2a): aggregate-tier incidence projection vs AIHW Table S1e.1.
build_table_s4 <- function(proj_annual, inc_agg_proj) {
  inc_agg_proj <- inc_agg_proj |> mutate(sex = tolower(sex))
  our_proj <- proj_annual |>
    filter(P >= 2026, P <= 2035) |>
    rename(year = P)

  val_hl <- our_proj |> filter(subtype == "hodgkin") |>
    select(year, sex, our_count = cases_mid, our_asr = asr_mid) |>
    left_join(
      inc_agg_proj |> filter(cancer_group == "Hodgkin lymphoma") |>
        select(year, sex, aihw_count = count, aihw_asr = asr_2001),
      by = c("year", "sex")
    ) |>
    mutate(count_diff_pct = (our_count - aihw_count) / aihw_count * 100,
           asr_diff_pct   = (our_asr   - aihw_asr)   / aihw_asr   * 100)

  val_nhl <- our_proj |> filter(subtype == "nhl") |>
    select(year, sex, our_count = cases_mid, our_asr = asr_mid) |>
    left_join(
      inc_agg_proj |> filter(cancer_group == "Non-Hodgkin lymphoma") |>
        select(year, sex, aihw_count = count, aihw_asr = asr_2001),
      by = c("year", "sex")
    ) |>
    mutate(count_diff_pct = (our_count - aihw_count) / aihw_count * 100,
           asr_diff_pct   = (our_asr   - aihw_asr)   / aihw_asr   * 100)

  bind_rows(
    val_hl  |> mutate(cancer_group = "Hodgkin lymphoma"),
    val_nhl |> mutate(cancer_group = "Non-Hodgkin lymphoma")
  ) |>
    group_by(cancer_group, sex) |>
    summarise(
      count_diff_mean = round(mean(count_diff_pct), 1),
      count_diff_min  = round(min(count_diff_pct), 1),
      count_diff_max  = round(max(count_diff_pct), 1),
      asr_diff_mean   = round(mean(asr_diff_pct), 1),
      asr_diff_min    = round(min(asr_diff_pct), 1),
      asr_diff_max    = round(max(asr_diff_pct), 1),
      .groups = "drop"
    ) |>
    arrange(cancer_group, sex)
}

# Table S8 (was S2b): 5- and 10-yr prevalence vs AIHW (per (subtype, sex)). The 4 modelled
# lymphomas validate against AIHW subtype prevalence (Book 11i); aggregate
# NHL validates against AIHW aggregate prevalence (Book 6) (QL #3).
build_table_s8 <- function(prev_5yr, prev_10yr, prev_subtype_aihw, prev_agg_aihw) {
  aihw_target <- function(dur) {
    sub <- prev_subtype_aihw |>
      filter(duration == dur) |>
      group_by(year, subtype, sex) |>
      summarise(prev_aihw = sum(prevalence, na.rm = TRUE), .groups = "drop")
    nhl <- prev_agg_aihw |>
      filter(duration == dur, cancer_group == "Non-Hodgkin lymphoma") |>
      group_by(year, sex) |>
      summarise(prev_aihw = sum(prevalence, na.rm = TRUE), .groups = "drop") |>
      mutate(subtype = "nhl")
    bind_rows(sub, nhl)
  }
  val <- function(prev_df, dur) {
    prev_df |>
      filter(duration == dur, year <= 2021) |>
      inner_join(aihw_target(dur), by = c("year", "subtype", "sex")) |>
      mutate(pct_diff = (prev_mid - prev_aihw) / prev_aihw * 100) |>
      group_by(subtype, sex) |>
      summarise(mean = round(mean(pct_diff), 1),
                min  = round(min(pct_diff), 1),
                max  = round(max(pct_diff), 1), .groups = "drop")
  }
  v5  <- val(prev_5yr, 5)   |> rename(dur5_mean = mean, dur5_min = min, dur5_max = max)
  v10 <- val(prev_10yr, 10) |> rename(dur10_mean = mean, dur10_min = min, dur10_max = max)
  inner_join(v5, v10, by = c("subtype", "sex")) |> arrange(subtype, sex)
}

# Table S9 (was S2c): 40-yr prevalence at 2021 vs AIHW.
build_table_s9 <- function(prev_40yr, prev_agg_aihw) {
  aihw_nhl_40yr <- prev_agg_aihw |>
    filter(cancer_group == "Non-Hodgkin lymphoma",
           year == 2021, duration == 40) |>
    group_by(sex) |>
    summarise(aihw_40yr = sum(prevalence, na.rm = TRUE), .groups = "drop")
  our_40yr_nhl <- prev_40yr |>
    filter(duration == 40, year == 2021, subtype %in% nhl_subtypes) |>
    group_by(sex) |>
    summarise(our_40yr = sum(prev_mid), .groups = "drop")
  nhl_comp <- our_40yr_nhl |>
    left_join(aihw_nhl_40yr, by = "sex") |>
    mutate(ratio_pct = our_40yr / aihw_40yr * 100)

  # Direct aggregate-NHL validation (QL #3): NHL computed from aggregate NHL
  # incidence + NHL survival, vs AIHW aggregate NHL 40-yr. Unlike the subtype
  # sum (which covers only ~60% of NHL), this should validate close to 100%.
  our_40yr_nhlagg <- prev_40yr |>
    filter(duration == 40, year == 2021, subtype == "nhl") |>
    group_by(sex) |>
    summarise(our_40yr = sum(prev_mid), .groups = "drop")
  nhlagg_comp <- our_40yr_nhlagg |>
    left_join(aihw_nhl_40yr, by = "sex") |>
    mutate(ratio_pct = our_40yr / aihw_40yr * 100)

  aihw_hl_40yr <- prev_agg_aihw |>
    filter(cancer_group == "Hodgkin lymphoma",
           year == 2021, duration == 40) |>
    group_by(sex) |>
    summarise(aihw_40yr = sum(prevalence, na.rm = TRUE), .groups = "drop")
  our_40yr_hl <- prev_40yr |>
    filter(duration == 40, year == 2021, subtype == "hodgkin") |>
    group_by(sex) |>
    summarise(our_40yr = sum(prev_mid), .groups = "drop")
  hl_comp <- our_40yr_hl |>
    left_join(aihw_hl_40yr, by = "sex") |>
    mutate(ratio_pct = our_40yr / aihw_40yr * 100)

  bind_rows(
    nhlagg_comp |> mutate(group = "NHL aggregate (direct)"),
    nhl_comp    |> mutate(group = "NHL subtypes (DLBCL+FL+MCL)"),
    hl_comp     |> mutate(group = "Hodgkin lymphoma")
  ) |> select(group, sex, our_40yr, aihw_40yr, ratio_pct)
}

# Table S5 (was S2d): NHL subtype sum / aggregate NHL ratio over time.
build_table_s5 <- function(proj_annual, hist_annual) {
  subtype_sum <- bind_rows(
    hist_annual |>
      filter(subtype %in% nhl_subtypes, year <= 2021) |>
      select(year, sex, subtype, cases = cases_obs),
    proj_annual |>
      filter(subtype %in% nhl_subtypes, P >= 2022) |>
      select(year = P, sex, subtype, cases = cases_mid)
  ) |>
    group_by(year, sex) |>
    summarise(subtype_total = sum(cases, na.rm = TRUE), .groups = "drop")

  agg_nhl <- bind_rows(
    hist_annual |> filter(subtype == "nhl", year <= 2021) |>
      select(year, sex, nhl_total = cases_obs),
    proj_annual |> filter(subtype == "nhl", P >= 2022) |>
      select(year = P, sex, nhl_total = cases_mid)
  )

  consistency <- subtype_sum |>
    inner_join(agg_nhl, by = c("year", "sex")) |>
    mutate(ratio_pct = subtype_total / nhl_total * 100)

  consistency |>
    filter(year %in% c(2003, 2010, 2021, 2030, 2040, 2045)) |>
    mutate(period    = if_else(year <= 2021, "Observed", "Projected"),
           ratio_pct = round(ratio_pct, 1)) |>
    select(year, period, sex, ratio_pct) |>
    pivot_wider(names_from = sex, values_from = ratio_pct,
                names_glue = "{sex}_pct") |>
    arrange(year)
}

# Table S6 (was S3a): suppressed cells per (subtype, sex).
build_table_s6 <- function(surv_obs) {
  surv_cleaning_summary <- surv_obs |>
    group_by(subtype, sex) |>
    summarise(total_cells  = n(),
              filled_cells = sum(was_filled),
              pct_filled   = round(100 * mean(was_filled), 1),
              .groups      = "drop")
  filled_ages <- surv_obs |>
    filter(period == "2017–2021", years_since_dx == 5, was_filled) |>
    group_by(subtype, sex) |>
    summarise(filled_age_groups = paste(age_group, collapse = "; "),
              .groups = "drop")
  surv_cleaning_summary |>
    left_join(filled_ages, by = c("subtype", "sex")) |>
    mutate(filled_age_groups = replace_na(filled_age_groups, "none"))
}

# Table S7 (was S3b): improvement rates summary (one row per subtype x sex x ysd) plus
# the full per-cell improvement_final.
build_table_s7 <- function(improvement_final) {
  list(
    summary = improvement_final |>
      group_by(subtype, sex, years_since_dx) |>
      summarise(mean_rate      = round(mean(annual_rate), 5),
                min_rate       = round(min(annual_rate), 5),
                max_rate       = round(max(annual_rate), 5),
                n_age_specific = sum(source == "age-specific"),
                n_fallback     = sum(source != "age-specific"),
                .groups        = "drop"),
    full = improvement_final
  )
}

# Table S3 (was S3c): back-estimated 2000 vs observed 2003 per subtype x sex.
build_table_s3 <- function(inc_back, inc_subtype) {
  map_dfr(sexes, function(sx) {
    back_2000 <- inc_back |>
      filter(year == 2000, sex == sx, subtype %in% nhl_subtypes) |>
      group_by(subtype) |>
      summarise(back_est_2000 = sum(count), .groups = "drop")
    obs_2003 <- inc_subtype |>
      filter(year == 2003, sex == sx, subtype %in% nhl_subtypes) |>
      group_by(subtype) |>
      summarise(observed_2003 = sum(count), .groups = "drop")
    back_2000 |>
      left_join(obs_2003, by = "subtype") |>
      mutate(sex   = sx,
             ratio = round(back_est_2000 / observed_2003, 3))
  }) |> arrange(sex, subtype)
}

# Table S12 (was S4): sensitivity scenarios at key years.
build_table_s12 <- function(prev_sensitivity) {
  # Total lymphoma = HL + aggregate NHL, taken directly from the draw-level
  # total (subtype "total", sex "persons") that prev_model.R now emits — the
  # same total used in table_3 and Figures 2/3. Do NOT re-sum the per-lymphoma
  # rows: prev_sensitivity carries both by-sex AND persons rows, so summing
  # over sex double-counts (persons = males + females), and DLBCL/FL/MCL are a
  # subset of NHL that must never enter the total.
  prev_sensitivity |>
    filter(year %in% c(2021, 2030, 2040, 2045),
           subtype == "total", sex == "persons") |>
    select(scenario, year, duration, total = prev_mid) |>
    pivot_wider(names_from = scenario, values_from = total) |>
    mutate(cons_vs_base_pct = (conservative - base) / base * 100,
           opt_vs_base_pct  = (optimistic  - base) / base * 100) |>
    arrange(duration, year)
}

# ----------------------
# Figure builders
# ----------------------

.theme_apc <- function() {
  theme_bw(base_size = 13) +
    theme(legend.position  = "none",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = "grey70"),
          strip.text       = element_text(face = "bold", size = 12),
          strip.text.y     = element_text(angle = 0, hjust = 0))
}

# Figures S1/S2/S3 (were S1a/b/c): APC age / period / cohort effects (long-format input).
.build_fig_s1 <- function(apc_effects, effect_name,
                          x_label, y_label, log_y = FALSE) {
  # All five series, headline tier first (HL, aggregate NHL) then the three
  # NHL subtypes — matching the two-tier convention in Figures 2/3, S7-S9, S12/S13.
  d <- apc_effects |>
    filter(effect == effect_name) |>
    mutate(subtype_label = subtype_labels[subtype])
  d$subtype_label <- factor(d$subtype_label,
                            levels = c("HL", "NHL", "DLBCL", "FL", "MCL"))

  d <- d |> mutate(sex_lab = factor(sex_labels[sex],
                                    levels = c("Females", "Males")))

  p <- ggplot(d, aes(x = value, y = rr, colour = sex_lab, fill = sex_lab)) +
    geom_ribbon(aes(ymin = rr_p025, ymax = rr_p975), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ subtype_label, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = sex_fig_colours) +
    scale_fill_manual(values = sex_fig_colours) +
    labs(x = x_label, y = y_label) +
    .theme_apc()

  if (effect_name %in% c("period", "cohort")) {
    p <- p + geom_hline(yintercept = 1.0, linetype = "dotted",
                        colour = "grey50", linewidth = 0.3)
  }
  if (log_y) {
    p <- p + scale_y_log10(labels = scales::label_number())
  } else {
    p <- p + scale_y_continuous(labels = scales::label_number())
  }
  # Five lymphomas in two columns leave a sixth cell empty; it carries the key,
  # as on Figures 1-3 (see _setup.R).
  entries <- c(sex_key(), list(key_band("95% credible interval")))
  if (effect_name %in% c("period", "cohort")) {
    entries <- c(entries, list(key_line("Rate ratio = 1", "grey50",
                                        lty = "dotted", lwd = 1.6)))
  }
  attach_figure_key(p, entries)
}

# Figure S4 (was S3a): survival curves for one age group (65-74) at dx_year = 2019. Uses
# point estimates (build_projected_surv_pt) since this figure shows the
# central trajectory by year-since-diagnosis.
build_figure_s4_curves <- function(surv_obs, improvement_final, build_fn) {
  # All five series (headline HL + aggregate NHL, then the three subtypes);
  # aggregate NHL survival (AIHW Book 11f1) drives the largest prevalence
  # component and must be shown.
  plot_series <- c("hodgkin", "nhl", nhl_subtypes)
  surv_curves <- map_dfr(plot_series, function(st) {
    map_dfr(sexes, function(sx) {
      build_fn(st, sx, dx_year = 2019, cap_year = improve_cap,
              surv_obs, improvement_final) |>
        filter(age_group == "65–74") |>
        mutate(subtype = st, sex = sx)
    })
  })
  surv_curves <- surv_curves |>
    mutate(subtype_label = subtype_labels[subtype],
           region        = if_else(years_since_dx <= 5,
                                   "Observed (AIHW)", "Extrapolated"))
  surv_curves$subtype_label <- factor(surv_curves$subtype_label,
                                      levels = c("HL", "NHL", "DLBCL", "FL", "MCL"))

  surv_curves_line <- bind_rows(
    surv_curves,
    surv_curves |> filter(years_since_dx == 5) |> mutate(region = "Extrapolated")
  ) |>
    mutate(region = factor(region, levels = c("Observed (AIHW)", "Extrapolated")))

  ggplot(surv_curves, aes(x = years_since_dx, colour = sex, fill = sex)) +
    geom_ribbon(aes(ymin = surv_lo, ymax = surv_hi),
                alpha = 0.15, colour = NA) +
    geom_line(data = surv_curves_line,
              aes(y = surv_mid, linetype = region,
                  group = interaction(sex, region)),
              linewidth = 0.7) +
    facet_wrap(~ subtype_label, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = sex_colours, labels = sex_labels) +
    scale_fill_manual(values = sex_colours, labels = sex_labels) +
    scale_linetype_manual(values = c("Observed (AIHW)" = "solid",
                                     "Extrapolated"    = "dashed")) +
    scale_x_continuous(breaks = c(1, 5, 10, 20, 30, 40)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Years since diagnosis",
         y = "Survival probability",
         colour = NULL, fill = NULL, linetype = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank())
}

# Figure S5 (was S3b): annotated worked example of a single modelled survival curve, making
# the counting-method survival assumptions explicit (co-author request):
# the sqrt(S(1)) half-year weight applied to the newly-diagnosed cohort, the
# observed AIHW window (years 1-5), and the constant-conditional-survival
# extrapolation to 40 years. Example: DLBCL, females, 65-74, diagnosed 2019.
build_figure_survival_example <- function(surv_obs, improvement_final, build_fn,
                                          subtype_name = "dlbcl", sex_name = "females",
                                          age_grp = "65–74", dx_year = 2019) {
  curve <- build_fn(subtype_name, sex_name, dx_year = dx_year, cap_year = improve_cap,
                    surv_obs, improvement_final) |>
    filter(age_group == age_grp) |>
    select(years_since_dx, surv_mid) |>
    arrange(years_since_dx)

  s1  <- curve$surv_mid[curve$years_since_dx == 1]
  s0h <- sqrt(s1)   # half-year weight for the year-of-diagnosis cohort

  obs      <- curve |> filter(years_since_dx <= 5)
  extr     <- bind_rows(curve |> filter(years_since_dx == 5),
                        curve |> filter(years_since_dx > 5))
  anchor   <- tibble(years_since_dx = 0, surv_mid = 1)
  obs_line <- bind_rows(anchor, obs)
  col      <- unname(line_colours["DLBCL"])

  ggplot() +
    annotate("rect", xmin = 0, xmax = 5,  ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.35) +
    annotate("rect", xmin = 5, xmax = 40, ymin = -Inf, ymax = Inf,
             fill = "grey95", alpha = 0.35) +
    geom_line(data = obs_line, aes(years_since_dx, surv_mid), colour = col, linewidth = 1) +
    geom_line(data = extr, aes(years_since_dx, surv_mid), colour = col,
              linewidth = 1, linetype = "dashed") +
    geom_point(data = obs, aes(years_since_dx, surv_mid), colour = col, size = 2.4) +
    geom_point(aes(x = 0, y = 1), colour = "grey30", size = 2.4) +
    geom_point(aes(x = 0.5, y = s0h), colour = "#b2182b", size = 3, shape = 18) +
    annotate("segment", x = 6, xend = 1.4, y = 0.95, yend = s0h + 0.005,
             colour = "#b2182b", linewidth = 0.4, arrow = arrow(length = unit(0.15, "cm"))) +
    annotate("text", x = 6.3, y = 0.965, hjust = 0, size = 3.1, colour = "#b2182b",
             label = "Half-year weight for the newly\ndiagnosed cohort: S(0) ≈ √S(1)") +
    annotate("text", x = 2.5, y = 0.22, hjust = 0.5, size = 3.1, colour = "grey30",
             label = "Observed AIHW survival\n(1–5 years)") +
    annotate("segment", x = 24, xend = 30, y = 0.52, yend = 0.34,
             colour = "grey45", linewidth = 0.4, arrow = arrow(length = unit(0.13, "cm"))) +
    annotate("text", x = 22, y = 0.62, hjust = 0.5, size = 3.1, colour = "grey30",
             label = "Extrapolated to 40 years:\nconstant conditional survival\n(geometric mean of years 3–5)") +
    scale_x_continuous(breaks = c(0, 1, 5, 10, 20, 30, 40)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Years since diagnosis", y = "Observed survival probability",
         title = "Worked example: modelled survival curve for the counting method",
         subtitle = "DLBCL, females, diagnosed 2019, aged 65–74 at diagnosis") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey40"))
}

# Figures S10 / S11 (were S3f / S3g): survival over time. Case-weighted (by recent incidence age
# distribution) survival at one horizon (5 or 10 years) by year of diagnosis,
# one line per lymphoma x sex. Split into two figures (5-year = S10, 10-year =
# S11) so each carries a single clean line per lymphoma. Headline tier (HL,
# aggregate NHL) solid, decomposition tier (DLBCL/FL/MCL) dashed. Supports the
# "prevalence growth is driven partly by improving survival" claim (QL).
build_figure_survival_trend <- function(surv_obs, improvement_final,
                                        inc_subtype, inc_agg, build_fn,
                                        horizon = 5) {
  # Weights per series: subtype counts, plus aggregate HL and aggregate NHL.
  agg_wt <- function(grp, st) inc_agg |>
    filter(cancer_group == grp, year %in% 2017:2021) |> mutate(subtype = st) |>
    group_by(subtype, sex, age_group) |> summarise(w = sum(count), .groups = "drop")
  wts <- bind_rows(
    inc_subtype |> filter(year %in% 2017:2021) |>
      group_by(subtype, sex, age_group) |> summarise(w = sum(count), .groups = "drop"),
    agg_wt("Hodgkin lymphoma",     "hodgkin"),
    agg_wt("Non-Hodgkin lymphoma", "nhl")
  )
  dx_years <- 2009:improve_cap
  # All five series (headline HL + aggregate NHL, then the subtypes).
  plot_series <- c("hodgkin", "nhl", nhl_subtypes)
  trend <- map_dfr(plot_series, function(st) map_dfr(sexes, function(sx) {
    w <- wts |> filter(subtype == st, sex == sx)
    map_dfr(dx_years, function(dy) {
      build_fn(st, sx, dx_year = dy, cap_year = improve_cap, surv_obs, improvement_final) |>
        filter(years_since_dx == horizon) |>
        left_join(w, by = "age_group") |>
        summarise(surv = sum(surv_mid * w) / sum(w), .groups = "drop") |>
        mutate(dx_year = dy, subtype = st, sex = sx)
    })
  }))
  # Panelled by lymphoma with the sexes separated by colour, matching Figures
  # 1-3 and S1-S3.
  trend <- trend |>
    mutate(subtype_label = factor(label_map[subtype],
                                  levels = c("HL", "NHL", "DLBCL", "FL", "MCL")),
           sex_lab = factor(sex_labels[sex], levels = c("Females", "Males")))

  p <- ggplot(trend, aes(dx_year, surv, colour = sex_lab)) +
    geom_vline(xintercept = hist_end + 0.5, linetype = "dotted",
               colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ subtype_label, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = sex_fig_colours) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Year of diagnosis", y = "Observed survival (case-weighted)") +
    theme_bw(base_size = 13) +
    theme(legend.position  = "none",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = "grey70"),
          strip.text       = element_text(face = "bold", size = 12))

  attach_figure_key(p, c(
    sex_key(),
    # The line marks the last observed year of diagnosis; improvement is
    # extrapolated beyond it to the cap year, NOT held constant there.
    list(key_line("Last observed year", "grey60",
                  lty = "dotted", lwd = 1.6))))
}

# Figures S12/S13 (were S4a/b): sensitivity prevalence figures.
build_figure_sensitivity <- function(prev_sensitivity, duration_value) {
  # Headline tier (HL, aggregate NHL) plus its decomposition (DLBCL/FL/MCL),
  # by sex; the leaked "persons" rows are excluded (this figure is by sex).
  draw_series <- c("hodgkin", "nhl", nhl_subtypes)
  d <- prev_sensitivity |>
    filter(duration == duration_value, subtype %in% draw_series, sex %in% sexes) |>
    mutate(subtype_label = subtype_labels[subtype],
           prev_mid_k    = prev_mid / 1000)
  d$subtype_label <- factor(d$subtype_label,
                            levels = c("HL", "NHL", "DLBCL", "FL", "MCL"))
  d$scenario <- factor(d$scenario,
                       levels = c("conservative", "base", "optimistic"))

  ggplot(d, aes(x = year, y = prev_mid_k,
                colour = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = 2021.5, linetype = "dotted",
               colour = "grey50", linewidth = 0.3) +
    facet_grid(subtype_label ~ sex, scales = "free_y",
               labeller = labeller(sex = sex_labels)) +
    scale_colour_manual(values = scenario_colours, labels = scenario_labels) +
    scale_linetype_manual(
      values = c(conservative = "dashed", base = "solid", optimistic = "dotted"),
      labels = scenario_labels
    ) +
    labs(x = "Year", y = "Prevalence (thousands)",
         colour = NULL, linetype = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank(),
          strip.text.y     = element_text(angle = 0, hjust = 0))
}

# ----------------------
# Figure S6: population pyramid, 2021 vs 2045 (ABS Series B)
# ----------------------
# Added round 3 [TEG]: illustrates the population-ageing driver. 2021 as
# filled bars, 2045 projection as outlines; males left, females right.

build_figure_population_pyramid <- function(pop_hist, pop_proj) {
  prep <- function(df, yr_val, lab) df |>
    dplyr::filter(year == yr_val, sex %in% c("males", "females")) |>
    dplyr::group_by(sex, age_group) |>
    dplyr::summarise(population = sum(population), .groups = "drop") |>
    dplyr::mutate(year = lab)
  d <- dplyr::bind_rows(prep(pop_hist, 2021, "2021"),
                        prep(pop_proj, 2045, "2045")) |>
    dplyr::mutate(pop_m  = population / 1e6,
                  signed = ifelse(sex == "males", -pop_m, pop_m))
  ord <- unique(d$age_group)
  ord <- ord[order(suppressWarnings(as.numeric(sub("[\u2013+].*", "", ord))))]
  d$age_group <- factor(d$age_group, levels = ord)
  ggplot() +
    geom_col(data = dplyr::filter(d, year == "2021"),
             aes(x = signed, y = age_group, fill = sex),
             width = 0.78, alpha = 0.85) +
    geom_col(data = dplyr::filter(d, year == "2045"),
             aes(x = signed, y = age_group,
                 colour = "2045 (ABS Series B projection)"),
             fill = NA, width = 0.78, linewidth = 0.55) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.3) +
    scale_fill_manual(values = c(males = "#9ecae1", females = "#fdc086"),
                      labels = c(males = "Males, 2021",
                                 females = "Females, 2021"), name = NULL) +
    scale_colour_manual(values = c("2045 (ABS Series B projection)" = "grey15"),
                        name = NULL) +
    scale_x_continuous(labels = function(x) abs(x)) +
    labs(x = "Population (millions)", y = "Age group") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

# ----------------------
# Driver
# ----------------------

run_supplement <- function(save_dir = "output", verbose = TRUE,
                           prev_model_path = "code/prev_model.R") {
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  # Source prev_model.R for build_surv_obs / build_improvement_final /
  # build_projected_surv_pt helpers (small functions, no side effects
  # because of the auto-run guard).
  source(prev_model_path, local = TRUE)

  # 1. Load inputs from disk
  inputs <- list(
    proj_annual        = read_csv(file.path(save_dir, "incidence_projections.csv"),
                                  show_col_types = FALSE),
    hist_annual        = read_csv(file.path(save_dir, "incidence_historical_asr.csv"),
                                  show_col_types = FALSE),
    apc_effects        = read_csv(file.path(save_dir, "apc_effects.csv"),
                                  show_col_types = FALSE),
    prev_full          = read_csv(file.path(save_dir, "prevalence_projections.csv"),
                                  show_col_types = FALSE),
    prev_sensitivity   = read_csv(file.path(save_dir, "prevalence_sensitivity.csv"),
                                  show_col_types = FALSE),
    prev_10yr          = read_csv(file.path(save_dir, "prevalence_validation_10yr.csv"),
                                  show_col_types = FALSE),

    surv_raw           = read_csv("data/survival.csv",            show_col_types = FALSE),
    inc_subtype        = read_csv("data/incidence_subtype.csv",   show_col_types = FALSE),
    inc_agg            = read_csv("data/incidence_agg.csv",       show_col_types = FALSE),
    inc_agg_proj       = read_csv("data/incidence_agg_proj.csv",  show_col_types = FALSE),
    inc_back           = read_csv("data/incidence_subtype_back.csv", show_col_types = FALSE),
    prev_subtype_aihw  = read_csv("data/prevalence_subtype.csv",  show_col_types = FALSE),
    prev_agg_aihw      = read_csv("data/prevalence_agg.csv",      show_col_types = FALSE)
  )

  pop_hist <- read_csv("data/pop_hist.csv", show_col_types = FALSE)
  pop_proj <- read_csv("data/pop_proj.csv", show_col_types = FALSE)

  # 2. Rebuild surv_obs and improvement_final (cheap; needed for Figures S4/S5
  # and Figure S2)
  surv_obs          <- build_surv_obs(inputs$surv_raw)
  improvement_final <- build_improvement_final(surv_obs)

  # 3. Build tables
  table_s4 <- build_table_s4(inputs$proj_annual, inputs$inc_agg_proj)
  table_s8 <- build_table_s8(inputs$prev_full, inputs$prev_10yr,
                               inputs$prev_subtype_aihw, inputs$prev_agg_aihw)
  table_s9 <- build_table_s9(inputs$prev_full, inputs$prev_agg_aihw)
  table_s5 <- build_table_s5(inputs$proj_annual, inputs$hist_annual)
  table_s6 <- build_table_s6(surv_obs)
  s3b       <- build_table_s7(improvement_final)
  table_s3 <- build_table_s3(inputs$inc_back, inputs$inc_subtype)
  table_s12  <- build_table_s12(inputs$prev_sensitivity)

  write_csv(table_s4, file.path(save_dir, "table_s4_incidence_aihw.csv"))
  write_csv(table_s8, file.path(save_dir, "table_s8_prevalence_validation.csv"))
  write_csv(table_s9, file.path(save_dir, "table_s9_40yr_prevalence.csv"))
  write_csv(table_s5, file.path(save_dir, "table_s5_consistency.csv"))
  write_csv(table_s6, file.path(save_dir, "table_s6_survival_cleaning.csv"))
  write_csv(s3b$summary, file.path(save_dir, "table_s7_improvement_rates.csv"))
  write_csv(s3b$full,    file.path(save_dir, "table_s7_improvement_rates_full.csv"))
  write_csv(table_s3, file.path(save_dir, "table_s3_back_estimation.csv"))
  write_csv(table_s12,  file.path(save_dir, "table_s12_sensitivity.csv"))
  # Table S1: disease definitions (AIHW Book 11h main framework + ICD-10
  # framework). Static reference data; source verified 2026-09-04 [MB].
  table_s1 <- tibble::tribble(
    ~lymphoma, ~icd10, ~histology_icdo3,
    "Hodgkin lymphoma",              "C81",     "9650\u20139655, 9659, 9661\u20139665, 9667",
    "Non Hodgkin lymphoma (aggregate)", "C82\u2013C86", "\u2013",
    "Diffuse large B-cell lymphoma", "\u2013", "9675, 9680, 9684, 9738",
    "Follicular lymphoma",           "\u2013", "9597, 9690, 9691, 9695, 9698",
    "Mantle cell lymphoma",          "\u2013", "9673")
  write_csv(table_s1, file.path(save_dir, "table_s1_disease_definitions.csv"))

  # Table S11: prevalence by sex at the additional durations (2, 3, 10 yr)
  # with draw-level change CrIs [QL round 3]. Depends on prev_model.R's
  # table_2_change_cri.csv export.
  table_s11 <- read_csv(file.path(save_dir, "table_2_change_cri.csv"),
                        show_col_types = FALSE) |>
    dplyr::filter(duration %in% c(2, 3, 10)) |>
    dplyr::arrange(duration, subtype, sex)
  write_csv(table_s11, file.path(save_dir, "table_s11_duration_by_sex.csv"))


  # 4. Build figures
  fig_s1 <- .build_fig_s1(inputs$apc_effects, "age",
                           "Age (midpoint)", "Rate per 100,000",
                           log_y = TRUE)
  fig_s2 <- .build_fig_s1(inputs$apc_effects, "period",
                           "Calendar year", "Period relative risk")
  fig_s3 <- .build_fig_s1(inputs$apc_effects, "cohort",
                           "Birth cohort (year)", "Cohort relative risk")
  fig_s4    <- build_figure_s4_curves(surv_obs, improvement_final, build_projected_surv_pt)
  fig_s5 <- build_figure_survival_example(surv_obs, improvement_final,
                                              build_projected_surv_pt)
  fig_s10    <- build_figure_survival_trend(surv_obs, improvement_final,
                                            inputs$inc_subtype, inputs$inc_agg,
                                            build_projected_surv_pt, horizon = 5)
  fig_s11    <- build_figure_survival_trend(surv_obs, improvement_final,
                                            inputs$inc_subtype, inputs$inc_agg,
                                            build_projected_surv_pt, horizon = 10)

  # Additional duration prevalence trajectories (2-year, 3-year, 10-year)
  # - same style as Figures 2 and 3 in the main paper. The 10-year
  # version overlays AIHW observed points; 2- and 3-year do not (AIHW
  # does not publish prevalence at those durations).
  fig_prev2  <- build_prev_fig(inputs$prev_full, duration_value = 2)
  fig_prev3  <- build_prev_fig(inputs$prev_full, duration_value = 3)
  fig_prev10 <- build_prev_fig(inputs$prev_full, duration_value = 10,
                               prev_subtype_aihw = inputs$prev_subtype_aihw,
                               prev_agg_aihw     = inputs$prev_agg_aihw)

  fig_s12 <- build_figure_sensitivity(inputs$prev_sensitivity, 5)
  fig_s13 <- build_figure_sensitivity(inputs$prev_sensitivity, 40)

  # One portrait 3x2 canvas for every re-panelled lymphoma figure (_setup.R).
  save_fig(fig_s1, file.path(save_dir, "figure_s1_age_effects"),
           width = fig_panel_w, height = fig_panel_h)
  save_fig(fig_s2, file.path(save_dir, "figure_s2_period_effects"),
           width = fig_panel_w, height = fig_panel_h)
  save_fig(fig_s3, file.path(save_dir, "figure_s3_cohort_effects"),
           width = fig_panel_w, height = fig_panel_h)
  save_fig(fig_s4, file.path(save_dir, "figure_s4_survival_curves"),
           width = 10, height = 11)   # 5 series -> 3 rows x 2 cols
  save_fig(fig_s5, file.path(save_dir, "figure_s5_survival_example"),
           width = 8.5, height = 5.2)
  save_fig(fig_prev2, file.path(save_dir, "figure_s7_prevalence_2yr"),
           width = fig_panel_w, height = fig_panel_h)
  save_fig(fig_prev3, file.path(save_dir, "figure_s8_prevalence_3yr"),
           width = fig_panel_w, height = fig_panel_h)
  save_fig(fig_prev10, file.path(save_dir, "figure_s9_prevalence_10yr"),
           width = fig_panel_w, height = fig_panel_h)
  save_fig(fig_s10, file.path(save_dir, "figure_s10_survival_trend_5yr"),
           width = fig_panel_w, height = fig_panel_h)
  save_fig(fig_s11, file.path(save_dir, "figure_s11_survival_trend_10yr"),
           width = fig_panel_w, height = fig_panel_h)
  save_fig(fig_s12, file.path(save_dir, "figure_s12_sensitivity_5yr"),
           width = 10, height = 10)
  save_fig(fig_s13, file.path(save_dir, "figure_s13_sensitivity_40yr"),
           width = 10, height = 10)

  fig_s6 <- build_figure_population_pyramid(pop_hist, pop_proj)
  save_fig(fig_s6, file.path(save_dir, "figure_s6_population_pyramid"),
           width = 9, height = 5.5)

  if (verbose) {
    cat("\n=== Supplement complete ===\n")
    cat(sprintf("Tables: 11 CSVs written to %s\n", save_dir))
    cat(sprintf("Figures: 13 figures (PDF + PNG = 26 files) written to %s\n", save_dir))
  }

  invisible(list(
    tables = list(s2a = table_s4, s2b = table_s8, s2c = table_s9,
                  s2d = table_s5, s3a = table_s6,
                  s3b_summary = s3b$summary, s3b_full = s3b$full,
                  s3c = table_s3, s4 = table_s12),
    figures = list(s1a = fig_s1, s1b = fig_s2, s1c = fig_s3,
                   s3a = fig_s4, s3b_ex = fig_s5,
                   s3f_5yr = fig_s10, s3g_10yr = fig_s11,
                   s3c = fig_prev2, s3d = fig_prev3, s3e = fig_prev10,
                   s4a = fig_s12, s4b = fig_s13)
  ))
}

# ----------------------
# Auto-run guard
# ----------------------

if (!interactive() && sys.nframe() == 0) {
  invisible(run_supplement())
}
