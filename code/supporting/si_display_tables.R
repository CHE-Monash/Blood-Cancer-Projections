# ================
# Blood Cancer Projections (lymphoma): Supporting Information display tables
# ================
# Builds the two Supporting Information tables that are not a plain copy of
# one existing export:
#
#   output/si_table_holdout.csv
#       Incidence holdout (temporal backtesting), one row per lymphoma and sex.
#       Source: output/table_holdout_validation.csv (code/supporting/holdout.R).
#
#   output/si_table_prevalence_validation.csv
#       Prevalence against AIHW published estimates: 5- and 10-year
#       (2012-2021, mean and range of annual differences) and 40-year (2021).
#       Sources: output/table_s8_prevalence_validation.csv and
#                output/table_s9_40yr_prevalence.csv (code/supplement.R).
#
# READ-ONLY with respect to the scientific outputs: this script reads three
# existing CSVs and writes two new ones. It fits nothing, samples nothing and
# does not source or run any model driver. Existing files in output/ are not
# modified.
#
#   Rscript code/supporting/si_display_tables.R
# ================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

lymphoma_levels <- c("hodgkin", "nhl", "dlbcl", "follicular", "mantle_cell")
lymphoma_labels <- c(hodgkin = "HL", nhl = "NHL", dlbcl = "DLBCL",
                     follicular = "FL", mantle_cell = "MCL")
sex_levels <- c("females", "males")
sex_lab    <- c(females = "Females", males = "Males")

build_si_table_holdout <- function(out_dir = "output") {
  src <- read_csv(file.path(out_dir, "table_holdout_validation.csv"),
                  show_col_types = FALSE)
  stopifnot(nrow(src) == 10, all(src$n_years == 10),
            setequal(src$subtype, lymphoma_levels),
            setequal(src$sex, sex_levels))
  # coverage_cases is exported as a percentage of the held-out years; the
  # table reports it as a count of years.
  years_in <- src$coverage_cases / 100 * src$n_years
  stopifnot(all(abs(years_in - round(years_in)) < 1e-9))
  src |>
    mutate(subtype = factor(subtype, levels = lymphoma_levels),
           sex     = factor(sex, levels = sex_levels)) |>
    arrange(subtype, sex) |>
    transmute(
      lymphoma              = unname(lymphoma_labels[as.character(subtype)]),
      sex                   = unname(sex_lab[as.character(sex)]),
      observed_2021         = cases_obs_2021,
      projected_2021        = round(cases_mid_2021),
      projected_2021_p025   = round(cases_p025_2021),
      projected_2021_p975   = round(cases_p975_2021),
      # signed error at the 2021 endpoint: (projected - observed) / observed
      error_2021_pct        = round(cases_pe_2021, 1),
      # mean absolute percentage error across the ten held-out years
      mape_diagnoses_pct    = round(cases_mape, 1),
      mape_asr_pct          = round(asr_mape, 1),
      years_within_cri      = as.integer(round(coverage_cases / 100 * n_years)),
      years_assessed        = as.integer(n_years)
    )
}

build_si_table_prevalence_validation <- function(out_dir = "output") {
  s8 <- read_csv(file.path(out_dir, "table_s8_prevalence_validation.csv"),
                 show_col_types = FALSE)
  s9 <- read_csv(file.path(out_dir, "table_s9_40yr_prevalence.csv"),
                 show_col_types = FALSE, na = c("", "NA"))
  stopifnot(nrow(s8) == 10, setequal(s8$subtype, lymphoma_levels))
  # The AIHW publishes 40-year prevalence for HL and aggregate NHL only, and
  # by sex. The subtype-sum rows of table_s9 are a coverage figure, not a
  # comparison with a published value, and are not carried into this table.
  group_map <- c("Hodgkin lymphoma" = "hodgkin",
                 "NHL aggregate (direct)" = "nhl")
  d40 <- s9 |>
    filter(group %in% names(group_map), sex %in% sex_levels) |>
    transmute(subtype = unname(group_map[group]), sex,
              dur40_2021 = round(ratio_pct - 100, 1))
  stopifnot(nrow(d40) == 4, !anyNA(d40$dur40_2021))
  s8 |>
    left_join(d40, by = c("subtype", "sex")) |>
    mutate(subtype = factor(subtype, levels = lymphoma_levels),
           sex     = factor(sex, levels = sex_levels)) |>
    arrange(subtype, sex) |>
    transmute(
      lymphoma = unname(lymphoma_labels[as.character(subtype)]),
      sex      = unname(sex_lab[as.character(sex)]),
      dur5_mean, dur5_min, dur5_max,
      dur10_mean, dur10_min, dur10_max,
      dur40_2021   # NA = no published AIHW comparator, not zero
    )
}

run_si_display_tables <- function(out_dir = "output") {
  h <- build_si_table_holdout(out_dir)
  v <- build_si_table_prevalence_validation(out_dir)
  write_csv(h, file.path(out_dir, "si_table_holdout.csv"))
  write_csv(v, file.path(out_dir, "si_table_prevalence_validation.csv"), na = "")
  cat(sprintf("Wrote si_table_holdout.csv (%d rows) and si_table_prevalence_validation.csv (%d rows) to %s\n",
              nrow(h), nrow(v), out_dir))
  invisible(list(holdout = h, prevalence_validation = v))
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_si_display_tables())
}
