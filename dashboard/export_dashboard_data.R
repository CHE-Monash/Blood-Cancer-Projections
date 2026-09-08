#!/usr/bin/env Rscript

# Build the deliberately public, browser-sized dataset used by the Quarto
# dashboard.
# Run from the repository root: Rscript dashboard/export_dashboard_data.R This script reads only committed aggregate/modelled outputs. It
# does not touch raw workbooks or fitted model objects.

read_required <- function(path, required_columns) {
  if (!file.exists(path)) {
    stop("Required dashboard input is missing: ", path, call. = FALSE)
  }

  value <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_columns, names(value))
  if (length(missing_columns) > 0L) {
    stop(
      "Dashboard input ", path, " is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  value
}

diagnosis_labels <- c(
  total = "Total lymphoma",
  hodgkin = "Hodgkin lymphoma (HL)",
  nhl = "Non-Hodgkin lymphoma (NHL)",
  dlbcl = "Diffuse large B-cell lymphoma (DLBCL)",
  follicular = "Follicular lymphoma (FL)",
  mantle_cell = "Mantle cell lymphoma (MCL)"
)

label_diagnosis <- function(code) {
  labels <- unname(diagnosis_labels[code])
  if (anyNA(labels)) {
    stop(
      "Unrecognised diagnosis code(s): ",
      paste(unique(code[is.na(labels)]), collapse = ", "),
      call. = FALSE
    )
  }
  labels
}

label_sex <- function(value) {
  labels <- c(persons = "Persons", females = "Females", males = "Males")
  result <- unname(labels[value])
  if (anyNA(result)) {
    stop(
      "Unrecognised sex value(s): ",
      paste(unique(value[is.na(result)]), collapse = ", "),
      call. = FALSE
    )
  }
  result
}

series_group <- function(code) {
  ifelse(code %in% c("total", "hodgkin", "nhl"), "Headline", "NHL subtype")
}

make_rows <- function(
  source,
  measure,
  year,
  diagnosis_code,
  sex,
  duration_years,
  period,
  estimate,
  lower,
  upper,
  unit,
  source_output
) {
  data.frame(
    measure = measure,
    diagnosis_code = diagnosis_code,
    diagnosis = label_diagnosis(diagnosis_code),
    series_group = series_group(diagnosis_code),
    sex = label_sex(sex),
    duration_years = duration_years,
    year = as.integer(year),
    period = period,
    estimate = as.numeric(estimate),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    unit = unit,
    source_output = source_output,
    stringsAsFactors = FALSE
  )
}

incidence_projection <- read_required(
  "output/incidence_projections.csv",
  c(
    "P", "cases_mid", "cases_p025", "cases_p975",
    "asr_mid", "asr_p025", "asr_p975", "subtype", "sex"
  )
)

incidence_history <- read_required(
  "output/incidence_historical_asr.csv",
  c("year", "cases_obs", "asr_obs", "subtype", "sex")
)

prevalence <- read_required(
  "output/prevalence_projections.csv",
  c(
    "year", "prev_mid", "prev_p025", "prev_p975",
    "duration", "subtype", "sex"
  )
)

# The historical output currently carries convenience rows after the final
# observed year. They are not observed data and must not enter the public
# dashboard as such.
incidence_history <- incidence_history[
  incidence_history$year <= 2021 & incidence_history$asr_obs > 0,
  ,
  drop = FALSE
]

dashboard_data <- rbind(
  make_rows(
    incidence_history,
    "Incidence (cases)",
    incidence_history$year,
    incidence_history$subtype,
    incidence_history$sex,
    NA_integer_,
    "Historical",
    incidence_history$cases_obs,
    NA_real_,
    NA_real_,
    "people diagnosed",
    "output/incidence_historical_asr.csv"
  ),
  make_rows(
    incidence_projection,
    "Incidence (cases)",
    incidence_projection$P,
    incidence_projection$subtype,
    incidence_projection$sex,
    NA_integer_,
    "Projected",
    incidence_projection$cases_mid,
    incidence_projection$cases_p025,
    incidence_projection$cases_p975,
    "people diagnosed",
    "output/incidence_projections.csv"
  ),
  make_rows(
    incidence_history,
    "Incidence (age-standardised rate)",
    incidence_history$year,
    incidence_history$subtype,
    incidence_history$sex,
    NA_integer_,
    "Historical",
    incidence_history$asr_obs,
    NA_real_,
    NA_real_,
    "per 100,000 population",
    "output/incidence_historical_asr.csv"
  ),
  make_rows(
    incidence_projection,
    "Incidence (age-standardised rate)",
    incidence_projection$P,
    incidence_projection$subtype,
    incidence_projection$sex,
    NA_integer_,
    "Projected",
    incidence_projection$asr_mid,
    incidence_projection$asr_p025,
    incidence_projection$asr_p975,
    "per 100,000 population",
    "output/incidence_projections.csv"
  ),
  make_rows(
    prevalence,
    "Prevalence (people)",
    prevalence$year,
    prevalence$subtype,
    prevalence$sex,
    prevalence$duration,
    ifelse(prevalence$year <= 2021, "Estimated", "Projected"),
    prevalence$prev_mid,
    prevalence$prev_p025,
    prevalence$prev_p975,
    "people living after diagnosis",
    "output/prevalence_projections.csv"
  )
)

measure_order <- c(
  "Incidence (cases)",
  "Incidence (age-standardised rate)",
  "Prevalence (people)"
)
diagnosis_order <- names(diagnosis_labels)
sex_order <- c("Persons", "Females", "Males")

dashboard_data <- dashboard_data[
  order(
    match(dashboard_data$measure, measure_order),
    match(dashboard_data$diagnosis_code, diagnosis_order),
    match(dashboard_data$sex, sex_order),
    dashboard_data$duration_years,
    dashboard_data$year
  ),
  ,
  drop = FALSE
]
row.names(dashboard_data) <- NULL

# Release checks: one row per selectable result, valid uncertainty ordering,
# and the expected temporal boundary between observed/estimated and projected.
key <- with(
  dashboard_data,
  paste(measure, diagnosis_code, sex, duration_years, year, sep = "|")
)
if (anyDuplicated(key)) {
  stop("Dashboard release contains duplicate result keys.", call. = FALSE)
}
if (anyNA(dashboard_data$estimate) || any(!is.finite(dashboard_data$estimate))) {
  stop("Dashboard release contains missing or non-finite estimates.", call. = FALSE)
}
with_interval <- !is.na(dashboard_data$lower) & !is.na(dashboard_data$upper)
if (!all(
  dashboard_data$lower[with_interval] <= dashboard_data$estimate[with_interval] &
    dashboard_data$estimate[with_interval] <= dashboard_data$upper[with_interval]
)) {
  stop("Dashboard release contains an invalid uncertainty interval.", call. = FALSE)
}
if (any(dashboard_data$period == "Historical" & dashboard_data$year > 2021)) {
  stop("Post-2021 rows have been labelled as historical.", call. = FALSE)
}
if (!all(measure_order %in% dashboard_data$measure)) {
  stop("One or more dashboard measures are missing.", call. = FALSE)
}
if (!all(c(2021L, 2045L) %in% dashboard_data$year)) {
  stop("The release does not contain both the 2021 baseline and 2045 horizon.", call. = FALSE)
}

dir.create("dashboard/data", recursive = TRUE, showWarnings = FALSE)
write.csv(
  dashboard_data,
  "dashboard/data/lymphoma-dashboard-data.csv",
  row.names = FALSE,
  na = ""
)

dictionary <- data.frame(
  variable = names(dashboard_data),
  definition = c(
    "Outcome displayed in the dashboard.",
    "Machine-readable lymphoma series code.",
    "Reader-facing lymphoma series name.",
    "Headline total/aggregate series or NHL subtype decomposition.",
    "Sex reported by the source result.",
    "Prevalence duration in years; blank for incidence measures.",
    "Calendar year.",
    "Historical observed incidence, estimated prevalence, or projected result.",
    "Point estimate. Modelled values are posterior medians.",
    "Lower limit of the 95% credible interval; blank for observed incidence.",
    "Upper limit of the 95% credible interval; blank for observed incidence.",
    "Unit attached to the selected measure.",
    "Canonical analysis output from which the release row was derived."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  dictionary,
  "dashboard/data/data-dictionary.csv",
  row.names = FALSE,
  na = ""
)

message(
  "Wrote ", nrow(dashboard_data), " dashboard rows to ",
  "dashboard/data/lymphoma-dashboard-data.csv"
)
