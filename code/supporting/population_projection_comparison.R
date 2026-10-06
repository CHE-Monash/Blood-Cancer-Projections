# ================
# Blood Cancer Projections (lymphoma): sensitivity of constant-rate NHL
# diagnoses to the ABS population projection release
# ================
# Luo et al. (Lancet Reg Health West Pac 2026;73:101946) project cancer
# prevalence on the ABS 2017 (base) to 2066 population projections (2018
# release). This study uses the 2022 (base) to 2071 release, Table B9,
# Series B. This script asks how much of the difference in projected NHL
# diagnosis growth between the two studies could be the population input
# alone: it holds age- and sex-specific NHL rates at their 2017-2021 mean
# and applies them to BOTH projections. Luo et al. do not state which series
# they used; Series B (medium) is assumed, matching this study.
#
# READ-ONLY with respect to the analysis. Nothing is fitted or sampled, and
# no committed input or output of the pipeline is changed.
#
# Inputs
#   data/raw/ABS/32220_Table_B9_2017base.xls   ABS 3222.0, 2017 (base) to
#       2066, Table B9 (Series B), 2018-release file, downloaded from the ABS site on
#       7 Oct 2026 (MD5 7ca0aceda4402b568badf43d88238247). data/raw/ is not
#       committed; if the file is absent the script falls back to the derived
#       CSV below, which is committed.
#   data/pop_proj_2017base.csv   derived from it in the ten age bands used
#       throughout (written by this script when the raw file is present)
#   data/pop_proj.csv, data/pop_hist.csv, data/incidence_agg.csv
#
# Writes output/table_population_projection_comparison.csv
#
#   Rscript code/supporting/population_projection_comparison.R
# ================

grDevices::pdf(NULL)
suppressPackageStartupMessages({
  source("code/apc_model.R")
  source("code/prev_model.R")
})

raw_2017 <- "data/raw/ABS/32220_Table_B9_2017base.xls"
csv_2017 <- "data/pop_proj_2017base.csv"

.band <- function(age) dplyr::case_when(
  age <= 14 ~ "0–14",  age <= 24 ~ "15–24", age <= 34 ~ "25–34", age <= 44 ~ "35–44",
  age <= 54 ~ "45–54", age <= 64 ~ "55–64", age <= 74 ~ "65–74", age <= 84 ~ "75–84",
  TRUE ~ "85+")

# Same parsing as import_data.R uses for the 2022-base Table B9: the sheet
# layout is identical, with the series description in row 1 and data from
# row 11. Males and females are kept; the Persons series are dropped.
build_pop_proj_2017base <- function(path = raw_2017) {
  all <- suppressMessages(readxl::read_excel(path, sheet = "Data1", col_names = FALSE))
  nm  <- as.character(unlist(all[1, ], use.names = FALSE)); nm[1] <- "date"
  dat <- all[-(1:10), ]
  names(dat) <- nm
  dat |>
    dplyr::mutate(year = as.integer(format(as.Date(as.numeric(date), origin = "1899-12-30"), "%Y"))) |>
    dplyr::select(-date) |>
    tidyr::pivot_longer(-year, names_to = "description", values_to = "population") |>
    dplyr::mutate(
      population = as.numeric(population),
      sex = dplyr::case_when(stringr::str_detect(description, "Male")   ~ "males",
                             stringr::str_detect(description, "Female") ~ "females",
                             TRUE ~ NA_character_),
      age = dplyr::if_else(stringr::str_detect(description, "100 and over"), "100",
                           stringr::str_extract(description, "(?<=;  )[0-9]+(?= ;)")),
      age = as.integer(age)) |>
    dplyr::filter(!is.na(sex), !is.na(age)) |>
    dplyr::mutate(age_group = .band(age)) |>
    dplyr::group_by(year, sex, age_group) |>
    dplyr::summarise(population = sum(population), .groups = "drop") |>
    dplyr::arrange(year, sex, age_group)
}

run_population_projection_comparison <- function(save_dir = "output",
                                                  years = c(2024, 2025, 2045, 2049, 2050)) {
  if (file.exists(raw_2017)) {
    pop_2017 <- build_pop_proj_2017base(raw_2017)
    stopifnot(nrow(pop_2017) == 50 * 2 * length(age_groups))
    readr::write_csv(pop_2017, csv_2017)
  } else {
    message("Raw 2017-base file absent; reading ", csv_2017)
    pop_2017 <- readr::read_csv(csv_2017, show_col_types = FALSE)
  }
  pop_2022 <- readr::read_csv("data/pop_proj.csv", show_col_types = FALSE)
  pop_hist <- readr::read_csv("data/pop_hist.csv", show_col_types = FALSE)
  inc_agg  <- readr::read_csv("data/incidence_agg.csv", show_col_types = FALSE)

  # Mean observed age- and sex-specific NHL rate, 2017-2021 (as in
  # prevalence_convention_comparison.R)
  rate <- purrr::map_dfr(sexes, function(sx) {
    d <- inc_agg |>
      dplyr::filter(cancer_group == "Non-Hodgkin lymphoma", sex == sx, year %in% 2017:2021) |>
      dplyr::group_by(age_group) |> dplyr::summarise(d = sum(count), .groups = "drop")
    y <- pop_hist |> dplyr::filter(sex == sx, year %in% 2017:2021) |>
      dplyr::group_by(age_group) |> dplyr::summarise(y = sum(population), .groups = "drop")
    dplyr::inner_join(d, y, by = "age_group") |> dplyr::mutate(sex = sx, rate = d / y)
  })
  stopifnot(nrow(rate) == 2 * length(age_groups))

  summarise_pop <- function(pop, label) {
    pop |>
      dplyr::filter(year %in% years, sex %in% sexes) |>
      dplyr::inner_join(rate[, c("sex", "age_group", "rate")], by = c("sex", "age_group")) |>
      dplyr::group_by(year) |>
      dplyr::summarise(
        population_total = sum(population),
        population_65plus = sum(population[age_group %in% c("65–74", "75–84", "85+")]),
        diagnoses_constant_rate = round(sum(rate * population)), .groups = "drop") |>
      dplyr::mutate(share_65plus_pct = round(100 * population_65plus / population_total, 1),
                    projection = label, .before = 1)
  }
  out <- dplyr::bind_rows(
    summarise_pop(pop_2022, "ABS 2022-base (this study)"),
    summarise_pop(pop_2017, "ABS 2017-base (Luo 2026)")) |>
    dplyr::mutate(subtype = "nhl", sex = "persons",
                  rates = "mean observed age- and sex-specific rate, 2017-2021",
                  series = "B") |>
    dplyr::select(subtype, sex, projection, series, year, population_total, population_65plus,
                  share_65plus_pct, diagnoses_constant_rate, rates)
  readr::write_csv(out, file.path(save_dir, "table_population_projection_comparison.csv"))
  print(as.data.frame(out[, c("projection", "year", "population_total", "share_65plus_pct",
                              "diagnoses_constant_rate")]))
  g <- function(lab, a, b) { x <- out[out$projection == lab, ]
    (x$diagnoses_constant_rate[x$year == b] / x$diagnoses_constant_rate[x$year == a] - 1) * 100 }
  for (lab in unique(out$projection))
    cat(sprintf("%-28s constant-rate growth in NHL diagnoses: 2025-2045 %+.1f%%; 2025-2050 %+.1f%%; 2024-2049 %+.1f%%\n",
                lab, g(lab, 2025, 2045), g(lab, 2025, 2050), g(lab, 2024, 2049)))
  lvl <- function(y) { a <- out[out$projection == "ABS 2017-base (Luo 2026)" & out$year == y, ]
    b <- out[out$projection == "ABS 2022-base (this study)" & out$year == y, ]
    sprintf("%d: population %+.1f%%, 65+ %+.1f%%, diagnoses %+.1f%%", y,
            100 * (a$population_total / b$population_total - 1),
            100 * (a$population_65plus / b$population_65plus - 1),
            100 * (a$diagnoses_constant_rate / b$diagnoses_constant_rate - 1)) }
  cat("2017-base relative to 2022-base --", lvl(2025), "|", lvl(2050), "\n")
  invisible(out)
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_population_projection_comparison())
  invisible(grDevices::dev.off())
}
