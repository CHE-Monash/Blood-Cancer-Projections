# =================================================================
# Tests for the presentation artefacts (Supporting Information display
# tables, figure data workbooks, dashboard additional results).
#
# Kept apart from code/_tests.R on purpose. _tests.R re-runs a pipeline phase
# when an output is older than its source file; these tests must be runnable
# without that, so this file sources no model driver through the harness and
# NEVER fits, samples or regenerates anything. It reads files and compares.
#
# Run the producers first:
#   Rscript code/supporting/si_display_tables.R
#   Rscript code/supporting/figure_data_export.R
#   Rscript dashboard/export_supplementary.R
# then:
#   Rscript code/_tests_presentation.R
# =================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr) })
stopifnot(requireNamespace("openxlsx"))
rd <- function(f) read_csv(f, show_col_types = FALSE, na = c("", "NA"))
lab <- c(hodgkin = "HL", nhl = "NHL", dlbcl = "DLBCL", follicular = "FL",
         mantle_cell = "MCL")
sx  <- c(females = "Females", males = "Males")

cat("\n========== presentation artefacts ==========\n")

# Test D1: holdout display table equals its source, cell by cell
h   <- rd("output/si_table_holdout.csv")
src <- rd("output/table_holdout_validation.csv") |>
  mutate(lymphoma = lab[subtype], sex = sx[sex])
stopifnot(nrow(h) == 10, !anyDuplicated(h[, c("lymphoma", "sex")]),
          identical(unique(h$lymphoma), unname(lab)))
j <- inner_join(h, src, by = c("lymphoma", "sex"))
stopifnot(nrow(j) == 10,
          all(j$observed_2021 == j$cases_obs_2021),
          all(j$projected_2021 == round(j$cases_mid_2021)),
          all(j$projected_2021_p025 == round(j$cases_p025_2021)),
          all(j$projected_2021_p975 == round(j$cases_p975_2021)),
          all(abs(j$error_2021_pct - j$cases_pe_2021) <= 0.05 + 1e-9),
          all(abs(j$mape_diagnoses_pct - j$cases_mape) <= 0.05 + 1e-9),
          all(abs(j$mape_asr_pct - j$asr_mape) <= 0.05 + 1e-9),
          all(j$years_within_cri * 10 == j$coverage_cases),
          all(j$years_assessed == 10),
          all(j$projected_2021_p025 <= j$projected_2021),
          all(j$projected_2021 <= j$projected_2021_p975))
# the signed error is (projected - observed) / observed
stopifnot(all(abs(j$cases_pe_2021 -
                  (j$cases_mid_2021 - j$cases_obs_2021) / j$cases_obs_2021 * 100) < 1e-6))
cat("  Test D1 (holdout display table equals table_holdout_validation.csv)  PASS  10 rows\n")

# Test D2: prevalence validation display table equals its two sources
v  <- rd("output/si_table_prevalence_validation.csv")
s8 <- rd("output/table_s8_prevalence_validation.csv") |>
  mutate(lymphoma = lab[subtype], sex = sx[sex])
s9 <- rd("output/table_s9_40yr_prevalence.csv")
stopifnot(nrow(v) == 10, !anyDuplicated(v[, c("lymphoma", "sex")]))
j <- inner_join(v, s8, by = c("lymphoma", "sex"), suffix = c("", ".src"))
for (k in c("dur5_mean", "dur5_min", "dur5_max", "dur10_mean", "dur10_min", "dur10_max"))
  stopifnot(all(j[[k]] == j[[paste0(k, ".src")]]))
want <- s9 |>
  filter(group %in% c("Hodgkin lymphoma", "NHL aggregate (direct)"),
         sex %in% c("females", "males")) |>
  transmute(lymphoma = if_else(group == "Hodgkin lymphoma", "HL", "NHL"),
            sex = sx[sex], want = round(ratio_pct - 100, 1))
j <- left_join(v, want, by = c("lymphoma", "sex"))
stopifnot(sum(!is.na(v$dur40_2021)) == 4,
          all(is.na(v$dur40_2021) == is.na(j$want)),
          all(v$dur40_2021[!is.na(v$dur40_2021)] == j$want[!is.na(j$want)]),
          all(is.na(v$dur40_2021[v$lymphoma %in% c("DLBCL", "FL", "MCL")])))
# negative control: a subtype-sum row must not have leaked in as a comparison
stopifnot(!any(abs(v$dur40_2021) > 20, na.rm = TRUE))
cat("  Test D2 (prevalence validation table equals table_s8 + table_s9; subtypes blank at 40 years)  PASS\n")

# Test D3: figure workbooks hold the values in the source outputs
xl <- function(f, s) openxlsx::read.xlsx(f, sheet = s)
ip <- rd("output/incidence_projections.csv") |>
  filter(subtype %in% names(lab), sex %in% names(sx)) |>
  transmute(lymphoma = lab[subtype], sex = sx[sex], year = P,
            asr_mid, asr_p025, asr_p975)
f1 <- xl("output/figure_1_data.xlsx", "Modelled rate")
j  <- inner_join(f1 |> filter(series == "Projected"), ip, by = c("lymphoma", "sex", "year"))
stopifnot(nrow(j) == 10 * 24, nrow(filter(f1, series == "Projected")) == 240,
          all(abs(j$asr_per_100000 - j$asr_mid) < 1e-9),
          all(abs(j$cri_lower - j$asr_p025) < 1e-9),
          all(abs(j$cri_upper - j$asr_p975) < 1e-9))
fit <- filter(f1, series == "Fitted")
stopifnot(nrow(fit) == 4 * 32 + 6 * 19, all(fit$year <= 2021),
          all(fit$asr_per_100000 > 0))
ho <- rd("output/incidence_historical_asr.csv") |>
  filter(subtype %in% names(lab), sex %in% names(sx)) |>
  transmute(lymphoma = lab[subtype], sex = sx[sex], year, asr_obs)
ob <- xl("output/figure_1_data.xlsx", "Observed incidence")
j  <- inner_join(ob, ho, by = c("lymphoma", "sex", "year"))
stopifnot(nrow(j) == nrow(ob), nrow(ob) == 4 * 32 + 6 * 19,
          all(abs(j$observed_asr_per_100000 - j$asr_obs) < 1e-9))
ai <- xl("output/figure_1_data.xlsx", "AIHW projected estimates")
stopifnot(nrow(ai) == 4 * 10, setequal(ai$lymphoma, c("HL", "NHL")),
          all(ai$year %in% 2026:2035))
pp <- rd("output/prevalence_projections.csv") |>
  filter(subtype %in% names(lab), sex %in% names(sx), year >= 2012) |>
  transmute(lymphoma = lab[subtype], sex = sx[sex], year, duration,
            prev_mid, prev_p025, prev_p975)
for (d in list(c(5, 2), c(40, 3))) {
  f <- xl(sprintf("output/figure_%d_data.xlsx", d[2]), "Modelled prevalence")
  j <- inner_join(f, filter(pp, duration == d[1]), by = c("lymphoma", "sex", "year"))
  stopifnot(nrow(f) == 10 * 34, nrow(j) == nrow(f),
            all(abs(j$prevalence_thousands * 1000 - j$prev_mid) < 1e-6),
            all(abs(j$cri_lower_thousands * 1000 - j$prev_p025) < 1e-6),
            all(abs(j$cri_upper_thousands * 1000 - j$prev_p975) < 1e-6))
}
a2 <- xl("output/figure_2_data.xlsx", "AIHW published estimates")
stopifnot(nrow(a2) == 10 * 10, all(a2$year %in% 2012:2021))
stopifnot(!"AIHW published estimates" %in% openxlsx::getSheetNames("output/figure_3_data.xlsx"))
cat("  Test D3 (figure workbooks equal the source outputs, layer by layer)  PASS\n")

# Test D4: dashboard copies are byte-identical to their sources
m <- rd("dashboard/supplementary/manifest.csv")
stopifnot(nrow(m) == 30, sum(m$kind == "figure") == 14, sum(m$kind == "table") == 16,
          all(file.exists(m$source)),
          all(file.exists(file.path("dashboard", m$published))),
          identical(unname(tools::md5sum(m$source)), m$md5),
          identical(unname(tools::md5sum(file.path("dashboard", m$published))), m$md5),
          setequal(unique(m$label[m$kind == "figure"]), sprintf("A%d", 1:14)),
          setequal(unique(m$label[m$kind == "table"]), sprintf("A%d", 1:9)))
dd <- rd("dashboard/supplementary/download-dictionary.csv")
for (f in unique(dd$file)) {
  hdr <- names(rd(file.path("dashboard/supplementary/tables", f)))
  stopifnot(identical(hdr, dd$column[dd$file == f]))
}
stopifnot(setequal(unique(dd$file), basename(m$published[m$kind == "table"])))
page <- readLines("dashboard/additional-results.qmd")
for (f in basename(m$published)) stopifnot(any(grepl(f, page, fixed = TRUE)))
cat("  Test D4 (30 dashboard files match sources; dictionary matches headers; page links every file)  PASS\n")

cat("\nAll presentation tests passed. Nothing was fitted or regenerated.\n")
