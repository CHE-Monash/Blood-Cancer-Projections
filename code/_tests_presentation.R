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
#   Rscript code/supporting/prevalence_convention_comparison.R
#   Rscript code/supporting/population_projection_comparison.R
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
cat("  Test D3 (figure workbooks: projected and observed layers equal the source outputs)  PASS\n")

# Test D3b: the remaining layers, value by value, from sources that do not
# pass through the figure builders.
#   AIHW comparator points come straight from data/.
#   The fitted history of Figure 1 exists in no CSV, so it is recomputed from
#   the stored fits with calc_fitted_asr(). Sourcing apc_model.R defines
#   functions only (its driver is behind a sys.nframe() guard); nothing is
#   fitted. Graphics are sent to a null device so no Rplots.pdf is written.
grDevices::pdf(NULL)
aihw_inc <- rd("data/incidence_agg_proj.csv") |>
  filter(cancer_group %in% c("Hodgkin lymphoma", "Non-Hodgkin lymphoma")) |>
  transmute(lymphoma = if_else(cancer_group == "Hodgkin lymphoma", "HL", "NHL"),
            sex = sx[tolower(sex)], year, asr_2001)
j <- inner_join(ai, aihw_inc, by = c("lymphoma", "sex", "year"))
stopifnot(nrow(j) == nrow(ai), nrow(aihw_inc) == nrow(ai),
          all(abs(j$aihw_projected_asr_per_100000 - j$asr_2001) < 1e-9))
agg <- rd("data/prevalence_agg.csv") |>
  filter(duration == 5, year %in% 2012:2021, sex %in% names(sx),
         cancer_group %in% c("Hodgkin lymphoma", "Non-Hodgkin lymphoma")) |>
  mutate(lymphoma = if_else(cancer_group == "Hodgkin lymphoma", "HL", "NHL"))
sub <- rd("data/prevalence_subtype.csv") |>
  filter(duration == 5, year %in% 2012:2021, sex %in% names(sx),
         subtype %in% c("dlbcl", "follicular", "mantle_cell")) |>
  mutate(lymphoma = lab[subtype])
aihw_prev <- bind_rows(agg, sub) |>
  group_by(lymphoma, sex = sx[sex], year) |>
  summarise(prev = sum(prevalence, na.rm = TRUE) / 1000, .groups = "drop")
j <- inner_join(a2, aihw_prev, by = c("lymphoma", "sex", "year"))
stopifnot(nrow(j) == nrow(a2),
          all(abs(j$aihw_prevalence_thousands - j$prev) < 1e-9))
suppressPackageStartupMessages(source("code/apc_model.R"))
res     <- readRDS("output/apc_results.rds")
std_pop <- build_std_pop(rd("data/pop_hist.csv"))
fitted  <- purrr::map_dfr(res, function(x)
  calc_fitted_asr(x$fit, x$sex, x$subtype,
                  if (x$tier == "agg") agg_start else subtype_start, std_pop)) |>
  dplyr::transmute(lymphoma = lab[subtype], sex = sx[sex], year = P, asr_fitted)
j <- dplyr::inner_join(fit, fitted, by = c("lymphoma", "sex", "year"))
stopifnot(nrow(j) == nrow(fit), nrow(fitted) == nrow(fit),
          all(abs(j$asr_per_100000 - j$asr_fitted) < 1e-9))
bd <- xl("output/figure_1_data.xlsx", "95% credible interval")
j  <- dplyr::inner_join(bd, ip, by = c("lymphoma", "sex", "year"))
stopifnot(nrow(bd) == 240, nrow(j) == 240, all(bd$year > 2021),
          all(abs(j$cri_lower - j$asr_p025) < 1e-9),
          all(abs(j$cri_upper - j$asr_p975) < 1e-9))
for (d in list(c(5, 2), c(40, 3))) {
  b <- xl(sprintf("output/figure_%d_data.xlsx", d[2]), "95% credible interval")
  j <- dplyr::inner_join(b, dplyr::filter(pp, duration == d[1]),
                         by = c("lymphoma", "sex", "year"))
  stopifnot(nrow(b) == 240, nrow(j) == 240, all(b$year > 2021),
            all(abs(j$cri_lower_thousands * 1000 - j$prev_p025) < 1e-6),
            all(abs(j$cri_upper_thousands * 1000 - j$prev_p975) < 1e-6))
}
invisible(grDevices::dev.off())
cat("  Test D3b (fitted history, AIHW comparators and interval bands equal their sources, value by value)  PASS\n")

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

# Named failure for the D5/D6 source checks, so a planted defect reports
# which responsibility caught it (see code/_tests_presentation_fixtures.R).
chk <- function(cond, code, detail = "") {
  if (!isTRUE(all(cond))) stop(sprintf("%s %s", code, detail), call. = FALSE)
  invisible(TRUE)
}
# Match two tables on a key set: every expected key present exactly once.
keyed <- function(got, expected, keys, code) {
  chk(!anyDuplicated(got[, keys, drop = FALSE]), paste0(code, "_KEYS"), "duplicate key rows")
  chk(nrow(got) == nrow(expected), paste0(code, "_KEYS"),
      sprintf("%d rows stored, %d expected", nrow(got), nrow(expected)))
  j <- inner_join(expected, got, by = keys, suffix = c(".exp", ".got"))
  chk(nrow(j) == nrow(expected), paste0(code, "_KEYS"), "stored keys do not match the expected set")
  j
}

# Test D5: prevalence-convention comparison (Discussion, comparison with
# Luo et al. 2026). Every stored value is recomputed here from the same
# unchanged sources with code written separately from the producer, and the
# producer's cohort loop is exercised on hand-derived controls. Needs
# output/apc_results.rds for the 2022+ diagnoses; without it the source
# recomputation of the prevalence columns is reported as not run.
cv <- rd("output/table_prevalence_convention.csv")
cr <- rd("output/table_constant_rate_diagnoses.csv")
pp <- rd("output/table_population_projection_comparison.csv")
chk(setequal(names(cv), c("subtype", "sex", "duration", "year", "end_of_year", "during_year", "difference_pct")),
    "D5_SCHEMA", "convention table columns")
chk(all(cv$end_of_year > 0) && all(cv$during_year > 0), "D5_SIGN", "non-positive prevalence")
chk(all(abs(cv$difference_pct - (cv$during_year / cv$end_of_year - 1) * 100) <= 0.06),
    "D5_DIFFERENCE_ARITHMETIC")

# D5a: the end-of-year column reproduces the published 5-year NHL prevalence
t3  <- rd("output/table_3_prevalence_combined.csv") |> filter(subtype == "nhl", duration == 5)
dd  <- rd("dashboard/data/lymphoma-dashboard-data.csv") |>
  filter(diagnosis_code == "nhl", sex == "Persons", duration_years == 5,
         measure == "Prevalence (people)")
pub <- tibble(duration = 5, year = c(2021, 2025, 2045),
              published = c(t3$prev_2021, dd$estimate[dd$year == 2025], t3$prev_2045))
j <- inner_join(cv, pub, by = c("duration", "year"))
chk(nrow(j) == 3, "D5_KEYS", "5-year rows for 2021, 2025, 2045 expected")
chk(abs(j$end_of_year / j$published - 1) < 0.002, "D5_END_PUBLISHED_MISMATCH",
    "end-of-year 5-year NHL prevalence differs from the published value by 0.2% or more")
chk(j$during_year < j$end_of_year, "D5_DIRECTION", "the during-year convention must count fewer people")

# D5b: both convention columns from source, independently assembled
grDevices::pdf(NULL)
suppressPackageStartupMessages(suppressWarnings(
  source("code/supporting/prevalence_convention_comparison.R")))
# the production functions .prevalence() and .nhl_cases() are now defined;
# their driver sits behind a sys.nframe() guard and did not run
prev_indep <- function(sx, Y, K, conv, cases, surv_obs, improvement) {
  curve <- function(dx) {
    S <- build_projected_surv_pt("nhl", sx, dx, improve_cap, surv_obs, improvement)
    function(t) { if (t == 0) return(setNames(rep(1, length(age_groups)), age_groups))
      v <- S$surv_mid[S$years_since_dx == t]; names(v) <- S$age_group[S$years_since_dx == t]; v[age_groups] }
  }
  N <- function(y) { n <- cases$count[cases$year == y]; names(n) <- cases$age_group[cases$year == y]
    chk(setequal(names(n), age_groups), "D5_SOURCE_CASES", sprintf("diagnoses for %d incomplete", y)); n[age_groups] }
  total <- 0
  if (conv == "end_of_year") {
    for (t in 0:(K - 1)) { S <- curve(Y - t)
      w <- if (t == 0) sqrt(S(1)) else S(t)
      total <- total + sum(N(Y - t) * w) }
  } else {
    for (k in 1:K) { S <- curve(Y - k)
      total <- total + sum(N(Y - k) * (S(k - 1) + S(k)) / 2) }
  }
  total
}
if (file.exists("output/apc_results.rds")) {
  inputs      <- load_prev_inputs("data")
  surv_obs    <- build_surv_obs(inputs$surv_raw)
  improvement <- build_improvement_final(surv_obs)
  apc_results <- readRDS("output/apc_results.rds")
  cases <- lapply(setNames(sexes, sexes), .nhl_cases, inputs = inputs, apc_results = apc_results)
  exp_cv <- cv |> distinct(duration, year) |> rowwise() |>
    mutate(end_exp = sum(sapply(sexes, function(sx) prev_indep(sx, year, duration, "end_of_year", cases[[sx]], surv_obs, improvement))),
           dur_exp = sum(sapply(sexes, function(sx) prev_indep(sx, year, duration, "during_year", cases[[sx]], surv_obs, improvement)))) |>
    ungroup()
  j <- keyed(cv, exp_cv, c("duration", "year"), "D5")
  chk(abs(j$end_of_year - j$end_exp) <= 0.51, "D5_END_SOURCE_MISMATCH",
      "end-of-year column differs from the source recomputation")
  chk(abs(j$during_year - j$dur_exp) <= 0.51, "D5_DURING_SOURCE_MISMATCH",
      "during-year column differs from the source recomputation")
  d5b <- "both convention columns recomputed from source"
} else {
  d5b <- "NOT RUN: output/apc_results.rds absent, convention columns not recomputed from source"
}

# D5c: hand-derived controls on the production cohort loop. One age group
# holds diagnoses 100 (year Y), 200 (Y-1), 300 (Y-2); the other groups hold
# none; survival is S(k) = 0.64^k for every group. Expected: end-of-year
# K=1 100*0.8 = 80, K=2 80 + 200*0.64 = 208; during-year K=1 200*(1+0.64)/2
# = 164, K=2 164 + 300*(0.64+0.4096)/2 = 321.44. The survival builder is
# replaced by a stub for the calls and restored afterwards.
.orig_surv <- build_projected_surv_pt
build_projected_surv_pt <- function(subtype_name, sex_name, dx_year, cap_year, surv_obs, improvement_final) {
  g <- expand.grid(age_group = age_groups, years_since_dx = 1:40, stringsAsFactors = FALSE)
  g$surv_mid <- 0.64^g$years_since_dx; g
}
ctrl_cases <- expand.grid(year = 2040:2045, age_group = age_groups, stringsAsFactors = FALSE)
ctrl_cases$count <- ifelse(ctrl_cases$age_group == age_groups[5],
                           c(`2045` = 100, `2044` = 200, `2043` = 300)[as.character(ctrl_cases$year)], 0)
ctrl_cases$count[is.na(ctrl_cases$count)] <- 0
ctrl <- c(.prevalence("males", 2045, 1, "end_of_year", ctrl_cases, NULL, NULL),
          .prevalence("males", 2045, 2, "end_of_year", ctrl_cases, NULL, NULL),
          .prevalence("males", 2045, 1, "during_year", ctrl_cases, NULL, NULL),
          .prevalence("males", 2045, 2, "during_year", ctrl_cases, NULL, NULL))
ctrl_i <- c(prev_indep("males", 2045, 1, "end_of_year", ctrl_cases, NULL, NULL),
            prev_indep("males", 2045, 2, "end_of_year", ctrl_cases, NULL, NULL),
            prev_indep("males", 2045, 1, "during_year", ctrl_cases, NULL, NULL),
            prev_indep("males", 2045, 2, "during_year", ctrl_cases, NULL, NULL))
build_projected_surv_pt <- .orig_surv
chk(abs(ctrl - c(80, 208, 164, 321.44)) < 1e-9, "D5_CONTROL_PRODUCER",
    sprintf("production loop gave %s", paste(round(ctrl, 4), collapse = ", ")))
chk(abs(ctrl_i - c(80, 208, 164, 321.44)) < 1e-9, "D5_CONTROL_COMPARATOR",
    sprintf("independent comparator gave %s", paste(round(ctrl_i, 4), collapse = ", ")))

# D5d / D6: constant-rate diagnoses and populations from source, for both
# ABS releases, matched by key
inc_agg  <- rd("data/incidence_agg.csv")
pop_hist <- rd("data/pop_hist.csv")
pop_2022 <- rd("data/pop_proj.csv")
pop_2017 <- rd("data/pop_proj_2017base.csv")
rate_exp <- inc_agg |>
  filter(cancer_group == "Non-Hodgkin lymphoma", year %in% 2017:2021) |>
  group_by(sex, age_group) |> summarise(d = sum(count), .groups = "drop") |>
  inner_join(pop_hist |> filter(year %in% 2017:2021) |>
               group_by(sex, age_group) |> summarise(y = sum(population), .groups = "drop"),
             by = c("sex", "age_group")) |>
  mutate(rate = d / y)
chk(nrow(rate_exp) == 2 * length(age_groups), "D5_SOURCE_RATES", "rate table incomplete")
pop_summary <- function(pop, label) pop |>
  filter(sex %in% c("males", "females")) |>
  inner_join(rate_exp[, c("sex", "age_group", "rate")], by = c("sex", "age_group")) |>
  group_by(year) |>
  summarise(total_exp = sum(population),
            old_exp   = sum(population[age_group %in% c("65–74", "75–84", "85+")]),
            diag_exp  = sum(rate * population), .groups = "drop") |>
  mutate(projection = label)
exp22 <- pop_summary(pop_2022, "ABS 2022-base (this study)")
exp17 <- pop_summary(pop_2017, "ABS 2017-base (Luo 2026)")
# D5d: constant-rate table
j <- keyed(cr, exp22 |> filter(year %in% c(2024, 2025, 2045, 2049, 2050)), "year", "D5")
chk(abs(j$diagnoses_constant_rate - j$diag_exp) <= 0.51, "D5_CONST_RATE_SOURCE_MISMATCH",
    "constant-rate diagnoses differ from rates x ABS 2022-base population")
# D6: population comparison table
exp_pp <- bind_rows(exp22, exp17) |> filter(year %in% c(2024, 2025, 2045, 2049, 2050))
j <- keyed(pp, exp_pp, c("projection", "year"), "D6")
chk(abs(j$population_total / j$total_exp - 1) <= 1e-6, "D6_POPULATION_SOURCE_MISMATCH",
    "population totals differ from the population files")
chk(abs(j$population_65plus / j$old_exp - 1) <= 1e-6, "D6_POPULATION_SOURCE_MISMATCH",
    "65+ counts differ from the population files")
chk(abs(j$share_65plus_pct - round(100 * j$old_exp / j$total_exp, 1)) <= 0.051, "D6_SHARE_ARITHMETIC")
chk(abs(j$diagnoses_constant_rate - j$diag_exp) <= 0.51, "D6_DIAGNOSES_SOURCE_MISMATCH",
    "constant-rate diagnoses differ from rates x population for at least one release")
ours <- j |> filter(projection == "ABS 2022-base (this study)") |> arrange(year)
luo  <- j |> filter(projection == "ABS 2017-base (Luo 2026)") |> arrange(year)
chk(identical(ours$diagnoses_constant_rate, cr$diagnoses_constant_rate[match(ours$year, cr$year)]),
    "D6_CONSISTENCY", "2022-base rows differ from table_constant_rate_diagnoses.csv")
chk(luo$population_total > ours$population_total & luo$population_65plus < ours$population_65plus &
      luo$diagnoses_constant_rate < ours$diagnoses_constant_rate, "D6_DIRECTION",
    "the 2017-base release should be larger, older-poorer and give fewer diagnoses")
# D6b: derived 2017-base population against the ABS workbook, value by value
tot17 <- pop_2017 |> group_by(year) |> summarise(population = sum(population), .groups = "drop")
chk(nrow(pop_2017) == 50 * 2 * length(age_groups) && setequal(pop_2017$age_group, age_groups),
    "D6_DERIVED_SHAPE")
chk(tot17$population[tot17$year == 2030] == 29931725 && tot17$population[tot17$year == 2050] == 37077210,
    "D6_ABS_TOTALS", "ABS 3222.0 2017-base Series B totals for 2030 and 2050")
raw17 <- "data/raw/ABS/32220_Table_B9_2017base.xls"
if (file.exists(raw17)) {
  all <- suppressMessages(readxl::read_excel(raw17, sheet = "Data1", col_names = FALSE))
  desc <- as.character(unlist(all[1, ], use.names = FALSE))
  parts <- strsplit(desc, ";", fixed = TRUE)
  sexv <- trimws(sapply(parts, function(x) if (length(x) >= 2) x[2] else NA))
  agev <- trimws(sapply(parts, function(x) if (length(x) >= 3) x[3] else NA))
  keep <- which(sexv %in% c("Male", "Female"))
  yrs  <- as.integer(format(as.Date(as.numeric(all[[1]][-(1:10)]), origin = "1899-12-30"), "%Y"))
  long <- bind_rows(lapply(keep, function(i) tibble(
    year = yrs, sex = if (sexv[i] == "Male") "males" else "females",
    age = if (agev[i] == "100 and over") 100L else as.integer(agev[i]),
    population = as.numeric(all[[i]][-(1:10)]))))
  chk(length(keep) == 202 && !anyNA(long$age), "D6_RAW_PARSE")
  brk <- c(0, 15, 25, 35, 45, 55, 65, 75, 85, Inf)
  long$age_group <- age_groups[findInterval(long$age, brk)]
  raw_agg <- long |> group_by(year, sex, age_group) |> summarise(population = sum(population), .groups = "drop")
  j <- keyed(pop_2017, raw_agg, c("year", "sex", "age_group"), "D6_RAW")
  chk(j$population.got == j$population.exp, "D6_RAW_DERIVED_MISMATCH",
      "derived 2017-base population differs from the ABS workbook")
  d6b <- sprintf("all %d derived values equal the ABS workbook", nrow(j))
} else {
  d6b <- "NOT RUN: raw ABS 2017-base workbook absent, derived CSV checked against ABS totals only"
}
cat(sprintf("  Test D5 (end-of-year column within 0.2%% of published 5-year NHL prevalence; %s; four hand-derived controls pass in producer and comparator; constant-rate diagnoses from source)  PASS\n", d5b))
cat(sprintf("  Test D6 (population comparison from source for both releases; %s; constant-rate diagnoses rise %.1f%% from 2025 to 2050 on the 2017-base against %.1f%% on the 2022-base)  PASS\n",
            d6b,
            (luo$diagnoses_constant_rate[luo$year == 2050] / luo$diagnoses_constant_rate[luo$year == 2025] - 1) * 100,
            (ours$diagnoses_constant_rate[ours$year == 2050] / ours$diagnoses_constant_rate[ours$year == 2025] - 1) * 100))

cat("\nAll presentation tests passed. Nothing was fitted or regenerated.\n")
