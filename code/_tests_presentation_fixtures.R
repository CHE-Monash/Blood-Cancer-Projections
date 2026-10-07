# =================================================================
# Planted-defect fixtures for code/_tests_presentation.R
#
# Builds scratch copies of the repository's inputs and outputs under
# tempdir(), plants one defect per fixture, runs the REAL test entrypoint
# there, and reports whether the planted defect was detected and by which
# named check. Nothing in the live repository is touched: data/ and output/
# are copied (the large saved fit is linked read-only), code/ and dashboard/
# are linked, and the two fixtures that regenerate outputs write only inside
# the scratch copy.
#
#   Rscript code/_tests_presentation_fixtures.R
#
# Exit status is nonzero if any fixture's outcome differs from expectation.
# =================================================================

suppressPackageStartupMessages({ library(readr); library(dplyr) })
rd <- function(f) read_csv(f, show_col_types = FALSE, na = c("", "NA"))

root <- normalizePath(".")
stopifnot(file.exists(file.path(root, "code/_tests_presentation.R")))

make_scratch <- function() {
  d <- tempfile("fixture_"); dir.create(d)
  file.symlink(file.path(root, "code"), file.path(d, "code"))
  file.symlink(file.path(root, "dashboard"), file.path(d, "dashboard"))
  dir.create(file.path(d, "data")); dir.create(file.path(d, "output"))
  for (f in list.files(file.path(root, "data"), full.names = TRUE, recursive = FALSE))
    if (!dir.exists(f)) file.copy(f, file.path(d, "data", basename(f)))
  if (dir.exists(file.path(root, "data/raw"))) file.symlink(file.path(root, "data/raw"), file.path(d, "data/raw"))
  for (f in list.files(file.path(root, "output"), full.names = TRUE)) {
    if (grepl("\\.rds$", f)) file.symlink(f, file.path(d, "output", basename(f)))
    else file.copy(f, file.path(d, "output", basename(f)))
  }
  d
}
run_tests <- function(d) {
  out <- suppressWarnings(system2("sh", c("-c", shQuote(sprintf("cd %s && Rscript code/_tests_presentation.R 2>&1", d))),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status"); if (is.null(status)) status <- 0L
  list(status = status, log = out)
}
regen <- function(d) {
  for (s in c("code/supporting/prevalence_convention_comparison.R",
              "code/supporting/population_projection_comparison.R"))
    stopifnot(system2("sh", c("-c", shQuote(sprintf("cd %s && Rscript %s >/dev/null 2>&1", d, s)))) == 0)
}
edit_csv <- function(path, fn) write_csv(fn(rd(path)), path)

fixtures <- list(
  list(name = "clean copy", expect = "pass", mutate = function(d) NULL),
  list(name = "during-year column set to 80% of end-of-year, difference recomputed", expect = "fail", code = "D5_DURING_SOURCE_MISMATCH",
       mutate = function(d) edit_csv(file.path(d, "output/table_prevalence_convention.csv"), function(x)
         mutate(x, during_year = round(0.8 * end_of_year), difference_pct = -20))),
  list(name = "2017-base 2050 constant-rate diagnoses lowered 10%", expect = "fail", code = "D6_DIAGNOSES_SOURCE_MISMATCH",
       mutate = function(d) edit_csv(file.path(d, "output/table_population_projection_comparison.csv"), function(x)
         mutate(x, diagnoses_constant_rate = ifelse(grepl("2017", projection) & year == 2050,
                                                    round(0.9 * diagnoses_constant_rate), diagnoses_constant_rate)))),
  list(name = "end-of-year column raised 10%, difference recomputed", expect = "fail", code = "D5_END_PUBLISHED_MISMATCH",
       mutate = function(d) edit_csv(file.path(d, "output/table_prevalence_convention.csv"), function(x)
         mutate(x, end_of_year = round(1.1 * end_of_year),
                difference_pct = round((during_year / end_of_year - 1) * 100, 1)))),
  list(name = "all 2017-base constant-rate diagnoses lowered 20%", expect = "fail", code = "D6_DIAGNOSES_SOURCE_MISMATCH",
       mutate = function(d) edit_csv(file.path(d, "output/table_population_projection_comparison.csv"), function(x)
         mutate(x, diagnoses_constant_rate = ifelse(grepl("2017", projection),
                                                    round(0.8 * diagnoses_constant_rate), diagnoses_constant_rate)))),
  list(name = "convention table year 2045 relabelled 2046 (shifted year)", expect = "fail", code = "D5_KEYS",
       mutate = function(d) edit_csv(file.path(d, "output/table_prevalence_convention.csv"), function(x)
         mutate(x, year = ifelse(year == 2045, 2046, year)))),
  list(name = "duplicate row appended to the population comparison table", expect = "fail", code = "D6_KEYS",
       mutate = function(d) edit_csv(file.path(d, "output/table_population_projection_comparison.csv"), function(x)
         bind_rows(x, x[1, ]))),
  list(name = "one row dropped from the constant-rate table", expect = "fail", code = "D5_KEYS",
       mutate = function(d) edit_csv(file.path(d, "output/table_constant_rate_diagnoses.csv"), function(x)
         filter(x, year != 2049))),
  list(name = "survival input altered (NHL males 65-74, 5 years, 2017-2021 period, -20 points), outputs stale", expect = "fail", code = "D5_END_SOURCE_MISMATCH",
       mutate = function(d) edit_csv(file.path(d, "data/survival.csv"), function(x) {
         i <- which(x$subtype == "nhl" & x$sex == "males" & x$age_group == "65–74" &
                    x$years_since_dx == 5 & x$period == "2017–2021" & !is.na(x$survival_pct))
         stopifnot(length(i) >= 1); x$survival_pct[i] <- x$survival_pct[i] - 20; x })),
  list(name = "historical NHL incidence altered (+200 diagnoses, males 65-74, 2019), outputs stale", expect = "fail", code = "D5_END_SOURCE_MISMATCH",
       mutate = function(d) edit_csv(file.path(d, "data/incidence_agg.csv"), function(x) {
         i <- which(x$cancer_group == "Non-Hodgkin lymphoma" & x$sex == "males" & x$age_group == "65–74" & x$year == 2019)
         stopifnot(length(i) == 1); x$count[i] <- x$count[i] + 200; x })),
  list(name = "same incidence alteration, outputs regenerated: still fails, because the end-of-year column must agree with the PUBLISHED prevalence, which the altered source no longer reproduces", expect = "fail", code = "D5_END_PUBLISHED_MISMATCH",
       mutate = function(d) {
         edit_csv(file.path(d, "data/incidence_agg.csv"), function(x) {
           i <- which(x$cancer_group == "Non-Hodgkin lymphoma" & x$sex == "males" & x$age_group == "65–74" & x$year == 2019)
           x$count[i] <- x$count[i] + 200; x })
         regen(d) }),
  list(name = "projected population altered (2022-base, males 75-84, 2045, -100,000), outputs stale", expect = "fail", code = "D5_CONST_RATE_SOURCE_MISMATCH",
       mutate = function(d) edit_csv(file.path(d, "data/pop_proj.csv"), function(x) {
         i <- which(x$sex == "males" & x$age_group == "75–84" & x$year == 2045)
         stopifnot(length(i) == 1); x$population[i] <- x$population[i] - 100000; x })),
  list(name = "projected population altered as above, outputs regenerated in the scratch copy", expect = "pass",
       mutate = function(d) {
         edit_csv(file.path(d, "data/pop_proj.csv"), function(x) {
           i <- which(x$sex == "males" & x$age_group == "75–84" & x$year == 2045)
           x$population[i] <- x$population[i] - 100000; x })
         regen(d) }),
  list(name = "derived 2017-base population altered (one cell -50,000), outputs stale", expect = "fail", code = "D6_POPULATION_SOURCE_MISMATCH",
       mutate = function(d) edit_csv(file.path(d, "data/pop_proj_2017base.csv"), function(x) {
         i <- which(x$sex == "females" & x$age_group == "85+" & x$year == 2050)
         stopifnot(length(i) == 1); x$population[i] <- x$population[i] - 50000; x })),
  list(name = "projected population perturbed by one person (2022-base, females 45-54, 2049), outputs stale", expect = "pass",
       mutate = function(d) edit_csv(file.path(d, "data/pop_proj.csv"), function(x) {
         i <- which(x$sex == "females" & x$age_group == "45–54" & x$year == 2049)
         x$population[i] <- x$population[i] + 1; x })),
  # Probe of the harness itself: an unrelated failure (D1's input removed)
  # must NOT be counted as a detected D5/D6 defect. The probe is correct
  # when the harness refuses to credit it.
  list(name = "PROBE: unrelated failure (si_table_holdout.csv removed) must not count as a D5/D6 detection",
       expect = "fail", code = "D5_NEVER_EMITTED", probe = TRUE,
       mutate = function(d) file.remove(file.path(d, "output/si_table_holdout.csv")))
)

results <- lapply(fixtures, function(f) {
  d <- make_scratch(); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  f$mutate(d)
  r <- run_tests(d)
  outcome <- if (r$status == 0) "pass" else "fail"
  caught <- if (r$status == 0) "" else {
    m <- regmatches(r$log, regexpr("D[56]_[A-Z_]+", r$log)); if (length(m)) m[1] else "unnamed" }
  err <- grep("^Error", r$log, value = TRUE); err <- if (length(err)) sub("^Error[^:]*: ?", "", err[1]) else ""
  # a defect counts as detected only when the named diagnostic the fixture
  # expects is the one raised; any other nonzero exit is a harness failure
  credited <- if (f$expect == "pass") outcome == "pass" else outcome == "fail" && caught %in% f$code
  probe <- isTRUE(f$probe)
  tibble(fixture = f$name, expected = if (probe) "not credited" else f$expect, outcome = outcome,
         as_expected = if (probe) !credited else credited, detected_by = caught, message = substr(err, 1, 90))
})
res <- bind_rows(results)
print(as.data.frame(res[, c("fixture", "expected", "outcome", "as_expected", "detected_by")]), right = FALSE)
cat(sprintf("\n%d of %d fixtures behaved as expected; %d planted defects detected by their expected diagnostic of %d; harness probe %s.\n",
            sum(res$as_expected), nrow(res), sum(res$expected == "fail" & res$as_expected), sum(res$expected == "fail"),
            if (all(res$as_expected[res$expected == "not credited"])) "refused the unrelated failure" else "CREDITED AN UNRELATED FAILURE"))
if (!all(res$as_expected)) quit(status = 1)
