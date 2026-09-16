# =================================================================
# reproduce.R - one-click reproduction of the Blood Cancer Projections
# (lymphoma) pipeline for an external reproducer.
# =================================================================
# HOW TO USE
#   1. Fresh download of the repository. Open it as the working directory
#      (the repository root, the folder that contains code/ and data/).
#      In RStudio: Session > Set Working Directory > To Project Directory,
#      or setwd("<path to the repository root>").
#   2. Session > Restart R, so that nothing from an earlier session leaks in.
#   3. Open this file and click Source (or run: source("code/reproduce.R")).
#
# WHY THIS FILE EXISTS
#   The pipeline scripts (code/apc_model.R, code/prev_model.R,
#   code/supplement.R) follow a library + driver pattern with the guard
#   `if (!interactive() && sys.nframe() == 0) run_*()` at the bottom. In
#   RStudio, clicking Source (or Run All) on one of them only defines its
#   functions: nothing runs. `Rscript code/apc_model.R` from a terminal
#   runs it. This script sources each file and then calls its driver, in
#   order, so that Source on THIS file runs the whole pipeline.
#
#   The output/ folder is committed, so a fresh download already contains
#   the published numbers. Checking the manuscript against a download that
#   has not been re-run proves nothing. This script therefore moves the
#   committed output/ (and copies data/) into _reproduce/committed/ before
#   running, and compares the fresh run against that snapshot afterwards.
#
# PINNED ENVIRONMENT
#   The published results were produced with R 4.5.2 and Epi 2.61 (seed
#   20260507, B = 1000 Monte Carlo draws). Other versions run the pipeline
#   but may refit the APC models slightly differently; see the
#   interpretation printed at the end.
#
# READING THE RESULTS (printed at the end of the run)
#   - every key output identical (md5 after line-ending normalisation)
#     -> exact reproduction.
#   - inputs identical (or not rebuilt) and the model deviances in
#     table_s2_apc_fit_stats.csv identical -> the same ten models were
#     fitted; remaining differences are Monte Carlo draws, which can differ
#     between R/Epi versions and platforms even at the same seed. With
#     B = 1000 the Monte Carlo standard error is roughly +/-4 on a median
#     and +/-8 on a credible-interval bound for the Table 1 totals, and
#     larger for MCL. The point rate ratios in apc_effects.csv (no Monte
#     Carlo) may still differ at the ends of the cohort range if the Epi
#     version allocates the drift between effects differently; that moves
#     the effect plots (Figures S1 to S3), not the projections.
#   - inputs identical but the deviances differ -> the fit itself differs;
#     install the pinned R and Epi versions and re-run.
#   - inputs differ -> look at the data before anything else.
#   The sentinel table prints the headline Table 1 and Table 3 values,
#   new against committed, so no manual lookup is needed. A copy is
#   written to _reproduce/sentinels.csv.
#
# Run time is about 90 seconds at B = 1000 (most of it in prev_model.R).
# =================================================================

rebuild_data <- FALSE   # TRUE rebuilds data/*.csv from the AIHW and ABS
                        # workbooks, which must be in data/raw/ (see README,
                        # "Files required to recreate the analysis from raw
                        # data"). FALSE uses the committed data/*.csv.

pinned_r   <- "4.5.2"
pinned_epi <- "2.61"

if (!file.exists("code/apc_model.R") || !dir.exists("data")) {
  stop("Run this script from the repository root (the folder containing code/ and data/).")
}

# -----------------------------------------------------------------
# 1. Snapshot the committed files
# -----------------------------------------------------------------
snap <- file.path("_reproduce", "committed")
if (!dir.exists(snap)) {
  dir.create(file.path(snap, "data"), recursive = TRUE, showWarnings = FALSE)
  for (f in list.files("data", pattern = "\\.csv$", full.names = TRUE)) {
    file.copy(f, file.path(snap, "data"))       # data/raw/ is not copied
  }
  if (dir.exists("output")) {
    ok <- file.rename("output", file.path(snap, "output"))
    if (!ok) stop("Could not move output/ to ", file.path(snap, "output"))
  } else {
    dir.create(file.path(snap, "output"), showWarnings = FALSE)
  }
  cat("Snapshot of the committed data/ and output/ taken in", snap, "\n")
} else {
  cat("Snapshot already exists in", snap, "and is kept; this run overwrites output/\n")
}
dir.create("output", showWarnings = FALSE)

# -----------------------------------------------------------------
# 2. Rebuild the data, if requested
# -----------------------------------------------------------------
if (isTRUE(rebuild_data)) {
  if (!dir.exists("data/raw")) {
    stop("rebuild_data is TRUE but data/raw/ is missing. Place the AIHW and ABS ",
         "workbooks there as listed in README.md, or set rebuild_data <- FALSE.")
  }
  cat("\nRebuilding data/*.csv from data/raw/ (import_data.R)\n")
  source("code/import_data.R")   # no auto-run guard: sourcing runs it
}

# -----------------------------------------------------------------
# 3. Run the pipeline
# -----------------------------------------------------------------
t0 <- Sys.time()
cat("\n=== apc_model.R (APC fits, Monte Carlo incidence, Figure 1, Tables 1 and S2,",
    "Table 1 change CrIs) ===\n")
source("code/apc_model.R");  invisible(run_apc_model())
cat("\n=== prev_model.R (prevalence, Figures 2 and 3, Table 2 and its change CrIs) ===\n")
source("code/prev_model.R"); invisible(run_prev_model())
cat("\n=== supplement.R (SI tables and figures) ===\n")
source("code/supplement.R"); invisible(run_supplement())
cat("\n=== _tests.R (structural test harness) ===\n")
source("code/_tests.R")
cat(sprintf("\nPipeline and tests complete in %.0f s\n",
            as.numeric(Sys.time() - t0, units = "secs")))

# -----------------------------------------------------------------
# 4. Report
# -----------------------------------------------------------------
# md5 of a file after normalising line endings (git on Windows may check the
# committed files out with CRLF while R writes LF; that is not a difference).
md5_of <- function(f) {
  if (!file.exists(f)) return(NA_character_)
  txt <- readChar(f, file.info(f)$size, useBytes = TRUE)
  txt <- gsub("\r\n", "\n", txt, fixed = TRUE)
  tmp <- tempfile(); on.exit(unlink(tmp)); writeBin(charToRaw(txt), tmp)
  unname(tools::md5sum(tmp))
}
read_plain <- function(f) utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)

# Numerical comparison of two CSVs: same shape (dims and names), then the
# maximum relative difference over numeric columns. Non-numeric columns must
# match exactly. Line endings cannot cause a false mismatch.
compare_csv <- function(new_f, old_f) {
  if (!file.exists(new_f) || !file.exists(old_f)) {
    return(list(same_shape = NA, max_rel_diff = NA_real_, note = "file missing"))
  }
  x <- read_plain(new_f); y <- read_plain(old_f)
  same_shape <- identical(dim(x), dim(y)) && identical(names(x), names(y))
  if (!same_shape) return(list(same_shape = FALSE, max_rel_diff = NA_real_, note = "dims or names differ"))
  num <- names(x)[vapply(x, is.numeric, logical(1)) & vapply(y, is.numeric, logical(1))]
  chr <- setdiff(names(x), num)
  chr_ok <- all(vapply(chr, function(v) identical(as.character(x[[v]]), as.character(y[[v]])), logical(1)))
  rel <- unlist(lapply(num, function(v) {
    d <- abs(x[[v]] - y[[v]]) / pmax(abs(y[[v]]), 1e-12)
    d[is.na(x[[v]]) & is.na(y[[v]])] <- 0
    d
  }))
  list(same_shape = TRUE, max_rel_diff = if (length(rel)) max(rel, na.rm = TRUE) else 0,
       note = if (chr_ok) "" else "non-numeric columns differ")
}

cat("\n=================================================================\n")
cat("REPRODUCTION REPORT\n")
cat("=================================================================\n")

# 4a. Versions
r_ver   <- as.character(getRversion())
epi_ver <- as.character(utils::packageVersion("Epi"))
cat(sprintf("R:   %s  (pinned %s)%s\n", R.version.string, pinned_r,
            if (r_ver == pinned_r) "" else "   <-- differs from the pinned version"))
cat(sprintf("Epi: %s  (pinned %s)%s\n", epi_ver, pinned_epi,
            if (epi_ver == pinned_epi) "" else "   <-- differs from the pinned version"))

# 4b. md5 identity of the key outputs
key_outputs <- c("table_1_incidence_summary.csv", "table_1_change_cri.csv",
                 "incidence_projections.csv", "apc_effects.csv",
                 "table_2_change_cri.csv", "table_3_prevalence_combined.csv",
                 "prevalence_projections.csv")
md5_same <- vapply(key_outputs, function(f) {
  a <- md5_of(file.path("output", f)); b <- md5_of(file.path(snap, "output", f))
  !is.na(a) && !is.na(b) && a == b
}, logical(1))
cat("\nKey outputs identical to the committed snapshot (md5 after line-ending normalisation):\n")
for (f in key_outputs) cat(sprintf("  %-36s %s\n", f, md5_same[[f]]))

# 4c. Inputs (only meaningful when the data were rebuilt)
inputs_identical <- NA
if (isTRUE(rebuild_data)) {
  cat("\nRebuilt data/*.csv against the committed snapshot (numerical comparison):\n")
  data_files <- list.files(file.path(snap, "data"), pattern = "\\.csv$")
  data_cmp <- lapply(data_files, function(f) compare_csv(file.path("data", f), file.path(snap, "data", f)))
  names(data_cmp) <- data_files
  for (f in data_files) {
    z <- data_cmp[[f]]
    cat(sprintf("  %-32s same shape: %-5s  max relative difference: %s  %s\n", f,
                z$same_shape, if (is.na(z$max_rel_diff)) "NA" else format(z$max_rel_diff, digits = 3), z$note))
  }
  inputs_identical <- all(vapply(data_cmp, function(z)
    isTRUE(z$same_shape) && z$note == "" && is.finite(z$max_rel_diff) && z$max_rel_diff < 1e-9, logical(1)))
} else {
  cat("\nInputs: data/*.csv were not rebuilt (rebuild_data = FALSE); the committed CSVs were used as they are.\n")
}

# 4d. Fitted models: deviance of the ten APC fits (no Monte Carlo involved).
#     Identical deviances mean the same models were fitted.
dev_cmp <- compare_csv("output/table_s2_apc_fit_stats.csv", file.path(snap, "output", "table_s2_apc_fit_stats.csv"))
dev_diff <- NA_real_
if (isTRUE(dev_cmp$same_shape)) {
  a <- read_plain("output/table_s2_apc_fit_stats.csv"); b <- read_plain(file.path(snap, "output", "table_s2_apc_fit_stats.csv"))
  dev_diff <- max(abs(a$deviance - b$deviance) / pmax(abs(b$deviance), 1e-12), na.rm = TRUE)
}
cat(sprintf("\nMaximum relative difference in model deviance, table_s2_apc_fit_stats.csv (the fit itself, no Monte Carlo): %s\n",
            if (is.na(dev_diff)) paste0("not comparable (", dev_cmp$note, ")") else format(dev_diff, digits = 3)))

# 4e. Point rate ratios (no Monte Carlo). These can differ at the ends of the
#     cohort range between Epi versions even when the fit is identical.
rr_cmp <- compare_csv("output/apc_effects.csv", file.path(snap, "output", "apc_effects.csv"))
rr_diff <- NA_real_
if (isTRUE(rr_cmp$same_shape)) {
  a <- read_plain("output/apc_effects.csv"); b <- read_plain(file.path(snap, "output", "apc_effects.csv"))
  rel <- abs(a$rr - b$rr) / pmax(abs(b$rr), 1e-12)
  rr_diff <- max(rel, na.rm = TRUE)
  cat(sprintf("Maximum relative difference in apc_effects.csv$rr (point rate ratios): %s; median %s\n",
              format(rr_diff, digits = 3), format(median(rel, na.rm = TRUE), digits = 3)))
  if (is.finite(rr_diff) && rr_diff > 1e-6) {
    by_fit <- tapply(rel, paste(a$subtype, a$sex), max)
    cat("  Largest by fit:", paste(sprintf("%s %.2e", names(by_fit), by_fit)[order(-by_fit)][1:3], collapse = "; "), "\n")
  }
} else {
  cat(sprintf("apc_effects.csv not comparable (%s)\n", rr_cmp$note))
}

# 4f. Sentinel values, new against committed
pick <- function(df, cond, cols) { r <- df[cond, , drop = FALSE]; if (nrow(r) != 1) rep(NA, length(cols)) else unname(unlist(r[1, cols])) }
sentinels <- function(dir) {
  t1 <- read_plain(file.path(dir, "table_1_change_cri.csv"))
  t3 <- read_plain(file.path(dir, "table_3_prevalence_combined.csv"))
  fmt_ci <- function(m, lo, hi, d = 0) sprintf("%s (%s to %s)", format(round(m, d), big.mark = ","),
                                               format(round(lo, d), big.mark = ","), format(round(hi, d), big.mark = ","))
  t1_row <- function(sub, sx) {
    v <- pick(t1, t1$subtype == sub & t1$sex == sx,
              c("cases_2045", "cases_2045_p025", "cases_2045_p975", "asr_2045", "asr_2045_p025", "asr_2045_p975",
                "pct_change_cases", "pct_change_cases_p025", "pct_change_cases_p975",
                "pct_change_asr", "pct_change_asr_p025", "pct_change_asr_p975"))
    c(cases_2045 = fmt_ci(v[1], v[2], v[3]), asr_2045 = fmt_ci(v[4], v[5], v[6], 1),
      change_cases_pct = fmt_ci(v[7], v[8], v[9], 1), change_asr_pct = fmt_ci(v[10], v[11], v[12], 1))
  }
  t3_row <- function(sub) {
    v <- pick(t3, t3$subtype == sub & t3$duration == 40,
              c("prev_2021", "prev_2045", "prev_2045_p025", "prev_2045_p975", "pct_change", "pct_change_p025", "pct_change_p975"))
    c(prev_2021 = format(round(v[1]), big.mark = ","), prev_2045 = fmt_ci(v[2], v[3], v[4]),
      change_pct = fmt_ci(v[5], v[6], v[7], 1))
  }
  a <- t1_row("total", "males"); b <- t1_row("mantle_cell", "females")
  d <- t3_row("total"); e <- t3_row("mantle_cell")
  c(setNames(a, paste("Table 1, Total, Males:", names(a))),
    setNames(b, paste("Table 1, MCL, Females:", names(b))),
    setNames(d, paste("Table 3, Total, 40-year:", names(d))),
    setNames(e[1:2], paste("Table 3, MCL, 40-year:", names(e)[1:2])))
}
sent <- data.frame(sentinel = names(sentinels("output")),
                   new = unname(sentinels("output")),
                   committed = unname(sentinels(file.path(snap, "output"))),
                   stringsAsFactors = FALSE)
sent$same <- sent$new == sent$committed
cat("\nSentinel values (rounded as displayed in the manuscript), new against committed:\n")
print(sent, row.names = FALSE, right = FALSE)
utils::write.csv(sent, file.path("_reproduce", "sentinels.csv"), row.names = FALSE)

# 4g. Interpretation
cat("\nINTERPRETATION\n")
inputs_ok <- is.na(inputs_identical) || isTRUE(inputs_identical)
if (all(md5_same)) {
  cat("  Every key output is identical to the committed snapshot: exact reproduction.\n")
} else if (!inputs_ok) {
  cat("  The rebuilt inputs differ from the committed data: investigate data/*.csv against the\n",
      " snapshot, and the workbook versions in data/raw/, before anything else.\n")
} else if (is.finite(dev_diff) && dev_diff < 1e-8) {
  cat("  Inputs", if (is.na(inputs_identical)) "not rebuilt (committed CSVs used)" else "identical",
      "and every model deviance identical: the same ten models were fitted.\n",
      " The remaining differences are Monte Carlo draws, which differ between R/Epi versions and\n",
      " platforms even at the same seed. Expect up to about +/-4 on medians and +/-8 on\n",
      " credible-interval bounds for the Table 1 totals, and more for MCL; compare the sentinel\n",
      " table above against those tolerances.\n")
  if (is.finite(rr_diff) && rr_diff > 1e-6) {
    cat("  The point rate ratios differ by up to", format(rr_diff, digits = 2),
        "at the ends of the cohort range in some fits: a difference in how this\n",
        " Epi version allocates the drift between the period and cohort effects, which moves the\n",
        " effect plots (Figures S1 to S3) but not the projections. Installing the pinned Epi 2.61\n",
        " should remove it.\n")
  }
} else if (is.finite(dev_diff)) {
  cat("  Inputs identical but the model deviances differ by up to", format(dev_diff, digits = 2),
      ": the fit itself differs.\n",
      " Install the pinned R 4.5.2 and Epi 2.61 and re-run before drawing any conclusion.\n")
} else {
  cat("  The fit statistics could not be compared (", dev_cmp$note, "); check that the run completed.\n")
}
cat("\nCommitted snapshot:", snap, "\n")
