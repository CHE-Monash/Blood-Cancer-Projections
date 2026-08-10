# =================================================================
# Structural test script for the Blood Cancer Projections (lymphoma) pipeline (MC-only).
#
# Run after any change to apc_model.R, prev_model.R, or supplement.R.
# Uses the cached outputs from a prior pipeline run if newer than the
# source files; otherwise re-runs the relevant pipeline phase. First
# failing stopifnot() halts execution with a readable message.
#
# Tests are structural and methodological only - file existence,
# reconstructed-component bitwise equality, distribution domain checks,
# CI ordering, and reproducibility under fixed seed. Cross-method
# (MC vs analytical) checks have been retired with the analytical
# pipeline.
# =================================================================

set.seed(20260507)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(MASS, exclude = "select")
})

source("code/apc_model.R")
source("code/prev_model.R")
source("code/supplement.R")

# Helper: run apc_model.R if its output is older than the source
.maybe_run_apc <- function() {
  src    <- file.info("code/apc_model.R")$mtime
  out_rds <- "output/apc_results.rds"
  if (!file.exists(out_rds) || file.info(out_rds)$mtime < src) {
    cat("  Running apc_model.R (output stale)\n")
    invisible(run_apc_model(verbose = FALSE))
  }
}
.maybe_run_prev <- function() {
  src     <- max(file.info("code/prev_model.R")$mtime,
                 file.info("code/apc_model.R")$mtime)
  out_csv <- "output/prevalence_projections.csv"
  if (!file.exists(out_csv) || file.info(out_csv)$mtime < src) {
    cat("  Running prev_model.R (output stale)\n")
    invisible(run_prev_model(verbose = FALSE))
  }
}
.maybe_run_supp <- function() {
  src     <- max(file.info("code/supplement.R")$mtime,
                 file.info("code/prev_model.R")$mtime)
  out_csv <- "output/table_s2a_incidence_aihw.csv"
  if (!file.exists(out_csv) || file.info(out_csv)$mtime < src) {
    cat("  Running supplement.R (output stale)\n")
    invisible(run_supplement(verbose = FALSE))
  }
}

# -----------------------------------------------------------------
# apc_model.R tests
# -----------------------------------------------------------------
cat("\n========== apc_model.R ==========\n")
.maybe_run_apc()

apc_results <- readRDS("output/apc_results.rds")
test_keys   <- names(apc_results)

# Test A1: parameter reconstruction (apc_components vs fit$Age/Per/Coh,
# bitwise across all 10 fits)
worst_recon <- 0
inputs_for_basis <- list(
  inc_subtype = read_csv("data/incidence_subtype.csv", show_col_types = FALSE),
  inc_agg     = read_csv("data/incidence_agg.csv",     show_col_types = FALSE),
  pop_hist    = read_csv("data/pop_hist.csv",          show_col_types = FALSE)
)
for (key in test_keys) {
  fit  <- apc_results[[key]]$fit
  tier <- apc_results[[key]]$tier
  df   <- prepare_fit_data(key, tier,
                           inputs_for_basis$inc_subtype,
                           inputs_for_basis$inc_agg,
                           inputs_for_basis$pop_hist)
  basis <- apc_basis(fit, df)
  comp  <- apc_components(coef(fit$Model), basis)
  err <- max(
    max(abs(comp$Age - as.numeric(fit$Age[, "Rate"])) /
        abs(as.numeric(fit$Age[, "Rate"]))),
    max(abs(comp$Per - as.numeric(fit$Per[, "P-RR"])) /
        abs(as.numeric(fit$Per[, "P-RR"]))),
    max(abs(comp$Coh - as.numeric(fit$Coh[, "C-RR"])) /
        abs(as.numeric(fit$Coh[, "C-RR"])))
  )
  worst_recon <- max(worst_recon, err)
  stopifnot(err < 1e-8)
}
cat(sprintf("  Test A1 (apc_components reconstruction)  PASS  worst rel err = %.2e\n",
            worst_recon))

# Test A2: cases_draws shape + non-negative + ordering
total_cells <- 0
for (key in test_keys) {
  cd <- apc_results[[key]]$cases_draws
  stopifnot(length(dim(cd)) == 3)
  stopifnot(dim(cd)[1] == 24, dim(cd)[2] == 9, dim(cd)[3] == 1000)
  stopifnot(all(cd >= 0))
  s <- apc_results[[key]]$cases_summary
  stopifnot(all(s$cases_p025 <= s$cases_mid + 1e-9))
  stopifnot(all(s$cases_mid  <= s$cases_p975 + 1e-9))
  stopifnot(all(s$cases_p025 >= 0))
  total_cells <- total_cells + nrow(s)
}
cat(sprintf("  Test A2 (cases_draws shape + p025<=mid<=p975)  PASS  %d cells\n",
            total_cells))

# Test A3: incidence_projections.csv structure
proj <- read_csv("output/incidence_projections.csv", show_col_types = FALSE)
stopifnot(all(c("P", "cases_mid", "cases_p025", "cases_p975",
                "asr_mid", "asr_p025", "asr_p975",
                "subtype", "sex", "tier") %in% names(proj)))
stopifnot(all(proj$cases_p025 <= proj$cases_mid + 1e-9))
stopifnot(all(proj$cases_mid  <= proj$cases_p975 + 1e-9))
stopifnot(all(proj$cases_p025 >= 0))
# 312 = 10 fits x 24 years + headline total (HL+NHL) x (2 sex + persons) x 24
stopifnot(nrow(proj) == 312)
stopifnot("total" %in% proj$subtype && "persons" %in% proj$sex)
# incidence total = HL + NHL at 2045 (draw-median ~= sum-of-medians), persons
# HL/NHL incidence are emitted by sex; the total is by sex and persons.
# Check total = HL + NHL where both exist (males, females).
i45 <- proj |> filter(P == 2045, sex %in% c("males", "females")) |>
  select(subtype, sex, cases_mid) |>
  tidyr::pivot_wider(names_from = subtype, values_from = cases_mid)
stopifnot(all(abs((i45$total - (i45$hodgkin + i45$nhl)) / i45$total) < 0.02))
cat(sprintf("  Test A3 (incidence_projections.csv structure + total)  PASS  %d rows\n",
            nrow(proj)))

# Test A4: apc_effects.csv structure (rr_p025 <= rr <= rr_p975)
eff <- read_csv("output/apc_effects.csv", show_col_types = FALSE)
stopifnot(all(c("effect", "value", "rr", "rr_p025", "rr_p975",
                "subtype", "sex", "tier") %in% names(eff)))
stopifnot(all(eff$rr_p025 <= eff$rr + 1e-9))
stopifnot(all(eff$rr      <= eff$rr_p975 + 1e-9))
stopifnot(all(eff$rr_p025 > 0))
cat(sprintf("  Test A4 (apc_effects.csv structure)  PASS  %d rows\n",
            nrow(eff)))

# Test A5: figure_1 + supporting tables exist (figure size > 1 KB; CSVs
# allowed smaller because Table S1 with 10 fits is genuinely tiny ~ 0.5 KB)
.exists_nonempty <- function(f, min_bytes = 200) {
  p <- file.path("output", f)
  stopifnot(file.exists(p))
  stopifnot(file.info(p)$size > min_bytes)
}
.exists_nonempty("figure_1_combined_projections.pdf", min_bytes = 1024)
.exists_nonempty("figure_1_combined_projections.png", min_bytes = 1024)
.exists_nonempty("table_1_incidence_summary.csv")
.exists_nonempty("table_s1_apc_fit_stats.csv")
.exists_nonempty("incidence_historical_asr.csv")
cat("  Test A5 (Figure 1 + Tables 1 / S1 / hist exist)  PASS\n")

# Test A6: project_one_draw runtime (under 50 ms/call)
fit  <- apc_results[["dlbcl_males"]]$fit
df   <- prepare_fit_data("dlbcl_males", "subtype",
                         inputs_for_basis$inc_subtype,
                         inputs_for_basis$inc_agg,
                         inputs_for_basis$pop_hist)
basis   <- apc_basis(fit, df)
pop_proj <- read_csv("data/pop_proj.csv", show_col_types = FALSE)
pop_mat <- build_pop_mat("males", pop_proj)
b0 <- coef(fit$Model)
n_iter <- 100
t0 <- Sys.time()
for (i in seq_len(n_iter)) project_one_draw(b0, basis, pop_mat)
elapsed <- as.numeric(Sys.time() - t0, units = "secs") / n_iter
stopifnot(elapsed < 0.05)
cat(sprintf("  Test A6 (project_one_draw <50 ms/call)  PASS  %.2f ms/call\n",
            elapsed * 1e3))

# -----------------------------------------------------------------
# prev_model.R tests
# -----------------------------------------------------------------
cat("\n========== prev_model.R ==========\n")
.maybe_run_prev()

inputs   <- load_prev_inputs()
surv_obs <- build_surv_obs(inputs$surv_raw)

# Test P1: surv_obs structure (donor_age_group set; no NA survival_pct)
stopifnot("donor_age_group" %in% names(surv_obs))
stopifnot(sum(is.na(surv_obs$survival_pct)) == 0)
filled <- surv_obs |> filter(was_filled)
stopifnot(all(!is.na(filled$donor_age_group)))
unfilled <- surv_obs |> filter(!was_filled)
stopifnot(all(unfilled$donor_age_group == unfilled$age_group))
cat(sprintf("  Test P1 (surv_obs + donor_age_group)  PASS  %d filled cells\n",
            nrow(filled)))

# Test P2: sample_surv_draws - in [0,1], filled cells share donor draws
draws <- sample_surv_draws(surv_obs, B = 200, seed = 20260507)
all_draws <- unlist(draws$surv_draws)
stopifnot(all(is.finite(all_draws)))
stopifnot(all(all_draws >= 0 & all_draws <= 1))
filled_idx <- which(draws$was_filled)
checked <- 0
for (i in head(filled_idx, 50)) {
  donor <- draws |>
    filter(subtype        == draws$subtype[i],
           sex            == draws$sex[i],
           period         == draws$period[i],
           years_since_dx == draws$years_since_dx[i],
           age_group      == draws$donor_age_group[i],
           !was_filled)
  stopifnot(nrow(donor) == 1)
  stopifnot(identical(draws$surv_draws[[i]], donor$surv_draws[[1]]))
  checked <- checked + 1
}
cat(sprintf("  Test P2 (surv draws domain + donor identity)  PASS  %d filled checked\n",
            checked))

# Test P3: sample_improvement_draws - finite, ordering of source classifications
imp_pt <- build_improvement_final(surv_obs)
imp_draws <- sample_improvement_draws(draws, imp_pt, B = 200)
stopifnot(all(is.finite(unlist(imp_draws$rate_draws))))
expected_sources <- c("age-specific", "fallback (filled cell)",
                      "fallback (no data)", "fallback (non-improving)")
stopifnot(all(imp_draws$source %in% expected_sources))
cat(sprintf("  Test P3 (improvement rate draws domain + source)  PASS  %d cells\n",
            nrow(imp_draws)))

# Test P4: build_surv_array - shape + in [0,1], stopifnot(!any(NA))
arr <- build_surv_array("dlbcl", "males",
                        dx_years = 1982:2045, cap_year = improve_cap,
                        draws, imp_draws, B = 200)
stopifnot(all(dim(arr) == c(64, 9, 40, 200)))
stopifnot(all(arr >= 0 & arr <= 1))
stopifnot(!any(is.na(arr)))
cat(sprintf("  Test P4 (surv_array shape + in [0,1])  PASS  %s\n",
            paste(dim(arr), collapse = "x")))

# Test P5: prevalence outputs - p025 <= mid <= p975, p025 >= 0
prev <- read_csv("output/prevalence_projections.csv", show_col_types = FALSE)
stopifnot(all(c("year", "prev_mid", "prev_p025", "prev_p975",
                "prev_mean", "duration", "subtype", "sex") %in% names(prev)))
stopifnot(all(prev$prev_p025 <= prev$prev_mid + 1e-9))
stopifnot(all(prev$prev_mid  <= prev$prev_p975 + 1e-9))
stopifnot(all(prev$prev_p025 >= 0))
stopifnot(setequal(unique(prev$duration), c(2, 3, 5, 10, 40)))
stopifnot(setequal(unique(prev$subtype),
                   c("hodgkin", "nhl", "dlbcl", "follicular", "mantle_cell", "total")))
stopifnot(setequal(unique(prev$sex), c("males", "females", "persons")))
# 3060 = [5 series x 2 sex] + [5 series persons] + [total x 2 sex] + [total persons]
#        all x 5 durations x 34 years = (10 + 5 + 2 + 1) x 5 x 34
stopifnot(nrow(prev) == 3060)
cat(sprintf("  Test P5 (prevalence_projections.csv structure)  PASS  %d rows\n",
            nrow(prev)))

# Test P5b: headline TOTAL = HL + NHL (draw-median ~= sum-of-medians), and the
# NHL subtypes are a SUBSET of NHL (sum <= NHL) and are NEVER the total.
pw <- prev |> select(subtype, sex, duration, year, prev_mid) |>
  tidyr::pivot_wider(names_from = subtype, values_from = prev_mid)
# total == HL + NHL (persons + by sex); median-of-sum vs sum-of-medians < 1.5%
rel_err <- abs(pw$total - (pw$hodgkin + pw$nhl)) / pw$total
stopifnot(all(rel_err < 0.015))
# subtypes (DLBCL+FL+MCL) <= aggregate NHL for every year, sex, duration
sub_sum <- pw$dlbcl + pw$follicular + pw$mantle_cell
stopifnot(all(sub_sum <= pw$nhl + 1e-6))
# and they cover 40-85% of NHL (a real subset, never ~100% and never the total)
share <- sub_sum / pw$nhl
stopifnot(all(share > 0.40 & share < 0.90))
stopifnot(all(sub_sum < pw$total))   # subtypes never equal/exceed the total
cat(sprintf("  Test P5b (total = HL+NHL; subtypes subset of NHL, %.0f-%.0f%%)  PASS\n",
            100 * min(share), 100 * max(share)))

# Test P6: sensitivity outputs structure
prev_sens <- read_csv("output/prevalence_sensitivity.csv", show_col_types = FALSE)
stopifnot(setequal(unique(prev_sens$duration), c(2, 3, 5, 10, 40)))
stopifnot(nrow(prev_sens) == 9180)   # 3060 x 3 scenarios
stopifnot(all(c("conservative", "base", "optimistic") %in% prev_sens$scenario))
stopifnot(all(prev_sens$prev_p025 <= prev_sens$prev_mid + 1e-9))
stopifnot(all(prev_sens$prev_mid  <= prev_sens$prev_p975 + 1e-9))
cat(sprintf("  Test P6 (prevalence_sensitivity.csv structure)  PASS  %d rows\n",
            nrow(prev_sens)))

# Test P6b: table_3_prevalence_combined.csv - combined-duration summary with
# the headline total; 6 series x 5 durations; CrI ordering + sign.
t3 <- read_csv("output/table_3_prevalence_combined.csv", show_col_types = FALSE)
stopifnot(nrow(t3) == 30)
stopifnot(setequal(unique(t3$subtype),
                   c("hodgkin", "nhl", "dlbcl", "follicular", "mantle_cell", "total")))
stopifnot(all(t3$subset_of_nhl[t3$subtype %in% c("dlbcl","follicular","mantle_cell")]))
stopifnot(all(!t3$subset_of_nhl[t3$subtype %in% c("hodgkin","nhl","total")]))
stopifnot(all(t3$prev_2021_p025 <= t3$prev_2021 & t3$prev_2021 <= t3$prev_2021_p975))
stopifnot(all(t3$prev_2045_p025 <= t3$prev_2045 & t3$prev_2045 <= t3$prev_2045_p975))
stopifnot(all(t3$pct_change_p025 <= t3$pct_change & t3$pct_change <= t3$pct_change_p975))
stopifnot(all(t3$prev_2045 > t3$prev_2021) && all(t3$pct_change > 0))  # all grow
cat(sprintf("  Test P6b (table_3 combined structure + CrIs)  PASS  %d rows\n", nrow(t3)))

# Test P7: 10-year validation, summary, figs 2 & 3 exist
.exists_nonempty("prevalence_validation_10yr.csv")
.exists_nonempty("table_2_prevalence_summary.csv")
.exists_nonempty("figure_2_prevalence_5yr.pdf",  min_bytes = 1024)
.exists_nonempty("figure_2_prevalence_5yr.png",  min_bytes = 1024)
.exists_nonempty("figure_3_prevalence_40yr.pdf", min_bytes = 1024)
.exists_nonempty("figure_3_prevalence_40yr.png", min_bytes = 1024)
cat("  Test P7 (validation_10yr + Table 2 + Figs 2 & 3 exist)  PASS\n")

# Test P8: back-estimation produced data/incidence_subtype_back.csv
stopifnot(file.exists("data/incidence_subtype_back.csv"))
inc_back <- read_csv("data/incidence_subtype_back.csv", show_col_types = FALSE)
stopifnot(nrow(inc_back) > 0)
stopifnot(all(c("year", "sex", "age_group", "subtype", "count") %in% names(inc_back)))
stopifnot(min(inc_back$year) <= 1982 && max(inc_back$year) <= 2002)
cat(sprintf("  Test P8 (back-estimation csv)  PASS  %d rows, years %d-%d\n",
            nrow(inc_back), min(inc_back$year), max(inc_back$year)))

# -----------------------------------------------------------------
# supplement.R tests
# -----------------------------------------------------------------
cat("\n========== supplement.R ==========\n")
.maybe_run_supp()

# Test S1: all 9 SI tables exist + parse + non-trivial
si_tables <- c(
  "table_s1_apc_fit_stats.csv",
  "table_s2a_incidence_aihw.csv",
  "table_s2b_prevalence_validation.csv",
  "table_s2c_40yr_prevalence.csv",
  "table_s2d_consistency.csv",
  "table_s3a_survival_cleaning.csv",
  "table_s3b_improvement_rates.csv",
  "table_s3b_improvement_rates_full.csv",
  "table_s3c_back_estimation.csv",
  "table_s4_sensitivity.csv"
)
for (f in si_tables) {
  p <- file.path("output", f)
  stopifnot(file.exists(p))
  d <- read_csv(p, show_col_types = FALSE)
  stopifnot(nrow(d) > 0)
}
cat(sprintf("  Test S1 (SI tables)  PASS  %d files\n", length(si_tables)))

# Test S2: all SI figures exist + non-empty (PDF and PNG)
si_figs <- c(
  "figure_s1a_age_effects",
  "figure_s1b_period_effects",
  "figure_s1c_cohort_effects",
  "figure_s3a_survival_curves",
  "figure_s3b_survival_example",
  "figure_s3c_prevalence_2yr",
  "figure_s3d_prevalence_3yr",
  "figure_s3e_prevalence_10yr",
  "figure_s3f_survival_trend",
  "figure_s3g_survival_trend",
  "figure_s4a_sensitivity_5yr",
  "figure_s4b_sensitivity_40yr"
)
for (f in si_figs) for (ext in c(".pdf", ".png")) {
  .exists_nonempty(paste0(f, ext), min_bytes = 1024)
}
cat(sprintf("  Test S2 (SI figures)  PASS  %d figs x 2 formats = %d files\n",
            length(si_figs), 2 * length(si_figs)))

# Test S3: sensitivity path (table_s4). These cover the S4 summarisation that
# the reframing regressed (a sum over sex that double-counted the persons
# rows); the total = HL+NHL check elsewhere did not reach this path.
s4 <- read_csv("output/table_s4_sensitivity.csv", show_col_types = FALSE)
# (a) base case == table_3 total at 2021 and 2045 (definitive doubling catch)
t3tot  <- t3 |> filter(subtype == "total") |> select(duration, t21 = prev_2021, t45 = prev_2045)
s4wide <- s4 |> filter(year %in% c(2021, 2045)) |> select(duration, year, base) |>
  tidyr::pivot_wider(names_from = year, values_from = base, names_prefix = "s")
s4cmp  <- dplyr::left_join(s4wide, t3tot, by = "duration")
stopifnot(all(abs(s4cmp$s2021 - s4cmp$t21) < 1), all(abs(s4cmp$s2045 - s4cmp$t45) < 1))
# (b) monotone non-decreasing in year, every duration x scenario (prevalence can't fall)
for (col in c("conservative", "base", "optimistic")) {
  for (d in unique(s4$duration)) {
    v <- s4 |> filter(duration == d) |> arrange(year) |> (\(x) x[[col]])()
    stopifnot(all(diff(v) >= -1e-6))
  }
}
# (c) double-count canary: the total (HL+NHL) is strictly less than
#     HL + NHL + subtypes (the subtypes are a subset of NHL, so extra mass)
allsum <- prev_sens |>
  filter(scenario == "base", sex == "persons",
         subtype %in% c("hodgkin", "nhl", "dlbcl", "follicular", "mantle_cell")) |>
  group_by(duration, year) |> summarise(s = sum(prev_mid), .groups = "drop")
tot_ps <- prev_sens |>
  filter(scenario == "base", sex == "persons", subtype == "total") |>
  select(duration, year, tot = prev_mid)
canary <- dplyr::left_join(tot_ps, allsum, by = c("duration", "year"))
stopifnot(all(canary$tot < canary$s - 1e-6))
# (d) conservative <= base <= optimistic at every duration x year
stopifnot(all(s4$conservative <= s4$base + 1e-6),
          all(s4$base <= s4$optimistic + 1e-6))
cat("  Test S3 (S4 sensitivity: base==table_3, monotone, no double-count, cons<=base<=opt)  PASS\n")

# Test S4: every SI figure that shows lymphoma as a dimension must carry all
# FIVE series (headline HL + aggregate NHL, and the three subtypes). The broken
# builders iterated the four-lymphoma `subtypes` set or filtered out "nhl";
# two prior re-runs each fixed only the figures they touched and nothing caught
# the rest. This calls each such builder and inspects its plotted data.
apc_eff <- read_csv("output/apc_effects.csv", show_col_types = FALSE)
inc_sub <- read_csv("data/incidence_subtype.csv", show_col_types = FALSE)
inc_ag  <- read_csv("data/incidence_agg.csv", show_col_types = FALSE)
five    <- c("hodgkin", "nhl", "dlbcl", "follicular", "mantle_cell")
.fig_series <- function(p) {
  d   <- p$data
  col <- intersect(c("subtype", "subtype_label", "label"), names(d))[1]
  sort(unique(as.character(d[[col]])))
}
lymphoma_figs <- list(
  "S1 APC effects"      = .build_fig_s1(apc_eff, "age", "x", "y"),
  "S3a survival curves" = build_figure_s3a(surv_obs, imp_pt, build_projected_surv_pt),
  "S3f survival trend"  = build_figure_survival_trend(surv_obs, imp_pt, inc_sub, inc_ag,
                                                      build_projected_surv_pt),
  "S4 sensitivity"      = build_figure_s4(prev_sens, 5)
)
for (nm in names(lymphoma_figs)) {
  got <- .fig_series(lymphoma_figs[[nm]])
  if (!setequal(got, five)) {
    stop("Figure ", nm, " is missing series: ",
         paste(setdiff(five, got), collapse = ", "),
         " (plots: ", paste(got, collapse = ", "), ")")
  }
}
cat(sprintf("  Test S4 (SI lymphoma figures carry all 5 series)  PASS  %d builders\n",
            length(lymphoma_figs)))

# -----------------------------------------------------------------
# Reproducibility test
# -----------------------------------------------------------------
cat("\n========== Reproducibility ==========\n")
# Sample beta draws twice with the same seed; verify identical
fit  <- apc_results[["dlbcl_males"]]$fit
set.seed(20260507)
b1 <- sample_beta_draws(fit, B = 100)
set.seed(20260507)
b2 <- sample_beta_draws(fit, B = 100)
stopifnot(identical(b1, b2))
# Same for survival draws
set.seed(20260507)
d1 <- sample_surv_draws(surv_obs, B = 50, seed = 20260507)
set.seed(20260507)
d2 <- sample_surv_draws(surv_obs, B = 50, seed = 20260507)
stopifnot(identical(d1$surv_draws, d2$surv_draws))
cat("  Test R1 (beta + surv draws bitwise reproducible)  PASS\n")

cat("\n=========================================\n")
cat("All structural tests passed.\n")
cat("=========================================\n")
