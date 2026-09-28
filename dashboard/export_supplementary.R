# ================
# Dashboard: additional results (static figures and downloads)
# ================
# Copies existing analysis outputs into dashboard/supplementary/ so that the
# "Additional results" page (additional-results.qmd) can show them, and writes
# a manifest recording where each published file came from and its MD5.
#
# READ-ONLY with respect to the analysis: files are copied byte for byte from
# output/. Nothing is fitted, sampled, redrawn or recalculated, and output/ is
# not written to. The copied files are committed, because the publishing
# workflow renders the site without running R.
#
# From the repository root:
#   Rscript dashboard/export_supplementary.R
# ================

out_dir  <- "output"
dest     <- file.path("dashboard", "supplementary")
fig_dest <- file.path(dest, "figures")
tab_dest <- file.path(dest, "tables")

figures <- data.frame(
  label = sprintf("A%d", 1:14),
  title = c(
    "Estimated age effects, by lymphoma and sex",
    "Estimated period effects, by lymphoma and sex",
    "Estimated cohort effects, by lymphoma and sex",
    "Incidence holdout: projections from models fitted to 2011 against observed incidence, 2012 to 2021",
    "Survival curves used in the prevalence computation, age group 65 to 74, diagnosed 2019",
    "Worked example of one modelled survival curve",
    "Australian population age structure, 2021 and 2045",
    "Projected 2-year prevalence, by lymphoma and sex",
    "Projected 3-year prevalence, by lymphoma and sex",
    "Projected 10-year prevalence, by lymphoma and sex",
    "Modelled case-weighted 5-year survival by year of diagnosis",
    "Modelled case-weighted 10-year survival by year of diagnosis",
    "Sensitivity of 5-year prevalence to the survival improvement cap year",
    "Sensitivity of 40-year prevalence to the survival improvement cap year"),
  source = c(
    "figure_s1_age_effects.png", "figure_s2_period_effects.png",
    "figure_s3_cohort_effects.png", "figure_holdout_validation.png",
    "figure_s4_survival_curves.png", "figure_s5_survival_example.png",
    "figure_s6_population_pyramid.png", "figure_s7_prevalence_2yr.png",
    "figure_s8_prevalence_3yr.png", "figure_s9_prevalence_10yr.png",
    "figure_s10_survival_trend_5yr.png", "figure_s11_survival_trend_10yr.png",
    "figure_s12_sensitivity_5yr.png", "figure_s13_sensitivity_40yr.png"),
  stringsAsFactors = FALSE
)

tables <- data.frame(
  label = c("A1", "A2", "A3", "A4", "A5", "A5", "A6", "A7", "A7", "A7",
            "A8", "A8", "A9", "A9", "A9", "A9"),
  title = c(
    "Fit statistics for the ten age-period-cohort models",
    "Back-estimated 2000 against observed 2003 subtype incidence",
    "Aggregate incidence projections against AIHW projections, 2026 to 2035",
    "Modelled NHL subtypes as a share of aggregate NHL incidence",
    "Annual survival improvement rates, summary",
    "Annual survival improvement rates, every age group, with the source of each rate",
    "40-year prevalence at 2021 against AIHW published values",
    "Incidence, 2021 and 2045, with the change and its 95% credible interval, by sex",
    "Prevalence at 2, 3, 5, 10 and 40 years, 2021 and 2045, with the change and its 95% credible interval, by sex",
    "Prevalence at 2, 3, 5, 10 and 40 years, 2021 and 2045, with the change and its 95% credible interval, persons",
    "Survival improvement scenarios: total lymphoma prevalence at key years",
    "Survival improvement scenarios: prevalence by lymphoma, sex, duration and year",
    "Knot selection grid", "Start-year sensitivity",
    "Damping sensitivity", "COVID-19 (2020 to 2021) exclusion sensitivity"),
  source = c(
    "table_s2_apc_fit_stats.csv", "table_s3_back_estimation.csv",
    "table_s4_incidence_aihw.csv", "table_s5_consistency.csv",
    "table_s7_improvement_rates.csv", "table_s7_improvement_rates_full.csv",
    "table_s9_40yr_prevalence.csv",
    "table_1_change_cri.csv", "table_2_change_cri.csv",
    "table_3_prevalence_combined.csv",
    "table_s12_sensitivity.csv", "prevalence_sensitivity.csv",
    "table_s_knot_selection.csv", "table_s_start_year_sensitivity.csv",
    "table_s_damping_sensitivity.csv", "table_s_covid_sensitivity.csv"),
  stringsAsFactors = FALSE
)

copy_set <- function(df, kind, to) {
  dir.create(to, showWarnings = FALSE, recursive = TRUE)
  from <- file.path(out_dir, df$source)
  stopifnot(all(file.exists(from)))
  ok <- file.copy(from, file.path(to, df$source), overwrite = TRUE,
                  copy.date = TRUE)
  stopifnot(all(ok))
  data.frame(
    kind = kind, label = df$label, title = df$title,
    source = file.path("output", df$source),
    published = file.path("supplementary", basename(to), df$source),
    md5 = unname(tools::md5sum(file.path(to, df$source))),
    source_md5 = unname(tools::md5sum(from)),
    stringsAsFactors = FALSE
  )
}

manifest <- rbind(copy_set(figures, "figure", fig_dest),
                  copy_set(tables, "table", tab_dest))
stopifnot(identical(manifest$md5, manifest$source_md5))
utils::write.csv(manifest[, c("kind", "label", "title", "source", "published", "md5")],
                 file.path(dest, "manifest.csv"), row.names = FALSE)
cat(sprintf("Copied %d figures and %d table files to %s; manifest written.\n",
            nrow(figures), nrow(tables), dest))
