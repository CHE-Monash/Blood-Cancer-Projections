# ================
# Blood Cancer Projections (lymphoma): data behind Figures 1-3
# ================
# Writes one Excel workbook per main-text figure, holding every plotted
# series, so that the figures can be redrawn from the numbers alone:
#
#   output/figure_1_data.xlsx   incidence: modelled rate (fitted history and
#                               projection with 95% CrI), observed incidence,
#                               AIHW projected estimates
#   output/figure_2_data.xlsx   5-year prevalence: modelled (with 95% CrI) and
#                               AIHW published estimates
#   output/figure_3_data.xlsx   40-year prevalence: modelled (with 95% CrI)
#
# HOW IDENTITY WITH THE FIGURES IS GUARANTEED. The script calls the same
# builders that drew the figures (build_figure_1, build_prev_fig) on the same
# inputs and exports the data frames held by the resulting plot's layers. It
# does not re-implement the plotting filters.
#
# NO REFIT. Sourcing apc_model.R and prev_model.R defines functions only
# (their drivers sit behind sys.nframe() == 0 guards). The fitted historical
# line of Figure 1 is not in any output CSV, so it is rebuilt from the stored
# fits in output/apc_results.rds with calc_fitted_asr(), the pure function the
# driver itself uses. Nothing is fitted, nothing is sampled, and no existing
# file in output/ or data/ is written.
#
#   Rscript code/supporting/figure_data_export.R
# ================

stopifnot(requireNamespace("openxlsx"))
source("code/apc_model.R")
source("code/prev_model.R")

.layer_data <- function(p) {
  lapply(p$layers, function(l) {
    d <- l$data
    if (inherits(d, "waiver")) p$data
    else if (is.function(d)) d(p$data)
    else d
  })
}
.geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
.one <- function(p, geom) {
  i <- which(.geoms(p) == geom)
  stopifnot(length(i) >= 1)
  .layer_data(p)[i]
}
.tidy <- function(d) {
  d |> dplyr::mutate(lymphoma = as.character(label),
                     sex = as.character(sex_lab)) |>
    dplyr::arrange(factor(lymphoma, levels = c("HL", "NHL", "DLBCL", "FL", "MCL")),
                   sex, year)
}
.write <- function(sheets, notes, path) {
  wb <- openxlsx::createWorkbook(creator = "")
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, as.data.frame(sheets[[nm]]))
  }
  openxlsx::addWorksheet(wb, "Notes")
  openxlsx::writeData(wb, "Notes", data.frame(note = notes))
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}

build_figure_1_data <- function(out_dir = "output") {
  apc_results   <- readRDS(file.path(out_dir, "apc_results.rds"))
  pop_hist      <- readr::read_csv("data/pop_hist.csv", show_col_types = FALSE)
  aihw_proj     <- readr::read_csv("data/incidence_agg_proj.csv", show_col_types = FALSE)
  proj_annual   <- readr::read_csv(file.path(out_dir, "incidence_projections.csv"),
                                   show_col_types = FALSE)
  hist_annual   <- readr::read_csv(file.path(out_dir, "incidence_historical_asr.csv"),
                                   show_col_types = FALSE)
  std_pop       <- build_std_pop(pop_hist)
  fitted_annual <- purrr::map_dfr(apc_results, function(x) {
    fit_start <- if (x$tier == "agg") agg_start else subtype_start
    calc_fitted_asr(x$fit, x$sex, x$subtype, fit_start, std_pop)
  })
  fig <- build_figure_1(proj_annual, fitted_annual, hist_annual, aihw_proj)
  p   <- attr(fig, "source_plot")
  stopifnot(identical(.geoms(p),
                      c("GeomRibbon", "GeomPoint", "GeomLine", "GeomPoint", "GeomVline")))
  ld <- .layer_data(p)
  line <- .tidy(ld[[3]]) |>
    dplyr::transmute(lymphoma, sex, year,
                     series = dplyr::if_else(year > hist_end, "Projected", "Fitted"),
                     asr_per_100000 = asr,
                     cri_lower = asr_p025, cri_upper = asr_p975)
  band <- .tidy(ld[[1]]) |>
    dplyr::transmute(lymphoma, sex, year, cri_lower = asr_p025, cri_upper = asr_p975)
  obs  <- .tidy(ld[[2]]) |>
    dplyr::transmute(lymphoma, sex, year, observed_asr_per_100000 = asr_obs)
  aihw <- .tidy(ld[[4]]) |>
    dplyr::transmute(lymphoma, sex, year, aihw_projected_asr_per_100000 = asr_2001)
  list(`Modelled rate` = line, `95% credible interval` = band,
       `Observed incidence` = obs, `AIHW projected estimates` = aihw)
}

build_prev_figure_data <- function(duration_value, with_aihw, out_dir = "output") {
  prev_base <- readr::read_csv(file.path(out_dir, "prevalence_projections.csv"),
                               show_col_types = FALSE)
  inputs <- if (with_aihw) list(
    prev_subtype = readr::read_csv("data/prevalence_subtype.csv", show_col_types = FALSE),
    prev_agg     = readr::read_csv("data/prevalence_agg.csv", show_col_types = FALSE)
  ) else list(prev_subtype = NULL, prev_agg = NULL)
  fig <- build_prev_fig(prev_base, duration_value = duration_value,
                        prev_subtype_aihw = inputs$prev_subtype,
                        prev_agg_aihw     = inputs$prev_agg)
  p  <- attr(fig, "source_plot")
  g  <- .geoms(p)
  stopifnot(identical(g[1:3], c("GeomRibbon", "GeomLine", "GeomVline")))
  ld <- .layer_data(p)
  line <- .tidy(ld[[2]]) |>
    dplyr::transmute(lymphoma, sex, year,
                     period = dplyr::if_else(year > hist_end, "Projected", "Estimated"),
                     prevalence_thousands = prev_mid,
                     cri_lower_thousands = prev_lo, cri_upper_thousands = prev_hi)
  band <- .tidy(ld[[1]]) |>
    dplyr::transmute(lymphoma, sex, year,
                     cri_lower_thousands = prev_lo, cri_upper_thousands = prev_hi)
  out <- list(`Modelled prevalence` = line, `95% credible interval` = band)
  if (with_aihw) {
    stopifnot(length(g) == 4, g[4] == "GeomPoint")
    out[["AIHW published estimates"]] <- .tidy(ld[[4]]) |>
      dplyr::transmute(lymphoma, sex, year, aihw_prevalence_thousands = prev_aihw)
  } else stopifnot(length(g) == 3)
  out
}

run_figure_data_export <- function(out_dir = "output") {
  common <- c(
    "Every sheet holds the data frame of one plotted layer, taken from the plot object built by the figure's own builder.",
    "The shaded band is drawn for projected years only (2022 to 2045); the line sheet also carries the interval limits for every modelled year.",
    "DLBCL, FL and MCL are subtypes of aggregate NHL and do not sum to it.")
  f1 <- build_figure_1_data(out_dir)
  .write(f1, c("Figure 1. Age-standardised incidence rate per 100,000 (2001 Australian Standard Population).",
               "Fitted: model-smoothed historical rate (point estimate; no interval is drawn). Projected: median of 1,000 Monte Carlo draws with 95% credible interval.",
               common), file.path(out_dir, "figure_1_data.xlsx"))
  f2 <- build_prev_figure_data(5, TRUE, out_dir)
  .write(f2, c("Figure 2. 5-year prevalence, thousands of people.",
               "AIHW published estimates are shown for 2012 to 2021.", common),
         file.path(out_dir, "figure_2_data.xlsx"))
  f3 <- build_prev_figure_data(40, FALSE, out_dir)
  .write(f3, c("Figure 3. 40-year prevalence, thousands of people.", common),
         file.path(out_dir, "figure_3_data.xlsx"))
  cat("Wrote figure_1_data.xlsx, figure_2_data.xlsx, figure_3_data.xlsx to", out_dir, "\n")
  invisible(list(figure_1 = f1, figure_2 = f2, figure_3 = f3))
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_figure_data_export())
}
