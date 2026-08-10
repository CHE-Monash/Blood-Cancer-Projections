# ================
# Blood Cancer Projections (lymphoma): fan chart for 5-year FL prevalence
# ================
# Standalone helper — not part of the main pipeline and not used by the
# manuscript or SI. Produced for a co-author who asked to see the full
# uncertainty distribution rather than only the 95% credible interval.
#
# Outputs:
#   output/fan_chart_fl_5yr.{pdf,png}
#   output/fan_chart_fl_5yr_percentiles.csv
#
# Usage (from project root):
#   source("code/fan_chart_fl_5yr.R")
#   run_fl_fan()
# or:
#   Rscript code/fan_chart_fl_5yr.R
#
# Requires apc_results.rds in output/ (run apc_model.R first).
# ================

set.seed(20260507)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

source("code/_setup.R")
source("code/prev_model.R")  # for load_prev_inputs, build_surv_obs,
# build_improvement_final, sample_surv_draws,
# sample_improvement_draws, run_back_estimation,
# build_cases_array, build_surv_array

# ----------------------
# Per-draw prevalence (variant of compute_prevalence that keeps the full matrix)
# ----------------------

compute_prev_draws <- function(cases_array, surv_array, calc_years, duration) {
  B <- dim(cases_array)[3]
  cases_years <- as.numeric(dimnames(cases_array)[[1]])
  surv_years  <- as.numeric(dimnames(surv_array)[[1]])
  min_dx_year <- max(min(cases_years), min(surv_years))
  
  prev_draws <- matrix(0, nrow = length(calc_years), ncol = B)
  rownames(prev_draws) <- as.character(calc_years)
  
  for (ci in seq_along(calc_years)) {
    cy <- calc_years[ci]
    for (t in 0:(duration - 1)) {
      dy <- cy - t
      if (dy < min_dx_year) next
      dy_str <- as.character(dy)
      cm <- cases_array[dy_str, , ]
      sm <- if (t == 0) sqrt(surv_array[dy_str, , "1", ])
      else surv_array[dy_str, , as.character(t), ]
      prev_draws[ci, ] <- prev_draws[ci, ] + colSums(cm * sm)
    }
  }
  prev_draws
}

# ----------------------
# Compute percentile bands for FL 5-year prevalence, both sexes
# ----------------------

build_fl_fan <- function(B = 1000, seed = 20260507,
                         duration = 5,
                         calc_years = 2012:proj_end) {
  set.seed(seed)
  
  apc_path <- "output/apc_results.rds"
  if (!file.exists(apc_path)) {
    stop(apc_path, " missing - run apc_model.R first")
  }
  apc_results <- readRDS(apc_path)
  inputs <- load_prev_inputs()
  
  # Re-create the survival / improvement draws (cheap; ~5 s)
  inc_back              <- run_back_estimation(inputs$inc_subtype, inputs$inc_agg)
  surv_obs              <- build_surv_obs(inputs$surv_raw)
  improvement_pt        <- build_improvement_final(surv_obs)
  surv_draws_tbl        <- sample_surv_draws(surv_obs, B = B, seed = seed)
  improvement_draws_tbl <- sample_improvement_draws(surv_draws_tbl,
                                                    improvement_pt, B = B)
  
  probs <- c(0.005, 0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975, 0.995)
  prob_names <- sprintf("p%03d", round(probs * 1000))
  
  fl_rows <- list()
  for (sx in sexes) {
    cases_array <- build_cases_array("follicular", sx, B,
                                     inc_subtype = inputs$inc_subtype,
                                     inc_agg     = inputs$inc_agg,
                                     inc_back    = inc_back,
                                     apc_results = apc_results)
    surv_array  <- build_surv_array("follicular", sx,
                                    dx_years = 1982:proj_end,
                                    cap_year = improve_cap,
                                    surv_draws_tbl,
                                    improvement_draws_tbl,
                                    B = B)
    prev_draws <- compute_prev_draws(cases_array, surv_array,
                                     calc_years = calc_years,
                                     duration   = duration)
    
    pct <- t(apply(prev_draws, 1, function(x) quantile(x, probs)))
    colnames(pct) <- prob_names
    df <- as.data.frame(pct)
    df$year <- as.numeric(rownames(prev_draws))
    df$sex  <- sx
    fl_rows[[sx]] <- df
    
    rm(cases_array, surv_array, prev_draws); gc(verbose = FALSE)
  }
  
  bind_rows(fl_rows) |>
    relocate(year, sex)
}

# ----------------------
# Plot — Figure 2 style, single subtype, multiple ribbons in projection period
# ----------------------

plot_fl_fan <- function(fan_df) {
  fl_colour <- unname(line_colours["FL"])
  
  # Convert percentile columns to thousands
  d <- fan_df |>
    mutate(across(starts_with("p"), ~ .x / 1000))
  
  proj <- d |> filter(year > hist_end)
  
  ggplot(d, aes(x = year)) +
    # 99% band (widest, palest)
    geom_ribbon(data = proj, aes(ymin = p005, ymax = p995),
                fill = fl_colour, alpha = 0.10, colour = NA) +
    # 95% band
    geom_ribbon(data = proj, aes(ymin = p025, ymax = p975),
                fill = fl_colour, alpha = 0.15, colour = NA) +
    # 80% band
    geom_ribbon(data = proj, aes(ymin = p100, ymax = p900),
                fill = fl_colour, alpha = 0.20, colour = NA) +
    # 50% band (interquartile)
    geom_ribbon(data = proj, aes(ymin = p250, ymax = p750),
                fill = fl_colour, alpha = 0.25, colour = NA) +
    # Median line across full period
    geom_line(aes(y = p500), colour = fl_colour, linewidth = 0.8) +
    geom_vline(xintercept = hist_end + 0.5, linetype = "dotted",
               colour = "grey50", linewidth = 0.4) +
    facet_wrap(~ sex, ncol = 2, labeller = labeller(sex = sex_labels)) +
    labs(x = "Year",
         y = "5-year prevalence (thousands)",
         subtitle = "Follicular lymphoma — 50%, 80%, 95%, 99% credible bands (light to dark from outside in)") +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.background  = element_rect(colour = NA, fill = NA))
}

# ----------------------
# Driver
# ----------------------

run_fl_fan <- function(B = 1000, save_dir = "output", verbose = TRUE) {
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  t0 <- Sys.time()
  if (verbose) cat("Building FL fan chart data...\n")
  fan_df <- build_fl_fan(B = B)
  if (verbose) cat(sprintf("  done in %.1f s\n",
                           as.numeric(Sys.time() - t0, units = "secs")))
  
  p <- plot_fl_fan(fan_df)
  
  save_fig(p, file.path(save_dir, "fan_chart_fl_5yr"),
           width = 12, height = 6)
  write_csv(fan_df, file.path(save_dir, "fan_chart_fl_5yr_percentiles.csv"))
  
  if (verbose) {
    cat(sprintf("Wrote fan_chart_fl_5yr.{pdf,png,csv} to %s/\n", save_dir))
  }
  invisible(list(data = fan_df, plot = p))
}

# ----------------------
# Auto-run guard
# ----------------------

if (!interactive() && sys.nframe() == 0) {
  invisible(run_fl_fan())
}