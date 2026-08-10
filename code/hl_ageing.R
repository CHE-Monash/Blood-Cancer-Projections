# ================
# Blood Cancer Projections (lymphoma): Hodgkin lymphoma ageing check
# ================
# Addresses the co-author question (TEG, EH, MB, JW): HL is (partly) a
# younger-onset disease, so why does population ageing drive its projected
# increase — is the older (second) bimodal peak growing, or is the median
# age at diagnosis shifting? Answers by computing, for HL by sex:
#   (i)  mean & median age at diagnosis per year (observed 1990-2021 +
#        projected 2022-2045), and the share of cases at 55+/65+;
#   (ii) the age-at-diagnosis distribution at 1990 / 2021 / 2045.
#
# Uses observed HL incidence (data/incidence_agg.csv) and the projected HL
# case draws (output/apc_results.rds). Diagnostic driver, run separately:
#   Rscript code/hl_ageing.R -> output/table_hl_ageing.csv + figure_hl_ageing.{pdf,png}
# ================

suppressPackageStartupMessages({
  library(tidyverse); library(cowplot); library(grid)
})
source("code/_setup.R")

.wmedian <- function(a, w) {
  o <- order(a); a <- a[o]; w <- w[o]
  a[which(cumsum(w) / sum(w) >= 0.5)[1]]
}

run_hl_ageing <- function(save_dir = "output", verbose = TRUE) {
  inc_agg <- readr::read_csv("data/incidence_agg.csv", show_col_types = FALSE)
  apc     <- readRDS(file.path(save_dir, "apc_results.rds"))

  obs <- inc_agg |>
    dplyr::filter(cancer_group == "Hodgkin lymphoma", year >= agg_start, year <= hist_end) |>
    dplyr::select(year, sex, age_group, count) |>
    dplyr::mutate(A = age_mid[age_group])

  proj <- purrr::map_dfr(sexes, function(sx) {
    cd <- apc[[paste0("hodgkin_", sx)]]$cases_draws       # (proj_years x age x B)
    m  <- apply(cd, c(1, 2), median)
    tibble::as_tibble(as.table(m)) |>
      setNames(c("year", "age_group", "count")) |>
      dplyr::mutate(year = as.integer(as.character(year)), sex = sx,
                    A = age_mid[age_group])
  })
  hl <- dplyr::bind_rows(obs, proj)

  age_trend <- hl |>
    dplyr::group_by(sex, year) |>
    dplyr::summarise(
      mean_age     = sum(A * count) / sum(count),
      median_age   = .wmedian(A, count),
      share_55plus = sum(count[A >= 55]) / sum(count) * 100,
      share_65plus = sum(count[A >= 65]) / sum(count) * 100,
      period       = if (first(year) <= hist_end) "Observed" else "Projected",
      .groups = "drop"
    ) |>
    dplyr::mutate(period = if_else(year <= hist_end, "Observed", "Projected"))
  readr::write_csv(age_trend, file.path(save_dir, "table_hl_ageing.csv"))

  # Figure: age-at-diagnosis distribution at 1990 / 2021 / 2045, by sex.
  dist <- hl |>
    dplyr::filter(year %in% c(1990, hist_end, proj_end)) |>
    dplyr::group_by(sex, year, age_group) |>
    dplyr::summarise(count = sum(count), .groups = "drop") |>
    dplyr::group_by(sex, year) |>
    dplyr::mutate(pct = count / sum(count) * 100) |>
    dplyr::ungroup() |>
    dplyr::mutate(A = age_mid[age_group], year = factor(year))

  p <- ggplot(dist, aes(x = A, y = pct, colour = year)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
    facet_wrap(~ sex, ncol = 2, labeller = labeller(sex = sex_labels)) +
    scale_colour_manual(values = c("1990" = "#4575b4", "2021" = "#984ea3",
                                   "2045" = "#d73027"), name = "Year") +
    scale_x_continuous(breaks = unname(age_mid), labels = age_groups) +
    labs(x = "Age at diagnosis", y = "Share of HL cases (%)",
         subtitle = "The younger-adult peak shrinks and the older peak grows as cases shift to older ages") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          legend.position = "bottom")
  save_fig(p, file.path(save_dir, "figure_hl_ageing"), width = 10, height = 5.5)

  if (verbose) {
    cat("=== HL age at diagnosis (mean) and older-age share, key years ===\n")
    print(as.data.frame(age_trend |>
      dplyr::filter(year %in% c(1990, 2010, hist_end, 2035, proj_end)) |>
      dplyr::select(sex, year, mean_age, median_age, share_55plus, share_65plus) |>
      dplyr::mutate(dplyr::across(where(is.double), ~round(., 1)))), row.names = FALSE)
    cat("\nWrote table_hl_ageing.csv + figure_hl_ageing.{pdf,png}\n")
  }
  invisible(age_trend)
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_hl_ageing())
}
