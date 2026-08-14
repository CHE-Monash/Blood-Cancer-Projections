# ================
# Blood Cancer Projections (lymphoma): shared setup
# ================
# Sourced at the top of apc_model.R, prev_model.R, and supplement.R.
# Defines constants, visual conventions, and figure helpers used across
# the pipeline. No side effects beyond defining values and functions in
# the calling environment; no auto-run guard.
#
# Sourcing this file does NOT load packages (each driver loads what it
# specifically needs, since requirements differ slightly between files
# and library() modifies the search path).
# ================

# ----------------------
# Years and projection parameters
# ----------------------
damping       <- 0.92    # period/cohort drift damping per year beyond hist_end
hist_end      <- 2021    # last year of AIHW incidence + survival data
proj_start    <- 2022    # first projection year
proj_end      <- 2045    # last projection year
agg_start     <- 1990    # aggregate-tier APC fit start (NHL, HL); see below
subtype_start <- 2003    # subtype-tier APC fit start (DLBCL, FL, MCL)
# agg_start selected by the start-year sensitivity (code/start_year_sensitivity.R;
# _notes.md "Start-year selection"). Observed incidence runs from 1982, but fitting
# the aggregate tier from 1990 markedly improves out-of-sample validation (NHL 2021
# holdout endpoint bias falls from ~+13% to ~+2%) by excluding the transient pre-1990
# NHL surge (HIV/AIDS-associated NHL + pre-WHO reclassification), which the APC model
# would otherwise extrapolate. Holdout error plateaus for any start >= 1990, so 1990
# is the earliest well-validated year — retaining the longest series (32 yr). Both
# prep_agg_data() and prep_subtype_data() filter P >= these starts.

# APC natural-cubic-spline knot counts (age, period, cohort), selected by
# out-of-sample holdout error + parsimony (see code/knot_selection.R and
# _notes.md "Knot selection"). The aggregate tier retains the Epi default
# (5,5,5), which the holdout confirmed as near-optimal; the shorter subtype
# series generalises out-of-sample with fewer AGE knots (4 vs 5) — this
# roughly halves subtype out-of-sample error, chiefly by stabilising MCL.
# Period/cohort knots kept at 5 (holdout was insensitive/non-monotone there).
agg_knots     <- c(A = 5, P = 5, C = 5)
subtype_knots <- c(A = 4, P = 5, C = 5)
improve_cap   <- 2031    # year beyond which survival improvement is held constant (base scenario)
prop_years    <- 2003:2007  # window for back-estimation subtype proportions

# ----------------------
# Domain constants
# ----------------------
age_mid <- c(
  "5–14"  = 9.5,  "15–24" = 19.5, "25–34" = 29.5,
  "35–44" = 39.5, "45–54" = 49.5, "55–64" = 59.5,
  "65–74" = 69.5, "75–84" = 79.5, "85+"   = 90
)
age_groups   <- names(age_mid)
sexes        <- c("males", "females")
nhl_subtypes <- c("dlbcl", "follicular", "mantle_cell")             # decomposition tier (subset of NHL)
subtypes     <- c("dlbcl", "follicular", "mantle_cell", "hodgkin")  # HL + the 3 NHL subtypes (per-lymphoma rows)
# Two-tier structure, consistent across incidence and prevalence:
#   Headline tier      - HL and aggregate NHL (disjoint; together = all lymphoma).
#                        Total lymphoma = HL + NHL, built at the DRAW level.
#   Decomposition tier - DLBCL/FL/MCL, a SUBSET of aggregate NHL (~53-73% of it),
#                        NEVER summed into the total (would double-count).
# Aggregate NHL is a HEADLINE result, not a validation-only series; the AIHW
# comparison is retained but reclassified as validation OF a headline result.
prev_series  <- c(subtypes, "nhl")                                   # HL, 3 NHL subtypes, aggregate NHL
agg_groups   <- c("nhl", "hodgkin")                                  # headline tier (2 aggregate APC fits)

period_levels    <- c("2007–2011", "2012–2016", "2017–2021")
period_midpoints <- c("2007–2011" = 2009, "2012–2016" = 2014, "2017–2021" = 2019)
period_map       <- c("2007–2011" = 1L, "2012–2016" = 2L, "2017–2021" = 3L)

# ----------------------
# Numerical edge-handling for logit-normal survival sampling
# ----------------------
SURV_LO    <- 0.001   # shrink survival to (LO, HI) before qlogis (avoids inf)
SURV_HI    <- 0.999
SURV_FLOOR <- 1e-6    # floor for S(3) when computing cond = sqrt(S(5)/S(3))

# ----------------------
# Visual constants (colours, labels, line types)
# ----------------------
# Lymphoma colours. HL purple and aggregate NHL navy; the three NHL subtypes
# form a green/teal ramp. NHL was moved off green (was #1b7837) because it read
# too close to DLBCL. Every series is drawn as a SOLID line of the SAME weight:
# the figures do not try to encode the aggregate/subtype relationship visually
# (that comes from the text), they only distinguish the five lymphomas.
line_colours <- c(
  "HL"    = "#762a83",  # purple
  "NHL"   = "#08306b",  # navy
  "DLBCL" = "#1b7837",  # dark green  - NHL subtype
  "FL"    = "#7fbc41",  # mid green   - NHL subtype
  "MCL"   = "#80cdc1"   # light teal  - NHL subtype
)
line_types <- c(
  "HL" = "solid", "NHL" = "solid",
  "DLBCL" = "solid", "FL" = "solid", "MCL" = "solid"
)
# Uniform line weight for all lymphomas. The tier names are retained only to
# control draw order (aggregate series drawn last, so they sit on top).
lymphoma_lwd <- c(headline = 0.8, decomposition = 0.8)
label_map <- c(
  hodgkin = "HL",  nhl = "NHL",
  dlbcl   = "DLBCL", follicular = "FL", mantle_cell = "MCL"
)
subtype_labels <- c(
  nhl = "NHL (aggregate)", hodgkin = "HL",
  dlbcl = "DLBCL", follicular = "FL", mantle_cell = "MCL"
)
sex_labels  <- c(males = "Males", females = "Females")
sex_colours <- c(males = "#2166ac", females = "#b2182b")
scenario_colours <- c(
  conservative = "#d95f02",
  base         = "#1b9e77",
  optimistic   = "#7570b3"
)
scenario_labels <- c(
  conservative = "Conservative (cap 2021)",
  base         = "Base case (cap 2031)",
  optimistic   = "Optimistic (cap 2036)"
)

# ----------------------
# Draw summary helper
# ----------------------
# Returns a list with `mid` (median), `p025`, `p975`, and `mean`,
# each computed by applying the relevant function across the requested
# margin of `arr`. Used to deduplicate the pattern that otherwise
# appears 4 places across the pipeline.
#
# Examples:
#   summarise_draws(cases_array, c(1, 2))   # 3D -> per-cell summary
#   summarise_draws(prev_draws, 1)          # 2D -> per-row (calc_year)
#   summarise_draws(M, 2)                   # 2D -> per-column (knot value)

summarise_draws <- function(arr, margin, probs = c(0.025, 0.975)) {
  list(
    mid  = apply(arr, margin, median),
    p025 = apply(arr, margin, function(x) unname(quantile(x, probs[1]))),
    p975 = apply(arr, margin, function(x) unname(quantile(x, probs[2]))),
    mean = apply(arr, margin, mean)
  )
}

# ----------------------
# Figure helpers: gtable legend grobs
# ----------------------
# Used by Figure 1 (apc_model.R) and Figures 2 & 3 (prev_model.R) to
# build the custom horizontal legend below the plot.

make_line_grob <- function(colour, lty = "solid", lwd = 2.5) {
  grid::segmentsGrob(x0 = 0.05, x1 = 0.95, y0 = 0.5, y1 = 0.5,
                     gp = grid::gpar(col = colour, lwd = lwd, lty = lty))
}
make_point_grob <- function(colour = "grey30", pch = 16, size = 0.5) {
  grid::pointsGrob(x = 0.5, y = 0.5, pch = pch,
                   size = grid::unit(size, "char"),
                   gp = grid::gpar(col = colour))
}
make_text_grob <- function(txt, fontsize = 12, fontface = "plain") {
  grid::textGrob(txt, x = 0, hjust = 0,
                 gp = grid::gpar(fontsize = fontsize, fontface = fontface))
}
make_spacer <- function() {
  grid::rectGrob(gp = grid::gpar(col = NA, fill = NA))
}

# ----------------------
# Shared lymphoma legend (Figures 1, 2, 3 and SI S3c-S3e)
# ----------------------
# One legend used by every lymphoma figure, so they render identically and the
# panels align. Lists ONLY the five lymphomas with their colour, drawn as solid
# lines of uniform weight. Deliberately has no surrounding box and no entry for
# the AIHW points or the credible-interval band - those are explained in the
# .docx figure captions.
# `attach_lymphoma_legend(p)` returns the plot with the legend beneath it.

.LYMPHOMA_LEGEND_ORDER <- c("HL", "NHL", "DLBCL", "FL", "MCL")

make_lymphoma_legend <- function(labels = .LYMPHOMA_LEGEND_ORDER) {
  seg_w     <- grid::unit(0.7, "cm")
  gap       <- grid::unit(0.1, "cm")
  big_gap   <- grid::unit(0.5, "cm")
  label_gap <- grid::unit(0.25, "cm")
  grobs <- list(); widths <- NULL
  for (i in seq_along(labels)) {
    lab <- labels[i]
    tier <- if (lab %in% c("HL", "NHL")) "headline" else "decomposition"
    if (i > 1) {
      sep <- if (lab == "DLBCL" || labels[i - 1] %in% c("HL", "NHL")) big_gap else label_gap
      grobs  <- c(grobs, list(make_spacer()))
      widths <- grid::unit.c(widths, sep)
    }
    grobs  <- c(grobs, list(
      make_line_grob(line_colours[[lab]], lty = "solid",
                     lwd = unname(lymphoma_lwd[tier]) * 2.6),
      make_spacer(), make_text_grob(lab)))
    w <- grid::unit.c(seg_w, gap, grid::stringWidth(lab))
    widths <- if (is.null(widths)) w else grid::unit.c(widths, w)
  }
  gtable::gtable_row("legend", grobs = grobs, widths = widths,
                     height = grid::unit(0.7, "cm"))
}

attach_lymphoma_legend <- function(p, rel_legend = 0.06) {
  legend_centred <- cowplot::ggdraw() +
    cowplot::draw_grob(make_lymphoma_legend(), x = 0.5, y = 0.5,
                       hjust = 0.5, vjust = 0.5)
  cowplot::plot_grid(p, legend_centred, ncol = 1,
                     rel_heights = c(1, rel_legend))
}

# ----------------------
# Figure helpers: save PDF + PNG together
# ----------------------
# Replaces the repeated `ggsave(...pdf...)` + `ggsave(...png...)` pair.
# Works for plain ggplot objects and for cowplot::plot_grid() outputs.

save_fig <- function(plot, path_stem,
                     width = 10, height = 7,
                     units = "in", dpi = 300) {
  ggplot2::ggsave(paste0(path_stem, ".pdf"), plot,
                  width = width, height = height, units = units)
  ggplot2::ggsave(paste0(path_stem, ".png"), plot,
                  width = width, height = height, units = units, dpi = dpi)
  invisible(c(pdf = paste0(path_stem, ".pdf"),
              png = paste0(path_stem, ".png")))
}
