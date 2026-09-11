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
# agg_start selected by the start-year sensitivity (code/supporting/start_year_sensitivity.R;
# _notes.md "Start-year selection"). Observed incidence runs from 1982, but fitting
# the aggregate tier from 1990 markedly improves out-of-sample backtest error (NHL 2021
# holdout endpoint bias falls from ~+13% to ~+2%) by excluding the transient pre-1990
# NHL surge (HIV/AIDS-associated NHL + pre-WHO reclassification), which the APC model
# would otherwise extrapolate. Holdout error plateaus for any start >= 1990, so 1990
# is the earliest start selected by this temporal backtesting — retaining the longest
# series (32 yr). The same 2012-2021 window scored the knot grid, so this is a
# selection criterion, not an independent test (R4). Both
# prep_agg_data() and prep_subtype_data() filter P >= these starts.

# APC natural-cubic-spline knot counts (age, period, cohort), selected by
# out-of-sample holdout error + parsimony (see code/supporting/knot_selection.R and
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
# Age coordinates are plain interval midpoints, as for every band. The
# youngest group spans 0-14, so its midpoint is 7.0 (R1,
# research-review-2026-09-09.md). This is a convention, not an empirically
# preferred representative age: the case-weighted mean age within 0-14 over
# 1982-2021 is 10.47 and 11.50 for male and female HL and 8.15 and 8.60 for
# male and female NHL, and those values, together with the pre-repair 9.5,
# are run as recorded sensitivities. The subtype tier has no five-year
# source detail, so a case-weighted coordinate cannot be derived for it and
# an aggregate coordinate must never be transferred to a subtype.
age_mid <- c(
  "0–14"  = 7.0,  "15–24" = 19.5, "25–34" = 29.5,
  "35–44" = 39.5, "45–54" = 49.5, "55–64" = 59.5,
  "65–74" = 69.5, "75–84" = 79.5, "85+"   = 90
)
# Alternative youngest-band coordinates for the R1 sensitivity. Set
# age_mid["0–14"] to one of these before sourcing apc_model.R.
age_mid_youngest_alt <- c(
  midpoint_7    = 7.0,     # adopted; interval midpoint, consistent with all other bands
  pre_repair_9_5 = 9.5,    # the coordinate used before the R1 repair
  cw_hodgkin_males   = 10.47, cw_hodgkin_females   = 11.50,
  cw_nhl_males       = 8.15,  cw_nhl_females       = 8.60
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
# Panel-strip labels for the SI figures. "(aggregate)" was dropped from NHL on
# 2026-09-07 (QL, round 1: misleading) so every panel strip in the manuscript
# and the SI reads the same as the tables.
subtype_labels <- c(
  nhl = "NHL", hodgkin = "HL",
  dlbcl = "DLBCL", follicular = "FL", mantle_cell = "MCL"
)
sex_labels  <- c(males = "Males", females = "Females")
sex_colours <- c(males = "#2166ac", females = "#b2182b")
# Keyed by the display labels used on the re-panelled figures, where the facet
# strip carries the lymphoma and colour carries the sex (see the figure-key
# section below).
sex_fig_colours <- c(Females = unname(sex_colours[["females"]]),
                     Males   = unname(sex_colours[["males"]]))
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
# Shared canvas for the re-panelled lymphoma figures
# ----------------------
# Portrait 3x2 (five lymphoma panels plus the key in the sixth cell). Fits
# the text width of both .docx files without a landscape section break, and
# matches the canvas S1-S3 already use.
fig_panel_w <- 7.5   # inches
fig_panel_h <- 9.5   # inches

# ----------------------
# Shared figure key (Figures 1-3 and SI S1-S3, S7-S11)
# ----------------------
# Every lymphoma figure panels by LYMPHOMA and separates the sexes by colour
# (round 3, QL: the per-cancer view of each trend matters more in the figures
# than a cross-cancer comparison of level, which the tables carry). Five
# lymphomas in a two-column facet leave an empty sixth cell, and that cell
# carries the whole key - sexes, credible band, observed and AIHW markers, the
# projection-start rule - so the .docx captions no longer have to define them.
#
# The key is built from grid grobs and placed directly into the empty panel cell
# of the ggplot gtable, so its position is derived from the layout rather than
# guessed in npc, and nothing depends on ggplot2's own legend placement (which
# changed to legend.position.inside at 3.5.0).
#
# `key_line()` / `key_band()` / `key_point()` build entries;
# `attach_figure_key(p, entries)` returns the plot with the key in place.
#
# Replaces make_lymphoma_legend() / attach_lymphoma_legend(), retired 2026-09-07
# with the re-panel: colour now encodes sex, and the lymphoma is named by the
# facet strip, so a five-lymphoma legend no longer has anything to say.

key_line <- function(label, colour, lty = "solid", lwd = 2.4) {
  list(label = label, grob = make_line_grob(colour, lty = lty, lwd = lwd))
}
key_point <- function(label, colour = "grey25", pch = 16, size = 0.6) {
  list(label = label, grob = make_point_grob(colour = colour, pch = pch,
                                             size = size))
}
key_band <- function(label, colour = "grey25", alpha = 0.15, lwd = 2.4) {
  list(label = label,
       grob  = grid::grobTree(
         grid::rectGrob(width = 0.9, height = 0.55,
                        gp = grid::gpar(col = NA, fill = colour, alpha = alpha)),
         make_line_grob(colour, lwd = lwd)))
}

# Keep labels short: the key is centred in one panel cell (about half the
# figure width), and clip = "off" means an over-long label runs off the
# canvas rather than wrapping. About 30 characters is the practical limit at
# fontsize 12 on the 7.5 in canvas.
make_figure_key <- function(entries, fontsize = 12, row_height = 0.62) {
  labs  <- vapply(entries, `[[`, character(1), "label")
  lab_w <- max(do.call(grid::unit.c, lapply(labs, grid::stringWidth)))
  n     <- length(entries)
  # The outer null rows/columns centre the key within the empty panel cell.
  g <- gtable::gtable(
    widths  = grid::unit.c(grid::unit(1, "null"), grid::unit(0.9, "cm"),
                           grid::unit(0.3, "cm"), lab_w, grid::unit(1, "null")),
    heights = grid::unit.c(grid::unit(1, "null"),
                           grid::unit(rep(row_height, n), "cm"),
                           grid::unit(1, "null")))
  for (i in seq_len(n)) {
    g <- gtable::gtable_add_grob(g, entries[[i]]$grob, t = i + 1L, l = 2L,
                                 name = paste0("key-sym-", i))
    g <- gtable::gtable_add_grob(g, make_text_grob(labs[i], fontsize = fontsize),
                                 t = i + 1L, l = 4L,
                                 name = paste0("key-lab-", i))
  }
  g
}

attach_figure_key <- function(p, entries) {
  g   <- ggplot2::ggplotGrob(p)
  idx <- which(grepl("^panel", g$layout$name))
  # An odd facet count leaves one panel cell empty. ggplot2 keeps that cell in
  # the layout with a zeroGrob in it (3.4.x); if a future version drops the cell
  # instead, fall back to the bottom-right position.
  empty <- idx[vapply(g$grobs[idx], inherits, logical(1), "zeroGrob")]
  if (length(empty) >= 1L) {
    t_key <- g$layout$t[empty[1]]
    l_key <- g$layout$l[empty[1]]
  } else {
    t_key <- max(g$layout$t[idx])
    l_key <- max(g$layout$l[idx])
    if (any(g$layout$t[idx] == t_key & g$layout$l[idx] == l_key)) {
      warning("attach_figure_key(): no empty panel cell, key not drawn")
      out <- cowplot::ggdraw(g)
      # The composed object is a new ggplot whose $data is empty; keep the
      # plot that was keyed reachable so tests can inspect what was drawn
      # (code/_tests.R Test S4). An attribute does not affect rendering.
      attr(out, "source_plot") <- p
      return(out)
    }
  }
  g <- gtable::gtable_add_grob(g, make_figure_key(entries),
                               t = t_key, l = l_key,
                               name = "figure-key", clip = "off")
  # An explicit white background. Figures 1-3 previously saved with a
  # transparent one (plot.background fill = NA through cowplot::ggdraw), which
  # Word hides but a journal production system may composite onto black; the SI
  # effect figures saved white, so the set was inconsistent. Set here so every
  # figure that carries the key is white by construction.
  out <- cowplot::ggdraw(g) +
    ggplot2::theme(plot.background =
                     ggplot2::element_rect(fill = "white", colour = NA))
  # See the early return above: expose the keyed plot for the test harness.
  attr(out, "source_plot") <- p
  out
}

# Standard key entries. `sex_key()` is the two-line sex key every re-panelled
# figure opens with; the rest are appended per figure as that figure needs them.
sex_key <- function() {
  list(key_line("Females", sex_colours[["females"]]),
       key_line("Males",   sex_colours[["males"]]))
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
