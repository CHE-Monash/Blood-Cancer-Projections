# ================
# table_1_change_cri.R - 2021 -> 2045 incidence change with draw-level CrIs
# ================
# Builds output/table_1_change_cri.csv: per (lymphoma, sex) plus the TOTAL
# (HL + aggregate NHL) by sex, the percentage change 2021 -> 2045 in annual
# cases and in the ASR, with 95% CrIs. Change draws divide the projected
# 2045 draw by the OBSERVED 2021 value (cases and ASR are observed in 2021,
# so the change CrI inherits its width from the 2045 projection). Uses the
# stored draws in output/apc_results.rds; no refitting, no RNG.
# Added round 3 [QL: CrI consistency across Tables 1-3].
#
#   Rscript code/table_1_change_cri.R
# ================

source("code/apc_model.R")   # functions only; auto-run guard not triggered
suppressPackageStartupMessages({ library(dplyr); library(readr) })

run_table1_change_cri <- function(save_dir = "output") {
  apc_results <- readRDS(file.path(save_dir, "apc_results.rds"))
  inc_agg  <- read_csv("data/incidence_agg.csv",  show_col_types = FALSE)
  pop_hist <- read_csv("data/pop_hist.csv",       show_col_types = FALSE)
  pop_proj <- read_csv("data/pop_proj.csv",       show_col_types = FALSE)
  std_pop  <- build_std_pop(pop_hist)
  pop_mat_by_sex <- list(
    males   = build_pop_mat("males",   pop_proj),
    females = build_pop_mat("females", pop_proj)
  )
  q <- function(x, p) unname(quantile(x, p))

  draws_2045 <- function(cases_age_b, pop_mat) {
    std_vec <- std_pop[dimnames(cases_age_b)[[2]]]
    sum_std <- sum(std_pop)
    yr <- match("2045", dimnames(cases_age_b)[[1]])
    cases_b <- apply(cases_age_b[yr, , , drop = FALSE], 3, sum)
    pop_row <- pop_mat[match("2045", rownames(pop_mat)), dimnames(cases_age_b)[[2]]]
    rates_b <- sweep(cases_age_b[yr, , ], 1, pop_row, "/") * 1e5   # (age x B)
    asr_b   <- apply(rates_b, 2, function(r) sum(r * std_vec) / sum_std)
    list(cases = cases_b, asr = asr_b)
  }
  obs_2021 <- function(grp_label, sx, source_df, label_col) {
    df <- source_df |> filter(.data[[label_col]] == grp_label, sex == sx, year == 2021)
    cases <- sum(df$count)
    asr <- calc_hist_asr(df |> select(year, age_group, count), sx, pop_hist, std_pop) |>
      filter(year == 2021)
    list(cases = cases, asr = asr$asr_obs)
  }

  rows <- list()
  add_row <- function(label, sx, cases_age_b, pop_mat, obs) {
    d <- draws_2045(cases_age_b, pop_mat)
    pc_cases <- (d$cases / obs$cases - 1) * 100
    pc_asr   <- (d$asr   / obs$asr   - 1) * 100
    rows[[length(rows) + 1]] <<- tibble::tibble(
      subtype = label, sex = sx,
      cases_2021 = obs$cases, asr_2021 = obs$asr,
      cases_2045 = median(d$cases), cases_2045_p025 = q(d$cases, .025), cases_2045_p975 = q(d$cases, .975),
      asr_2045 = median(d$asr), asr_2045_p025 = q(d$asr, .025), asr_2045_p975 = q(d$asr, .975),
      pct_change_cases = median(pc_cases),
      pct_change_cases_p025 = q(pc_cases, .025), pct_change_cases_p975 = q(pc_cases, .975),
      pct_change_asr = median(pc_asr),
      pct_change_asr_p025 = q(pc_asr, .025), pct_change_asr_p975 = q(pc_asr, .975)
    )
  }

  inc_subtype <- read_csv("data/incidence_subtype.csv", show_col_types = FALSE)
  for (key in names(apc_results)) {
    e <- apc_results[[key]]
    obs <- if (e$tier == "agg") {
      obs_2021(switch(e$subtype, nhl = "Non-Hodgkin lymphoma",
                                 hodgkin = "Hodgkin lymphoma"), e$sex,
               inc_agg |> rename(label = cancer_group), "label")
    } else {
      obs_2021(e$subtype, e$sex, inc_subtype |> rename(label = subtype), "label")
    }
    add_row(e$subtype, e$sex, e$cases_draws, pop_mat_by_sex[[e$sex]], obs)
  }
  cd <- function(grp, sx) apc_results[[paste(grp, sx, sep = "_")]]$cases_draws
  for (sx in sexes) {
    tot_b <- cd("hodgkin", sx) + cd("nhl", sx)
    df21 <- inc_agg |> dplyr::filter(sex == !!sx, year == 2021) |>
      group_by(year, age_group) |> summarise(count = sum(count), .groups = "drop")
    obs <- list(
      cases = sum(df21$count),
      asr = calc_hist_asr(df21, sx, pop_hist, std_pop) |> filter(year == 2021) |> pull(asr_obs)
    )
    add_row("total", sx, tot_b, pop_mat_by_sex[[sx]], obs)
  }
  out <- bind_rows(rows) |>
    mutate(across(c(cases_2045, cases_2045_p025, cases_2045_p975), ~ round(.x)),
           across(starts_with("asr_"),  ~ round(.x, 2)),
           across(starts_with("pct_"),  ~ round(.x, 1)))
  write_csv(out, file.path(save_dir, "table_1_change_cri.csv"))
  invisible(out)
}

if (!interactive() && sys.nframe() == 0) {
  invisible(run_table1_change_cri())
}
