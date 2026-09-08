# Co-author question (round 3): aggregate NHL rises for males while DLBCL and
# FL look stable over roughly 2005-2020. What is driving it?
#
# One-off diagnostic; writes no output. The answer is recorded in the project
# notes (Lymphoma Projections _notes.md, "male NHL question"). Run from the
# repository root: Rscript code/diagnostics/male_nhl_decomposition.R
#
# Decomposes observed male NHL into the three modelled subtypes and the
# unmodelled residual, on the same 2001-standardised age bands, and separates
# rate change from population ageing.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
})
source("code/_setup.R")

agg  <- read_csv("data/incidence_agg.csv", show_col_types = FALSE)
sub  <- read_csv("data/incidence_subtype.csv", show_col_types = FALSE)
pop  <- read_csv("data/pop_hist.csv", show_col_types = FALSE)

# Standard population exactly as the pipeline builds it: sourcing apc_model.R
# defines its functions without running the pipeline (library + driver pattern).
source("code/apc_model.R")
std <- build_std_pop(pop)

asr <- function(d) {
  d |>
    left_join(pop |> rename(pop = population), by = c("year", "sex", "age_group")) |>
    mutate(rate = count / pop * 1e5,
           w = std[age_group] / sum(std)) |>
    group_by(year) |>
    summarise(asr = sum(rate * w, na.rm = TRUE),
              cases = sum(count, na.rm = TRUE), .groups = "drop")
}

SEX <- "males"
nhl <- agg |> filter(sex == SEX, cancer_group == "Non-Hodgkin lymphoma") |>
  select(year, sex, age_group, count) |> asr() |> rename(nhl_asr = asr, nhl_n = cases)

subs <- sub |> filter(sex == SEX, subtype %in% c("dlbcl", "follicular", "mantle_cell")) |>
  group_by(year, sex, age_group) |> summarise(count = sum(count), .groups = "drop") |>
  asr() |> rename(sub_asr = asr, sub_n = cases)

d <- inner_join(nhl, subs, by = "year") |>
  mutate(res_asr = nhl_asr - sub_asr, res_n = nhl_n - sub_n,
         res_share = res_n / nhl_n)

cat("=== Male NHL: modelled subtypes vs unmodelled residual (ASR per 100,000, 2001 ASP)\n\n")
print(as.data.frame(d |> filter(year %in% c(2003, 2005, 2010, 2015, 2020, 2021)) |>
  mutate(across(where(is.numeric), ~round(.x, 2)))))

f <- function(y) d |> filter(year == y)
a3 <- f(2003); a21 <- f(2021)
cat(sprintf("\n2003 -> 2021 change in ASR:\n"))
cat(sprintf("  aggregate NHL      %+.2f  (%.2f -> %.2f, %+.1f%%)\n",
            a21$nhl_asr-a3$nhl_asr, a3$nhl_asr, a21$nhl_asr,
            100*(a21$nhl_asr/a3$nhl_asr-1)))
cat(sprintf("  3 modelled subtypes%+.2f  (%.2f -> %.2f, %+.1f%%)\n",
            a21$sub_asr-a3$sub_asr, a3$sub_asr, a21$sub_asr,
            100*(a21$sub_asr/a3$sub_asr-1)))
cat(sprintf("  unmodelled residual%+.2f  (%.2f -> %.2f, %+.1f%%)\n",
            a21$res_asr-a3$res_asr, a3$res_asr, a21$res_asr,
            100*(a21$res_asr/a3$res_asr-1)))
tot <- a21$nhl_asr - a3$nhl_asr
cat(sprintf("\n  share of the aggregate rise from the residual: %.0f%%\n",
            100*(a21$res_asr-a3$res_asr)/tot))
cat(sprintf("  residual as a share of male NHL cases: %.0f%% (2003) -> %.0f%% (2021)\n",
            100*a3$res_share, 100*a21$res_share))

# per-subtype ASR change, to test "DLBCL and FL are stable"
cat("\n=== Each modelled subtype, male ASR\n")
for (st in c("dlbcl", "follicular", "mantle_cell")) {
  x <- sub |> filter(sex == SEX, subtype == st) |>
    select(year, sex, age_group, count) |> asr()
  x3 <- x$asr[x$year == 2003]; x21 <- x$asr[x$year == 2021]
  cat(sprintf("  %-12s %.2f -> %.2f  (%+.2f, %+.1f%%)\n", st, x3, x21, x21-x3,
              100*(x21/x3-1)))
}

# counts vs rates: how much of the CASE growth is ageing rather than rate?
cat("\n=== Cases vs age-standardised rate, male aggregate NHL\n")
cat(sprintf("  cases  %d -> %d  (%+.1f%%)\n", a3$nhl_n, a21$nhl_n,
            100*(a21$nhl_n/a3$nhl_n-1)))
cat(sprintf("  ASR    %.2f -> %.2f  (%+.1f%%)\n", a3$nhl_asr, a21$nhl_asr,
            100*(a21$nhl_asr/a3$nhl_asr-1)))
cat("  -> the gap between the two is population growth and ageing.\n")

# the 2005-2020 window the co-author actually named
b5 <- f(2005); b20 <- f(2020)
cat("\n=== The 2005-2020 window the co-author named\n")
cat(sprintf("  aggregate NHL ASR  %.2f -> %.2f (%+.1f%%)\n", b5$nhl_asr, b20$nhl_asr,
            100*(b20$nhl_asr/b5$nhl_asr-1)))
cat(sprintf("  3 subtypes ASR     %.2f -> %.2f (%+.1f%%)\n", b5$sub_asr, b20$sub_asr,
            100*(b20$sub_asr/b5$sub_asr-1)))
cat(sprintf("  residual ASR       %.2f -> %.2f (%+.1f%%)\n", b5$res_asr, b20$res_asr,
            100*(b20$res_asr/b5$res_asr-1)))
