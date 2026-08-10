# ================
# Blood Cancer Projections (lymphoma): Import Data
# ================
# Transforms raw AIHW (CDiA) and ABS workbooks in data/raw/ into the
# processed CSVs in data/.
#
# NOT RUNNABLE FOR EXTERNAL USERS: data/raw/ is excluded from this
# repository via .gitignore. To reproduce the analysis from the
# included CSVs, start at apc_model.R. This script is provided for
# transparency, documenting how the CSVs in data/ were derived from
# the public AIHW and ABS source workbooks (see README.md for
# download links). Run from the repository root.
# ================
library(tidyverse)
library(readxl)

aihw_path <- "data/raw/AIHW"

# Target subtypes (match start of I.D. string)
target_subtypes <- c(
  "B2.02.01.07",  # Follicular lymphoma
  "B2.02.01.08",  # Mantle cell lymphoma
  "B2.02.01.09",  # DLBCL
  "B2.02.03"      # Hodgkin lymphoma
)

subtype_labels <- c(
  "B2.02.01.07" = "follicular",
  "B2.02.01.08" = "mantle_cell",
  "B2.02.01.09" = "dlbcl",
  "B2.02.03"    = "hodgkin"
)

# Helper function to match subtypes
match_subtype <- function(id_string) {
  code <- str_extract(id_string, "^[A-Z0-9.]+")
  if (!is.na(code) && code %in% target_subtypes) return(code)
  return(NA_character_)
}

# ==========
# Aggregate incidence (1982-2021)
# ==========
# Source: AIHW Cancer Data in Australia, Table S1a.1
# Incidence by cancer group, sex, and 5-year age group
# Three uses:
#   1. NHL series: back-estimation of DLBCL, FL, MCL incidence for
#      1982–2002 using proportional allocation from 2003–2007 shares.
#   2. HL series: used directly as observed HL incidence for 1982–2002
#      (HL maps 1:1 to the aggregate category, no allocation needed).
#   3. Aggregate NHL + HL APC models for AIHW incidence validation.

incidence_agg_raw <- read_excel(
  file.path(aihw_path, "CDiA-2025-Book-1a-Cancer-incidence-age-standardised-rates-5-year-age-groups.xlsx"),
  sheet = "Table S1a.1",
  skip = 5
)

# Check column names
names(incidence_agg_raw)

# Check cancer group values to find exact strings
incidence_agg_raw |>
  count(`Cancer group/site`) |>
  filter(str_detect(`Cancer group/site`, regex("hodgkin|lymphoma", ignore_case = TRUE))) |>
  print()

incidence_agg <- incidence_agg_raw |>
  rename(
    data_type    = `Data type`,
    cancer_group = `Cancer group/site`,
    year         = Year,
    sex          = Sex,
    age_group    = `Age group (years)`,
    count        = Count
  ) |>
  filter(cancer_group %in% c("Non-Hodgkin lymphoma", "Hodgkin lymphoma")) |>
  filter(!is.na(year)) |>
  mutate(
    sex = tolower(sex),
    count = as.numeric(count)
  ) |>
  filter(sex != "persons") |>
  filter(!str_starts(age_group, "All ages")) |>
  select(year, sex, age_group, cancer_group, count)

# Check
incidence_agg |> count(cancer_group)
incidence_agg |> count(year) |> print(n = 50)
incidence_agg |> count(age_group) |> print(n = 30)
incidence_agg |> count(sex)

# Aggregate into 10-year age bands to match subtype data
incidence_agg <- incidence_agg |>
  mutate(
    age_num = as.integer(str_extract(age_group, "^[0-9]+")),
    age_band = case_when(
      age_num >= 0  & age_num < 15  ~ "5–14",
      age_num >= 15 & age_num < 25  ~ "15–24",
      age_num >= 25 & age_num < 35  ~ "25–34",
      age_num >= 35 & age_num < 45  ~ "35–44",
      age_num >= 45 & age_num < 55  ~ "45–54",
      age_num >= 55 & age_num < 65  ~ "55–64",
      age_num >= 65 & age_num < 75  ~ "65–74",
      age_num >= 75 & age_num < 85  ~ "75–84",
      age_num >= 85                 ~ "85+",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(age_band)) |>
  group_by(year, sex, age_band, cancer_group) |>
  summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
  rename(age_group = age_band)

# Check
incidence_agg |> count(age_group)
incidence_agg |> filter(year == 2021, cancer_group == "Non-Hodgkin lymphoma") |> print(n = 20)

# Save
write_csv(incidence_agg, "data/incidence_agg.csv")

# ==========
# Subtype incidence (2003–2021)
# ==========
# Source: AIHW Cancer Data in Australia, Table S11a.1
# Subtype-specific incidence by sex and 10-year age group.
# Primary input for APC model fitting (DLBCL, FL, MCL) and observed
# incidence in the prevalence counting method. Also used to compute
# subtype proportions within NHL for back-estimation to 1982.
# HL is available here from 2003 but can also be sourced from aggregate
# tables back to 1982 (see below).

incidence_subtype_raw <- read_excel(
  file.path(aihw_path, "CDiA-2025-Book-11a-Blood-cancer-incidence-by-histology-main-framework.xlsx"),
  sheet = "Table S11a.1",
  skip = 4
)

incidence_subtype <- incidence_subtype_raw |>
  rename(
    year = Year,
    sex = Sex,
    age_group = `Age group (years)`,
    cancer_type = `Cancer type`,
    count = Count, 
    rate = `Rate (per 100,000)`,
    id_full = `I.D. and cancer type`
  ) |>
  select(year, sex, age_group, cancer_type, count, rate, id_full) |>
  filter(!is.na(id_full)) |>
  mutate(
    subtype_code = map_chr(id_full, match_subtype),
    sex = tolower(sex),
    age_group = if_else(age_group == "0–14", "5–14", age_group)
  ) |>
  filter(!is.na(subtype_code)) |>
  filter(!str_starts(age_group, "All ages")) |>
  filter(sex != "persons") |>
  mutate(subtype = subtype_labels[subtype_code]) |>
  select(year, sex, age_group, subtype, count, rate)

# Check
incidence_subtype |> count(subtype)
incidence_subtype |> count(age_group)
incidence_subtype |> count(year)

# Save 
write_csv(incidence_subtype, "data/incidence_subtype.csv")

# ==========
# Aggregate prevalence (2000–2021)
# ==========
# Source: AIHW Cancer Data in Australia, Table S6.2
# Used for validating 30-year prevalence estimates (ceiling check)

prevalence_agg_raw <- read_excel(
  file.path(aihw_path, "CDiA-2025-Book-6-Cancer-prevalence.xlsx"),
  sheet = "Table S6.2",
  skip = 5
)

names(prevalence_agg_raw)

# Check cancer group strings
prevalence_agg_raw |>
  count(`Cancer group/site`) |>
  filter(str_detect(`Cancer group/site`, regex("hodgkin|lymphoma", ignore_case = TRUE))) |>
  print()

prevalence_agg <- prevalence_agg_raw |>
  rename(
    year         = `Year (as at 31 Dec)`,
    cancer_group = `Cancer group/site`,
    duration     = `Duration (years)`,
    sex          = Sex,
    age_group    = `Age group (years)`,
    prevalence   = Prevalence
  ) |>
  filter(cancer_group %in% c("Non-Hodgkin lymphoma", "Hodgkin lymphoma")) |>
  filter(!is.na(year)) |>
  mutate(
    sex        = tolower(sex),
    prevalence = as.numeric(prevalence),
    duration   = as.numeric(duration)
  ) |>
  filter(sex != "persons") |>
  filter(!str_starts(age_group, "All ages")) |>
  select(year, cancer_group, duration, sex, age_group, prevalence)

# Check
prevalence_agg |> count(cancer_group)
prevalence_agg |> count(duration)
prevalence_agg |> count(year) |> print(n = 30)
prevalence_agg |> count(age_group) |> print(n = 20)
prevalence_agg |> count(sex)

# Aggregate 5-year age bands into 10-year to match our data
prevalence_agg <- prevalence_agg |>
  mutate(
    age_num = as.integer(str_extract(age_group, "^[0-9]+")),
    age_band = case_when(
      age_num >= 0  & age_num < 15  ~ "5–14",
      age_num >= 15 & age_num < 25  ~ "15–24",
      age_num >= 25 & age_num < 35  ~ "25–34",
      age_num >= 35 & age_num < 45  ~ "35–44",
      age_num >= 45 & age_num < 55  ~ "45–54",
      age_num >= 55 & age_num < 65  ~ "55–64",
      age_num >= 65 & age_num < 75  ~ "65–74",
      age_num >= 75 & age_num < 85  ~ "75–84",
      age_num >= 85                 ~ "85+",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(age_band)) |>
  group_by(year, cancer_group, duration, sex, age_band) |>
  summarise(prevalence = sum(prevalence, na.rm = TRUE), .groups = "drop") |>
  rename(age_group = age_band)

# Check
prevalence_agg |> count(age_group)

# Save
write_csv(prevalence_agg, "data/prevalence_agg.csv")

# ==========
# Subtype prevalence (2012–2021)
# ==========
# Source: AIHW Cancer Data in Australia, Table S11i.1
# Subtype-specific prevalence by sex, age group, and duration
# Used as the validation benchmark for computed prevalence from the
# modified counting method (diagnostics checks 4 and 5).

prevalence_subtype_raw <- read_excel (
  file.path(aihw_path, "CDiA-2025-Book-11i-Blood-cancer-prevalence-main-framework.xlsx"),
  sheet = "Table S11i.1",
  skip = 4
)

prevalence_subtype = prevalence_subtype_raw |>
  rename(
    year = `Year (as at 31 Dec)`,
    duration = `Duration (years)`,
    sex = Sex,
    age_group = `Age group (years)`,
    cancer_type = `Cancer type`,
    prevalence = Prevalence,
    rate = `Rate\r\n(per 100,000)`,
    id_full = `I.D. and cancer type`
  ) |>
  select(year, duration, sex, age_group, cancer_type, prevalence, rate, id_full) |>
  filter(!is.na(id_full)) |>
  mutate(
    subtype_code = map_chr(id_full, match_subtype),
    sex = tolower(sex),
    age_group = if_else(age_group == "0–14", "5–14", age_group)
  ) |>
  filter(!is.na(subtype_code)) |>
  filter(!str_starts(age_group, "All ages")) |>
  filter(sex != "persons") |>
  mutate(subtype = subtype_labels[subtype_code]) |>
  select(year, duration, sex, age_group, subtype, prevalence, rate)
  
# Check
prevalence_subtype |> count(subtype)
prevalence_subtype |> count(duration)
prevalence_subtype |> count(year)

# Save
write_csv(prevalence_subtype, "data/prevalence_subtype.csv")

# ==========
# Subtype survival (2007–2021)
# ==========
# Source: AIHW Cancer Data in Australia, Table S11e.1
# Observed survival by subtype, sex, age group, period, and years since
# diagnosis. Used in the prevalence counting method: observed survival to
# 5 years post-diagnosis, with annual improvement rates extrapolated to
# project survival beyond the last observed period. Suppressed (n.p.)
# cells are filled from the nearest available age group.
#
# Covers all 4 prevalence subtypes: DLBCL, follicular, mantle cell, Hodgkin.
# Hodgkin maps 1:1 between the blood cancer framework and ICD-10, so Book 11e
# provides the same data that would come from an aggregate-category survival
# table, but with the advantage of all three periods for improvement modelling.
#
# Aggregate NHL survival (Book 3a) is NOT imported — NHL is excluded from
# prevalence projections as it is a classification umbrella, not a treatment
# entity. Aggregate NHL incidence projections are retained in apc_model.R for
# AIHW validation only.

survival_raw <- read_excel(
  file.path(aihw_path, "CDiA-2025-Book-11e-Blood-cancer-survival-and-age-group-proportions-by-histology-main-framework.xlsx"),
  sheet = "Table S11e.1",
  skip = 4
)

survival <- survival_raw |>
  rename(
    survival_type = `Survival type`,
    period = `Period (years)`,
    sex = Sex,
    age_group = `Age group (years)`,
    years_since_dx = `Years after diagnosis`,
    cancer_type = `Cancer type`,
    survival_pct = `Survival (%)`,
    ci_lower = `95% CI lower bound (%)`,
    ci_upper = `95% CI upper bound (%)`,
    id_full = `I.D. and cancer type`
  ) |>
  select(survival_type, period, sex, age_group, years_since_dx,
         cancer_type, survival_pct, ci_lower, ci_upper, id_full) |>
  filter(!is.na(id_full)) |>
  mutate(
    subtype_code = map_chr(id_full, match_subtype),
    sex = tolower(sex),
    age_group = if_else(age_group == "0–14", "5–14", age_group),
    # Convert n.p. → NA and coerce to numeric
    survival_pct = as.numeric(survival_pct),
    ci_lower     = as.numeric(ci_lower),
    ci_upper     = as.numeric(ci_upper)
  ) |>
  filter(!is.na(subtype_code)) |>
  filter(!str_starts(age_group, "All ages")) |>
  mutate(subtype = subtype_labels[subtype_code]) |>
  select(survival_type, period, sex, age_group, years_since_dx,
         subtype, survival_pct, ci_lower, ci_upper)

# ---- Aggregate NHL survival (ICD-10 C82–C86) from Book 11f1 ----
# Book 11e (histology framework) has no aggregate NHL row, so aggregate NHL
# survival is taken from Book 11f1 (ICD-10 framework), which carries
# "Non-Hodgkin lymphoma" as a group with the same age/sex/period/ysd
# structure. Needed to compute and validate aggregate-NHL prevalence (QL #3;
# see prev_model.R prev_series). Same cleaning as the subtype survival above.

nhl_surv_raw <- bind_rows(
  read_excel(file.path(aihw_path, "CDiA-2025-Book-11f1-Part1-Blood-cancer-survival-by-histology-ICD10-framework-10-year-age-groups.xlsx"),
             sheet = "Table S11f.1", skip = 4),
  read_excel(file.path(aihw_path, "CDiA-2025-Book-11f1-Part2-Blood-cancer-survival-by-histology-ICD10-framework-10-year-age-groups.xlsx"),
             sheet = "Table S11f.1", skip = 4)
)
names(nhl_surv_raw)[1:10] <- c("survival_type", "period", "sex", "age_group",
                               "years_since_dx", "cancer_group", "cancer_type",
                               "survival_pct", "ci_lower", "ci_upper")

nhl_survival <- nhl_surv_raw |>
  filter(cancer_type == "Non-Hodgkin lymphoma") |>
  transmute(
    survival_type, period, sex = tolower(sex),
    age_group      = if_else(age_group == "0–14", "5–14", age_group),
    years_since_dx = as.integer(years_since_dx),
    subtype        = "nhl",
    survival_pct   = as.numeric(survival_pct),
    ci_lower       = as.numeric(ci_lower),
    ci_upper       = as.numeric(ci_upper)
  ) |>
  filter(sex != "persons", !str_starts(age_group, "All ages")) |>
  select(survival_type, period, sex, age_group, years_since_dx,
         subtype, survival_pct, ci_lower, ci_upper)

survival <- bind_rows(survival, nhl_survival)

# Check
survival |> count(subtype)
survival |> count(survival_type)
survival |> count(period)

# Save
write_csv(survival, "data/survival.csv")

# ==========
# Aggregate incidence projections (2026–2035)
# ==========
# Source: AIHW Cancer Data in Australia, Table S1e.1
# Used for plausibility validation of the incidence projections
# - HL: direct comparison (counts, ASR)
# - NHL: subtype sum (DLBCL + FL + MCL) must be < NHL total

agg_proj_path <- "data/raw/AIHW/CDiA-2025-Book-1e-Long-term-cancer-incidence-projections.xlsx"

incidence_agg_proj_raw <- read_excel(
  agg_proj_path,
  sheet = "Table S1e.1",
  skip = 4
)

names(incidence_agg_proj_raw)

# Check available cancer groups
incidence_agg_proj_raw |>
  count(`Cancer group/site`) |>
  filter(str_detect(`Cancer group/site`, 
                    regex("hodgkin|lymphoma", ignore_case = TRUE))) |>
  print()

incidence_agg_proj <- incidence_agg_proj_raw |>
  rename(
    cancer_group = `Cancer group/site`,
    year         = Year,
    sex          = Sex,
    count        = Count,
    crude_rate   = `Crude rate\r\n(per 100,000)`,
    asr_2001     = `Age-standardised rate\r\n2001 Australian Standard Population \r\n(per 100,000)`,
    asr_2025     = `Age-standardised rate\r\n2025 Australian population \r\n(per 100,000)`
  ) |>
  filter(cancer_group %in% c("Non-Hodgkin lymphoma", "Hodgkin lymphoma")) |>
  filter(!is.na(year)) |>
  mutate(
    sex   = tolower(sex),
    count = as.numeric(count),
    asr_2001 = as.numeric(asr_2001),
    asr_2025 = as.numeric(asr_2025)
  ) |>
  filter(sex != "persons") |>
  select(year, sex, cancer_group, count, crude_rate, asr_2001, asr_2025)

# Check
incidence_agg_proj |> count(cancer_group)
incidence_agg_proj |> count(year)
incidence_agg_proj |> count(sex)
incidence_agg_proj |> print(n = 20)

# Save
write_csv(incidence_agg_proj, "data/incidence_agg_proj.csv")

# ==========
# Population Projections (2022–2050)
# ==========
# Source: ABS 3222.0, Table B9 (Series B, medium growth)
# Projected population by sex and single year of age,
# aggregated to 10-year bands.
# Two uses:
#   1. Denominator (person-years) for converting projected age-specific
#      rates to case counts in apc_model.R.
#   2. Weights for computing projected ASR (via the 2001 ASP).

abs_proj_path <- "data/raw/ABS/3222_Table_B9.xlsx"

# Get column names from first read
col_names <- names(read_excel(abs_proj_path, sheet = "Data1", skip = 0))

# Read data rows, apply those names
pop_raw <- read_excel(abs_proj_path, sheet = "Data1", skip = 9, col_names = FALSE)
names(pop_raw) <- col_names

# Clean up
pop_raw <- pop_raw |>
  rename(date = 1) |>
  mutate(
    date = as.Date(as.numeric(date), origin = "1899-12-30"),
    year = year(date)
  ) |>
  filter(year >= 2022, year <= 2050)

# Pivot to long format
pop_long <- pop_raw |>
  select(-date) |>
  pivot_longer(
    cols = -year,
    names_to = "description",
    values_to = "population"
  ) |>
  mutate(
    population = as.numeric(population),
    sex = case_when(
      str_detect(description, "Male") ~ "males",
      str_detect(description, "Female") ~ "females",
      TRUE ~ NA_character_
    ),
    # Extract age
    age = str_extract(description, "(?<=;  )[0-9]+(?= ;)"),
    age = if_else(str_detect(description, "100 and over"), "100", age),
    age = as.integer(age)
  ) |>
  # Keep only male/female, drop "Persons"
  filter(!is.na(sex)) |>
  select(year, sex, age, population)

# Check
pop_long |> count(sex)
pop_long |> count(age) |> print(n = 30)

# Aggregate into 10-year age bands
# Ages 0–4 are retained as a separate band: lymphoma incidence in 0–4 is
# essentially zero and the AIHW data does not provide it, so 0–4 is not
# used as a numerator in any model. The 0–4 population is needed in the
# 2001 Australian Standard Population denominator so that ASRs match the
# AIHW convention (which uses the full 2001 ASP).
pop_proj <- pop_long |>
  mutate(
    age_group = case_when(
      age <= 4              ~ "0–4",
      age >= 5 & age <= 14 ~ "5–14",
      age >= 15 & age <= 24 ~ "15–24",
      age <= 34 ~ "25–34",
      age <= 44 ~ "35–44",
      age <= 54 ~ "45–54",
      age <= 64 ~ "55–64",
      age <= 74 ~ "65–74",
      age <= 84 ~ "75–84",
      age >= 85 ~ "85+",
      TRUE ~ NA_character_
    )
  ) |>
  group_by(year, sex, age_group) |>
  summarise(population = sum(population), .groups = "drop")

# Check
pop_proj |> count(age_group)
pop_proj |> filter(year == 2025)

# Save
write_csv(pop_proj, "data/pop_proj.csv")

# ==========
# Historical Population (1982–2021)
# ==========
# Source: ABS 3101.0, Table 59
# Estimated resident population by sex and single year of age,
# aggregated to 10-year bands. The series must start in 1982 so the
# aggregate-tier APC models (NHL, HL) fit the full 1982–2021 incidence
# window; the raw workbook covers 1971–2025 (single-year ERP at 30 June).
# Two uses:
#   1. Denominator (person-years) for APC model fitting in apc_model.R.
#   2. Source of 2001 Australian Standard Population (ASP) weights for
#      age-standardised rates (ERP at 30 June 2001, both sexes).

abs_hist_path <- "data/raw/ABS/3101059.xlsx"

col_names_hist <- names(read_excel(abs_hist_path, sheet = "Data1", skip = 0))

pop_hist_raw <- read_excel(abs_hist_path, sheet = "Data1", skip = 9, col_names = FALSE)
names(pop_hist_raw) <- col_names_hist

pop_hist <- pop_hist_raw |>
  rename(date = 1) |>
  mutate(
    date =as.Date(as.numeric(date), origin = "1899-12-30"),
    year = year(date)
  ) |>
  filter(year >= 1982, year <= 2021) |>
  select(-date) |>
  pivot_longer(
    cols = -year,
    names_to = "description",
    values_to = "population"
  )|>
  mutate(
    population = as.numeric(population),
    sex = case_when(
      str_detect(description, "Male") ~ "males",
      str_detect(description, "Female") ~ "females",
      TRUE ~ NA_character_
    ),
    age = str_extract(description, "(?<=;  )[0-9]+(?= ;)"),
    age = if_else(str_detect(description, "100 and over"), "100", age),
    age = as.integer(age)
  ) |>
  filter(!is.na(sex)) |>
  mutate(
    age_group = case_when(
      age <= 4              ~ "0–4",
      age >= 5 & age <= 14 ~ "5–14",
      age >= 15 & age <= 24 ~ "15–24",
      age <= 34 ~ "25–34",
      age <= 44 ~ "35–44",
      age <= 54 ~ "45–54",
      age <= 64 ~ "55–64",
      age <= 74 ~ "65–74",
      age <= 84 ~ "75–84",
      age >= 85 ~ "85+",
      TRUE ~ NA_character_
    )
  ) |>
  group_by(year, sex, age_group) |>
  summarise(population = sum(population), .groups = "drop")

# Check
pop_hist |> count(year)
pop_hist |> count(sex)
pop_hist |> count(age_group)

# Save
write_csv(pop_hist, "data/pop_hist.csv")

cat("\n=== import_data.R complete ===\n")
