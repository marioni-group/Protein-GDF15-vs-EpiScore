
## Script5_Results_Tables.R
## Results Tables

### Table 1: Baseline characteristics
rm(list = ls())
options(stringsAsFactors = FALSE)

packages <- c(
  "dplyr", "readr", "stringr", "tidyr",
  "glue", "openxlsx", "janitor"
)

invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))


#  File paths
## Data files are not included in this repository.

# Required inputs:
## data/raw/original raw covariate and biomarker files
## results/cleaned_data/analysis_df_final_v2.csv
## results/cox_models/cox_results_main.csv

# Outputs:
## results/tables/

DIR_RAW <- file.path("data", "raw")
DIR_CLEAN <- file.path("results", "cleaned_data")
DIR_COX <- file.path("results", "cox_models")
DIR_TABLES <- file.path("results", "tables")

dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)

cat("Using DIR_RAW:\n", DIR_RAW, "\n")
cat("Using DIR_CLEAN:\n", DIR_CLEAN, "\n")
cat("Using DIR_COX:\n", DIR_COX, "\n")
cat("Saving tables to:\n", DIR_TABLES, "\n")

#Load 
analysis_file <- file.path(DIR_CLEAN, "analysis_df_final_v2.csv")
df <- read_csv(analysis_file, show_col_types = FALSE)
cat("\nLoaded:\n", analysis_file, "\n")

covar_raw <- read_csv(
  file.path(DIR_RAW, "2026-01-26_covariates.csv"),
  show_col_types = FALSE
) %>%
  mutate(id = as.character(id))

# Collapse to participant-level cohort
df_table1 <- df %>%
  mutate(id = as.character(id)) %>%
  arrange(id, disease) %>%
  distinct(id, .keep_all = TRUE)

covar_raw <- covar_raw %>%
  mutate(id = as.character(id))

df_table1 <- df_table1 %>%
  left_join(
    covar_raw %>%
      select(id, qualification),
    by = "id"
  )

print(
  tibble(
    n_rows = nrow(df_table1),
    n_ids = n_distinct(df_table1$id),
    duplicated_ids = sum(duplicated(df_table1$id)),
    n_deaths = if ("death" %in% names(df_table1)) sum(df_table1$death == 1, na.rm = TRUE) else NA_integer_
  )
)

### Helper functions
fmt_mean_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  glue("{round(mean(x, na.rm = TRUE), 2)} ({round(sd(x, na.rm = TRUE), 2)})")
}

fmt_median_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  glue("{round(q[2], 2)} ({round(q[1], 2)}, {round(q[3], 2)})")
}

fmt_n_pct <- function(x) {
  n <- sum(x, na.rm = TRUE)
  pct <- 100 * n / length(x)
  glue("{n} ({round(pct, 1)}%)")
}

add_cont_row <- function(data, variable, label, summary_type = c("mean_sd", "median_iqr")) {
  summary_type <- match.arg(summary_type)
  
  if (!variable %in% names(data)) {
    return(tibble(
      Characteristic = label,
      Level = "",
      Value = "Variable not found"
    ))
  }
  
  value <- if (summary_type == "mean_sd") {
    fmt_mean_sd(data[[variable]])
  } else {
    fmt_median_iqr(data[[variable]])
  }
  
  tibble(
    Characteristic = label,
    Level = "",
    Value = as.character(value)
  )
}

add_cat_rows <- function(data, variable, label) {
  if (!variable %in% names(data)) {
    return(tibble(
      Characteristic = label,
      Level = "",
      Value = "Variable not found"
    ))
  }
  
  data %>%
    filter(!is.na(.data[[variable]])) %>%
    mutate(Level = as.character(.data[[variable]])) %>%
    count(Level) %>%
    mutate(
      pct = 100 * n / sum(n),
      Characteristic = label,
      Value = glue("{n} ({round(pct, 1)}%)")
    ) %>%
    select(Characteristic, Level, Value)
}

death_row <- function(data) {
  if (!"death" %in% names(data)) {
    return(tibble(
      Characteristic = "Deaths during available follow-up",
      Level = "",
      Value = "Variable not found"
    ))
  }
  
  dead <- data$death == 1
  
  tibble(
    Characteristic = "Deaths during available follow-up",
    Level = "",
    Value = as.character(fmt_n_pct(dead))
  )
}


# SIMD quintile 
df_table1 <- df_table1 %>%
  mutate(
    simd_quintile = if ("rank" %in% names(.)) {
      ntile(rank, 5)
    } else {
      NA_integer_
    },
    simd_quintile = factor(
      simd_quintile,
      levels = 1:5,
      labels = c("Q1", "Q2", "Q3", "Q4", "Q5")
    )
  )


### Build Table 1 
table1 <- bind_rows(
  tibble(
    Characteristic = "Analytic cohort, N",
    Level = "",
    Value = as.character(nrow(df_table1))
  ),
  
  add_cont_row(df_table1, "age", "Age at baseline, years", "mean_sd"),
  add_cat_rows(df_table1, "sex", "Sex"),
  add_cont_row(df_table1, "bmi", "Body mass index, kg/m²", "mean_sd"),
  add_cat_rows(df_table1, "smoking", "Smoking status"),
  add_cont_row(df_table1, "pack_years_num", "Smoking pack-years", "median_iqr"),
  add_cont_row(df_table1, "alcohol_units_wins", "Alcohol units per week", "median_iqr"),
  add_cat_rows(df_table1, "qualification", "Educational qualification code"),
  add_cont_row(df_table1, "rank", "SIMD rank", "median_iqr"),
  add_cont_row(df_table1, "gdf15", "Measured GDF15 protein", "median_iqr"),
  add_cont_row(df_table1, "gdf15_log2", "Log2 measured GDF15 protein", "mean_sd"),
  add_cont_row(df_table1, "DNAmGDF15.1", "DNAm GDF15 EpiScore", "median_iqr"),
  add_cont_row(df_table1, "t_censor_years", "Available follow-up time, years", "median_iqr"),
  death_row(df_table1)
)

print(table1, n = Inf)

# Save
csv_out <- file.path(DIR_TABLES, "results_table1_baseline_characteristics.csv")
xlsx_out <- file.path(DIR_TABLES, "results_table1_baseline_characteristics.xlsx")

write_csv(table1, csv_out)

openxlsx::write.xlsx(
  table1,
  file = xlsx_out,
  overwrite = TRUE
)

cat(csv_out, "\n")
cat(xlsx_out, "\n")

############################################################
### Table 2: Concordant and discordant disease associations
############################################################

# Load 
top_concordant_file <- file.path(DIR_TABLES, "scriptG_DE_top20_concordant.csv")
top_discordant_file <- file.path(DIR_TABLES, "scriptG_top20_discordant.csv")

top_concordant <- read_csv(top_concordant_file, show_col_types = FALSE)
top_discordant <- read_csv(top_discordant_file, show_col_types = FALSE)

cat(top_concordant_file, "\n")
cat(top_discordant_file, "\n")
  
print(names(top_concordant))
print(names(top_discordant))


# Table 2
table2_concordant <- top_concordant %>%
  arrange(rank_concordant) %>%
  slice_head(n = 8) %>%
  transmute(
    Section = "Top concordant associations",
    Disease = disease_label,
    `Protein log(HR)` = round(logHR_protein, 3),
    `EpiScore log(HR)` = round(logHR_episcore, 3),
    `Absolute difference` = round(abs_logHR_diff, 3),
    `Bonferroni concordance category` = bonf_concordance_group
  )

table2_discordant <- top_discordant %>%
  arrange(rank_abs_diff) %>%
  slice_head(n = 8) %>%
  transmute(
    Section = "Top discordant associations",
    Disease = disease_label,
    `Protein log(HR)` = round(logHR_protein, 3),
    `EpiScore log(HR)` = round(logHR_episcore, 3),
    `Absolute difference` = round(abs_logHR_diff, 3),
    `Bonferroni concordance category` = bonf_concordance_group
  )

table2 <- bind_rows(
  table2_concordant,
  table2_discordant
)

print(table2, n = Inf)

# Save 
csv_out_table2 <- file.path(DIR_TABLES, "results_table2_concordant_discordant_associations.csv")
xlsx_out_table2 <- file.path(DIR_TABLES, "results_table2_concordant_discordant_associations.xlsx")

write_csv(table2, csv_out_table2)

openxlsx::write.xlsx(
  table2,
  file = xlsx_out_table2,
  overwrite = TRUE
)

cat(csv_out_table2, "\n")
cat(xlsx_out_table2, "\n")
