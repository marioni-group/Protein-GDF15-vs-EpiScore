### =========================================================
### Script F: Full looping Cox analysis across 173 diseases
### =========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

packages <- c(
  "dplyr", "readr", "stringr", "survival",
  "broom", "purrr", "tibble", "tidyr"
)

invisible(lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) install.packages(x)
  library(x, character.only = TRUE)
}))


# File paths
## Data files are not included in this repository

#Required inputs:
## results/cleaned_data/analysis_df_final_v2.csv
## results/metadata/df_metadata_final.csv

#Outputs:
## results/cox_models/

DIR_CLEANED_DATA <- file.path("results", "cleaned_data")
DIR_METADATA <- file.path("results", "metadata")
DIR_COX_RESULTS <- file.path("results", "cox_models")

dir.create(DIR_COX_RESULTS, recursive = TRUE, showWarnings = FALSE)

FILE_ANALYSIS <- file.path(DIR_CLEANED_DATA, "analysis_df_final_v2.csv")
FILE_METADATA <- file.path(DIR_METADATA, "df_metadata_final.csv")

analysis_df <- read_csv(FILE_ANALYSIS, show_col_types = FALSE)
metadata_df <- read_csv(FILE_METADATA, show_col_types = FALSE)

### ====================================
### 10-year survival variables
### ====================================

analysis_df <- analysis_df %>%
  mutate(
    event_10y = case_when(
      event == 1 & !is.na(time_to_event_years) & time_to_event_years <= 10 ~ 1L,
      TRUE ~ 0L
    ),
    time_10y = case_when(
      event_10y == 1L ~ time_to_event_years,
      event_10y == 0L ~ pmin(t_censor_years, 10),
      TRUE ~ NA_real_
    )
  )

print(table(analysis_df$event, useNA = "ifany"))
print(table(analysis_df$event_10y, useNA = "ifany"))

print(table(
  original_event = analysis_df$event,
  event_10y = analysis_df$event_10y,
  useNA = "ifany"
))

print(summary(analysis_df$time_to_event_years))
print(summary(analysis_df$t_censor_years))
print(summary(analysis_df$time_10y))
print(sum(is.na(analysis_df$time_10y)))
print(sum(analysis_df$time_10y < 0, na.rm = TRUE))
print(sum(analysis_df$time_10y > 10, na.rm = TRUE))

print(
  analysis_df %>%
    filter(event == 1, !is.na(time_to_event_years), time_to_event_years > 10) %>%
    summarise(n = n())
)

print(
  analysis_df %>%
    filter(event_10y == 1) %>%
    summarise(
      n_missing_time_10y = sum(is.na(time_10y)),
      n_nonpositive_time_10y = sum(time_10y <= 0, na.rm = TRUE)
    )
)

print(
  analysis_df %>%
    filter(event_10y == 0) %>%
    summarise(
      n_missing_t_censor_years = sum(is.na(t_censor_years)),
      n_missing_time_10y = sum(is.na(time_10y))
    )
)

print(
  analysis_df %>%
    select(id, disease, prevalent, event, time_to_event_years, t_censor_years, event_10y, time_10y) %>%
    slice_head(n = 10)
)

# =========================================================
# Define disease list for looping
# =========================================================

disease_list <- metadata_df %>%
  filter(include_model == TRUE) %>%
  arrange(disease) %>%
  pull(disease)

print(length(disease_list))
print(sum(disease_list %in% analysis_df$disease))

# =========================================================
# Test run on disease-specific dataset (AF)
# =========================================================

test_disease <- "AF"

df_test <- analysis_df %>%
  filter(disease == test_disease)

print(nrow(df_test))
print(n_distinct(df_test$id))
print(table(df_test$event_10y, useNA = "ifany"))

df_test <- df_test %>%
  filter(prevalent == 0)


print(nrow(df_test))
print(n_distinct(df_test$id))
print(table(df_test$event_10y, useNA = "ifany"))
print(summary(df_test$time_10y))
print(sum(is.na(df_test$time_10y)))
print(sum(df_test$time_10y > 10, na.rm = TRUE))
print(sum(df_test$time_10y < 0, na.rm = TRUE))

# =========================================================
# Function to run Cox models (1 disease) 
# =========================================================

run_cox_models <- function(df, disease_name, metadata_df) {
  
  meta_row <- metadata_df %>%
    filter(disease == disease_name) %>%
    slice(1)
  
  if (nrow(meta_row) == 0) {
    return(tibble())
  }
  
  disease_label <- meta_row$disease_label[[1]]
  category <- meta_row$category[[1]]
  
  if (nrow(df) == 0) {
    return(tibble())
  }
  
  required_cols_local <- c(
    "time_10y", "event_10y", "sex", "age",
    "bmi", "pack_years_num", "education_cont",
    "alcohol_units_wins", "rank",
    "gdf15_std", "gdf15_episcore_std"
  )
  
  if (!all(required_cols_local %in% names(df))) {
    return(tibble())
  }
  
  n_events_pre <- sum(df$event_10y, na.rm = TRUE)
  if (n_events_pre == 0) {
    return(tibble())
  }
  
  # Sex-specific rule - set at 90%
  female_event_n <- sum(df$event_10y == 1 & df$sex == "Female", na.rm = TRUE)
  total_event_n  <- sum(df$event_10y == 1, na.rm = TRUE)
  
  prop_f <- female_event_n / total_event_n
  
  sex_rule <- case_when(
    prop_f <= 0.1 ~ "male_only",
    prop_f >= 0.9 ~ "female_only",
    TRUE ~ "mixed"
  )
  
  if (sex_rule == "male_only") {
    df_model <- df %>% filter(sex == "Male")
    sex_in_model <- FALSE
  } else if (sex_rule == "female_only") {
    df_model <- df %>% filter(sex == "Female")
    sex_in_model <- FALSE
  } else {
    df_model <- df
    sex_in_model <- TRUE
  }
  
  n_total_post_sex <- nrow(df_model)
  n_events_post_sex <- sum(df_model$event_10y, na.rm = TRUE)
  
  if (n_total_post_sex == 0 || n_events_post_sex == 0) {
    return(tibble())
  }
  
  if (sex_in_model) {
    formula_list <- list(
      A = Surv(time_10y, event_10y) ~ gdf15_std + age + sex,
      B = Surv(time_10y, event_10y) ~ gdf15_episcore_std + age + sex,
      C = Surv(time_10y, event_10y) ~ gdf15_std + gdf15_episcore_std + age + sex,
      D = Surv(time_10y, event_10y) ~ gdf15_std + age + sex + bmi + pack_years_num + education_cont + alcohol_units_wins + rank,
      E = Surv(time_10y, event_10y) ~ gdf15_episcore_std + age + sex + bmi + pack_years_num + education_cont + alcohol_units_wins + rank,
      F = Surv(time_10y, event_10y) ~ gdf15_std + gdf15_episcore_std + age + sex + bmi + pack_years_num + education_cont + alcohol_units_wins + rank
    )
  } else {
    formula_list <- list(
      A = Surv(time_10y, event_10y) ~ gdf15_std + age,
      B = Surv(time_10y, event_10y) ~ gdf15_episcore_std + age,
      C = Surv(time_10y, event_10y) ~ gdf15_std + gdf15_episcore_std + age,
      D = Surv(time_10y, event_10y) ~ gdf15_std + age + bmi + pack_years_num + education_cont + alcohol_units_wins + rank,
      E = Surv(time_10y, event_10y) ~ gdf15_episcore_std + age + bmi + pack_years_num + education_cont + alcohol_units_wins + rank,
      F = Surv(time_10y, event_10y) ~ gdf15_std + gdf15_episcore_std + age + bmi + pack_years_num + education_cont + alcohol_units_wins + rank
    )
  }
  
  results_list <- list()
  
  for (model_name in names(formula_list)) {
    
    fit <- tryCatch(
      coxph(formula_list[[model_name]], data = df_model),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    tidy_fit <- tryCatch(
      broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(tidy_fit)) next
    
    tidy_fit <- tidy_fit %>%
      filter(term %in% c("gdf15_std", "gdf15_episcore_std"))
    
    if (nrow(tidy_fit) == 0) next
    
    fit_n <- fit$n
    fit_nevent <- fit$nevent
    
    tidy_fit <- tidy_fit %>%
      mutate(
        disease = disease_name,
        disease_label = disease_label,
        category = category,
        model = model_name,
        exposure = term,
        n = fit_n,
        n_events = fit_nevent,
        n_pre_sex_filter = nrow(df),
        n_events_pre_sex_filter = n_events_pre,
        n_post_sex_filter = n_total_post_sex,
        n_events_post_sex_filter = n_events_post_sex,
        HR = estimate,
        lower_CI = conf.low,
        upper_CI = conf.high,
        p_value = p.value,
        prop_f = prop_f,
        sex_rule = sex_rule,
        sex_in_model = sex_in_model
      ) %>%
      select(
        disease,
        disease_label,
        category,
        model,
        exposure,
        n,
        n_events,
        n_pre_sex_filter,
        n_events_pre_sex_filter,
        n_post_sex_filter,
        n_events_post_sex_filter,
        HR,
        lower_CI,
        upper_CI,
        p_value,
        prop_f,
        sex_rule,
        sex_in_model
      )
    
    results_list[[model_name]] <- tidy_fit
  }
  
  bind_rows(results_list)
}

# ===============================
# Loop analysis
# ===============================

all_results_list <- list()
failed_log <- list()

for (d in disease_list) {
  
  cat("\nRunning disease:", d)
  
  df_d <- analysis_df %>%
    filter(disease == d) %>%
    filter(prevalent == 0)
  
  if (nrow(df_d) == 0) {
    failed_log[[d]] <- tibble(
      disease = d,
      reason = "no_rows_after_prevalent_filter"
    )
    next
  }
  
  if (sum(df_d$event_10y, na.rm = TRUE) == 0) {
    failed_log[[d]] <- tibble(
      disease = d,
      reason = "no_10y_events_after_prevalent_filter"
    )
    next
  }
  
  res <- tryCatch(
    run_cox_models(
      df = df_d,
      disease_name = d,
      metadata_df = metadata_df
    ),
    error = function(e) {
      failed_log[[d]] <<- tibble(
        disease = d,
        reason = paste0("model_error: ", conditionMessage(e))
      )
      NULL
    }
  )
  
  if (is.null(res) || nrow(res) == 0) {
    if (is.null(failed_log[[d]])) {
      failed_log[[d]] <- tibble(
        disease = d,
        reason = "no_model_output"
      )
    }
    next
  }
  
  all_results_list[[d]] <- res
}

final_results <- bind_rows(all_results_list)
failed_models_log <- bind_rows(failed_log)

# ======================
# Check results
# ======================

print(length(disease_list))
print(length(all_results_list))
print(nrow(failed_models_log))
print(failed_models_log)
print(nrow(final_results))
print(n_distinct(final_results$disease))
print(table(final_results$model))
print(table(final_results$exposure))
print(table(final_results$sex_rule))
print(table(final_results$sex_in_model))

print(
  final_results %>%
    distinct(disease, sex_rule) %>%
    count(sex_rule)
)

print(
  final_results %>%
    summarise(across(everything(), ~ sum(is.na(.))))
)

# =========================================================
# Save
# =========================================================

DIR_COX <- file.path(DIR_RESULTS, "Cox analysis")

if (!dir.exists(DIR_COX)) {
  dir.create(DIR_COX, recursive = TRUE)
}

write.csv(
  final_results,
  file = file.path(DIR_COX, "cox_results_main.csv"),
  row.names = FALSE
)

write.csv(
  failed_models_log,
  file = file.path(DIR_COX, "cox_failed_models_log.csv"),
  row.names = FALSE
)

disease_summary <- final_results %>%
  distinct(
    disease, disease_label, category,
    sex_rule, prop_f,
    n_pre_sex_filter, n_events_pre_sex_filter,
    n_post_sex_filter, n_events_post_sex_filter
  )

write.csv(
  disease_summary,
  file = file.path(DIR_COX, "cox_disease_summary.csv"),
  row.names = FALSE
)

saveRDS(
  final_results,
  file = file.path(DIR_COX, "cox_results_main.rds")
)

saveRDS(
  disease_summary,
  file = file.path(DIR_COX, "cox_disease_summary.rds")
)

cat(file.path(DIR_COX, "cox_results_main.csv"), "\n")
cat(file.path(DIR_COX, "cox_failed_models_log.csv"), "\n")
cat(file.path(DIR_COX, "cox_disease_summary.csv"), "\n")
cat(file.path(DIR_COX, "cox_results_main.rds"), "\n")
cat(file.path(DIR_COX, "cox_disease_summary.rds"), "\n")
