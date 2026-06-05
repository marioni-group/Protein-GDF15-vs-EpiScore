### =========================================================
### Script B: Rebuild analysis_df_final_v2
### =========================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

packages <- c(
  "dplyr", "readr", "stringr", "lubridate",
  "tidyr", "forcats", "tibble"
)

invisible(lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) install.packages(x)
  library(x, character.only = TRUE)
}))

# File paths
raw_path_options <- c(
  "C:/Users/lowyi/OneDrive/Documents/Edinburgh Masters/Year 3/2ND RUN OF EVERYTHING/Raw Data csv",   # laptop
  "C:/Users/YI MEI/OneDrive/Documents/Edinburgh Masters/Year 3/2ND RUN OF EVERYTHING/Raw Data csv"    # PC
)

results_path_options <- c(
  "C:/Users/lowyi/OneDrive/Documents/Edinburgh Masters/Year 3/2ND RUN OF EVERYTHING/Results/Cleaned Data",
  "C:/Users/YI MEI/OneDrive/Documents/Edinburgh Masters/Year 3/2ND RUN OF EVERYTHING/Results/Cleaned Data"
)

DIR_RAW <- raw_path_options[file.exists(raw_path_options)][1]
DIR_RESULTS <- results_path_options[file.exists(results_path_options)][1]

cat("Using DIR_RAW:\n", DIR_RAW, "\n")
cat("Using DIR_RESULTS:\n", DIR_RESULTS, "\n")

covar    <- read_csv(file.path(DIR_RAW, "2026-01-26_covariates.csv"), show_col_types = FALSE)
ehr      <- read_csv(file.path(DIR_RAW, "2025-09-11_EHR_diseases.csv"), show_col_types = FALSE)
deaths   <- read_csv(file.path(DIR_RAW, "2025-09-11_deaths.csv"), show_col_types = FALSE)
episcore <- read_csv(file.path(DIR_RAW, "2025-09-11_GDF15_episcore.csv"), show_col_types = FALSE)
protein  <- read_csv(file.path(DIR_RAW, "2025-09-11_measured_GDF15.csv"), show_col_types = FALSE)
metadata <- read_csv(file.path(DIR_RESULTS, "df_metadata_final.csv"), show_col_types = FALSE)


# Standardise IDs
covar    <- covar %>% mutate(id = as.character(id))
ehr      <- ehr %>% mutate(id = as.character(id))
deaths   <- deaths %>% mutate(id = as.character(id))
episcore <- episcore %>% mutate(id = as.character(id))
protein  <- protein %>% mutate(id = as.character(id))
metadata <- metadata %>% mutate(disease = as.character(disease))

cat("\n--- Raw dimensions ---\n")
cat("covar    :", nrow(covar), "x", ncol(covar), "\n")
cat("ehr      :", nrow(ehr), "x", ncol(ehr), "\n")
cat("deaths   :", nrow(deaths), "x", ncol(deaths), "\n")
cat("episcore :", nrow(episcore), "x", ncol(episcore), "\n")
cat("protein  :", nrow(protein), "x", ncol(protein), "\n")
cat("metadata :", nrow(metadata), "x", ncol(metadata), "\n")

# Helpers
safe_ym <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  lubridate::ym(x)
}

winsorise <- function(x, probs = c(0.01, 0.99)) {
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  x <- pmax(x, qs[1], na.rm = FALSE)
  x <- pmin(x, qs[2], na.rm = FALSE)
  x
}

# =========================
# Disease list from metadata
# =========================
disease_list <- metadata %>%
  filter(include_model == TRUE) %>%
  distinct(disease) %>%
  arrange(disease) %>%
  pull(disease)

print(length(disease_list))

# =========================
# Exposure-complete IDs
# =========================
protein_clean <- protein %>%
  mutate(
    gdf15_log2 = if_else(!is.na(gdf15) & gdf15 > 0, log2(gdf15), NA_real_)
  )

exposure_ids <- episcore %>%
  select(id, DNAmGDF15.1) %>%
  inner_join(
    protein_clean %>% select(id, gdf15, gdf15_log2),
    by = "id"
  ) %>%
  filter(!is.na(DNAmGDF15.1), !is.na(gdf15_log2)) %>%
  distinct(id)

print(n_distinct(exposure_ids$id))

dup_check <- function(df, df_name) {
  n_dup <- df %>%
    count(id) %>%
    filter(n > 1) %>%
    nrow()
  
  cat("\n--- Duplicate ID check:", df_name, "---\n")
  cat("Duplicate IDs:", n_dup, "\n")
  
  if (n_dup > 0) {
    stop(paste("Duplicate IDs found in", df_name, "- fix before proceeding."))
  }
}

dup_check(covar, "covar")
dup_check(episcore, "episcore")
dup_check(protein, "protein")

# =========================
# Restrict covariate cohort
# =========================
covar_cohort <- covar %>%
  semi_join(exposure_ids, by = "id") %>%
  filter(!is.na(rank))

cat("\n--- Covariate cohort after rank restriction ---\n")
cat("Rows      :", nrow(covar_cohort), "\n")
cat("Unique IDs:", n_distinct(covar_cohort$id), "\n")
cat("Missing rank:", sum(is.na(covar_cohort$rank)), "\n")

# Restrict EHR and deaths to same cohort IDs
ehr_cohort <- ehr %>%
  semi_join(covar_cohort %>% distinct(id), by = "id")

deaths_cohort <- deaths %>%
  semi_join(covar_cohort %>% distinct(id), by = "id")

cat("\n--- EHR cohort after restriction ---\n")
cat("Rows      :", nrow(ehr_cohort), "\n")
cat("Unique IDs:", n_distinct(ehr_cohort$id), "\n")

cat("\n--- Deaths cohort after restriction ---\n")
cat("Rows      :", nrow(deaths_cohort), "\n")
cat("Unique IDs:", n_distinct(deaths_cohort$id), "\n")

# =========================
# Clean dates
# =========================
ehr_clean <- ehr_cohort %>%
  mutate(
    dt1_ym_date = safe_ym(dt1_ym),
    gs_appt_date = safe_ym(gs_appt)
  )

deaths_clean <- deaths_cohort %>%
  mutate(
    dod_ym_date = safe_ym(dod_ym)
  )

cat("\n--- Missing date QC ---\n")
cat("EHR missing dt1_ym_date :", sum(is.na(ehr_clean$dt1_ym_date)), "\n")
cat("EHR missing gs_appt_date:", sum(is.na(ehr_clean$gs_appt_date)), "\n")
cat("Deaths missing dod_ym_date:", sum(is.na(deaths_clean$dod_ym_date)), "\n")

# =========================
# Build ID-level baseline/follow-up table
# =========================
baseline_df <- covar_cohort %>%
  mutate(
    baseline_appt = safe_ym(appt),
    rank = suppressWarnings(as.numeric(rank))
  ) %>%
  distinct(id, .keep_all = TRUE)

baseline_df2 <- baseline_df %>%
  left_join(
    deaths_clean %>%
      distinct(id, dod_ym_date),
    by = "id"
  ) %>%
  mutate(
    death_date_final = dod_ym_date,
    admin_censor_date = as.Date("2023-10-01"),
    last_seen = case_when(
      !is.na(death_date_final) ~ death_date_final,
      TRUE ~ admin_censor_date
    ),
    last_seen = pmin(last_seen, admin_censor_date, na.rm = TRUE)
  )

baseline_df2 %>%
  summarise(
    n_rows = n(),
    n_ids = n_distinct(id),
    missing_baseline_appt = sum(is.na(baseline_appt)),
    missing_rank = sum(is.na(rank)),
    missing_last_seen = sum(is.na(last_seen)),
    negative_followup = sum(last_seen < baseline_appt, na.rm = TRUE)
  ) %>%
  print()


# Keep only included diseases in EHR
ehr_clean_inc <- ehr_clean %>%
  filter(disease %in% disease_list)

cat("\n--- Included-disease EHR QC ---\n")
cat("Rows      :", nrow(ehr_clean_inc), "\n")
cat("Unique IDs:", n_distinct(ehr_clean_inc$id), "\n")
cat("Diseases  :", n_distinct(ehr_clean_inc$disease), "\n")

# =========================
# Build full ID x disease grid
# =========================
full_grid <- tidyr::crossing(
  id = baseline_df2$id,
  disease = disease_list
)

cat("\n--- Full grid QC ---\n")
cat("Rows expected:", length(unique(baseline_df2$id)) * length(disease_list), "\n")
cat("Rows actual  :", nrow(full_grid), "\n")
cat("Unique IDs   :", n_distinct(full_grid$id), "\n")
cat("Diseases     :", n_distinct(full_grid$disease), "\n")


# Collapse EHR to ID x disease summary
ehr_summary <- ehr_clean_inc %>%
  left_join(
    baseline_df2 %>% select(id, baseline_appt),
    by = "id"
  ) %>%
  mutate(
    prevalent_flag = if_else(!is.na(dt1_ym_date) & dt1_ym_date <= baseline_appt, 1L, 0L),
    incident_flag  = if_else(!is.na(dt1_ym_date) & dt1_ym_date > baseline_appt, 1L, 0L)
  ) %>%
  group_by(id, disease) %>%
  summarise(
    prevalent = as.integer(any(prevalent_flag == 1L, na.rm = TRUE)),
    event = as.integer(any(incident_flag == 1L, na.rm = TRUE)),
    event_date = suppressWarnings(min(dt1_ym_date[incident_flag == 1L], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    event_date = if_else(is.infinite(event_date), as.Date(NA), event_date)
  )

ehr_summary %>%
  summarise(
    n_rows = n(),
    n_ids = n_distinct(id),
    n_diseases = n_distinct(disease),
    n_prevalent = sum(prevalent, na.rm = TRUE),
    n_event = sum(event, na.rm = TRUE),
    n_prev_and_event = sum(prevalent == 1 & event == 1, na.rm = TRUE)
  ) %>%
  print()


# Join back to full grid
analysis_base <- full_grid %>%
  left_join(ehr_summary, by = c("id", "disease")) %>%
  mutate(
    prevalent = if_else(is.na(prevalent), 0L, prevalent),
    event = if_else(is.na(event), 0L, event)
  )

cat("\n--- Analysis base QC ---\n")
analysis_base %>%
  summarise(
    n_rows = n(),
    n_ids = n_distinct(id),
    n_diseases = n_distinct(disease),
    n_prevalent = sum(prevalent, na.rm = TRUE),
    n_event = sum(event, na.rm = TRUE),
    n_prev_and_event = sum(prevalent == 1 & event == 1, na.rm = TRUE)
  ) %>%
  print()

cat("\nCheck one row per ID x disease:\n")
print(
  analysis_base %>%
    count(id, disease) %>%
    summarise(max_n = max(n), n_duplicates = sum(n > 1))
)

# =========================
# Add follow-up and time variables
# =========================
analysis_time <- analysis_base %>%
  left_join(
    baseline_df2 %>%
      select(id, baseline_appt, death_date_final, last_seen),
    by = "id"
  ) %>%
  mutate(
    end_of_followup = case_when(
      event == 1 & !is.na(event_date) ~ event_date,
      TRUE ~ last_seen
    ),
    time_to_event_days = case_when(
      event == 1 & !is.na(event_date) ~ as.numeric(event_date - baseline_appt),
      TRUE ~ as.numeric(end_of_followup - baseline_appt)
    ),
    time_to_event_years = time_to_event_days / 365.25,
    t_censor_years = as.numeric(last_seen - baseline_appt) / 365.25,
    death = if_else(!is.na(death_date_final), 1L, 0L),
    date_death = death_date_final
  )

analysis_time %>%
  summarise(
    n_rows = n(),
    n_ids = n_distinct(id),
    n_diseases = n_distinct(disease),
    missing_baseline = sum(is.na(baseline_appt)),
    missing_last_seen = sum(is.na(last_seen)),
    missing_end_of_followup = sum(is.na(end_of_followup)),
    missing_time_to_event_years = sum(is.na(time_to_event_years)),
    missing_t_censor_years = sum(is.na(t_censor_years)),
    negative_time_to_event = sum(time_to_event_years < 0, na.rm = TRUE),
    negative_t_censor = sum(t_censor_years < 0, na.rm = TRUE)
  ) %>%
  print()

# =========================
# Clean and derive covariates
# =========================
has_sex        <- "sex" %in% names(covar_cohort)
has_bmi        <- "bmi" %in% names(covar_cohort)
has_smk        <- "ever_smoke" %in% names(covar_cohort)
has_pack       <- "pack_years" %in% names(covar_cohort)
has_qual       <- "qualification" %in% names(covar_cohort)
has_rank       <- "rank" %in% names(covar_cohort)
has_alc_units  <- "units" %in% names(covar_cohort)
has_drink_stat <- "drink_status" %in% names(covar_cohort)

covar_clean <- covar_cohort %>%
  mutate(
    sex = if (has_sex) {
      case_when(
        toupper(as.character(sex)) %in% c("M", "MALE", "1") ~ "Male",
        toupper(as.character(sex)) %in% c("F", "FEMALE", "2") ~ "Female",
        TRUE ~ NA_character_
      )
    } else NA_character_,
    sex = factor(sex, levels = c("Female", "Male")),
    
    bmi = if (has_bmi) suppressWarnings(as.numeric(bmi)) else NA_real_,
    
    ever_smoke_num = if (has_smk) suppressWarnings(as.integer(ever_smoke)) else NA_integer_,
    
    smoking_4cat = case_when(
      ever_smoke_num == 1 ~ "Current",
      ever_smoke_num == 2 ~ "Quit<1y",
      ever_smoke_num == 3 ~ "Quit≥1y",
      ever_smoke_num == 4 ~ "Never",
      TRUE ~ NA_character_
    ),
    smoking_4cat = factor(
      smoking_4cat,
      levels = c("Never", "Quit≥1y", "Quit<1y", "Current")
    ),
    
    smoking = case_when(
      ever_smoke_num == 1 ~ "Current",
      ever_smoke_num %in% c(2, 3) ~ "Former",
      ever_smoke_num == 4 ~ "Never",
      TRUE ~ NA_character_
    ),
    smoking = factor(
      smoking,
      levels = c("Never", "Former", "Current")
    ),
    
    pack_years_num = if (has_pack) suppressWarnings(as.numeric(pack_years)) else NA_real_,
    pack_years_num = case_when(
      is.na(pack_years_num) ~ NA_real_,
      pack_years_num < 0 ~ NA_real_,
      TRUE ~ pack_years_num
    ),
    pack_years_num = ifelse(ever_smoke_num == 4, 0, pack_years_num),
    
    qual_num = if (has_qual) suppressWarnings(as.numeric(qualification)) else NA_real_,
    education_cont = qual_num,
    
    alcohol_units_num = if (has_alc_units) suppressWarnings(as.numeric(units)) else NA_real_,
    alcohol_units_num = case_when(
      is.na(alcohol_units_num) ~ NA_real_,
      alcohol_units_num < 0 ~ NA_real_,
      TRUE ~ alcohol_units_num
    ),
    alcohol_units_num = if (has_drink_stat) case_when(
      trimws(tolower(as.character(drink_status))) %in% c("3", "0", "none", "non", "never", "no") ~ 0,
      TRUE ~ alcohol_units_num
    ) else alcohol_units_num,
    
    rank = suppressWarnings(as.numeric(rank))
  )

covar_clean %>%
  summarise(
    n_rows = n(),
    n_ids = n_distinct(id),
    missing_age = sum(is.na(age)),
    missing_sex = sum(is.na(sex)),
    missing_bmi = sum(is.na(bmi)),
    missing_pack_years = sum(is.na(pack_years_num)),
    missing_education = sum(is.na(education_cont)),
    missing_alcohol = sum(is.na(alcohol_units_num)),
    missing_rank = sum(is.na(rank))
  ) %>%
  print()

# Winsorise continuous covariates
cont_trim_vars <- c("bmi", "pack_years_num", "education_cont", "alcohol_units_num")
cont_trim_vars <- cont_trim_vars[cont_trim_vars %in% names(covar_clean)]

covar_clean <- covar_clean %>%
  mutate(
    across(all_of(cont_trim_vars), ~ winsorise(.x, probs = c(0.01, 0.99))),
    alcohol_units_wins = alcohol_units_num
  )

# Mean imputation for continuous covariates only
set.seed(130196)

cont_vars <- c("bmi", "pack_years_num", "education_cont", "alcohol_units_wins")
cont_vars <- cont_vars[cont_vars %in% names(covar_clean)]

covar_clean_imp <- covar_clean %>%
  mutate(
    across(
      .cols = all_of(cont_vars),
      .fns  = ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)
    )
  )

covariates_panel <- covar_clean_imp %>%
  select(
    id,
    age,
    sex,
    bmi,
    smoking,
    smoking_4cat,
    pack_years_num,
    education_cont,
    alcohol_units_wins,
    rank
  )

covariates_panel %>%
  summarise(
    n_rows = n(),
    n_ids = n_distinct(id),
    missing_age = sum(is.na(age)),
    missing_sex = sum(is.na(sex)),
    missing_bmi = sum(is.na(bmi)),
    missing_pack_years = sum(is.na(pack_years_num)),
    missing_education = sum(is.na(education_cont)),
    missing_alcohol = sum(is.na(alcohol_units_wins)),
    missing_rank = sum(is.na(rank))
  ) %>%
  print()

# =========================
# Exposure panel
# =========================
exposure_panel <- episcore %>%
  select(id, DNAmGDF15.1) %>%
  inner_join(
    protein_clean %>% select(id, gdf15, gdf15_log2),
    by = "id"
  ) %>%
  semi_join(covar_cohort %>% distinct(id), by = "id") %>%
  filter(!is.na(DNAmGDF15.1), !is.na(gdf15_log2)) %>%
  mutate(
    gdf15_log2_wins = winsorise(gdf15_log2, probs = c(0.01, 0.99)),
    gdf15_episcore_wins = winsorise(DNAmGDF15.1, probs = c(0.01, 0.99))
  ) %>%
  mutate(
    gdf15_std = as.numeric(scale(gdf15_log2_wins)),
    gdf15_episcore_std = as.numeric(scale(gdf15_episcore_wins))
  )

exposure_panel %>%
  summarise(
    n_rows = n(),
    n_ids = n_distinct(id),
    missing_gdf15 = sum(is.na(gdf15)),
    missing_gdf15_log2 = sum(is.na(gdf15_log2)),
    missing_gdf15_log2_wins = sum(is.na(gdf15_log2_wins)),
    missing_epi = sum(is.na(DNAmGDF15.1)),
    missing_epi_wins = sum(is.na(gdf15_episcore_wins)),
    missing_gdf15_std = sum(is.na(gdf15_std)),
    missing_epi_std = sum(is.na(gdf15_episcore_std)),
    mean_gdf15_std = mean(gdf15_std, na.rm = TRUE),
    sd_gdf15_std = sd(gdf15_std, na.rm = TRUE),
    mean_epi_std = mean(gdf15_episcore_std, na.rm = TRUE),
    sd_epi_std = sd(gdf15_episcore_std, na.rm = TRUE)
  ) %>%
  print()

# =========================
# Final analysis dataset
# =========================
analysis_df_final <- analysis_time %>%
  left_join(covariates_panel, by = "id") %>%
  left_join(exposure_panel, by = "id")

analysis_df_final %>%
  summarise(
    n_rows = n(),
    n_ids = n_distinct(id),
    n_diseases = n_distinct(disease),
    missing_age = sum(is.na(age)),
    missing_sex = sum(is.na(sex)),
    missing_bmi = sum(is.na(bmi)),
    missing_pack_years = sum(is.na(pack_years_num)),
    missing_education = sum(is.na(education_cont)),
    missing_alcohol = sum(is.na(alcohol_units_wins)),
    missing_rank = sum(is.na(rank)),
    missing_gdf15_std = sum(is.na(gdf15_std)),
    missing_epi_std = sum(is.na(gdf15_episcore_std))
  ) %>%
  print()

print(
  analysis_df_final %>%
    count(id, disease) %>%
    summarise(max_n = max(n), n_duplicates = sum(n > 1))
)

analysis_df_final %>%
  select(
    id, disease, baseline_appt, prevalent, event,
    event_date, time_to_event_years, t_censor_years,
    age, sex, bmi, pack_years_num, education_cont,
    alcohol_units_wins, rank, gdf15_std, gdf15_episcore_std
  ) %>%
  slice_head(n = 10) %>%
  print(width = Inf)

# =========================
# Save 
# =========================
write_csv(
  analysis_df_final,
  file.path(DIR_RESULTS, "analysis_df_final_v2.csv")
)

saveRDS(
  analysis_df_final,
  file.path(DIR_RESULTS, "analysis_df_final_v2.rds")
)

write_csv(
  tibble(id = sort(unique(covar_cohort$id))),
  file.path(DIR_RESULTS, "valid_ids_exposure_complete.csv")
)

cat("\nSaved files:\n")
cat(file.path(DIR_RESULTS, "analysis_df_final_v2.csv"), "\n")
cat(file.path(DIR_RESULTS, "analysis_df_final_v2.rds"), "\n")
cat(file.path(DIR_RESULTS, "valid_ids_exposure_complete.csv"), "\n")
