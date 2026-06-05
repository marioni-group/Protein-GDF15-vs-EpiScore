### Script2_Disease_Metadata.R
### Disease Metadata and Eligibility Assessment

# Disease eligibility rules for full loop analysis:
# 1. Person-level modelling restricted to prevalent == 0
# 2. Disease-level inclusion requires >= 30 incident cases
# 3. Final incident/event definition will be assessed within a 10-year follow-up window

packages <- c("readr", "dplyr", "stringr", "janitor", "tibble")

invisible(lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) install.packages(x)
  library(x, character.only = TRUE)
}))

# File paths
## Data files are not included in this repository.

# Required inputs:
## data/raw/2025-09-11_EHR_diseases.csv
## results/cleaned_data/analysis_df_final_v2.csv

#Outputs:
## results/metadata/

DIR_RAW <- file.path("data", "raw")
DIR_RESULTS <- file.path("results", "cleaned_data")
DIR_METADATA <- file.path("results", "metadata")

dir.create(DIR_METADATA, recursive = TRUE, showWarnings = FALSE)

FILE_ANALYSIS <- file.path(DIR_RESULTS, "analysis_df_final_v2.csv")
FILE_EHR_RAW <- file.path(DIR_RAW, "2025-09-11_EHR_diseases.csv")

analysis_df <- read_csv(FILE_ANALYSIS, show_col_types = FALSE) %>%
clean_names()

ehr_raw <- read_csv(FILE_EHR_RAW, show_col_types = FALSE) %>%
clean_names()

# Unique IDs
qc_ids <- tibble(
  dataset = c("analysis_df", "ehr_raw"),
  n_unique_ids = c(
    dplyr::n_distinct(analysis_df$id),
    dplyr::n_distinct(ehr_raw$id)
  )
)

print(qc_ids)

# Unique diseases
qc_diseases <- tibble(
  dataset = c("analysis_df", "ehr_raw"),
  n_unique_diseases = c(
    dplyr::n_distinct(analysis_df$disease),
    dplyr::n_distinct(ehr_raw$disease)
  )
)

print(qc_diseases)
print(names(analysis_df))
print(names(ehr_raw))

print(head(analysis_df, 10))
print(head(ehr_raw, 10))

### =========================================================
### Disease-level summary 
### =========================================================
df_disease_summary <- analysis_df %>%
  group_by(disease) %>%
  summarise(
    n_total = n(),
    n_ids = n_distinct(id),
    n_prevalent = sum(prevalent == 1),
    n_incident = sum(event == 1),
    n_non_cases = sum(event == 0),
    .groups = "drop"
  ) %>%
  arrange(desc(n_incident))

print(df_disease_summary, n = 50)

### =========================================================
### Define 10-year event variables
### =========================================================

analysis_df <- analysis_df %>%
  mutate(
    # Cap follow-up at 10 years
    time_10y = pmin(time_to_event_years, 10),
    # Event occurs only if within 10 years
    event_10y = case_when(
      event == 1 & time_to_event_years <= 10 ~ 1,
      TRUE ~ 0
    )
  )

qc_event_10y <- analysis_df %>%
  summarise(
    n_event_original = sum(event == 1),
    n_event_10y = sum(event_10y == 1)
  )

print(qc_event_10y)

### =========================================================
### Disease summary (10-year events)
### =========================================================

df_disease_summary_10y <- analysis_df %>%
  group_by(disease) %>%
  summarise(
    n_total = n(),
    n_ids = n_distinct(id),
    n_prevalent = sum(prevalent == 1),
    n_incident = sum(event == 1),
    n_incident_10y = sum(event_10y == 1),
    n_non_cases_10y = sum(event_10y == 0),
    .groups = "drop"
  ) %>%
  arrange(desc(n_incident_10y))

print(df_disease_summary_10y, n = 50)

### =========================================================
### Disease metadata table
### =========================================================

df_metadata <- df_disease_summary_10y %>%
  mutate(
    include_model = n_incident_10y >= 30,
    # Exclusion incident <30
    exclude_reason = case_when(
      n_incident_10y < 30 ~ "Fewer than 30 incident cases within 10 years",
      TRUE ~ NA_character_
    ),
    
    disease_label = disease,
    category = NA_character_
  ) %>%
  select(
    disease,
    disease_label,
    category,
    n_total,
    n_prevalent,
    n_incident,
    n_incident_10y,
    include_model,
    exclude_reason
  )

metadata_summary <- df_metadata %>%
  summarise(
    n_total_diseases = n(),
    n_included = sum(include_model),
    n_excluded = sum(!include_model)
  )

print(metadata_summary)

# Excluded diseases
df_metadata %>%
  filter(include_model == FALSE) %>%
  arrange(n_incident_10y) %>%
  print(n = 50)

### =========================================================
### Load disease readme / mapping file
### =========================================================

# Input relevant FILE_DISEASE_MAP 

disease_map <- read_tsv(FILE_DISEASE_MAP, show_col_types = FALSE) %>%
  clean_names()

print(names(disease_map))
print(head(disease_map, 10))

qc_disease_map <- disease_map %>%
  summarise(
    n_rows = n(),
    n_unique_name = n_distinct(name, na.rm = TRUE),
    n_missing_name = sum(is.na(name) | name == ""),
    n_missing_phenotype = sum(is.na(phenotype) | phenotype == ""),
    n_missing_group = sum(is.na(group) | group == "")
  )

print(qc_disease_map)

### =========================================================
### Construct metadata table with labels + categories
### =========================================================

df_metadata <- df_disease_summary_10y %>%
  left_join(
    disease_map %>%
      transmute(
        disease = name,
        disease_label = phenotype,
        category = group
      ),
    by = "disease"
  ) %>%
  mutate(
    include_model = n_incident_10y >= 30,
    exclude_reason = case_when(
      n_incident_10y < 30 ~ "Fewer than 30 incident cases within 10 years",
      TRUE ~ NA_character_
    ),
    disease_label = if_else(is.na(disease_label) | disease_label == "", disease, disease_label)
  ) %>%
  select(
    disease,
    disease_label,
    category,
    n_total,
    n_prevalent,
    n_incident,
    n_incident_10y,
    include_model,
    exclude_reason
  )

metadata_summary <- df_metadata %>%
  summarise(
    n_total_diseases = n(),
    n_included = sum(include_model),
    n_excluded = sum(!include_model),
    n_missing_category = sum(is.na(category) | category == ""),
    n_missing_label = sum(is.na(disease_label) | disease_label == "")
  )

print(metadata_summary)

### =========================================================
### Save - analysis metadata table 
### =========================================================

FILE_METADATA_CSV  <- file.path(DIR_RESULTS, "df_metadata_final.csv")
FILE_METADATA_XLSX <- file.path(DIR_RESULTS, "df_metadata_final.xlsx")

write_csv(df_metadata, FILE_METADATA_CSV)

if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}
library(writexl)

write_xlsx(df_metadata, FILE_METADATA_XLSX)

### =========================================================
### Presentation metadata table (included diseases only)
### =========================================================

included_diseases <- df_metadata %>%
  filter(include_model == TRUE) %>%
  select(disease, disease_label, category)

df_presentation_metadata <- analysis_df %>%
  inner_join(included_diseases, by = "disease") %>%
  filter(prevalent == 0) %>%
  mutate(
    case_10y = event_10y == 1,
    control_10y = event_10y == 0,
    age_at_diagnosis = if_else(case_10y, age + time_to_event_years, NA_real_),
    female = case_when(
      sex %in% c("Female", "F", "female") ~ 1,
      sex %in% c("Male", "M", "male") ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(category, disease_label, disease) %>%
  summarise(
    ncase = sum(case_10y, na.rm = TRUE),
    ncontrol = sum(control_10y, na.rm = TRUE),
    
    age_cases_blood_mean = mean(age[case_10y], na.rm = TRUE),
    age_cases_blood_sd   = sd(age[case_10y], na.rm = TRUE),
    
    age_cases_dx_mean = mean(age_at_diagnosis[case_10y], na.rm = TRUE),
    age_cases_dx_sd   = sd(age_at_diagnosis[case_10y], na.rm = TRUE),
    
    time_to_dx_mean = mean(time_to_event_years[case_10y], na.rm = TRUE),
    time_to_dx_sd   = sd(time_to_event_years[case_10y], na.rm = TRUE),
    
    age_controls_blood_mean = mean(age[control_10y], na.rm = TRUE),
    age_controls_blood_sd   = sd(age[control_10y], na.rm = TRUE),
    
    female_cases_pct = 100 * mean(female[case_10y], na.rm = TRUE),
    female_controls_pct = 100 * mean(female[control_10y], na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(category, disease_label)

print(df_presentation_metadata, n = 50)

df_presentation_metadata_export <- df_presentation_metadata %>%
  transmute(
    Group = category,
    Disease = disease_label,
    ncase = ncase,
    ncontrol = ncontrol,
    
    `Age of cases at blood draw (years) - Mean` = round(age_cases_blood_mean, 1),
    `Age of cases at blood draw (years) - SD`   = round(age_cases_blood_sd, 1),
    
    `Age of cases at diagnosis (years) - Mean` = round(age_cases_dx_mean, 1),
    `Age of cases at diagnosis (years) - SD`   = round(age_cases_dx_sd, 1),
    
    `Time to diagnosis (years) - Mean` = round(time_to_dx_mean, 1),
    `Time to diagnosis (years) - SD`   = round(time_to_dx_sd, 1),
    
    `Age of controls at blood draw - Mean` = round(age_controls_blood_mean, 1),
    `Age of controls at blood draw - SD`   = round(age_controls_blood_sd, 1),
    
    `Female Cases %`    = round(female_cases_pct, 1),
    `Female controls %` = round(female_controls_pct, 1)
  )

print(df_presentation_metadata_export, n = 30)

FILE_PRESENT_META_CSV  <- file.path(DIR_RESULTS, "df_presentation_metadata.csv")
FILE_PRESENT_META_XLSX <- file.path(DIR_RESULTS, "df_presentation_metadata.xlsx")

write_csv(df_presentation_metadata_export, FILE_PRESENT_META_CSV)
write_xlsx(df_presentation_metadata_export, FILE_PRESENT_META_XLSX)
