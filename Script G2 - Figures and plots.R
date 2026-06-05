### =========================
### Script G: Final post-Cox analysis figures 
### =========================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

### Packages
packages <- c(
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "ggplot2",
  "forcats",
  "purrr",
  "tibble",
  "ggrepel"
)

invisible(lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) install.packages(x)
  library(x, character.only = TRUE)
}))

### =========================
### File paths
### =========================

results_path_options <- c(
  "C:/Users/lowyi/OneDrive/Documents/Edinburgh Masters/Year 3/2ND RUN OF EVERYTHING/Results",
  "C:/Users/YI MEI/OneDrive/Documents/Edinburgh Masters/Year 3/2ND RUN OF EVERYTHING/Results"
)

DIR_RESULTS <- results_path_options[file.exists(results_path_options)][1]

DIR_COX <- file.path(DIR_RESULTS, "Cox analysis")
DIR_FIG <- file.path(DIR_COX, "figures")
DIR_TAB <- file.path(DIR_COX, "tables")

dir.create(DIR_COX, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TAB, showWarnings = FALSE, recursive = TRUE)

cox_results_file <- file.path(DIR_COX, "cox_results_main.csv")

cox_results <- read_csv(cox_results_file, show_col_types = FALSE)

### =========================
### Prepare Model D/E paired dataset
### =========================

de_long <- cox_results %>%
  filter(model %in% c("D", "E"))

de_exposure_check <- de_long %>%
  distinct(model, exposure) %>%
  arrange(model, exposure)

print(de_exposure_check)

de_wide <- de_long %>%
  mutate(
    source = case_when(
      model == "D" ~ "protein",
      model == "E" ~ "episcore",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    disease,
    disease_label,
    category,
    sex_rule,
    sex_in_model,
    source,
    n,
    n_events,
    HR,
    lower_CI,
    upper_CI,
    p_value,
    prop_f
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = c(n, n_events, HR, lower_CI, upper_CI, p_value, prop_f),
    names_sep = "_"
  ) %>%
  arrange(disease)

de_comp <- de_wide %>%
  mutate(
    logHR_protein = log(HR_protein),
    logHR_episcore = log(HR_episcore),
    logHR_diff = logHR_protein - logHR_episcore,
    abs_logHR_diff = abs(logHR_diff),
    direction_protein = case_when(
      logHR_protein > 0 ~ "positive",
      logHR_protein < 0 ~ "negative",
      TRUE ~ "null"
    ),
    direction_episcore = case_when(
      logHR_episcore > 0 ~ "positive",
      logHR_episcore < 0 ~ "negative",
      TRUE ~ "null"
    ),
    opposite_direction = case_when(
      direction_protein == "positive" & direction_episcore == "negative" ~ TRUE,
      direction_protein == "negative" & direction_episcore == "positive" ~ TRUE,
      TRUE ~ FALSE
    )
  )

de_concordance <- de_comp %>%
  mutate(
    sig_protein = p_value_protein < 0.05,
    sig_episcore = p_value_episcore < 0.05,
    sig_concordance_group = case_when(
      sig_protein & sig_episcore ~ "both",
      sig_protein & !sig_episcore ~ "protein_only",
      !sig_protein & sig_episcore ~ "episcore_only",
      TRUE ~ "neither"
    )
  )

n_diseases <- n_distinct(de_concordance$disease)
bonf_threshold <- 0.05 / n_diseases

de_concordance <- de_concordance %>%
  mutate(
    sig_protein_bonf = p_value_protein < bonf_threshold,
    sig_episcore_bonf = p_value_episcore < bonf_threshold,
    bonf_concordance_group = case_when(
      sig_protein_bonf & sig_episcore_bonf ~ "both",
      sig_protein_bonf & !sig_episcore_bonf ~ "protein_only",
      !sig_protein_bonf & sig_episcore_bonf ~ "episcore_only",
      TRUE ~ "neither"
    )
  )

discordant_ranked <- de_concordance %>%
  arrange(desc(abs_logHR_diff)) %>%
  mutate(rank_abs_diff = row_number())

concordant_ranked <- de_concordance %>%
  filter(sig_concordance_group == "both") %>%
  arrange(abs_logHR_diff) %>%
  mutate(rank_concordant = row_number())

write_csv(de_concordance, file.path(DIR_TAB, "scriptG_modelDE_with_concordance.csv"))
write_csv(discordant_ranked, file.path(DIR_TAB, "scriptG_all_ranked_discordance.csv"))
write_csv(concordant_ranked, file.path(DIR_TAB, "scriptG_DE_ranked_concordant.csv"))

### =========================
### FIGURE 1: Master scatter plot
### =========================

pearson_test <- cor.test(
  de_concordance$logHR_protein,
  de_concordance$logHR_episcore,
  method = "pearson"
)

pearson_label <- paste0(
  "Pearson r = ", round(unname(pearson_test$estimate), 2),
  "\nP < 0.001"
)

master_scatter_df <- de_concordance %>%
  mutate(
    bonf_protein = p_value_protein < bonf_threshold,
    bonf_episcore = p_value_episcore < bonf_threshold,
    bonf_group = case_when(
      bonf_protein & bonf_episcore ~ "Both",
      bonf_protein & !bonf_episcore ~ "Protein only",
      !bonf_protein & bonf_episcore ~ "EpiScore only",
      TRUE ~ "Neither"
    ),
    bonf_group = factor(
      bonf_group,
      levels = c("Both", "Protein only", "EpiScore only", "Neither")
    )
  )

top_concordant_labels <- concordant_ranked %>%
  slice(1:5) %>%
  mutate(label_type = "Top concordant")

top_discordant_labels <- discordant_ranked %>%
  slice(1:5) %>%
  mutate(label_type = "Top discordant")

label_df <- bind_rows(
  top_concordant_labels,
  top_discordant_labels
) %>%
  distinct(disease, .keep_all = TRUE)

fig_master_scatter <- ggplot(
  master_scatter_df,
  aes(x = logHR_protein, y = logHR_episcore)
) +
  geom_point(
    aes(colour = bonf_group),
    size = 2.2,
    alpha = 0.72
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.5
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = disease_label),
    size = 3.8,
    max.overlaps = 30,
    box.padding = 0.65,
    point.padding = 0.45,
    segment.color = "grey55",
    seed = 130196
  ) +
  annotate(
    "label",
    x = -0.9,
    y = 1.12,
    label = pearson_label,
    hjust = 0,
    vjust = 1,
    size = 3.5,
    fill = "white",
    colour = "black"
  ) +
  scale_colour_manual(
    values = c(
      "Both" = "#1B7837",
      "Protein only" = "#2166AC",
      "EpiScore only" = "#D8A48F",
      "Neither" = "grey75"
    ),
    breaks = c("Both", "Protein only", "EpiScore only", "Neither"),
    drop = FALSE
  ) +
  coord_cartesian(
    xlim = c(-1, 1.3),
    ylim = c(-1, 1.3)
  ) +
  labs(
    title = "Disease-specific associations",
    subtitle = paste0(
      "Measured protein vs DNAm EpiScore\n",
      "Fully adjusted models (Model D vs Model E); n = ",
      n_diseases,
      " diseases"
    ),
    x = "log(HR) per SD increase — Measured GDF15",
    y = "log(HR) per SD increase — GDF15 EpiScore",
    colour = "Bonferroni significance"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11.5, hjust = 0.5, lineheight = 1.1),
    legend.position = "top",
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

print(fig_master_scatter)

ggsave(
  filename = file.path(DIR_FIG, "Figure_MASTER_scatter_Bonferroni_v2.png"),
  plot = fig_master_scatter,
  width = 9.5,
  height = 8.5,
  dpi = 300
)

ggsave(
  filename = file.path(DIR_FIG, "Figure_MASTER_scatter_Bonferroni_v2.pdf"),
  plot = fig_master_scatter,
  width = 9.5,
  height = 8.5
)

### =========================
### FIGURE 2: Violin plot
### =========================

plot_violin <- de_concordance %>%
  select(disease, disease_label, logHR_protein, logHR_episcore) %>%
  pivot_longer(
    cols = c(logHR_protein, logHR_episcore),
    names_to = "biomarker",
    values_to = "logHR"
  ) %>%
  mutate(
    biomarker = recode(
      biomarker,
      logHR_protein = "Measured GDF15 protein",
      logHR_episcore = "DNAm GDF15 EpiScore"
    ),
    biomarker = factor(
      biomarker,
      levels = c("Measured GDF15 protein", "DNAm GDF15 EpiScore")
    )
  )

median_violin <- plot_violin %>%
  group_by(biomarker) %>%
  summarise(
    median_logHR = median(logHR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    median_label = paste0("Median = ", round(median_logHR, 2))
  )

fig_violin_refined <- ggplot(
  plot_violin,
  aes(x = biomarker, y = logHR, fill = biomarker)
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.55,
    colour = "grey35",
    linewidth = 0.4
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    alpha = 0.85,
    colour = "grey20",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.5
  ) +
  geom_text(
    data = median_violin,
    aes(
      x = biomarker,
      y = median_logHR + 0.12,
      label = median_label
    ),
    inherit.aes = FALSE,
    size = 3.7,
    colour = "grey20"
  ) +
  scale_fill_manual(
    values = c(
      "Measured GDF15 protein" = "#8FB9A8",
      "DNAm GDF15 EpiScore" = "#D8A48F"
    )
  ) +
  coord_cartesian(
    ylim = c(-0.8, 1.25)
  ) +
  labs(
    title = "Distribution of disease-specific GDF15 associations",
    subtitle = paste0(
      "Fully adjusted models (Model D vs Model E); n = ",
      n_diseases,
      " diseases"
    ),
    x = NULL,
    y = "log(HR) per SD increase"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11.5, hjust = 0.5),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11.5),
    axis.text.y = element_text(size = 10.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

print(fig_violin_refined)

ggsave(
  filename = file.path(DIR_FIG, "Figure_MASTER_violin_v2.png"),
  plot = fig_violin_refined,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(DIR_FIG, "Figure_MASTER_violin_v2.pdf"),
  plot = fig_violin_refined,
  width = 7,
  height = 5
)

### =========================
### FIGURE 3: Combined concordant/discordant dumbbell plot
### =========================

top_concordant_plot <- concordant_ranked %>%
  slice(1:8) %>%
  arrange(abs_logHR_diff) %>%
  mutate(plot_group = "Top concordant")

top_discordant_plot <- discordant_ranked %>%
  slice(1:8) %>%
  arrange(desc(abs_logHR_diff)) %>%
  mutate(plot_group = "Top discordant")

combined_dumbbell <- bind_rows(
  top_concordant_plot,
  top_discordant_plot
) %>%
  mutate(
    disease_label_clean = stringr::str_replace_all(disease_label, "_", " "),
    bonf_protein = p_value_protein < bonf_threshold,
    bonf_episcore = p_value_episcore < bonf_threshold,
    bonf_label = case_when(
      bonf_protein & bonf_episcore ~ "Both",
      bonf_protein & !bonf_episcore ~ "Protein",
      !bonf_protein & bonf_episcore ~ "EpiScore",
      TRUE ~ "Neither"
    )
  )

desired_order <- combined_dumbbell %>%
  arrange(
    plot_group,
    if_else(plot_group == "Top concordant", abs_logHR_diff, -abs_logHR_diff)
  ) %>%
  pull(disease_label_clean)

combined_dumbbell <- combined_dumbbell %>%
  mutate(
    disease_label_clean = factor(
      disease_label_clean,
      levels = rev(desired_order)
    )
  )

combined_points <- combined_dumbbell %>%
  select(
    disease_label_clean,
    plot_group,
    logHR_protein,
    logHR_episcore
  ) %>%
  pivot_longer(
    cols = c(logHR_protein, logHR_episcore),
    names_to = "biomarker",
    values_to = "logHR"
  ) %>%
  mutate(
    biomarker = recode(
      biomarker,
      logHR_protein = "Measured GDF15 protein",
      logHR_episcore = "DNAm GDF15 EpiScore"
    ),
    biomarker = factor(
      biomarker,
      levels = c("Measured GDF15 protein", "DNAm GDF15 EpiScore")
    )
  )

x_sig_label <- max(
  combined_dumbbell$logHR_protein,
  combined_dumbbell$logHR_episcore,
  na.rm = TRUE
) + 0.08

fig_dumbbell_combined <- ggplot(combined_dumbbell) +
  geom_segment(
    aes(
      x = logHR_protein,
      xend = logHR_episcore,
      y = disease_label_clean,
      yend = disease_label_clean
    ),
    colour = "grey65",
    linewidth = 0.8
  ) +
  geom_point(
    data = combined_points,
    aes(
      x = logHR,
      y = disease_label_clean,
      colour = biomarker
    ),
    size = 3.2,
    alpha = 0.95
  ) +
  geom_text(
    aes(
      x = x_sig_label,
      y = disease_label_clean,
      label = bonf_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.1,
    colour = "grey30"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.5
  ) +
  scale_colour_manual(
    values = c(
      "Measured GDF15 protein" = "#8FB9A8",
      "DNAm GDF15 EpiScore" = "#D8A48F"
    )
  ) +
  coord_cartesian(
    xlim = c(-0.65, x_sig_label + 0.18),
    clip = "off"
  ) +
  labs(
    title = "Concordant and discordant disease associations",
    subtitle = "Top 8 concordant/discordant diseases; fully adjusted models",
    x = "log(HR) per SD increase",
    y = NULL,
    colour = NULL,
    caption = "Right-side labels indicate Bonferroni significance category."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11.5, hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(size = 9, hjust = 0),
    plot.margin = margin(10, 35, 10, 10)
  )

print(fig_dumbbell_combined)

ggsave(
  filename = file.path(DIR_FIG, "Figure_DUMBBELL_concordant_discordant_v2.png"),
  plot = fig_dumbbell_combined,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  filename = file.path(DIR_FIG, "Figure_DUMBBELL_concordant_discordant_v2.pdf"),
  plot = fig_dumbbell_combined,
  width = 8.5,
  height = 6.5
)

### =========================
### FIGURE 4: Category concordance plot
### =========================

category_concordance <- de_concordance %>%
  group_by(category, sig_concordance_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(category) %>%
  mutate(
    total = sum(n),
    prop = n / total
  ) %>%
  ungroup()

category_totals <- category_concordance %>%
  distinct(category, total) %>%
  mutate(category_label = paste0(category, " (n=", total, ")"))

category_concordance_plot <- category_concordance %>%
  left_join(category_totals, by = c("category", "total"))

category_order <- category_concordance_plot %>%
  filter(sig_concordance_group == "both") %>%
  select(category_label, prop) %>%
  right_join(
    category_concordance_plot %>% distinct(category_label),
    by = "category_label"
  ) %>%
  mutate(prop = ifelse(is.na(prop), 0, prop)) %>%
  arrange(prop) %>%
  pull(category_label)

category_concordance_plot <- category_concordance_plot %>%
  mutate(
    category_label = factor(category_label, levels = category_order),
    sig_concordance_group = factor(
      sig_concordance_group,
      levels = c("both", "protein_only", "episcore_only", "neither"),
      labels = c("Both", "Protein only", "EpiScore only", "Neither")
    )
  )

fig_category_concordance <- ggplot(
  category_concordance_plot,
  aes(x = category_label, y = prop, fill = sig_concordance_group)
) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
      "Both" = "#7FB77E",
      "Protein only" = "#8FB9A8",
      "EpiScore only" = "#D8A48F",
      "Neither" = "grey75"
    )
  ) +
  coord_flip() +
  labs(
    title = "Concordance patterns across disease categories",
    subtitle = "Fully adjusted models (Model D vs Model E)",
    x = NULL,
    y = "Proportion of diseases",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11.5, hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9.5)
  )

print(fig_category_concordance)

ggsave(
  filename = file.path(DIR_FIG, "Figure_MASTER_concordance_both_ordered_v2.png"),
  plot = fig_category_concordance,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

ggsave(
  filename = file.path(DIR_FIG, "Figure_MASTER_concordance_both_ordered_v2.pdf"),
  plot = fig_category_concordance,
  width = 8.5,
  height = 6.5
)