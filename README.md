# Protein-GDF15-vs-EpiScore
MSc (YM) - GDF15 vs EpiScore in 173 Diseases in GS

# Protein GDF15 versus its Epigenetic Surrogate: A 173-Outcome Disease Prediction Study in Generation Scotland

## Overview
This repository contains the R scripts used for the analyses presented in the MSc dissertation:

**Dr Low Yi Mei**
*MSc Data Science for Health and Social Care*
*University of Edinburgh*

The study evaluates the extent to which a DNA methylation-derived GDF15 EpiScore recapitulates disease-association patterns observed for measured circulating GDF15 across a phenome-wide set of disease outcomes in Generation Scotland.

A total of 173 incident disease outcomes were analysed using Cox proportional hazards models, with comparisons performed between measured GDF15 and its epigenetic surrogate.

---

## Analytical Workflow

```text
Raw Generation Scotland Data
        ↓
Script1_Data_Prep.R
        ↓
Analysis Dataset
        ↓
Script2_Disease_Metadata.R
        ↓
Disease Metadata Tables
        ↓
Script3_Cox_Disease_Analysis.R
        ↓
Phenome-wide Cox Results
        ↓
 ┌───────────────┬───────────────┐
 ↓                               ↓
Script4_Figures_and_Plots.R   Script5_Results_Tables.R
 ↓                               ↓
Figures                       Tables
```

---

## Scripts

### Script1_Data_Prep.R

Prepares the analytical dataset used for analyses.
- Merges participant-level datasets
- Excludes prevalent disease cases
- Generates disease-specific follow-up times
- Applies a 10-year follow-up window
- Performs censoring at death or administrative end of follow-up
- Cleans and derives covariates
- Produces the final analytical dataset

---

### Script2_Disease_Metadata.R

Generates disease-level metadata and descriptive summaries.
* Supplementary Table S1: Disease Metadata and Eligibility Assessment
* Supplementary Table S2: Descriptive Characteristics of Disease Outcomes

- Calculates disease frequencies
- Identifies prevalent and incident cases
- Assesses outcome eligibility for analysis
- Records exclusion reasons
- Produces descriptive statistics for included outcomes

---

### Script3_Cox_Disease_Analysis.R

Phenome-wide Cox proportional hazards analyses.
- Analyses 173 disease outcomes
- Requires ≥30 incident cases for inclusion
- Applies sex-specific disease restrictions where appropriate
- Fits six Cox regression models

Models:

* Model A: Protein GDF15 (age and sex adjusted)
* Model B: DNAm GDF15 EpiScore (age and sex adjusted)
* Model C: Joint model (age and sex adjusted)
* Model D: Protein GDF15 (fully adjusted)
* Model E: DNAm GDF15 EpiScore (fully adjusted)
* Model F: Joint model (fully adjusted)

---

### Script4_Figures_and_Plots.R

Generates all figures presented in the dissertation.
- Master scatter plot
- Concordance and discordance analyses
- Dumbbell plots
- Violin plots
- Disease-group summary figures

---

### Script5_Results_Tables.R

Generates dissertation tables
- Table 1: Baseline Characteristics
- Table 2: Association Concordance Categories
- Table 3: Disease Group Summaries
- Additional summary tables used in the dissertation

---

## Software
Analyses were performed using R (version 4.4.x).
Key packages include:
* survival
* tidyverse
* dplyr
* ggplot2
* cowplot
* readr
* stringr

---

## Data Availability
Generation Scotland data are available through application to Generation Scotland

---
## Author

**Dr Low Yi Mei**
Senior Hospital Clinician
Alexandra Hospital, National University Health System (NUHS), Singapore

MSc Data Science for Health and Social Care
University of Edinburgh
