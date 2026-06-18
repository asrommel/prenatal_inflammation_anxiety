# Prenatal inflammation and postpartum anxiety analysis

This repository contains the R code used to produce the analyses for the manuscript:

**"Inflammatory Predictors of Postpartum Anxiety in a Prospective Pregnancy Cohort in New York City."**

## Overview

This analysis examines associations between third-trimester inflammatory markers (IL-6, IL-17A, IL-1β, CRP) and postpartum anxiety symptoms (GAD-7) using quantile regression models in a prospective pregnancy cohort of 237 participants.

## Contents

```
├── README.md                          # This file
├── data_cleaning_and_analysis.R       # Main analysis script
├── data/                              # Data directory (not included; see below)
│   ├── README_DATA.md                 # Data dictionary
│   └── raw/                           # Raw data files
├── output/                            # Generated tables and figures
└── .gitignore
```

## Key Results

- **Sample size:** n = 237 (from 670 eligible participants with GAD-7 data)
- **Primary outcome:** Postpartum GAD-7 score (range: 0-21)
- **Exposures:** Log-transformed inflammatory markers from third-trimester serum samples
- **Analysis method:** Quantile regression at 50th, 75th, and 90th percentiles

**Main outputs:**
- Figure 1: Study sample flowchart
- Table 1: Participant characteristics by GAD-7 percentile
- Figure 2: Forest plot of adjusted associations
- Supplemental Tables S1-S3: Full regression results and sensitivity analyses

## Software Requirements

### R Version
- R ≥ 4.0.0 (tested with R 4.3.1)

### Required Packages
See Section 0 of `data_cleaning_and_analysis.R` for complete package list:

```r
library(Hmisc)           # Data labels
library(labelled)        # Label management
library(dplyr)           # Data manipulation
library(tidyverse)       # Data wrangling
library(janitor)         # Data cleaning
library(mice)            # Multiple imputation (50 datasets)
library(quantreg)        # Quantile regression
library(purrr)           # Functional programming
library(tidyr)           # Pivot operations
library(gt)              # Table formatting
library(gtsummary)       # Summary tables
library(officer)         # Word document output
library(flextable)       # Flexible tables
library(ggplot2)         # Plotting
library(DiagrammeR)      # Flowchart
library(DiagrammeRsvg)   # SVG export
library(rsvg)            # PDF export
```

### Installation

```r
# Install required packages
packages <- c("Hmisc", "labelled", "dplyr", "tidyverse", "janitor", "mice", 
              "quantreg", "purrr", "tidyr", "gt", "gtsummary", "officer", 
              "flextable", "ggplot2", "DiagrammeR", "DiagrammeRsvg", "rsvg")

install.packages(packages)
```

**Note:** `rsvg` may require system-level dependencies on Linux/Mac. See [rsvg installation guide](https://github.com/jeroen/rsvg).

## Data Availability

The datasets generated and/or analyzed during the current study are **not publicly available** because they contain potentially identifiable participant information.

### Accessing the Data

Data may be available from the corresponding author upon reasonable request, subject to:
- Institutional Review Board (IRB) approval
- Data Use Agreement (DUA) execution
- Confirmation of authorized use

**Requests should be directed to:** [corresponding author contact info]

### Data Requirements

The analysis requires three CSV files (locations specified at lines 20-22 of the script):

1. **Mental health history data** (`Mental health hx 2024-6-6_all variables.csv`)
   - Subject ID, depression history, anxiety history, etc.

2. **Main study dataset** (`GenC20-CarlyGADAndImmuneAct_DATA_2024-07-15_1505_FM.csv`)
   - Demographic variables, inflammatory markers (IL-6, IL-17A, IL-1β, CRP)
   - GAD-7 symptom items and derived scores
   - Pregnancy and delivery outcomes

3. **Sample logbook** (`Logbook 4.20.21_ML2 1.20.23 Frederieke_[55].xlsx - Sheet 1...`)
   - Sample collection and processing dates/times

### Data Structure

Key variables in the analysis dataset:

**Demographics:**
- `subject_id`: Unique participant identifier
- `maternalage`: Age at enrollment (years)
- `raceethnicitycombined`: Race/ethnicity (coded 1-8)
- `para`: Parity (number of prior pregnancies)
- `prepregnancybmi`: Pre-pregnancy BMI (kg/m²)

**Inflammatory markers** (collected at 6 timepoints):
- `il6_1` through `il6_6`: IL-6 (pg/mL)
- `il1b_1` through `il1b_6`: IL-1β (pg/mL)
- `il17a_1` through `il17a_6`: IL-17A (pg/mL)
- `crp_1` through `crp_6`: C-reactive protein (pg/mL)
- `date_sample_1` through `date_sample_6`: Sample collection dates

**Anxiety assessment:**
- `p3_2w_nerv`, `p3_2w_uncontrworry`, etc.: Individual GAD-7 items (0-3 scale)
- `p3_2w_GAD.score`: Total GAD-7 score (0-21)
- `p3_timestamp`: Date of GAD-7 completion

## Reproduction Guide

### Step 1: Set Up File Structure

```bash
git clone <repository-url>
cd prenatal-inflammation-anxiety
mkdir -p data/raw output
```

### Step 2: Obtain Data

Request access from corresponding author and place the three CSV files in `data/raw/`:

```
data/raw/
├── Mental health hx 2024-6-6_all variables.csv
├── GenC20-CarlyGADAndImmuneAct_DATA_2024-07-15_1505_FM.csv
└── Logbook 4.20.21_ML2 1.20.23 Frederieke_[55].xlsx - Sheet 1 - ...
```

### Step 3: Update File Paths

Edit lines 20-22 in `data_cleaning_and_analysis.R` to point to your data directory:

```r
# Current (will not work on your machine):
setwd("/Users/rommea01/Dropbox/Carly/")
mental_health_raw <- read.csv("/Users/rommea01/Dropbox/Carly/...")

# Change to:
data_dir <- "data/raw"
mental_health_raw <- read.csv(file.path(data_dir, "Mental health hx 2024-6-6_all variables.csv"))
# ... etc
```

### Step 4: Run the Analysis

```r
# Open RStudio or R console
setwd("~/path/to/prenatal-inflammation-anxiety")
source("data_cleaning_and_analysis.R")
```

**Expected runtime:** ~5-10 minutes (depending on system; multiple imputation is computationally intensive)

### Step 5: Check Outputs

Generated files will be saved to the working directory:

```
Figure_1_Flowchart_of_Study_Sample.pdf
Figure_2.forest_plot_gad7.pdf
Table_1_Participant_Characteristics_by_GAD7_Percentile_Group.html
Table_2.docx
Supplemental_Figure_S1_Timing_Distributions.pdf
Supplemental_Figure_S2_Collinearity.pdf
Supplemental_Table_S1_Full_Multivariable_Quantile_Regression_Results.docx
Supplemental_Table_S2_Complete_Case_Quantile_Regression_Results.docx
Supplementary_Table_S3_Imputed_Quantile_Regression_12_Weeks.docx
```

## Analysis Details

### Study Sample

**Inclusion criteria:**
- Enrolled in the Generation C prospective pregnancy cohort
- Completed postpartum survey
- Completed GAD-7 assessment ≤24 weeks postpartum
- Had ≥1 inflammatory marker measurement

**Exclusions:**
- GAD-7 completed during pregnancy (n=13)
- No third-trimester sample (≥28 weeks gestation, >7 days before delivery) (n=205)
- Missing all cytokine values (n=21)

**Final sample:** n = 237 unique participants with 237 observations

### Primary Analysis

**Exposure:** Log2-transformed third-trimester inflammatory markers (IL-6, IL-17A, IL-1β, CRP)

**Outcome:** Postpartum GAD-7 score (continuous, 0-21)

**Covariates:**
- Maternal age (continuous)
- Race/ethnicity (categorical)
- Parity (binary: nulliparous vs. multiparous)
- Education (3 categories)
- SARS-CoV-2 infection during pregnancy (binary)
- Days postpartum at GAD-7 completion (continuous)
- Days since pandemic onset (continuous, COVID-19 context)
- Pre-pregnancy BMI (continuous)
- History of anxiety/depression (binary)

**Method:** Quantile regression at 50th, 75th, and 90th percentiles of GAD-7 scores

**Imputation:** Multiple imputation by chained equations (MICE) with 50 datasets using predictive mean matching (PMM)

**Multiple testing correction:** Benjamini-Hochberg FDR adjustment applied to 12 primary exposure tests (4 markers × 3 percentiles)

### Sensitivity Analyses

1. **Complete-case analysis** (n=220): Excludes participants with missing covariates
2. **Restricted timing** (n=187): GAD-7 completed ≤12 weeks postpartum

## Key Decisions and Assumptions

- **Missing data handling:** MICE with PMM was chosen due to missing education data (3.8% of analytic sample)
- **Log transformation:** Applied to all inflammatory markers due to right skewness
- **Out-of-range values:** Replaced with assay limit of detection (IL-6: 0.01, IL-1β: 0.03, IL-17A: 0.64 pg/mL)
- **CRP missing values:** Imputed to maximum measured value (195.69 pg/mL) when other markers present
- **Random seed:** Set to 42 for reproducibility of imputation
- **Quantile selection:** 50th, 75th, 90th percentiles chosen to examine associations across anxiety symptom distribution
- **Third-trimester timing:** Last sample ≥28 weeks gestation, >7 days before delivery to capture late pregnancy immune state

## Interpretation Notes

- **Positive β:** Higher inflammatory marker level associated with higher GAD-7 score
- **Effect size:** β represents change in GAD-7 score per log2 increase in marker
  - Example: IL-6 β=0.50 at 90th percentile means for each doubling of IL-6, GAD-7 increases by 0.50 points
- **Quantile-specific:** Associations may differ across anxiety symptom distribution; stronger effects at higher percentiles would suggest inflammation more predictive of severe anxiety

## Citation

If you use this code or reproduce these analyses, please cite:

[Insert full manuscript citation once published]

## Contact

**Corresponding author:** [Name, email, institution]

**Questions about the code:** [Repository maintainer contact info]

## License

[Specify license: MIT, CC-BY-4.0, etc.]

## Funding & Acknowledgments

[Funding sources, IRB approval number, acknowledgments]

---

## Troubleshooting

### Common Issues

**Q: "File not found" error when running the script**
- A: Update the file paths in lines 20-22 to match your data directory location

**Q: Error installing `rsvg` package**
- A: `rsvg` requires system dependencies. See [installation guide](https://github.com/jeroen/rsvg#installation) for your OS

**Q: Multiple imputation produces different results each time I run it**
- A: This is normal due to randomness in MICE. The script sets `set.seed(42)` on lines 783 and 1467 to ensure reproducibility; results should match exactly with this seed

**Q: Analysis takes >30 minutes to run**
- A: Multiple imputation (50 datasets) × 4 markers × 3 quantiles = 600 quantile regressions. This is computationally intensive. Consider reducing to fewer imputed datasets (line 783: `m = 20`) for testing

**Q: Tables saved to wrong directory**
- A: Update `setwd()` on line 9 to your desired output directory, or modify output file paths to include directory name

## Version History

- **2024-07:** Initial release (analysis completed July 2024)
- **Date:** Analysis completed [specific date]
- **Last updated:** [date]

---

**Last updated:** [Current date]
