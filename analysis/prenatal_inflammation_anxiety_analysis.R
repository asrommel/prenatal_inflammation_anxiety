# ==============================================================================
# Inflammatory Predictors of Postpartum Anxiety
# Clean analysis script aligned with manuscript
#
# Main manuscript outputs:
#   Figure 1. Flowchart of study sample
#   Table 1. Participant characteristics by GAD-7 score percentile group
#   Figure 2. Adjusted associations with postpartum GAD-7 scores
#
# Supplementary outputs:
#   Supplemental Figure S1. Cytokine sample timing across pregnancy
#   Supplemental Figure S2. Collinearity among continuous covariates
#   Supplemental Table S1. Full multivariable quantile regression results
#   Supplemental Table S2. Complete-case quantile regression results
#   Supplemental Table S3. Quantile regression results restricted
#                         to participants completing the GAD-7 within 12 weeks postpartum
#
# Analytic sample in manuscript:
#   n = 237
#   complete-case sensitivity sample: n = 220
#   <=12 week sensitivity sample: n = 187
# ==============================================================================

rm(list = ls())
graphics.off()
setwd("/Users/rommea01/Dropbox/Carly/")

# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(Hmisc)
library(labelled)
library(dplyr)
library(tidyverse)
library(janitor)
library(mice)
library(quantreg)
library(purrr)
library(tidyr)
library(gt)
library(gtsummary)
library(officer)
library(flextable)
library(ggplot2)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

mental_health_raw <- read.csv("/Users/rommea01/Dropbox/Carly/Mental health hx 2024-6-6_all variables.csv")
data <- read.csv("/Users/rommea01/Dropbox/Carly/00 Updated Analysis/00 Feb 2026/Re_ Paper suggestions/GenC20-CarlyGADAndImmuneAct_DATA_2024-07-15_1505_FM.csv")
logbook <- read.csv("/Users/rommea01/Dropbox/Carly/Merging Data with Logbook/Logbook 4.20.21_ML2 1.20.23 Frederieke_[55].xlsx - Sheet 1 - Logbook 4.20.21_ML2 1 (1).csv")

# ==============================================================================
# 2. VARIABLE LABELS
# ==============================================================================

label(data$subject_id)            <- "Subject ID"
label(data$maternalage)           <- "Maternal Age (years) at enrollment"
label(data$raceethnicitycombined) <- "Race/Ethnicity"
label(data$insurance_cat)         <- "Insurance Category"
label(data$para)                  <- "Para"
label(data$preg_pos)              <- "SARS-CoV-2 infection during pregnancy"
label(data$vaccine_timing)        <- "Vaccination initiation timing (Dose 1)"
label(data$birthdate)             <- "Birth Date"
label(data$birthweight_grams)     <- "Birthweight (grams)"
label(data$gestationalagedays)    <- "Gestational Age at Delivery (days)"
label(data$prepregnancybmi)       <- "Pre-pregnancy BMI"

for (tp in c("p3", "s2", "s3")) {
  label(data[[paste0(tp, "_2w_nerv")]])         <- "Feeling nervous, anxious or on edge"
  label(data[[paste0(tp, "_2w_uncontrworry")]]) <- "Not being able to stop or control worrying"
  label(data[[paste0(tp, "_2w_pleasur")]])      <- "Little interest or pleasure in doing things"
  label(data[[paste0(tp, "_2w_depr")]])         <- "Feeling down, depressed, or hopeless"
  label(data[[paste0(tp, "_2w_toomuchworry")]]) <- "Worrying too much about different things"
  label(data[[paste0(tp, "_2w_relax")]])        <- "Trouble relaxing"
  label(data[[paste0(tp, "_2w_restless")]])     <- "Being so restless that it's hard to sit still"
  label(data[[paste0(tp, "_2w_annoy")]])        <- "Becoming easily annoyed or irritable"
  label(data[[paste0(tp, "_2w_afraid")]])       <- "Feeling afraid as if something awful might happen"
}

for (i in 1:6) {
  label(data[[paste0("il6_", i)]])         <- "IL-6"
  label(data[[paste0("il17a_", i)]])       <- "IL-17A"
  label(data[[paste0("il1b_", i)]])        <- "IL-1b"
  label(data[[paste0("crp_", i)]])         <- "CRP"
  label(data[[paste0("date_sample_", i)]]) <- paste("Date sample", i)
  label(data[[paste0("processtime_", i)]]) <- paste("Sample", i, "processing time")
}

for (tp in c("p1", "s1", "s3")) {
  label(data[[paste0(tp, "_income")]])   <- "Total household income (last year, before taxes)"
  label(data[[paste0(tp, "_mom_educ")]]) <- "Highest level of education completed"
}

# ==============================================================================
# 3. FACTOR VARIABLES
# ==============================================================================

data$raceethnicitycombined.factor <- factor(
  data$raceethnicitycombined,
  levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
  labels = c(
    "American Indian or Alaska Native",
    "Asian",
    "Black",
    "Hispanic",
    "Native Hawaiian or Other Pacific Islander",
    "White",
    "Other",
    "Unknown"
  )
)

data$insurance_cat.factor <- factor(
  data$insurance_cat,
  levels = c("1", "2", "3"),
  labels = c("Private", "Public", "Self-pay")
)

data$preg_pos.factor <- factor(
  data$preg_pos,
  levels = c("1", "0"),
  labels = c("Yes", "No")
)

data$vaccine_timing.factor <- factor(
  data$vaccine_timing,
  levels = c("1", "2", "3"),
  labels = c("Before pregnancy", "During pregnancy", "After pregnancy")
)

data$duplicate_enrollment.factor <- factor(
  data$duplicate_enrollment,
  levels = c("1", "0"),
  labels = c("Yes", "No")
)

income_levels <- c(
  "Less than $25,000",
  "$25,001 - $50,000",
  "$50,001 - $75,000",
  "$75,001 - $100,000",
  "$100,001 - $125,000",
  "More than $125,000"
)

educ_levels <- c(
  "Less than 9th grade",
  "Some high school; no diploma",
  "High school diploma/GED",
  "Trade/Technical/Vocational Training",
  "Some college",
  "Bachelors degree",
  "Post-graduate degree"
)

for (tp in c("p1", "s1", "s3")) {
  data[[paste0(tp, "_income.factor")]] <-
    factor(data[[paste0(tp, "_income")]], levels = 1:6, labels = income_levels)
  
  data[[paste0(tp, "_mom_educ.factor")]] <-
    factor(data[[paste0(tp, "_mom_educ")]], levels = 1:7, labels = educ_levels)
}

# Treat NA preg_pos as "No"
data$preg_pos.factor[is.na(data$preg_pos.factor)] <- "No"
data$preg_pos[is.na(data$preg_pos)] <- 0

# ==============================================================================
# 4. CORRECT s2 SURVEY SCORING
# ==============================================================================

s2_items <- paste0(
  "s2_2w_",
  c("nerv", "uncontrworry", "pleasur", "depr", "toomuchworry", "relax", "restless", "annoy", "afraid")
)

data[, s2_items] <- data[, s2_items] - 1

# ==============================================================================
# 5. GAD-7 SCORES
# ==============================================================================

gad_items <- c("nerv", "uncontrworry", "toomuchworry", "relax", "restless", "annoy", "afraid")
phq_items <- c("depr", "pleasur")

for (tp in c("p3", "s2", "s3")) {
  data[[paste0(tp, "_2w_GAD.score")]] <-
    rowSums(data[, paste0(tp, "_2w_", gad_items)], na.rm = FALSE)
  
  data[[paste0(tp, "_2w_PHQ.score")]] <-
    rowSums(data[, paste0(tp, "_2w_", phq_items)], na.rm = FALSE)
  
  data[[paste0(tp, "_2w_GAD.score")]] <-
    unclass(remove_labels(data[[paste0(tp, "_2w_GAD.score")]], user_na_to_na = TRUE))
}

data$GAD.score <- coalesce(data$p3_2w_GAD.score, data$s2_2w_GAD.score, data$s3_2w_GAD.score)
data$PHQ.score <- coalesce(data$p3_2w_PHQ.score, data$s2_2w_PHQ.score, data$s3_2w_PHQ.score)

# Remove duplicate GAD completion
data$GAD.score[data$subject_id == 1559] <- NA

data <- data %>%
  mutate(
    anxiety.factor = case_when(
      GAD.score < 5 ~ "Minimal Anxiety",
      GAD.score >= 5 & GAD.score < 10 ~ "Mild Anxiety",
      GAD.score >= 10 & GAD.score < 15 ~ "Moderate Anxiety",
      GAD.score >= 15 ~ "Severe Anxiety"
    ),
    para_factor = case_when(
      para == 0 ~ "0",
      para == 1 ~ "1",
      para == 2 ~ "2",
      para == 3 ~ "3",
      para >= 4 ~ "4+"
    ),
    para_binary = factor(
      ifelse(para == 0, "Nulliparous", "Multiparous"),
      levels = c("Nulliparous", "Multiparous")
    )
  )

# ==============================================================================
# 6. HANDLE OUT-OF-RANGE CYTOKINE VALUES
# ==============================================================================

for (i in 1:6) {
  data[[paste0("sample.empty_", i)]] <-
    is.na(as.numeric(data[[paste0("il6_", i)]])) &
    is.na(as.numeric(data[[paste0("il1b_", i)]])) &
    is.na(as.numeric(data[[paste0("il17a_", i)]])) &
    is.na(as.numeric(data[[paste0("crp_", i)]]))
}

replace_oor <- function(value, all_empty, limit) {
  ifelse(!all_empty & value %in% c("OOR", "0", "OOR <"), limit, value)
}

detect_limits <- list(il6 = "0.01", il1b = "0.03", il17a = "0.64")

for (marker in names(detect_limits)) {
  for (i in 1:6) {
    data[[paste0(marker, "_", i)]] <-
      mapply(
        replace_oor,
        data[[paste0(marker, "_", i)]],
        data[[paste0("sample.empty_", i)]],
        detect_limits[[marker]]
      )
  }
}

for (i in 1:6) {
  crp_col   <- paste0("crp_", i)
  empty_col <- paste0("sample.empty_", i)
  data[[crp_col]] <- as.numeric(data[[crp_col]])
  data[[crp_col]][!data[[empty_col]] & is.na(data[[crp_col]])] <- 195.69
}

# ==============================================================================
# 7. CONVERT CYTOKINES AND DATES
# ==============================================================================

cyto_cols <- c(
  paste0("il6_", 1:6),
  paste0("il17a_", 1:6),
  paste0("il1b_", 1:6),
  paste0("crp_", 1:6)
)

data[, cyto_cols] <- lapply(data[, cyto_cols], as.numeric)

data$birthdate <- as.Date(data$birthdate)

for (i in 1:6) {
  data[[paste0("date_sample_", i)]] <- as.Date(data[[paste0("date_sample_", i)]])
}

data$lmp_estimate <- data$birthdate - as.numeric(data$gestationalagedays)

# ==============================================================================
# 8. KEEP PARTICIPANTS WITH GAD + AT LEAST ONE CYTOKINE
# ==============================================================================

has_cyto <- Reduce(`|`, lapply(cyto_cols, function(col) !is.na(data[[col]])))
subset <- data[!is.na(data$GAD.score) & has_cyto, ]

# ==============================================================================
# 9. SURVEY TIMESTAMPS AND GAD TIMING
# ==============================================================================

subset$s3_timestamp[subset$subject_id == 476] <- "2020-12-08"

subset$mom_educ.factor <- coalesce(
  subset$p1_mom_educ.factor,
  subset$s1_mom_educ.factor,
  subset$s3_mom_educ.factor
)

subset$p3_timestamp <- as.Date(subset$p3_timestamp)
subset$s2_timestamp <- as.Date(subset$s2_timestamp)
subset$s3_timestamp <- as.Date(subset$s3_timestamp)

subset$timestamp_GAD <- coalesce(subset$p3_timestamp, subset$s2_timestamp, subset$s3_timestamp)
subset$GAD_days_postpart <- subset$timestamp_GAD - as.Date(subset$birthdate)
subset$days_since_pandemic <- subset$timestamp_GAD - as.Date("2020-03-01")

# Remove participant with GAD completed during pregnancy
subset <- subset[!(subset$GAD_days_postpart < 1 & !is.na(subset$GAD_days_postpart)), ]

for (i in 1:6) {
  subset[[paste0("sample_", i, "_gest_age")]] <-
    subset[[paste0("date_sample_", i)]] - subset$lmp_estimate
}

# ==============================================================================
# 10. COLLAPSE RACE/ETHNICITY AND EDUCATION
# ==============================================================================

subset <- subset %>%
  mutate(
    raceethnicitycombined.factor = fct_collapse(
      raceethnicitycombined.factor,
      "Other" = c(
        "American Indian or Alaska Native",
        "Native Hawaiian or Other Pacific Islander",
        "Other"
      )
    ),
    mom_educ.factor = fct_collapse(
      mom_educ.factor,
      "< High school diploma" = c("Less than 9th grade", "Some high school; no diploma"),
      "≤ Some college" = c("High school diploma/GED", "Trade/Technical/Vocational Training", "Some college"),
      "≥ Bachelors" = c("Bachelors degree", "Post-graduate degree")
    )
  )

# ==============================================================================
# 11. MERGE MENTAL HEALTH HISTORY
# ==============================================================================

subset$subject_id <- as.integer(remove_labels(subset$subject_id))

mhhx <- mental_health_raw %>%
  select(subject_id, depression, anxiety, depranx, affective, anymental, depressionandppd) %>%
  mutate(subject_id = as.integer(subject_id))

subset <- left_join(subset, mhhx, by = "subject_id")

if ("X" %in% names(subset)) subset <- select(subset, -X)

# ==============================================================================
# 12. LOG2-TRANSFORM CYTOKINES
# ==============================================================================

for (marker in c("il6", "il17a", "il1b", "crp")) {
  for (i in 1:6) {
    raw <- paste0(marker, "_", i)
    subset[[paste0(raw, ".log")]] <- log(subset[[raw]], 2)
    subset[[raw]] <- remove_labels(subset[[raw]])
  }
}

# ==============================================================================
# 13. RESTRICT GAD TO <=24 WEEKS POSTPARTUM
# ==============================================================================

subset$GAD_wks_postpart <- as.numeric(floor(subset$GAD_days_postpart / 7))
subset <- subset[subset$GAD_wks_postpart <= 24, ]

cat("N after GAD timing restriction:", nrow(subset), "\n")

# ==============================================================================
# 14. PIVOT TO LONG FORMAT
# ==============================================================================

for (i in 1:6) {
  subset[[paste0("processtime_", i)]] <- as.character(subset[[paste0("processtime_", i)]])
}

subset.rename <- subset %>% select(-starts_with("old."))

for (i in 1:6) {
  names(subset.rename)[names(subset.rename) == paste0("sample_", i, "_gest_age")] <-
    paste0("gest_age_sample_", i)
}

for (marker in c("il6", "il17a", "il1b", "crp")) {
  for (i in 1:6) {
    names(subset.rename)[names(subset.rename) == paste0(marker, "_", i, ".log")] <-
      paste0(marker, ".log_", i)
  }
}

lmmData <- subset.rename %>%
  pivot_longer(
    matches("^(.*)_([0-9])$"),
    names_to = c(".value"),
    names_pattern = "^(.*)_[0-9]$"
  ) %>%
  rename(gest_age_sample_day = gest_age_sample) %>%
  mutate(gest_age_sample_wk = as.numeric(floor(gest_age_sample_day / 7)))

# ==============================================================================
# 15. MERGE LOGBOOK
#     KEEP THIS SECTION: logbook is used for sample-level timing information
# ==============================================================================

logbook <- logbook %>%
  clean_names() %>%
  select(-starts_with("x")) %>%
  mutate(
    pick.date  = as.Date(pickup_date, format = "%m/%d/%y"),
    drop.date  = as.Date(dropoff_date, format = "%m/%d/%y"),
    subject_id = suppressWarnings(as.numeric(subject_id))
  ) %>%
  filter(!is.na(subject_id))

lmmData$date_sample <- as.Date(lmmData$date_sample, format = "%Y-%m-%d")
lmmData$subject_id  <- as.numeric(lmmData$subject_id)

# Remove known duplicate
logbook <- logbook[-1153, ]

join <- left_join(
  lmmData,
  logbook,
  by = join_by(subject_id, closest(date_sample <= drop.date))
)

lmmData <- join %>%
  relocate(subject_id, date_sample, drop.date, pick.date, il6, il1b, il17a, crp, notes)

for (marker in c("il6", "il17a", "il1b", "crp")) {
  lmmData[[paste0(marker, ".log")]] <- log(lmmData[[marker]], 2)
}

# ==============================================================================
# 16. DEFINE PRIMARY ANALYTIC SAMPLE
#     Keep last third-trimester sample (>=28 weeks, >7 days before delivery)
# ==============================================================================

lmmData$days_before_delivery <- as.numeric(lmmData$birthdate - lmmData$date_sample)

lmmData.T3 <- lmmData %>%
  filter(
    gest_age_sample_wk >= 28,
    days_before_delivery > 7,
    if_any(c("il1b", "il6", "il17a", "crp"), is.finite)
  ) %>%
  group_by(subject_id) %>%
  slice_max(order_by = gest_age_sample_wk, n = 1, with_ties = FALSE) %>%
  ungroup()

lmmData.T3 <- lmmData.T3 %>%
  select(-any_of(c("depression", "anxiety", "depranx", "affective", "anymental", "depressionandppd"))) %>%
  left_join(mhhx, by = "subject_id")

lmmData.T3$depranx <- factor(lmmData.T3$depranx, levels = c(0, 1), labels = c("No", "Yes"))
lmmData.T3$GAD.binary <- as.integer(lmmData.T3$GAD.score >= 10)

cat("N participants with qualifying T3 sample:", n_distinct(lmmData.T3$subject_id), "\n")
cat("Total rows in analytic sample:", nrow(lmmData.T3), "\n")

# ==============================================================================
# 17. DERIVED VARIABLES FOR ANALYTIC SAMPLE
# ==============================================================================

lmmData.T3 <- lmmData.T3 %>%
  mutate(
    GAD_days_postpart = as.numeric(GAD_days_postpart),
    days_since_pandemic = as.numeric(days_since_pandemic),
    GAD_wks_postpart = GAD_days_postpart / 7
  )

# Gestational age at sampling in weeks
lmmData.T3$date_sample <- as.Date(lmmData.T3$date_sample)
lmmData.T3$birthdate <- as.Date(lmmData.T3$birthdate)
lmmData.T3$days_before_delivery <- as.numeric(lmmData.T3$birthdate - lmmData.T3$date_sample)
lmmData.T3$ga_sample_days <- lmmData.T3$gestationalagedays - lmmData.T3$days_before_delivery
lmmData.T3$ga_sample_weeks <- lmmData.T3$ga_sample_days / 7

# ==============================================================================
# 18. FLOWCHART NUMBERS FOR FIGURE 1
# ==============================================================================

total_enrolled <- nrow(data)
n_postpartum_survey <- 959
n_complete_GAD_plus_any_marker_within_24wks <- 670
n_final <- nrow(lmmData.T3)

cat("\n=== FLOWCHART COUNTS ===\n")
cat("Total enrolled:", total_enrolled, "\n")
cat("Completed postpartum survey:", n_postpartum_survey, "\n")
cat("Complete GAD-7 within 24 weeks postpartum + inflammatory marker data available:", n_complete_GAD_plus_any_marker_within_24wks, "\n")
cat("Final analytic sample:", n_final, "\n")

# ==============================================================================
# 19. MISSINGNESS CHECK
# ==============================================================================

cat("\n--- Missingness in analytic covariates ---\n")

covariate_cols <- c(
  "GAD.score", "il6.log", "il17a.log", "il1b.log", "crp.log",
  "maternalage", "raceethnicitycombined.factor", "para_binary",
  "mom_educ.factor", "preg_pos.factor", "GAD_days_postpart",
  "days_since_pandemic", "prepregnancybmi", "depranx"
)

miss_summary <- sapply(covariate_cols, function(v) {
  x <- lmmData.T3[[v]]
  n_miss <- sum(is.na(x))
  pct <- round(100 * n_miss / nrow(lmmData.T3), 1)
  c(n_missing = n_miss, pct_missing = pct)
})

print(t(miss_summary))

# ==============================================================================
# 20. DESCRIPTIVE STATISTICS FOR MANUSCRIPT TEXT
# ==============================================================================

gad_iqr <- quantile(lmmData.T3$GAD.score, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
gad_pcts <- quantile(lmmData.T3$GAD.score, probs = c(0.75, 0.90), na.rm = TRUE)
gad_range <- range(lmmData.T3$GAD.score, na.rm = TRUE)

gad_time_iqr <- quantile(lmmData.T3$GAD_wks_postpart, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
age_iqr <- quantile(lmmData.T3$maternalage, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
ga_sample_iqr <- quantile(lmmData.T3$ga_sample_weeks, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

race_tab <- prop.table(table(lmmData.T3$raceethnicitycombined.factor)) * 100
race_tab <- round(race_tab, 1)

cat("\n--- Manuscript-ready descriptive values ---\n")
cat(
  "The median postpartum GAD-7 score was", gad_iqr[2],
  "(IQR,", gad_iqr[1], "–", gad_iqr[3], "), with the 75th and 90th percentiles at",
  gad_pcts[1], "and", gad_pcts[2], ", respectively. GAD-7 scores ranged from",
  gad_range[1], "to", gad_range[2], ", with a mean of",
  round(mean(lmmData.T3$GAD.score, na.rm = TRUE), 1),
  "(SD", round(sd(lmmData.T3$GAD.score, na.rm = TRUE), 1), ").",
  "The median timing of GAD-7 completion was", round(gad_time_iqr[2], 1),
  "weeks postpartum (IQR,", round(gad_time_iqr[1], 1), "–", round(gad_time_iqr[3], 1), ").",
  "The median maternal age at enrollment was", round(age_iqr[2], 0),
  "years (IQR,", round(age_iqr[1], 0), "–", round(age_iqr[3], 0), ").",
  "Participants were predominantly White (", race_tab["White"], "%) and Hispanic (", race_tab["Hispanic"],
  "%), with smaller proportions identifying as Asian (", race_tab["Asian"],
  "%), Black (", race_tab["Black"], "%), or Other (", race_tab["Other"], "%).",
  "\n"
)

# ==============================================================================
# 21. SUPPLEMENTAL FIGURE S1. Timing distributions for cytokine sampling and postpartum anxiety assessment
# ==============================================================================

png("Supplemental_Figure_S1_Timing_Distributions.png",
    width = 2400, height = 800, res = 300)

par(mfrow = c(1,3), mar = c(4,4,3,1))

hist(
  lmmData$gest_age_sample_wk,
  breaks = 20,
  col = "#acdffbff",
  xlab = "Gestational age (weeks)",
  main = "A. Cytokine sample timing",
  xlim = c(0,42),
  xaxp = c(0,42,7)
)

hist(
  lmmData.T3$GAD_wks_postpart,
  breaks = seq(-0.5, 24.5, by = 1),
  xlim = c(-0.5, 24.5),
  xaxt = "n",
  col = "#acdffbff",
  border = "black",
  xlab = "Weeks postpartum",
  main = "B. GAD-7 completion timing"
)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

hist(
  as.numeric(lmmData.T3$days_since_pandemic),
  breaks = 20,
  col = "#acdffbff",
  xlab = "Days since March 1, 2020",
  main = "C. Timing relative to pandemic onset"
)

dev.off()

pdf("Supplemental_Figure_S1_Timing_Distributions.pdf",
    width = 12, height = 4)

par(mfrow = c(1,3), mar = c(4,4,3,1))

hist(
  lmmData$gest_age_sample_wk,
  breaks = 20,
  col = "#acdffbff",
  xlab = "Gestational age (weeks)",
  main = "A. Cytokine sample timing",
  xlim = c(0,42),
  xaxp = c(0,42,7)
)

hist(
  lmmData.T3$GAD_wks_postpart,
  breaks = seq(-0.5, 24.5, by = 1),
  xlim = c(-0.5, 24.5),
  xaxt = "n",
  col = "#acdffbff",
  border = "black",
  xlab = "Weeks postpartum",
  main = "B. GAD-7 completion timing"
)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

hist(
  as.numeric(lmmData.T3$days_since_pandemic),
  breaks = 20,
  col = "#acdffbff",
  xlab = "Days since March 1, 2020",
  main = "C. Timing relative to pandemic onset"
)

dev.off()

# ==============================================================================
# 22. SUPPLEMENTAL FIGURE S2. COLLINEARITY CHECK
# ==============================================================================

cont_vars <- c("maternalage", "prepregnancybmi", "GAD_days_postpart", "days_since_pandemic")

cont_mat <- lmmData.T3 %>%
  select(all_of(cont_vars)) %>%
  mutate(across(everything(), as.numeric))

cor_matrix <- cor(cont_mat, use = "pairwise.complete.obs")
print(round(cor_matrix, 3))

pdf("Supplemental_Figure_S2_Collinearity.pdf", width = 7, height = 7)
pairs(
  cont_mat,
  main = "Supplemental Figure S2. Collinearity among continuous covariates",
  pch = 16,
  col = "#acdffbff",
  cex = 0.5,
  labels = c("Maternal age", "Pre-pregnancy BMI", "Days postpartum", "Days since pandemic")
)
dev.off()

png("Supplemental_Figure_S2_Collinearity.png", width = 1800, height = 1800, res = 300)
pairs(
  cont_mat,
  main = "Supplemental Figure S2. Collinearity among continuous covariates",
  pch = 16,
  col = "#acdffbff",
  cex = 0.5,
  labels = c("Maternal age", "Pre-pregnancy BMI", "Days postpartum", "Days since pandemic")
)
dev.off()

# ==============================================================================
# 23. MULTIPLE IMPUTATION
# ==============================================================================

cat("\n--- Missing education before imputation ---\n")
cat("N missing mom_educ.factor:", sum(is.na(lmmData.T3$mom_educ.factor)), "\n")

imp_vars <- c(
  "GAD.score", "mom_educ.factor", "maternalage", "raceethnicitycombined.factor",
  "para_binary", "preg_pos.factor", "prepregnancybmi", "depranx",
  "GAD_days_postpart", "days_since_pandemic",
  "il6.log", "il17a.log", "il1b.log", "crp.log"
)

imp_data <- lmmData.T3 %>%
  select(all_of(imp_vars)) %>%
  mutate(
    GAD_days_postpart = as.numeric(GAD_days_postpart),
    days_since_pandemic = as.numeric(days_since_pandemic)
  )

set.seed(42)
imp <- mice(imp_data, m = 50, method = "pmm", printFlag = FALSE)

# ==============================================================================
# 24. PRIMARY QUANTILE REGRESSION
# ==============================================================================

covariates <- paste(
  "maternalage + raceethnicitycombined.factor + para_binary + mom_educ.factor",
  "+ preg_pos.factor + GAD_days_postpart + days_since_pandemic",
  "+ prepregnancybmi + depranx"
)

taus <- c(0.5, 0.75, 0.90)

pool_qr <- function(marker, tau_val, covariates, imp) {
  fits <- lapply(1:imp$m, function(i) {
    d <- complete(imp, i)
    rq(
      as.formula(paste("GAD.score ~", marker, "+", covariates)),
      tau = tau_val,
      data = d,
      model = TRUE
    )
  })
  
  coef_mat <- do.call(cbind, lapply(fits, function(f) {
    cf <- coef(f)
    matrix(cf, ncol = 1, dimnames = list(names(cf), NULL))
  }))
  
  var_mat <- do.call(cbind, lapply(fits, function(f) {
    s <- summary(f, se = "iid")$coefficients
    matrix(s[, "Std. Error"]^2, ncol = 1, dimnames = list(rownames(s), NULL))
  }))
  
  Q_bar <- rowMeans(coef_mat)
  U_bar <- rowMeans(var_mat)
  B <- apply(coef_mat, 1, var)
  T_var <- U_bar + (1 + 1 / imp$m) * B
  se <- sqrt(T_var)
  z <- Q_bar / se
  p <- 2 * pnorm(-abs(z))
  
  data.frame(
    term = names(Q_bar),
    estimate = Q_bar,
    se = se,
    z = z,
    p = p,
    row.names = NULL
  )
}

primary_qr_results <- bind_rows(
  lapply(c("il1b.log", "il6.log", "il17a.log", "crp.log"), function(marker) {
    bind_rows(lapply(taus, function(tau_val) {
      res <- pool_qr(marker, tau_val, covariates, imp)
      res$marker <- marker
      res$tau <- tau_val
      res
    }))
  })
)

# Apply BH FDR across the 12 primary exposure tests only
primary_qr_results <- primary_qr_results %>%
  mutate(p_fdr = NA_real_)

exposure_idx <- primary_qr_results$term == primary_qr_results$marker

primary_qr_results$p_fdr[exposure_idx] <- p.adjust(
  primary_qr_results$p[exposure_idx],
  method = "BH"
)

# ==============================================================================
# 25. TABLE 1. PARTICIPANT CHARACTERISTICS BY GAD-7 SCORE PERCENTILE GROUP
# ==============================================================================

tab1_base <- lmmData.T3 %>%
  mutate(
    gad_timing_wks = as.numeric(GAD_days_postpart / 7),
    days_since_pandemic = as.numeric(days_since_pandemic),
    parity_binary = factor(
      case_when(para == 0 ~ "Nulliparous", para >= 1 ~ "Parous"),
      levels = c("Nulliparous", "Parous")
    ),
    gad_pct = factor(
      case_when(
        GAD.score <= quantile(GAD.score, 0.50) ~ "50th percentile",
        GAD.score <= quantile(GAD.score, 0.75) ~ "75th percentile",
        GAD.score <= quantile(GAD.score, 0.90) ~ "90th percentile",
        GAD.score > quantile(GAD.score, 0.90) ~ "Above 90th percentile"
      ),
      levels = c("50th percentile", "75th percentile", "90th percentile", "Above 90th percentile")
    ),
    preg_pos.factor = as.factor(preg_pos.factor),
    raceethnicitycombined.factor = fct_relevel(
      droplevels(raceethnicitycombined.factor),
      "Asian", "Black", "Hispanic", "White", "Other"
    )
  )

tab1_vars <- c(
  "il1b", "il6", "il17a", "crp",
  "maternalage", "raceethnicitycombined.factor", "mom_educ.factor",
  "parity_binary", "prepregnancybmi", "preg_pos.factor", "depranx",
  "gest_age_sample_wk", "gad_timing_wks", "days_since_pandemic"
)

tab1_labels <- list(
  il1b ~ "IL-1β (pg/mL)",
  il6 ~ "IL-6 (pg/mL)",
  il17a ~ "IL-17A (pg/mL)",
  crp ~ "CRP (pg/mL)",
  maternalage ~ "Maternal age (years)",
  raceethnicitycombined.factor ~ "Race/ethnicity",
  parity_binary ~ "Parity",
  mom_educ.factor ~ "Education",
  preg_pos.factor ~ "SARS-CoV-2 during pregnancy",
  gest_age_sample_wk ~ "Gestational age at sample collection (weeks)",
  gad_timing_wks ~ "GAD-7 completion (weeks postpartum)",
  days_since_pandemic ~ "Days since pandemic onset",
  prepregnancybmi ~ "Pre-pregnancy BMI (kg/m²)",
  depranx ~ "History of anxiety or depression"
)

table1_ft <- tab1_base %>%
  select(all_of(tab1_vars), gad_pct) %>%
  tbl_summary(
    by = gad_pct,
    type = all_dichotomous() ~ "categorical",
    label = tab1_labels,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)",
      all_continuous() ~ "{median} ({p25}, {p75})"
    )
  ) %>%
  add_overall() %>%
  add_p(test = list(raceethnicitycombined.factor ~ "chisq.test")) %>%
  bold_p() %>%
  modify_caption("Table 1. Participant characteristics by GAD-7 score percentile group") %>%
  modify_header(
    label ~ "**Characteristic**",
    stat_0 ~ "**Overall**  \nN = {N}",
    all_stat_cols(stat_0 = FALSE) ~ "**{level}**  \nn = {n}"
  ) %>%
  as_flex_table()

print(table1_ft)
save_as_html(table1_ft, path = "Table_1_Participant_Characteristics_by_GAD7_Percentile_Group.html")
save_as_docx("Table 1" = table1_ft, path = "Table_1.docx")

# ==============================================================================
# 26. TABLE 2. PRIMARY ADJUSTED ASSOCIATIONS
# ==============================================================================

table2_dat <- primary_qr_results %>%
  filter(term == marker) %>%
  mutate(
    Marker = recode(
      marker,
      "il1b.log" = "IL-1β",
      "il6.log" = "IL-6",
      "il17a.log" = "IL-17A",
      "crp.log" = "CRP"
    ),
    Percentile = recode(
      as.character(tau),
      "0.5" = "50th percentile",
      "0.75" = "75th percentile",
      "0.9" = "90th percentile"
    ),
    `Estimate (FDR-adjusted p-value)` = sprintf("%.2f (p=%.3f)", estimate, p_fdr)
  ) %>%
  select(Marker, Percentile, `Estimate (FDR-adjusted p-value)`)

table2_wide <- table2_dat %>%
  mutate(
    Marker = factor(Marker, levels = c("IL-1β", "IL-6", "IL-17A", "CRP")),
    Percentile = factor(Percentile, levels = c("50th percentile", "75th percentile", "90th percentile"))
  ) %>%
  arrange(Marker, Percentile) %>%
  pivot_wider(
    names_from = Percentile,
    values_from = `Estimate (FDR-adjusted p-value)`
  )

table2_gt <- table2_wide %>%
  gt() %>%
  cols_label(
    Marker = "Inflammatory marker",
    `50th percentile` = "50th percentile",
    `75th percentile` = "75th percentile",
    `90th percentile` = "90th percentile"
  ) %>%
  tab_header(
    title = "Table 2. Adjusted associations between third-trimester inflammatory markers and postpartum GAD-7 scores across the symptom distribution"
  ) %>%
  cols_align(align = "left", columns = Marker) %>%
  cols_align(align = "center", columns = c(`50th percentile`, `75th percentile`, `90th percentile`)) %>%
  tab_options(
    table.font.size = 11,
    data_row.padding = px(4),
    table.border.top.width = px(1),
    table.border.bottom.width = px(1),
    heading.border.bottom.width = px(1)
  ) %>%
  tab_source_note(
    source_note = md("Values are pooled β coefficients from adjusted quantile regression models fit across 50 imputed datasets and reported as estimate with FDR-adjusted p-value. Each model included one inflammatory marker at a time and adjusted for maternal age, race/ethnicity, parity, education, SARS-CoV-2 infection during pregnancy, days postpartum at GAD-7 completion, days since pandemic onset, pre-pregnancy BMI, and history of anxiety or depression.")
  )

print(table2_gt)
gtsave(table2_gt, "Table_2.docx")

# ==============================================================================
# 27. FIGURE 1. FLOWCHART OF STUDY SAMPLE
# ==============================================================================

flowchart <- grViz("
digraph flowchart {
  graph [
    layout = dot,
    rankdir = TB,
    splines = ortho,
    nodesep = 0.45,
    ranksep = 0.6
  ]
  node [
    shape = box,
    style = 'filled',
    fillcolor = white,
    color = black,
    fontname = Helvetica,
    fontsize = 11,
    margin = '0.2,0.12'
  ]
  edge [
    color = black,
    arrowsize = 0.7
  ]

  main1 [label = 'Enrolled in the Generation C Study\\n(n = 3,157)']
  main2 [label = 'Completed postpartum survey\\n(n = 959)']
  main3 [label = 'Complete GAD-7 within 24 weeks postpartum\\n+ inflammatory marker data available\\n(n = 670)']
  main4 [label = 'Final analytic sample\\n(n = 237)']

  ex1 [label = 'Excluded (n = 2,198):\\nLost to follow-up (208)\\nDelivered outside MSHS (117)\\nPregnancy loss (98)\\nMultiple gestation (46)\\nMissing postpartum survey (1,729)']

  ex2 [label = 'Excluded (n = 289)*:\\nIncomplete GAD-7 (235)\\nNo inflammatory markers (45)\\nGAD-7 >24 weeks postpartum (13)']

  ex3 [label = 'Excluded (n = 433)*:\\nNo third-trimester sample ≥28 weeks (205)\\nSample collected ≤7 days before delivery (263)\\nAll cytokine values missing (21)']

  main1 -> main2 -> main3 -> main4

  main1 -> ex1 [dir = none]
  main2 -> ex2 [dir = none]
  main3 -> ex3 [dir = none]

  {rank = same; main1; ex1}
  {rank = same; main2; ex2}
  {rank = same; main3; ex3}
}
")

flowchart %>%
  export_svg() %>%
  charToRaw() %>%
  rsvg_pdf("Figure_1_Flowchart_of_Study_Sample.pdf")

flowchart %>%
  export_svg() %>%
  charToRaw() %>%
  rsvg_png("Figure_1_Flowchart_of_Study_Sample.png", width = 2400, height = 1800)

flowchart

# ==============================================================================
# 28. FIGURE 2. ADJUSTED ASSOCIATIONS WITH POSTPARTUM GAD-7 SCORES
#     Includes inflammatory markers + psychiatric history
# ==============================================================================

forest_rows <- list()
psych_rows_collected <- character(0)

for (marker in c("il6.log", "il17a.log", "il1b.log", "crp.log")) {
  for (tau_val in taus) {
    
    res <- tryCatch(
      pool_qr(marker, tau_val, covariates, imp),
      error = function(e) NULL
    )
    if (is.null(res)) next
    
    cyto_row <- res %>%
      filter(term == marker) %>%
      mutate(marker = marker, tau = tau_val)
    
    if (nrow(cyto_row) == 1) {
      forest_rows[[length(forest_rows) + 1]] <- cyto_row
    }
    
    key <- as.character(tau_val)
    if (!key %in% psych_rows_collected) {
      psych_row <- res %>%
        filter(term == "depranxYes") %>%
        mutate(marker = "depranx", tau = tau_val)
      
      if (nrow(psych_row) == 1) {
        forest_rows[[length(forest_rows) + 1]] <- psych_row
        psych_rows_collected <- c(psych_rows_collected, key)
      }
    }
  }
}

forest_df <- bind_rows(forest_rows) %>%
  distinct(marker, tau, .keep_all = TRUE)

if (nrow(forest_df) == 0) {
  stop("forest_df is empty. Check that pool_qr() is running correctly.")
}

exposure_fdr <- primary_qr_results %>%
  filter(term == marker) %>%
  select(marker, tau, p_fdr)

forest_df <- forest_df %>%
  left_join(exposure_fdr, by = c("marker", "tau")) %>%
  mutate(
    significant = ifelse(!is.na(p_fdr), p_fdr < 0.05, FALSE),
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    marker_label = factor(
      marker,
      levels = c("depranx", "crp.log", "il17a.log", "il6.log", "il1b.log"),
      labels = c("History of anxiety/depression", "CRP", "IL-17A", "IL-6", "IL-1\u03b2")
    ),
    tau_label = factor(
      tau,
      levels = c(0.50, 0.75, 0.90),
      labels = c("50th percentile", "75th percentile", "90th percentile")
    )
  )

pal <- c(
  "50th percentile" = "#1b7837",
  "75th percentile" = "#762a83",
  "90th percentile" = "#e08214"
)

shapes <- c(
  "50th percentile" = 16,
  "75th percentile" = 17,
  "90th percentile" = 15
)

pd <- position_dodge(width = 0.6)

forest_plot <- ggplot(
  forest_df,
  aes(
    x = estimate,
    y = marker_label,
    xmin = ci_lo,
    xmax = ci_hi,
    colour = tau_label,
    shape = tau_label,
    group = tau_label
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(
    width = 0.18,
    linewidth = 0.7,
    position = pd
  ) +
  geom_point(
    aes(size = significant),
    position = pd,
    stroke = 0
  ) +
  scale_colour_manual(name = "GAD-7 percentile", values = pal) +
  scale_shape_manual(name = "GAD-7 percentile", values = shapes) +
  scale_size_manual(
    values = c("TRUE" = 4.5, "FALSE" = 2.5),
    guide = "none"
  ) +
  labs(
    x = "Quantile regression coefficient (change in GAD-7 score)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 11),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  "Figure 2.forest_plot_gad7.pdf",
  plot = forest_plot,
  width = 7,
  height = 5,
  device = cairo_pdf
)

ggsave(
  "Figure 2.forest_plot_gad7.png",
  plot = forest_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# ==============================================================================
# 29. SUPPLEMENTAL TABLE S1. FULL MULTIVARIABLE QUANTILE REGRESSION RESULTS
# ==============================================================================
pretty_term_names <- c(
  "(Intercept)" = "Intercept",
  "il1b.log" = "IL-1β",
  "il6.log" = "IL-6",
  "il17a.log" = "IL-17A",
  "crp.log" = "CRP",
  "maternalage" = "Maternal age",
  "raceethnicitycombined.factorAsian" = "Asian",
  "raceethnicitycombined.factorBlack" = "Black",
  "raceethnicitycombined.factorHispanic" = "Hispanic",
  "raceethnicitycombined.factorWhite" = "White",
  "raceethnicitycombined.factorOther" = "Other",
  "para_binaryMultiparous" = "Multiparous",
  "mom_educ.factor≤ Some college" = "≤ Some college",
  "mom_educ.factor≥ Bachelors" = "≥ Bachelor's degree",
  "preg_pos.factorNo" = "No SARS-CoV-2 infection during pregnancy",
  "preg_pos.factorYes" = "SARS-CoV-2 infection during pregnancy",
  "GAD_days_postpart" = "Days postpartum at GAD-7",
  "days_since_pandemic" = "Days since pandemic onset",
  "prepregnancybmi" = "Pre-pregnancy BMI",
  "depranxYes" = "History of anxiety/depression"
)

term_order <- c(
  "Intercept",
  "IL-1β",
  "IL-6",
  "IL-17A",
  "CRP",
  "Maternal age",
  "Asian",
  "Black",
  "Hispanic",
  "White",
  "Other",
  "Multiparous",
  "≤ Some college",
  "≥ Bachelor's degree",
  "No SARS-CoV-2 infection during pregnancy",
  "SARS-CoV-2 infection during pregnancy",
  "Days postpartum at GAD-7",
  "Days since pandemic onset",
  "Pre-pregnancy BMI",
  "History of anxiety/depression"
)

marker_order <- c("IL-1β", "IL-6", "IL-17A", "CRP")
percentile_order <- c("50th percentile", "75th percentile", "90th percentile")

marker_map <- c(
  "il1b.log" = "IL-1β",
  "il6.log" = "IL-6",
  "il17a.log" = "IL-17A",
  "crp.log" = "CRP"
)

percentile_map <- c(
  "0.5" = "50th percentile",
  "0.75" = "75th percentile",
  "0.9" = "90th percentile"
)

all_results_s1 <- map_dfr(c("il1b.log", "il6.log", "il17a.log", "crp.log"), function(marker) {
  map_dfr(c(0.5, 0.75, 0.9), function(tau_val) {
    pool_qr(marker, tau_val, covariates, imp) %>%
      mutate(
        marker = marker,
        tau = tau_val
      )
  })
})

all_results_s1 <- all_results_s1 %>%
  left_join(
    primary_qr_results %>%
      filter(term == marker) %>%
      select(marker, tau, term, p_fdr),
    by = c("marker", "tau", "term")
  ) %>%
  mutate(
    Marker = unname(marker_map[marker]),
    Percentile = unname(percentile_map[as.character(tau)]),
    Variable = dplyr::recode(term, !!!pretty_term_names),
    ci_lo = estimate - 1.96 * se,
    ci_hi = estimate + 1.96 * se,
    `β (95% CI)` = sprintf("%.2f (%.2f, %.2f)", estimate, ci_lo, ci_hi),
    `p-value` = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)),
    `FDR-adjusted p-value` = case_when(
      term == marker & is.na(p_fdr) ~ "",
      term == marker & p_fdr < 0.001 ~ "<0.001",
      term == marker ~ sprintf("%.3f", p_fdr),
      TRUE ~ ""
    )
  ) %>%
  filter(!is.na(Variable)) %>%
  filter(Variable %in% term_order) %>%
  mutate(
    Marker = factor(Marker, levels = marker_order),
    Percentile = factor(Percentile, levels = percentile_order),
    Variable = factor(Variable, levels = term_order)
  ) %>%
  arrange(Marker, Percentile, Variable) %>%
  group_by(Marker) %>%
  mutate(Marker_display = if_else(row_number() == 1, as.character(Marker), "")) %>%
  group_by(Marker, Percentile) %>%
  mutate(Percentile_display = if_else(row_number() == 1, as.character(Percentile), "")) %>%
  ungroup()

supp_table_s1 <- all_results_s1 %>%
  select(
    Marker = Marker_display,
    Percentile = Percentile_display,
    Variable,
    `β (95% CI)`,
    `p-value`,
    `FDR-adjusted p-value`
  ) %>%
  gt() %>%
  cols_label(
    Marker = "Marker",
    Percentile = "Percentile",
    Variable = "Variable",
    `β (95% CI)` = "β (95% CI)",
    `p-value` = "p-value",
    `FDR-adjusted p-value` = "FDR-adjusted p-value"
  ) %>%
  tab_header(
    title = "Supplemental Table S1. Full multivariable quantile regression results"
  ) %>%
  cols_align(
    align = "left",
    columns = c(Marker, Percentile, Variable)
  ) %>%
  cols_align(
    align = "center",
    columns = c(`β (95% CI)`, `p-value`, `FDR-adjusted p-value`)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Marker, rows = Marker != "")
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Percentile, rows = Percentile != "")
  ) %>%
  tab_options(
    table.font.size = 11,
    data_row.padding = px(4),
    table.border.top.width = px(1),
    table.border.bottom.width = px(1),
    heading.border.bottom.width = px(1)
  )

supp_table_s1

gtsave(
  supp_table_s1,
  "Supplemental_Table_S1_Full_Multivariable_Quantile_Regression_Results.docx"
)

write.csv(
  all_results_s1 %>%
    select(
      Marker,
      Percentile,
      Variable,
      `β (95% CI)`,
      `p-value`,
      `FDR-adjusted p-value`
    ),
  "Supplemental_Table_S1_Full_Multivariable_Quantile_Regression_Results.csv",
  row.names = FALSE
)
# ==============================================================================
# 30. SENSITIVITY ANALYSIS 1. COMPLETE-CASE QUANTILE REGRESSION
# ==============================================================================

cc_vars <- c(
  "GAD.score", "maternalage", "raceethnicitycombined.factor", "para_binary",
  "mom_educ.factor", "preg_pos.factor", "GAD_days_postpart",
  "days_since_pandemic", "prepregnancybmi", "depranx",
  "il6.log", "il17a.log", "il1b.log", "crp.log"
)

lmmData.T3.cc <- lmmData.T3 %>%
  mutate(
    GAD_days_postpart = as.numeric(GAD_days_postpart),
    days_since_pandemic = as.numeric(days_since_pandemic),
    mom_educ.factor = factor(mom_educ.factor),
    raceethnicitycombined.factor = factor(raceethnicitycombined.factor),
    para_binary = factor(para_binary),
    preg_pos.factor = factor(preg_pos.factor),
    depranx = factor(depranx)
  ) %>%
  filter(if_all(all_of(cc_vars), ~ !is.na(.)))

cat("\n--- Sensitivity analysis 1: complete-case sample ---\n")
cat("N complete-case:", nrow(lmmData.T3.cc), "\n")

extract_qr_journal <- function(data_in, marker, marker_label) {
  percentiles <- c(
    "0.5" = "50th percentile",
    "0.75" = "75th percentile",
    "0.9" = "90th percentile"
  )
  
  map_dfr(names(percentiles), function(tau_chr) {
    tau_val <- as.numeric(tau_chr)
    
    fml <- reformulate(
      c(
        marker,
        "maternalage",
        "raceethnicitycombined.factor",
        "para_binary",
        "mom_educ.factor",
        "preg_pos.factor",
        "GAD_days_postpart",
        "days_since_pandemic",
        "prepregnancybmi",
        "depranx"
      ),
      response = "GAD.score"
    )
    
    fit <- rq(fml, tau = tau_val, data = data_in, model = TRUE)
    summ <- summary(fit, se = "boot")$coefficients
    
    tibble(
      Marker = marker_label,
      Variable = rownames(summ),
      Percentile = percentiles[[tau_chr]],
      beta = summ[, 1],
      se = summ[, 2],
      ci_lo = summ[, 1] - 1.96 * summ[, 2],
      ci_hi = summ[, 1] + 1.96 * summ[, 2],
      p = summ[, 4]
    )
  })
}

pretty_variable_names_cc <- c(
  "(Intercept)" = "Intercept",
  "il1b.log" = "IL-1β",
  "il6.log" = "IL-6",
  "il17a.log" = "IL-17A",
  "crp.log" = "CRP",
  "maternalage" = "Maternal age",
  "raceethnicitycombined.factorAsian" = "Asian",
  "raceethnicitycombined.factorBlack" = "Black",
  "raceethnicitycombined.factorHispanic" = "Hispanic",
  "raceethnicitycombined.factorWhite" = "White",
  "para_binaryMultiparous" = "Multiparous",
  "mom_educ.factor≤ Some college" = "≤ Some college",
  "mom_educ.factor≥ Bachelors" = "≥ Bachelor's degree",
  "preg_pos.factorNo" = "No SARS-CoV-2 infection during pregnancy",
  "GAD_days_postpart" = "Days postpartum at GAD-7",
  "days_since_pandemic" = "Days since pandemic onset",
  "prepregnancybmi" = "Pre-pregnancy BMI",
  "depranxYes" = "History of anxiety or depression"
)

marker_order_cc <- c("IL-1β", "IL-6", "IL-17A", "CRP")
variable_order_cc <- c(
  "Intercept", "IL-1β", "IL-6", "IL-17A", "CRP", "Maternal age",
  "Asian", "Black", "Hispanic", "White", "Multiparous",
  "≤ Some college", "≥ Bachelor's degree",
  "No SARS-CoV-2 infection during pregnancy",
  "Days postpartum at GAD-7", "Days since pandemic onset",
  "Pre-pregnancy BMI", "History of anxiety or depression"
)

make_journal_qr_table <- function(data_in, table_title = NULL) {
  results_raw <- bind_rows(
    extract_qr_journal(data_in, "il1b.log", "IL-1β"),
    extract_qr_journal(data_in, "il6.log", "IL-6"),
    extract_qr_journal(data_in, "il17a.log", "IL-17A"),
    extract_qr_journal(data_in, "crp.log", "CRP")
  )
  
  results_fmt <- results_raw %>%
    mutate(
      Variable = recode(Variable, !!!pretty_variable_names_cc),
      `β (95% CI)` = sprintf("%.2f (%.2f, %.2f)", beta, ci_lo, ci_hi)
    ) %>%
    select(Marker, Variable, Percentile, `β (95% CI)`) %>%
    pivot_wider(names_from = Percentile, values_from = `β (95% CI)`) %>%
    mutate(
      Marker = factor(Marker, levels = marker_order_cc),
      Variable = factor(Variable, levels = variable_order_cc)
    ) %>%
    arrange(Marker, Variable) %>%
    mutate(Marker = ifelse(duplicated(Marker), "", as.character(Marker)))
  
  gt_tbl <- results_fmt %>%
    gt() %>%
    cols_label(
      Marker = "Marker",
      Variable = "Variable",
      `50th percentile` = "50th percentile",
      `75th percentile` = "75th percentile",
      `90th percentile` = "90th percentile"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(columns = Marker, rows = Marker != "")
    ) %>%
    cols_align(align = "left", columns = c(Marker, Variable)) %>%
    cols_align(
      align = "center",
      columns = c(`50th percentile`, `75th percentile`, `90th percentile`)
    ) %>%
    tab_options(
      table.font.size = 11,
      data_row.padding = px(4),
      table.border.top.width = px(1),
      table.border.bottom.width = px(1),
      heading.border.bottom.width = px(1)
    )
  
  if (!is.null(table_title)) {
    gt_tbl <- gt_tbl %>% tab_header(title = table_title)
  }
  
  gt_tbl
}

supp_table_s2 <- make_journal_qr_table(
  lmmData.T3.cc,
  table_title = "Supplemental Table S2. Complete-case quantile regression results"
)

supp_table_s2
gtsave(supp_table_s2, "Supplemental_Table_S2_Complete_Case_Quantile_Regression_Results.docx")

# ==============================================================================
# 31. SENSITIVITY ANALYSIS 2. <=12 WEEKS POSTPARTUM
# ==============================================================================

lmmData.T3.12wk <- lmmData.T3 %>%
  mutate(
    GAD_days_postpart   = as.numeric(GAD_days_postpart),
    days_since_pandemic = as.numeric(days_since_pandemic),
    mom_educ.factor = factor(mom_educ.factor),
    raceethnicitycombined.factor = factor(raceethnicitycombined.factor),
    para_binary = factor(para_binary),
    preg_pos.factor = factor(preg_pos.factor),
    depranx = factor(depranx)
  ) %>%
  filter(GAD_wks_postpart <= 12)

cat("\n--- Sensitivity analysis 2: <=12 weeks postpartum before imputation ---\n")
cat("N <=12 weeks:", nrow(lmmData.T3.12wk), "\n")

imp_vars_12wk <- c(
  "GAD.score",
  "mom_educ.factor",
  "maternalage",
  "raceethnicitycombined.factor",
  "para_binary",
  "preg_pos.factor",
  "prepregnancybmi",
  "depranx",
  "GAD_days_postpart",
  "days_since_pandemic",
  "il6.log",
  "il17a.log",
  "il1b.log",
  "crp.log"
)

imp_data_12wk <- lmmData.T3.12wk %>%
  select(all_of(imp_vars_12wk))

cat("N missing education in <=12 week sample:",
    sum(is.na(imp_data_12wk$mom_educ.factor)), "\n")

set.seed(42)
imp_12wk <- mice(
  imp_data_12wk,
  m = 50,
  method = "pmm",
  printFlag = FALSE
)

cat("Imputed sample size:", nrow(complete(imp_12wk, 1)), "\n")

covariates <- paste(
  "maternalage + raceethnicitycombined.factor + para_binary + mom_educ.factor",
  "+ preg_pos.factor + GAD_days_postpart + days_since_pandemic",
  "+ prepregnancybmi + depranx"
)

pool_qr_12wk <- function(marker, tau_val, covariates, imp_obj) {
  fits <- lapply(1:imp_obj$m, function(i) {
    d <- complete(imp_obj, i)
    rq(
      as.formula(paste("GAD.score ~", marker, "+", covariates)),
      tau = tau_val,
      data = d,
      model = TRUE
    )
  })
  
  coef_mat <- do.call(cbind, lapply(fits, function(f) {
    cf <- coef(f)
    matrix(cf, ncol = 1, dimnames = list(names(cf), NULL))
  }))
  
  var_mat <- do.call(cbind, lapply(fits, function(f) {
    s <- summary(f, se = "iid")$coefficients
    matrix(
      s[, "Std. Error"]^2,
      ncol = 1,
      dimnames = list(rownames(s), NULL)
    )
  }))
  
  Q_bar <- rowMeans(coef_mat)
  U_bar <- rowMeans(var_mat)
  B     <- apply(coef_mat, 1, var)
  T_var <- U_bar + (1 + 1 / imp_obj$m) * B
  se    <- sqrt(T_var)
  z     <- Q_bar / se
  p     <- 2 * pnorm(-abs(z))
  
  data.frame(
    term = names(Q_bar),
    estimate = Q_bar,
    se = se,
    ci_lo = Q_bar - 1.96 * se,
    ci_hi = Q_bar + 1.96 * se,
    z = z,
    p = p,
    row.names = NULL
  )
}

make_journal_qr_table_imputed <- function(imp_obj, table_title = NULL) {
  
  percentiles <- c(
    "0.5"  = "50th percentile",
    "0.75" = "75th percentile",
    "0.9"  = "90th percentile"
  )
  
  marker_map <- c(
    "il1b.log"  = "IL-1β",
    "il6.log"   = "IL-6",
    "il17a.log" = "IL-17A",
    "crp.log"   = "CRP"
  )
  
  pretty_variable_names <- c(
    "(Intercept)" = "Intercept",
    "il1b.log" = "IL-1β",
    "il6.log" = "IL-6",
    "il17a.log" = "IL-17A",
    "crp.log" = "CRP",
    "maternalage" = "Maternal age",
    "raceethnicitycombined.factorAsian" = "Asian",
    "raceethnicitycombined.factorBlack" = "Black",
    "raceethnicitycombined.factorHispanic" = "Hispanic",
    "raceethnicitycombined.factorWhite" = "White",
    "para_binaryMultiparous" = "Multiparous",
    "mom_educ.factor≤ Some college" = "≤ Some college",
    "mom_educ.factor≥ Bachelors" = "≥ Bachelor's degree",
    "preg_pos.factorNo" = "No SARS-CoV-2 infection during pregnancy",
    "GAD_days_postpart" = "Days postpartum at GAD-7",
    "days_since_pandemic" = "Days since pandemic onset",
    "prepregnancybmi" = "Pre-pregnancy BMI",
    "depranxYes" = "History of anxiety or depression"
  )
  
  marker_order <- c("IL-1β", "IL-6", "IL-17A", "CRP")
  
  variable_order <- c(
    "Intercept",
    "IL-1β",
    "IL-6",
    "IL-17A",
    "CRP",
    "Maternal age",
    "Asian",
    "Black",
    "Hispanic",
    "White",
    "Multiparous",
    "≤ Some college",
    "≥ Bachelor's degree",
    "No SARS-CoV-2 infection during pregnancy",
    "Days postpartum at GAD-7",
    "Days since pandemic onset",
    "Pre-pregnancy BMI",
    "History of anxiety or depression"
  )
  
  results_raw <- bind_rows(
    lapply(names(marker_map), function(marker) {
      bind_rows(
        lapply(names(percentiles), function(tau_chr) {
          tau_val <- as.numeric(tau_chr)
          
          res <- pool_qr_12wk(marker, tau_val, covariates, imp_obj)
          
          res %>%
            mutate(
              Marker = marker_map[[marker]],
              Variable = recode(term, !!!pretty_variable_names),
              Percentile = percentiles[[tau_chr]],
              `β (95% CI)` = sprintf("%.2f (%.2f, %.2f)", estimate, ci_lo, ci_hi)
            ) %>%
            select(Marker, Variable, Percentile, `β (95% CI)`)
        })
      )
    })
  )
  
  results_fmt <- results_raw %>%
    pivot_wider(
      names_from = Percentile,
      values_from = `β (95% CI)`
    ) %>%
    mutate(
      Marker = factor(Marker, levels = marker_order),
      Variable = factor(Variable, levels = variable_order)
    ) %>%
    arrange(Marker, Variable) %>%
    mutate(
      Marker = ifelse(duplicated(Marker), "", as.character(Marker))
    )
  
  gt_tbl <- results_fmt %>%
    gt() %>%
    cols_label(
      Marker = "Marker",
      Variable = "Variable",
      `50th percentile` = "50th percentile",
      `75th percentile` = "75th percentile",
      `90th percentile` = "90th percentile"
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = Marker,
        rows = Marker != ""
      )
    ) %>%
    cols_align(
      align = "left",
      columns = c(Marker, Variable)
    ) %>%
    cols_align(
      align = "center",
      columns = c(`50th percentile`, `75th percentile`, `90th percentile`)
    ) %>%
    tab_options(
      table.font.size = 11,
      data_row.padding = px(4),
      table.border.top.width = px(1),
      table.border.bottom.width = px(1),
      heading.border.bottom.width = px(1)
    )
  
  if (!is.null(table_title)) {
    gt_tbl <- gt_tbl %>%
      tab_header(title = table_title)
  }
  
  gt_tbl
}

supp_table_s3 <- make_journal_qr_table_imputed(
  imp_obj = imp_12wk,
  table_title = "Supplementary Table S3. Quantile regression results restricted to participants completing the GAD-7 within 12 weeks postpartum"
)

supp_table_s3

gtsave(
  supp_table_s3,
  "Supplementary_Table_S3_Imputed_Quantile_Regression_12_Weeks.docx"
)

# ==============================================================================
# END
# ==============================================================================

sessionInfo()
