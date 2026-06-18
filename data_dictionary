# Data Dictionary

## Overview
This document describes the key variables used in the analysis. The main dataset is a prospective pregnancy cohort with demographic, clinical, and inflammatory marker data.

## Main Study Dataset

### Participant Identifiers
| Variable | Type | Description | Range/Categories |
|----------|------|-------------|-------------------|
| `subject_id` | Integer | Unique participant identifier | 1-3157 |

### Demographics (Assessed at Enrollment, ~20 weeks gestation)
| Variable | Type | Description | Range/Categories |
|----------|------|-------------|-------------------|
| `maternalage` | Numeric | Maternal age at enrollment | 18-49 years |
| `raceethnicitycombined` | Categorical | Self-reported race/ethnicity | 1=AI/AN, 2=Asian, 3=Black, 4=Hispanic, 5=NHOPI, 6=White, 7=Other, 8=Unknown |
| `insurance_cat` | Categorical | Primary health insurance type | 1=Private, 2=Public, 3=Self-pay |
| `para` | Integer | Number of prior live births (parity) | 0-5+ |
| `prepregnancybmi` | Numeric | BMI calculated from pre-pregnancy weight/height | 16.5-52.3 kg/m² |

### Education & Income (May come from multiple visits; p1/s1/s3)
| Variable | Type | Description | Categories |
|----------|------|-------------|------------|
| `p1_mom_educ` | Categorical | Highest education level completed | 1=<9th grade, 2=Some HS, 3=HS/GED, 4=Trade/Tech, 5=Some college, 6=Bachelors, 7=Post-grad |
| `s1_mom_educ` | Categorical | (same as above) | 1-7 |
| `s3_mom_educ` | Categorical | (same as above) | 1-7 |
| `p1_income` | Categorical | Total household income (last year, before taxes) | 1=<$25k, 2=$25-50k, 3=$50-75k, 4=$75-100k, 5=$100-125k, 6=>$125k |

### Obstetric History
| Variable | Type | Description | Range/Categories |
|----------|------|-------------|-------------------|
| `gestationalagedays` | Integer | Gestational age at delivery | 168-294 days (24-42 weeks) |
| `birthweight_grams` | Integer | Infant birth weight | 450-5000 grams |
| `birthdate` | Date | Infant date of birth | YYYY-MM-DD |

### Infection & Vaccination
| Variable | Type | Description | Categories |
|----------|------|-------------|------------|
| `preg_pos` | Binary | SARS-CoV-2 infection during pregnancy | 0=No, 1=Yes, NA=Unknown |
| `vaccine_timing` | Categorical | COVID-19 vaccine timing relative to pregnancy | 1=Before pregnancy, 2=During pregnancy, 3=After pregnancy |

### Psychiatric History (Merged from separate dataset)
| Variable | Type | Description | Categories |
|----------|------|-------------|------------|
| `depression` | Binary | History of depression | 0=No, 1=Yes |
| `anxiety` | Binary | History of anxiety disorder | 0=No, 1=Yes |
| `depranx` | Binary | History of anxiety **or** depression | 0=No, 1=Yes |
| `affective` | Binary | Any affective disorder history | 0=No, 1=Yes |
| `anymental` | Binary | Any mental health diagnosis history | 0=No, 1=Yes |
| `depressionandppd` | Binary | Depression or postpartum depression history | 0=No, 1=Yes |

## Inflammatory Markers

### Collection Schedule
Serum samples collected at up to 6 timepoints across pregnancy and postpartum:
- Timepoint 1-3: During pregnancy (typically 1st, 2nd, 3rd trimester)
- Timepoint 4-6: Postpartum (typically 0-6 months postpartum)

### Marker Variables (6 timepoints each: _1 through _6)
| Variable | Type | Units | Detection Limit | Notes |
|----------|------|-------|-----------------|-------|
| `il6_i` | Numeric | pg/mL | 0.01 | Interleukin-6 |
| `il1b_i` | Numeric | pg/mL | 0.03 | Interleukin-1 beta |
| `il17a_i` | Numeric | pg/mL | 0.64 | Interleukin-17A |
| `crp_i` | Numeric | pg/mL | Varies | C-reactive protein |

**Special handling in analysis:**
- Out-of-range values ("OOR", "0", "OOR <") replaced with assay limit of detection
- CRP missing values imputed to 195.69 pg/mL (maximum measured) when other markers present
- All markers log2-transformed for analysis

### Sample Timing Variables
| Variable | Type | Description |
|----------|------|-------------|
| `date_sample_i` | Date | Collection date of sample i (YYYY-MM-DD) |
| `processtime_i` | Character | Time from collection to processing for sample i |
| `sample_i_gest_age` | Numeric | Gestational age (days) at sample collection i |
| `gest_age_sample_wk` | Numeric | Gestational age (weeks) at sample collection (derived) |

## Postpartum Anxiety Assessment (GAD-7)

### Administration
| Variable | Type | Description |
|----------|------|-------------|
| `p3_timestamp` | Date | Date of GAD-7 completion at postpartum visit 3 (YYYY-MM-DD) |
| `s2_timestamp` | Date | Date of GAD-7 completion at study visit 2 (YYYY-MM-DD) |
| `s3_timestamp` | Date | Date of GAD-7 completion at study visit 3 (YYYY-MM-DD) |

### Individual Items (0-3 response scale)
Questions asked: "Over the last 2 weeks, how often have you been bothered by..."

| Variable | Description | Scoring |
|----------|-------------|---------|
| `{p3,s2,s3}_2w_nerv` | Feeling nervous, anxious, or on edge | 0=Not at all, 1=Several days, 2=More than half, 3=Nearly every day |
| `{p3,s2,s3}_2w_uncontrworry` | Not being able to stop or control worrying | (same as above) |
| `{p3,s2,s3}_2w_toomuchworry` | Worrying too much about different things | (same as above) |
| `{p3,s2,s3}_2w_relax` | Trouble relaxing | (same as above) |
| `{p3,s2,s3}_2w_restless` | Being so restless that it's hard to sit still | (same as above) |
| `{p3,s2,s3}_2w_annoy` | Becoming easily annoyed or irritable | (same as above) |
| `{p3,s2,s3}_2w_afraid` | Feeling afraid as if something awful might happen | (same as above) |

### Derived Variables
| Variable | Type | Description | Range |
|----------|------|-------------|-------|
| `{p3,s2,s3}_2w_GAD.score` | Integer | Total GAD-7 score from 7 items | 0-21 |
| `{p3,s2,s3}_2w_PHQ.score` | Integer | PHQ-2 score from 2 depression items | 0-6 |
| `GAD.score` | Integer | Best available GAD-7 score (coalesced priority: p3 > s2 > s3) | 0-21 |
| `PHQ.score` | Integer | Best available PHQ-2 score (coalesced priority: p3 > s2 > s3) | 0-6 |
| `anxiety.factor` | Categorical | Categorical anxiety severity | Minimal (<5), Mild (5-9), Moderate (10-14), Severe (≥15) |

### Timing Variables
| Variable | Type | Description |
|----------|------|-------------|
| `GAD_days_postpart` | Integer | Days postpartum when GAD-7 completed | 0-168 days |
| `GAD_wks_postpart` | Numeric | Weeks postpartum when GAD-7 completed | 0-24 weeks |

## Contextual Variables

| Variable | Type | Description | Notes |
|----------|------|-------------|-------|
| `days_since_pandemic` | Numeric | Days since March 1, 2020 (COVID-19 pandemic onset) | Used to control for pandemic timing |
| `duplicate_enrollment` | Binary | Whether participant enrolled twice | 1=Yes (excluded), 0=No |

## Data Transformations Applied in Analysis

### Log-Transformation
All inflammatory markers were log2-transformed:
- `il6.log = log(il6, base=2)`
- `il1b.log = log(il1b, base=2)`
- `il17a.log = log(il17a, base=2)`
- `crp.log = log(crp, base=2)`

### Factor Recoding
- **Race/Ethnicity:** Collapsed to 5 categories (White, Hispanic, Black, Asian, Other)
- **Education:** Collapsed to 3 categories (<HS, HS-some college, ≥Bachelors)
- **Parity:** Binary (Nulliparous [para=0] vs. Multiparous [para≥1])
- **Psychiatric history:** Binary (depranx = 0 or 1)

## Missing Data

### Reasons for Missingness
- **Inflammatory markers:** Not all participants sampled at all timepoints; some samples failed processing or had out-of-range values
- **Education:** 8 participants (3.8%) missing education data—imputed using MICE
- **Pre-pregnancy BMI:** Calculated from self-reported pre-pregnancy weight/height; missing if not reported
- **Psychiatric history:** Not assessed in all participants; filled from available records

### Imputation Details
- **Method:** Multiple imputation by chained equations (MICE) with predictive mean matching (PMM)
- **Number of imputations:** 50 datasets
- **Variables included:** GAD-7 score, education, maternal age, race/ethnicity, parity, COVID-19 status, BMI, psychiatric history, timing variables, and all inflammatory markers
- **Seed:** 42 (for reproducibility)

## Analytic Sample Definitions

### Primary Analytic Sample (n=237)
Inclusion criteria:
- Complete GAD-7 assessment ≤24 weeks postpartum
- ≥1 inflammatory marker measured
- Third-trimester sample (≥28 weeks gestation, >7 days before delivery)
- Sample not completed during pregnancy

### Complete-Case Sensitivity Sample (n=220)
- Primary sample with no missing covariates
- Used for robustness checking

### ≤12 Weeks Postpartum Sensitivity Sample (n=187)
- Primary sample restricted to GAD-7 completion ≤12 weeks postpartum
- Used to examine whether earlier anxiety assessment changes findings

---

## Notes for Data Users

1. **Subject ID assignment:** Unique identifiers assigned at enrollment; some gaps due to screening failures or duplicate enrollments (marked in `duplicate_enrollment`)

2. **Multiple measurements:** Some variables measured at multiple study visits (denoted with prefix p1, s1, s3 for visits); analysis uses coalesced/prioritized versions

3. **Date formatting:** All dates in YYYY-MM-DD format. Some dates may be imputed if exact date not recorded but month/year known

4. **Unit consistency:** All cytokines in pg/mL; weights in grams; ages in years; dates in days where appropriate

5. **Confidentiality:** All data de-identified; subject IDs are sequential but non-personally identifiable

---

Last updated: [Date]
Data dictionary version: 1.0
