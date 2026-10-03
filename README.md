# Drug X vs Placebo in Type 2 Diabetes: SDTM → ADaM → TLF Pipeline and ANCOVA Analysis in SAS

An end-to-end clinical SAS programming and biostatistics project. Raw trial data are converted to **CDISC SDTM**, derived into **ADaM** analysis datasets, and analysed to produce **Tables, Listings and Figures (TLF)**. The primary endpoint is the **change from baseline in fasting glucose (mg/dL) at Week 12**, analysed with a two-sample t-test and an **ANCOVA** adjusted for baseline glucose.

**Pipeline:** Raw CSV → SDTM (DM, EX, LB) → ADaM (ADSL, ADLB) → Analysis (t-test, ANCOVA) → TLF (RTF table, CSV listing, PNG figure)

## Key results at a glance

| Statistic | Result |
|---|---|
| Primary endpoint | Change in fasting glucose from baseline to Week 12 |
| Drug X mean change | −31.2 mg/dL (SD 9.7, n = 35) |
| Placebo mean change | −7.9 mg/dL (SD 11.4, n = 39) |
| Unadjusted difference | −23.30 mg/dL, 95% CI −28.23 to −18.36 |
| **ANCOVA adjusted difference** | **−23.27 mg/dL, 95% CI −28.23 to −18.30, p < 0.0001** |
| Effect size (Cohen's d) | 2.19, 95% CI 1.61 to 2.77 |
| Responder rate (≥ 20 mg/dL reduction) | 88.6% vs 15.4% |
| Model fit | R² = 0.552, root MSE = 10.70 mg/dL |

## Study design

| Item | Detail |
|---|---|
| Design | Two-arm, parallel-group, 1:1 allocation |
| Treatments | Drug X 10 mg (n = 40) vs placebo (n = 40) |
| Population | 80 adults with Type 2 Diabetes across 4 sites |
| Visits | Screening, Baseline, Week 4, Week 8, Week 12 |
| Records | 80 subjects, 394 glucose records |
| Completion | 74 of 80 (92.5%) had a Week 12 value: Drug X 35/40 (87.5%), placebo 39/40 (97.5%) |
| Hypotheses | H0: μ(Drug X) − μ(Placebo) = 0 vs H1: ≠ 0, two-sided α = 0.05 |

## Workflow

1. **Raw data:** `DM.csv`, `EX.csv`, `LB.csv` imported with `PROC IMPORT` and checked with `PROC CONTENTS`.
2. **SDTM:** `DM`, `EX`, `LB` domains with `USUBJID`, `DOMAIN`, ISO 8601 `--DTC` dates and standardised numeric result `LBSTRESN`.
3. **ADaM:**
   - `ADSL` (80 rows): `TRT01P`, `TRT01A`, treatment dates, duration, `SAFFL`.
   - `ADLB` (394 rows): `AVAL`, `BASE`, `CHG = AVAL − BASE`, `PCHG = 100 × CHG / BASE`.
4. **Analysis:** `PROC MEANS`, `PROC FREQ`, `PROC TTEST`, `PROC GLM` (ANCOVA with `LSMEANS`).
5. **TLF:** demographics table, Week 12 efficacy table (RTF), subject listing (CSV), mean-change-by-visit figure (`PROC SGPLOT`).

## Statistical methods

**Two-sample t-test** (pooled variance, with Satterthwaite as a sensitivity check):

t = (M₁ − M₂) / ( s_p · √(1/n₁ + 1/n₂) ),  s_p² = [(n₁−1)s₁² + (n₂−1)s₂²] / (n₁ + n₂ − 2), df = n₁ + n₂ − 2

**ANCOVA** (`PROC GLM`):

CHG_i = β₀ + β₁·Treatment_i + β₂·Baseline_i + ε_i,  ε_i ~ N(0, σ²)

β₁ is the treatment effect adjusted for baseline. Adjusted (least-squares) means are evaluated at the mean baseline glucose of 165.8 mg/dL. ANCOVA is the standard primary analysis for change-from-baseline endpoints because it corrects for chance baseline imbalance and improves precision.

**Effect size:** Cohen's d = (M₁ − M₂) / s_p.

## Results

### 1. Baseline characteristics

| Characteristic | Drug X (n = 40) | Placebo (n = 40) | Total (n = 80) |
|---|---|---|---|
| Age, mean (SD) | 55.7 (9.3) | 53.5 (11.3) | 54.6 |
| Age range | 33 – 75 | 32 – 75 | 32 – 75 |
| Female, n (%) | 19 (47.5) | 18 (45.0) | 37 (46.3) |
| Male, n (%) | 21 (52.5) | 22 (55.0) | 43 (53.8) |
| Asian, n (%) | 18 (45.0) | 18 (45.0) | 36 (45.0) |
| Black or African American, n (%) | 7 (17.5) | 6 (15.0) | 13 (16.3) |
| White, n (%) | 15 (37.5) | 16 (40.0) | 31 (38.8) |
| Baseline glucose, mean (SD) | 166.4 (19.3) | 165.2 (18.0) | 165.8 |

The arms are balanced on age, sex, race and baseline glucose, so differences in outcome are unlikely to reflect baseline imbalance.

### 2. Week 12 glucose (mg/dL)

| Treatment | Variable | n | Mean | SD | Median | Min | Max |
|---|---|---|---|---|---|---|---|
| Drug X | Baseline | 35 | 166.4 | 19.3 | 167.0 | 123.3 | 202.0 |
| Drug X | Week 12 | 35 | 135.2 | 21.5 | 135.4 | 99.0 | 177.1 |
| Drug X | **Change** | 35 | **−31.2** | 9.7 | −30.9 | −51.3 | −12.5 |
| Placebo | Baseline | 39 | 165.2 | 18.0 | 169.5 | 130.5 | 198.2 |
| Placebo | Week 12 | 39 | 157.4 | 20.5 | 158.3 | 106.3 | 200.1 |
| Placebo | **Change** | 39 | **−7.9** | 11.4 | −6.5 | −33.8 | 8.9 |

Percent change from baseline at Week 12: Drug X −19.0% (SD 6.4) vs placebo −4.7% (SD 7.1). Every Drug X patient had a lower glucose at Week 12 than at baseline (largest change −12.5 mg/dL).

### 3. Response over time (change from baseline, mg/dL)

| Visit | Drug X n | Drug X mean (SD) | Drug X SE | Placebo n | Placebo mean (SD) | Placebo SE |
|---|---|---|---|---|---|---|
| Week 4 | 40 | −14.2 (10.4) | 1.64 | 40 | −2.7 (8.8) | 1.39 |
| Week 8 | 40 | −24.4 (12.1) | 1.91 | 40 | −2.9 (12.7) | 2.00 |
| Week 12 | 35 | −31.2 (9.7) | 1.64 | 39 | −7.9 (11.4) | 1.83 |

The separation between arms appears by Week 4 and widens at each visit. The Drug X effect builds steadily, while placebo stays close to zero until Week 12.

### 4. Two-sample t-test on Week 12 change

| Group / contrast | n | Mean | SD | SE | 95% CI for mean |
|---|---|---|---|---|---|
| Drug X | 35 | −31.17 | 9.68 | 1.64 | (−34.49, −27.84) |
| Placebo | 39 | −7.87 | 11.42 | 1.83 | (−11.57, −4.17) |
| Difference (pooled) | — | −23.30 | 10.64 | 2.48 | (−28.23, −18.36) |
| Difference (Satterthwaite) | — | −23.30 | — | 2.45 | (−28.19, −18.40) |

- Pooled: t(72) = −9.41, p < 0.0001.
- Satterthwaite: t(71.8) = −9.49, p < 0.0001.
- Equality of variances (folded F): F(38, 34) = 1.39, p = 0.33.

### 5. ANCOVA adjusted for baseline glucose

**Analysis of variance**

| Source | DF | Sum of squares | Mean square | F | Pr > F |
|---|---|---|---|---|---|
| Model | 2 | 10,030.03 | 5,015.02 | 43.82 | < 0.0001 |
| Error | 71 | 8,126.27 | 114.45 | — | — |
| Corrected total | 73 | 18,156.30 | — | — | — |

**Parameter estimates (Type III)**

| Effect | Type III F | Estimate | SE | t | p |
|---|---|---|---|---|---|
| Treatment (Drug X vs placebo) | 87.15 | −23.27 | 2.49 | −9.34 | < 0.0001 |
| Baseline glucose (slope) | 0.17 | −0.028 | 0.068 | −0.41 | 0.685 |
| Intercept | — | −3.32 | 11.33 | −0.29 | 0.771 |

**Least-squares means**

| Treatment | Adjusted mean | SE | 95% CI |
|---|---|---|---|
| Drug X | −31.15 | 1.81 | (−34.76, −27.54) |
| Placebo | −7.88 | 1.71 | (−11.30, −4.47) |
| Difference | −23.27 | 2.49 | (−28.23, −18.30) |

R² = 0.552, root MSE = 10.70, mean change = −18.89 mg/dL. Treatment accounts for almost all explained variation (Type I SS 10,011 of 10,030). Baseline glucose is not a significant covariate (p = 0.69), so the adjusted and unadjusted estimates agree within 0.03 mg/dL.

### 6. Model diagnostics

| Assumption | Method | Result | Conclusion |
|---|---|---|---|
| Normal residuals | Shapiro–Wilk | W = 0.978, p = 0.225 | Satisfied |
| Equal variances | Folded F / Levene | p = 0.33 / p = 0.43 | Satisfied |
| Parallel slopes | Treatment × baseline interaction | p = 0.70 | Satisfied |
| Outcome distribution | Skewness of change | Drug X 0.01, placebo −0.34 | Approximately symmetric |

### 7. Effect size and robustness

| Measure | Estimate | 95% CI | p-value |
|---|---|---|---|
| Cohen's d | 2.19 | 1.61 to 2.77 | — |
| Hedges' g (small-sample corrected) | 2.17 | — | — |
| Mann–Whitney U (non-parametric check) | U = 93 | — | < 0.0001 |
| Unadjusted mean difference | −23.30 mg/dL | −28.23 to −18.36 | < 0.0001 |
| ANCOVA adjusted difference | −23.27 mg/dL | −28.23 to −18.30 | < 0.0001 |

The t-test, ANCOVA and rank-based test all reach the same conclusion, so the result does not depend on the choice of method or on a normality assumption.

### 8. Exploratory responder analysis (reduction ≥ 20 mg/dL at Week 12)

| Measure | Value | 95% CI |
|---|---|---|
| Responders, Drug X | 31 / 35 (88.6%) | — |
| Responders, placebo | 6 / 39 (15.4%) | — |
| Risk difference | 73.2 percentage points | 57.7 to 88.7 |
| Relative risk | 5.76 | 2.73 to 12.13 |
| Odds ratio | 42.6 | 11.0 to 165.6 |
| Fisher exact test | — | p < 0.0001 |
| Number needed to treat | 1.4 | — |

The responder threshold is descriptive and was not pre-specified. Subgroup check by sex (exploratory): Drug X minus placebo difference was −21.7 mg/dL in males and −25.2 mg/dL in females, both in the same direction as the overall result.

### 9. Missing data

Six patients (5 Drug X, 1 placebo) have no Week 12 value; they were not imputed. Their baseline glucose averaged 171.9 mg/dL (Drug X) and 168.1 mg/dL (placebo), close to the overall mean of 165.8, and the dropout imbalance between arms was not statistically significant (Fisher exact p = 0.20). A completers analysis is therefore a reasonable primary analysis here, with MMRM as the recommended sensitivity analysis.

Statistics in sections 7 and 8 other than the t-test and ANCOVA values were calculated as supplementary checks outside the SAS program. SAS equivalents are `PROC NPAR1WAY` (Mann–Whitney) and `PROC FREQ` with the `RISKDIFF`, `RELRISK` and `FISHER` options (responder analysis).

## SAS procedures and skills demonstrated

`PROC IMPORT` · `DATA` step · `MERGE` / `PROC SORT` · `PROC CONTENTS` · `PROC MEANS` · `PROC FREQ` · `PROC TTEST` · `PROC GLM` · `PROC PRINT` · `PROC SGPLOT` · `PROC EXPORT` · `ODS RTF` · SAS macro variables and libraries
## How to run

1. Open SAS Studio (SAS OnDemand for Academics is free).
2. Upload the three raw CSV files into a `raw` folder in your home directory.
3. Open `programs/drugx_t2d_pipeline.sas` and change the path in `%let root = ...;`.
4. Run the program top to bottom and check the Log for errors.
5. Outputs are written to the `sdtm`, `adam` and `output` folders.

## Data note

The datasets are simulated mock-study data created to demonstrate the clinical programming workflow and analysis methods. They contain no real patient information, so results illustrate the methodology and are not evidence about a real compound.

## Limitations and next steps

- Add an **MMRM** analysis (`PROC MIXED`) using all visits as a sensitivity analysis for missing data.
- Define a full analysis set based on all randomised patients.
- Validate ADaM datasets with `PROC COMPARE` against an independent program.
- Add an adverse-event (AE) domain and safety tables.
- Add `define.xml` metadata and an ADaM specification document.

## Author

**Prakhar Gupta**, M.Sc. Statistics and Computing
