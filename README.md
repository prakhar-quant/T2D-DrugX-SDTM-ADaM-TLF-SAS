# Drug X vs Placebo in Type 2 Diabetes: SDTM → ADaM → TLF Pipeline and ANCOVA Analysis in SAS

An end-to-end clinical SAS programming project: raw trial data are converted to **CDISC SDTM**, derived into **ADaM** analysis datasets, and analysed to produce **Tables, Listings and Figures (TLF)**. The primary endpoint is the **change from baseline in fasting glucose (mg/dL) at Week 12**, analysed with a two-sample t-test and an **ANCOVA** adjusted for baseline glucose.

![Pipeline](figures/01_pipeline.png)

## Headline result

| | Drug X (n = 35) | Placebo (n = 39) |
|---|---|---|
| Baseline glucose, mean (SD) | 166.4 (19.3) | 165.2 (18.0) |
| Week 12 glucose, mean (SD) | 135.2 (21.5) | 157.4 (20.5) |
| Change from baseline, mean (SD) | **−31.2 (9.7)** | **−7.9 (11.4)** |
| ANCOVA adjusted mean change (SE) | −31.15 (1.81) | −7.88 (1.71) |

**Treatment difference (ANCOVA, adjusted for baseline): −23.27 mg/dL, 95% CI −28.23 to −18.30, p < 0.0001.**
Unadjusted t-test: −23.30 mg/dL, t(72) = −9.41, p < 0.0001. Cohen's d = 2.19.
Responder rate (reduction ≥ 20 mg/dL, exploratory): 88.6% vs 15.4%.

![Mean change by visit](figures/02_mean_change_by_visit.png)

## Study design

| Item | Detail |
|---|---|
| Design | Two-arm, parallel-group, 1:1 allocation |
| Treatments | Drug X 10 mg (n = 40) vs placebo (n = 40) |
| Population | 80 adults with Type 2 Diabetes, 4 sites |
| Visits | Screening, Baseline, Week 4, Week 8, Week 12 |
| Primary endpoint | Change from baseline in fasting glucose at Week 12 |
| Completion | 74 / 80 (92.5%) had a Week 12 value (5 Drug X, 1 placebo missing) |

## Workflow

1. **Raw data** (`data/raw`): `DM.csv`, `EX.csv`, `LB.csv` imported with `PROC IMPORT`, checked with `PROC CONTENTS`.
2. **SDTM**: `DM`, `EX`, `LB` domains with `USUBJID`, `DOMAIN`, ISO 8601 `--DTC` dates and standardised result `LBSTRESN`.
3. **ADaM**:
   - `ADSL` (80 rows): treatment (`TRT01P`, `TRT01A`), dates, duration, `SAFFL`.
   - `ADLB` (394 rows): `AVAL`, `BASE`, `CHG = AVAL − BASE`, `PCHG`.
4. **Analysis**: `PROC MEANS`, `PROC FREQ`, `PROC TTEST`, `PROC GLM` (ANCOVA with `LSMEANS`).
5. **TLF**: demographics table, Week 12 efficacy table (RTF), subject listing (CSV), mean-change figure (`PROC SGPLOT`).

## Statistical methods

- **Two-sample t-test** (pooled, with Satterthwaite sensitivity check) on Week 12 change.
- **ANCOVA**: `CHG = β0 + β1·Treatment + β2·Baseline + ε`, adjusted means evaluated at mean baseline glucose (165.8 mg/dL).
- **Assumption checks**: residual normality (Shapiro–Wilk p = 0.23, Q–Q plot), equal variances (folded F p = 0.33, Levene p = 0.43), parallel slopes (treatment × baseline p = 0.70).
- Baseline was not a significant covariate (p = 0.69), so adjusted and unadjusted estimates agree within 0.03 mg/dL.

![ANCOVA fit](figures/05_ancova_fit.png)
![Adjusted means and difference](figures/06_adjusted_means_difference.png)

## SAS procedures and skills demonstrated

`PROC IMPORT` · `DATA` step · `MERGE` / `PROC SORT` · `PROC CONTENTS` · `PROC MEANS` · `PROC FREQ` · `PROC TTEST` · `PROC GLM` · `PROC PRINT` · `PROC SGPLOT` · `PROC EXPORT` · `ODS RTF` · SAS macro variables and libraries

## How to run

1. Open **SAS Studio** (SAS OnDemand for Academics is free).
2. Upload `data/raw/*.csv` into a `raw` folder in your home directory.
3. Open `programs/drugx_t2d_pipeline.sas` and change the path in `%let root = ...;`.
4. Run the program top to bottom and check the Log for errors.
5. Outputs are written to the `sdtm`, `adam` and `output` folders.

## Data note

The datasets are simulated mock-study data created to demonstrate the clinical programming workflow and analysis methods. They contain no real patient information, so results illustrate the methodology and are not evidence about a real compound.

## Limitations and next steps

- Week 12 analysis uses observed values (completers); a full analysis set with **MMRM** (`PROC MIXED`) or multiple imputation would be the preferred sensitivity analysis.
- Add dataset validation with `PROC COMPARE` against an independent program.
- Add an adverse-event (AE) domain and a safety table.
- Add `define.xml` metadata and ADaM specification document.

## Author

**Prakhar Gupta**, M.Sc. Statistics and Computing
