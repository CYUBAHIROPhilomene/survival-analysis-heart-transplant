# Stanford Heart Transplant Survival Analysis

## Overview
This repository contains the code and presentation materials for a
parametric survival analysis of the **Stanford Heart Transplant Study**
using the `stanford2` dataset from R's `survival` package.

This project was completed as part of the **Survival Analysis course**
at the **African Institute for Mathematical Sciences (AIMS), Rwanda**.

---

## Group 6 Members
- Philomene CYUBAHIRO
- Peter Onimisi ATTAIGU
- Foumban Marah Khadij NOURIA
- Yosamu MUHANZI
- Zo Lalaina Andrianina ANDRIANANTENAINA

---

## Dataset
- **Name:** `stanford2`
- **Source:** R `survival` package (Escobar & Meeker, 1992)
- **Description:** Survival data from the Stanford Heart Transplant
  Programme (1967--1980)
- **Observations:** 184 patients
- **Events (deaths):** 113
- **Censored:** 71

### Key Variables
| Variable | Description |
|----------|-------------|
| `time`   | Follow-up time in days |
| `status` | 1 = death, 0 = censored |
| `age`    | Age at transplant (years) |
| `t5`     | T5 mismatch score |
| `id`     | Patient ID |

---

## Analysis Workflow

### Step 1 — Data Exploration
- Load and inspect the `stanford2` dataset
- Compute sample size, events, and censored counts
- Summarise survival time and age variables
- Check for missing values

### Step 2 — Data Preparation
- Confirm event coding (1 = death, 0 = censored)
- Remove rows with missing `time` or `status`
- Create survival object: `Surv(time, status)`

### Step 3 — Model Fitting
Six parametric models fitted using `flexsurvreg()`:
- Log-Normal
- Generalised Gamma
- Log-Logistic
- Weibull
- Gamma
- Exponential


### Step 4 — Model Selection
- Models compared using **AIC** (Akaike Information Criterion)
- **Log-Normal** selected as the best-fitting model (AIC = 1741.611)

### Step 5 — Life Functions
Computed for the Log-Normal model:
- Survival function S(t)
- Hazard function h(t)
- Density function f(t)
- Cumulative distribution function F(t)
- Mean and variance of survival time

### Step 6 — MLE Parameter Estimates
| Parameter | Estimate | SE | 95% CI |
|-----------|----------|----|--------|
| μ (meanlog) | 6.280 | 0.202 | (5.885, 6.675) |
| σ (sdlog) | 2.431 | 0.172 | (2.116, 2.792) |

### Step 7 — Interpretation
- Median survival ≈ **534 days**
- Mean survival ≈ **10,248 days**
- Non-monotone hazard: risk peaks early then declines
- High variance reflects substantial heterogeneity in outcomes

---


## Requirements

### R Packages
```r
install.packages("survival")
install.packages("flexsurv")
```

### LaTeX Packages
- `beamer`
- `biblatex`
- `booktabs`
- `colortbl`
- `palatino`

---

## How to Run
1. Clone the repository
2. Open `group6_SA1.R` in RStudio
3. Run the script from top to bottom
5. Compile `group6_SA1.tex` in Overleaf or locally


