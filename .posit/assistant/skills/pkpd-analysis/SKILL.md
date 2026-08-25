---
name: pkpd-analysis
description: Guides PK/PD analysis in this workshop's R/Quarto document — which package to use for each analysis step, data conventions for pk_data, and when NCA vs population PK modeling applies. Use whenever writing R code for pharmacokinetic/pharmacodynamic analysis, non-compartmental analysis (NCA), population PK modeling, dose simulation, or diagnostics in this project.
---

# PK/PD Analysis Methodology

This project analyzes concentration-time (PK) and effect (PD) data loaded into `pk_data` by `generate_pk_datasets.R`. Follow this methodology when writing analysis code.

## Package for each task

Use exactly one package per task below — don't mix alternatives unless the primary choice is unavailable.

| Task | Package |
|---|---|
| Data wrangling, plotting | `tidyverse` (dplyr, ggplot2) |
| Combining plots | `patchwork` |
| Tables | `gt` or `gtsummary` |
| Non-compartmental analysis (Cmax, Tmax, AUC, half-life) | `PKNCA` |
| Population PK modeling | `nlmixr2` (SAEM estimation) |
| Goodness-of-fit / Visual Predictive Check | `vpc` |
| Dose-scenario simulation | `mrgsolve` |
| Data quality validation | `pointblank` |
| Dashboard | `bslib` (Quarto Dashboard) or `shiny`/`shinydashboard` |

## `pk_data` conventions

`pk_data` (produced by `generate_pk_datasets.R`) has these columns: `SUBJID`, `ARM`/`COHORT`, `DOSE`, demographics (`AGE`, `SEX`, `WEIGHT`), `TIME_DAYS`, `CONC`, `BLQ_FLAG`, and a PD endpoint (`CRP` for the mAb dataset). The Warfarin dataset (`nlmixr2data::warfarin`) uses `id`, `time`, `dv`, `dvid` instead — check which one is loaded before assuming column names.

**BLQ handling:** `CONC` is already `NA` for below-limit-of-quantification samples. Never zero-fill or drop these rows silently — treat `NA` as the correct representation, and mention the BLQ count when summarizing.

## NCA vs. population PK modeling

- **NCA (`PKNCA`) is the primary method** for long-half-life biologics (e.g. the mAb dataset) and whenever sampling is sparse relative to the drug's half-life. It's model-independent and always applicable.
- **Population PK modeling (`nlmixr2`) is secondary/optional.** It requires starting values scaled to the drug class:
  - Oral/IV small molecule (e.g. Warfarin): `CL` ~5-15, `V` ~30-100, `Ka` ~0.5-1.5 (typical units, adjust to the data's time scale).
  - Subcutaneous monoclonal antibody: `CL` ~0.2-0.5, `V` ~3-7, `Ka` ~0.2-0.5 — roughly 10-30x smaller than oral small-molecule defaults. Using oral-scale starting values on mAb data is the most common cause of SAEM convergence failure — it's a starting-value problem, not evidence the data is unmodelable.
  - If modeling doesn't converge cleanly after trying appropriate starting values, that's a valid outcome: report it and fall back to NCA rather than forcing a fit.

## General principles

- **Validate before modeling.** Check for negative concentrations, consistent units, and reasonable demographic ranges before fitting anything.
- **Never install packages from `cloud.r-project.org` or other public CRAN mirrors** — always let R use its configured default repository (already set up in this environment).
- **State assumptions.** If you simulate demographics, choose starting values, or make a modeling decision, say so in the text, not just the code.
- **Sanity-check before presenting.** Before finalizing any table or plot, check it against the data it came from — units consistent, numbers in a plausible range, no BLQ artifacts.
