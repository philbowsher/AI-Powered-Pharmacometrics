---
name: pkpd-analysis
description: Provides exact data schema, package choices, and a fix for a known nlmixr2 convergence failure for this workshop's pk_data. Use whenever writing R code for PK/PD analysis, NCA, population PK modeling, dose simulation, or diagnostics in this project.
---

# PK/PD Analysis Methodology

## `pk_data` schema — use exactly these column names, don't guess

**mAb dataset** (`generate_pk_datasets.R`, `dataset_choice <- "mab"`):
`SUBJID`, `ARM`, `COHORT`, `DOSE`, `AGE`, `SEX`, `WEIGHT`, `TIME_DAYS`, `CONC` (NA = BLQ), `BLQ_FLAG`, `CRP` (PD endpoint).

**Warfarin dataset** (`dataset_choice <- "warfarin"`, from `nlmixr2data::warfarin`):
`id`, `time`, `dv`, `dvid` (`"cp"` = concentration, `"pca"` = PD effect), `amt`, `wt`. Lowercase, different shape than the mAb set — check which is loaded before writing code, don't assume mAb's column names apply.

`CONC`/`dv` is already `NA` for BLQ samples in the mAb set. Don't zero-fill it.

## Package choice — non-default packages to use, not substitute

| Task | Use | Not |
|---|---|---|
| NCA (Cmax, Tmax, AUC, half-life) | `PKNCA` | hand-rolled trapezoidal code |
| Population PK modeling | `nlmixr2`, SAEM | `nlme`, `nls`, or other fitting packages |
| VPC / GOF diagnostics | `vpc` | manual percentile plots, unless `vpc` fails |
| Dose-scenario simulation | `mrgsolve` | manual ODE solving |

## Known issue: nlmixr2 fails on the mAb dataset with default starting values

Generic oral-drug starting values (`CL~10`, `V~80`, `Ka~0.8`) are 10-30x too large for the mAb dataset and are the most common cause of SAEM non-convergence here — not evidence the data can't be modeled. Use mAb-scale values instead:

```r
model_mab <- function() {
  ini({
    tCL <- log(0.3); tV <- log(5); tKa <- log(0.35)
    eta.CL ~ 0.1; eta.V ~ 0.1; eta.Ka ~ 0.1
    prop.err <- 0.2
  })
  model({
    CL <- exp(tCL + eta.CL); V <- exp(tV + eta.V); Ka <- exp(tKa + eta.Ka)
    ke <- CL / V
    d/dt(depot)   <- -Ka * depot
    d/dt(central) <- Ka * depot - ke * central
    cp <- central / V
    cp ~ prop(prop.err)
  })
}
```

If it still won't converge after this, that's a valid result — report it and use NCA instead, which is the primary method for this drug class regardless.

## Package installs

This project's default repo is already configured for this workshop's environment (`dev.workshop.posit.team`). Never pass `repos=` to `install.packages()` or override `options(repos=...)` — that would bypass it.
