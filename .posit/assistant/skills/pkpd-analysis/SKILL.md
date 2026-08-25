---
name: pkpd-analysis
description: Provides exact data schema, required reshaping steps, package API patterns, and a fix for a known nlmixr2 convergence failure for this workshop's pk_data. Use whenever writing R code for PK/PD analysis, NCA, population PK modeling, dose simulation, or diagnostics in this project.
---

# PK/PD Analysis Methodology

## `pk_data` schema

**mAb dataset** (`dataset_choice <- "mab"`): `SUBJID`, `ARM`, `COHORT`, `DOSE`, `AGE`, `SEX`, `WEIGHT`, `TIME_DAYS`, `CONC` (NA = BLQ), `BLQ_FLAG`, `CRP` (PD endpoint). One row per subject per timepoint, `DOSE` repeated on every row — this is **not** modeling-ready, see below.

**Warfarin dataset** (`dataset_choice <- "warfarin"`, `nlmixr2data::warfarin`): `id`, `time`, `dv`, `dvid` (`"cp"` = concentration, `"pca"` = anticoagulant effect — **mixed in the same long dataset**), `amt`, `wt`, `evid`, `cmt`. Already event-structured for nlmixr2.

## Reshaping `pk_data` for `nlmixr2` (mAb dataset only — Warfarin already has this)

`nlmixr2` needs one dosing row (`EVID=1`, `AMT`=dose, `CMT`=1/depot) and one row per observation (`EVID=0`, `DV`=concentration, `CMT`=2/central) per subject — not a `DOSE` column repeated on every row. Build it like this before fitting:

```r
dosing_rows <- pk_data |>
  distinct(SUBJID, DOSE, WEIGHT) |>
  mutate(TIME_DAYS = 0, AMT = DOSE, EVID = 1, DV = NA_real_, CMT = 1)

obs_rows <- pk_data |>
  mutate(AMT = 0, EVID = 0, DV = CONC, CMT = 2) |>
  filter(!is.na(DV))  # drop BLQ rows, don't feed NA concentrations to nlmixr2

nlmix_data <- bind_rows(dosing_rows, obs_rows) |>
  arrange(SUBJID, TIME_DAYS) |>
  rename(ID = SUBJID, TIME = TIME_DAYS, WT = WEIGHT)
```

## Warfarin: filter to PK before modeling concentration

```r
warfarin_pk <- nlmixr2data::warfarin |> filter(dvid == "cp")
```
Don't fit a PK model on the unfiltered dataset — `dvid == "pca"` rows are the PD effect, not concentration, and will corrupt the fit.

## `nlmixr2` starting values — use drug-class-appropriate values

Generic oral starting values (`CL~10`, `V~80`, `Ka~0.8`) are scaled for small molecules, not biologics. If a mAb model fails to converge, mismatched starting values are the first thing to check — but SAEM often converges fine even from imperfect starting points, so don't assume failure. Use mAb-scale values as the default for that dataset regardless:

```r
# mAb (use nlmix_data built above)
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

```r
# Warfarin (use warfarin_pk, already event-structured)
model_warfarin <- function() {
  ini({
    tCL <- log(10); tV <- log(80); tKa <- log(0.8)
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

If it doesn't converge, that's a valid result — report it and use NCA instead, which is the primary method for this drug class regardless of whether modeling succeeds.

## `PKNCA` — object-based API, not plain dplyr

```r
library(PKNCA)
conc_obj <- PKNCAconc(obs_rows, DV ~ TIME_DAYS | SUBJID)
dose_obj <- PKNCAdose(dosing_rows, AMT ~ TIME_DAYS | SUBJID)
data_obj <- PKNCAdata(conc_obj, dose_obj)
results  <- pk.nca(data_obj)
summary(results)
```
`pk.nca()` returns Cmax, Tmax, AUClast, half-life, etc. per subject automatically — don't hand-write trapezoidal AUC or half-life regression.

## `vpc` — needs explicit column mapping, can fail silently

```r
vpc_result <- tryCatch(
  vpc::vpc(sim = fitted_model, obs = nlmix_data,
           obs_cols = list(dv = "DV", idv = "TIME", id = "ID"),
           sim_cols = list(dv = "DV", idv = "TIME", id = "ID"),
           pi = c(0.05, 0.95)),
  error = function(e) NULL
)
```
If `vpc_result` is `NULL`, fall back to a manual percentile-band plot: simulate several replicates with `rxode2::rxSolve()`, compute 5th/50th/95th percentiles per time bin with dplyr, plot with `geom_ribbon()`.

## `mrgsolve` — model code is a text block (a DSL), not an R function

```r
mod <- mrgsolve::mcode("pk_model", '
$PARAM CL = 0.3, V = 5, KA = 0.35
$CMT DEPOT CENTRAL
$ODE
dxdt_DEPOT = -KA * DEPOT;
dxdt_CENTRAL = KA * DEPOT - (CL/V) * CENTRAL;
$TABLE
capture CP = CENTRAL / V;
')
mod |> mrgsolve::ev(amt = 300, ii = 24, addl = 0) |>
  mrgsolve::idata_set(data.frame(ID = 1:100)) |>
  mrgsolve::mrgsim(end = 56, delta = 1) |> as_tibble()
```
Use the fitted `nlmixr2` (or NCA-derived) parameter estimates as `$PARAM` values — not arbitrary numbers.

## `pointblank` — agent/interrogate pattern, not manual if-checks

```r
agent <- pointblank::create_agent(pk_data) |>
  pointblank::col_vals_gte(columns = "CONC", value = 0, na_pass = TRUE) |>
  pointblank::col_vals_not_null(columns = "SUBJID") |>
  pointblank::interrogate()
agent  # prints a pass/fail report
```

## Package installs

This workshop's environment (`dev.workshop.posit.team`) is already configured to use:
- CRAN: `https://packagemanager.posit.co/cran/__linux__/jammy/latest`
- Bioconductor: `https://p3m.dev/bioconductor`

Never pass `repos=` to `install.packages()` or call `options(repos=...)` — that overrides these and can pull from a public CRAN mirror instead.
