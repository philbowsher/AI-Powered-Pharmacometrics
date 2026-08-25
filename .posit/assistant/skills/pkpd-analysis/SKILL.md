---
name: pkpd-analysis
description: Provides exact data schema, required reshaping steps, package API patterns, and a fix for a known nlmixr2 convergence failure for this workshop's pk_data. Use whenever writing R code for PK/PD analysis, NCA, population PK modeling, dose simulation, or diagnostics in this project.
---

# PK/PD Analysis Methodology

## Workshop structure

This workshop is split across four Quarto documents (created by `create_workshop_docs.R`), not one: `01_data_setup.qmd`, `02_nca.qmd`, `03_modeling.qmd`, `04_simulation.qmd`. Each has its section headers already in place — write code into the existing sections rather than proposing a new outline. Data validation/QA is handled by ready-made scripts in `validation/` (run directly, not built via prompting) — don't suggest rebuilding that logic inside the four documents.

## `pk_data` schema

**mAb dataset** (`dataset_choice <- "mab"`): `SUBJID`, `ARM`, `COHORT`, `DOSE`, `AGE`, `SEX`, `WEIGHT`, `TIME_DAYS`, `CONC` (NA = BLQ), `BLQ_FLAG`, `CRP` (PD endpoint). One row per subject per timepoint, `DOSE` repeated on every row — this is **not** modeling-ready, see below.

**Warfarin dataset** (`dataset_choice <- "warfarin"`, `nlmixr2data::warfarin`): `id`, `time`, `amt`, `dv`, `dvid` (`"cp"` = concentration, `"pca"` = anticoagulant effect — **mixed in the same long dataset**), `evid`, `wt`, `age`, `sex`. No `cmt` column — add one (1 for all rows is fine for a 1-compartment oral model) if a package requires it. Already has `amt`/`evid` dosing-event structure, unlike the mAb dataset.

## Reshaping `pk_data` for `nlmixr2` (mAb dataset only — Warfarin already has this)

`nlmixr2` needs one dosing row (`EVID=1`, `AMT`=dose, `CMT`=1/depot) and one row per observation (`EVID=0`, `DV`=concentration, `CMT`=2/central) per subject — not a `DOSE` column repeated on every row. Build it like this before fitting:

```r
dosing_rows <- pk_data |>
  distinct(SUBJID, DOSE, WEIGHT) |>
  mutate(TIME_DAYS = 0, AMT = DOSE * 1000, EVID = 1, DV = NA_real_, CMT = 1)  # see gotcha below

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

**Gotcha (tested, confirmed on the mAb dataset):** `generate_pk_datasets.R` scales dose by 1000 internally when simulating `CONC` (`F_sc * DOSE * 1000 * Ka / V / ...`). If you set `AMT = DOSE` (raw mg) instead of `AMT = DOSE * 1000` when reshaping for `nlmixr2`, the fit will still "succeed" numerically but `CL` and `V` come back ~700x too small while `Ka` stays correct — a dose/concentration unit mismatch masquerading as a fit result, not a convergence failure. If `Ka` looks right but `CL`/`V` look absurd, check dose scaling before anything else.

**Also expected:** `nlmixr2` may throw `Error calculating covariance via linearization` even on a good fit — this is common with small/sparse data and doesn't invalidate the point estimates (`fixef(fit)` still works). Wrap the covariance step, not the whole fit, if you want to handle it gracefully:
```r
fit <- nlmixr(model_mab, nlmix_data, est = "saem", control = saemControl(print = 0))
if (inherits(try(fit$cov, silent = TRUE), "try-error")) {
  message("Covariance failed (common with sparse data) - point estimates below are still valid")
}
exp(fixef(fit))  # CL, V, Ka on the natural scale
```

## `PKNCA` — object-based API, not plain dplyr

```r
library(PKNCA)
conc_data <- pk_data |> filter(!is.na(CONC))  # drop the BLQ/NA pre-dose row, see below
dose_data <- pk_data |> distinct(SUBJID, DOSE) |> mutate(TIME_DAYS = 0)

conc_obj <- PKNCAconc(conc_data, CONC ~ TIME_DAYS | SUBJID)
dose_obj <- PKNCAdose(dose_data, DOSE ~ TIME_DAYS | SUBJID)

# Explicit intervals starting at the first REAL measurement time (see gotcha below)
my_intervals <- data.frame(start = 1, end = Inf, auclast = TRUE, cmax = TRUE,
                            tmax = TRUE, half.life = TRUE, aucinf.obs = TRUE)
data_obj <- PKNCAdata(conc_obj, dose_obj, intervals = my_intervals)
results  <- pk.nca(data_obj)
summary(results)
```
`pk.nca()` returns Cmax, Tmax, AUClast, half-life, etc. per subject automatically — don't hand-write trapezoidal AUC or half-life regression.

**Gotcha (tested, confirmed on the mAb dataset):** the pre-dose `TIME_DAYS=0` row has `CONC=NA` (BLQ). If you let `PKNCAdata()` auto-derive intervals, it anchors the AUC start to the dose time (0), which precedes the first real measurement (day 1) — this throws dozens of "before the first measurement" warnings and `aucinf.obs` comes back `NC` (not calculated). Fix: filter out `NA` concentration rows before building `PKNCAconc`, and pass an explicit `intervals` data frame with `start` set to the first real observation time, not 0.

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
Use the fitted `nlmixr2` (or NCA-derived) parameter estimates as `$PARAM` values — not arbitrary numbers. If those parameters came from fitting the mAb dataset, the same `AMT = DOSE * 1000` scaling gotcha above applies here too — simulate doses in the same scaled units the model was fit on, or Cmax/AUC will come out ~1000x wrong.

## `pointblank` — agent/interrogate pattern, not manual if-checks

```r
agent <- pointblank::create_agent(pk_data) |>
  pointblank::col_vals_gte(columns = "CONC", value = 0, na_pass = TRUE) |>
  pointblank::col_vals_not_null(columns = "SUBJID") |>
  pointblank::interrogate()
agent  # prints a pass/fail report
```

## Dose-response summaries — compare to each subject's own baseline

**Gotcha (tested):** comparing a cohort's min PD value to that cohort's *own* max (e.g. `min(CRP) / max(CRP)` within a dose group) mixes between-subject baseline variability into the result and can produce a flat or nonsensical dose-response curve (identical values across doses). Compute suppression relative to each *subject's own* pre-dose baseline first, then average by cohort:

```r
pk_data |>
  group_by(SUBJID, ARM, DOSE) |>
  mutate(baseline = CRP[TIME_DAYS == 0]) |>
  summarise(pct_suppression = 100 * (1 - min(CRP, na.rm = TRUE) / baseline[1]), .groups = "drop") |>
  group_by(ARM, DOSE) |>
  summarise(mean_pct_suppression = mean(pct_suppression), .groups = "drop")
```
On the mAb dataset this correctly shows suppression rising slightly with dose (~91% → ~93%) rather than a flat/backwards line — expected, since simulated concentrations are far above the PD model's IC50 at every tested dose.

## Plotting `pk_data` — irregular sampling days break default axis breaks

**Gotcha (tested):** `TIME_DAYS` values are irregular (0, 1, 3, 7, 14, 21, 28, 42, 56). Default `ggplot2` continuous x-axis breaks produce overlapping, unreadable tick labels. Fix with explicit breaks and rotated text:
```r
scale_x_continuous(breaks = unique(pk_data$TIME_DAYS)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))
```
When combining a PK and PD plot side by side with `patchwork`, also collect one shared legend instead of duplicating it in both panels — otherwise the legend eats most of the panel width:
```r
(p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
```

## Package installs

This workshop's environment (`dev.workshop.posit.team`) is already configured to use:
- CRAN: `https://packagemanager.posit.co/cran/__linux__/jammy/latest`
- Bioconductor: `https://p3m.dev/bioconductor`

Never pass `repos=` to `install.packages()` or call `options(repos=...)` — that overrides these and can pull from a public CRAN mirror instead.
