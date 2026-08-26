---
name: pkpd-analysis
description: Provides exact data schema, required reshaping steps, package API patterns, and a fix for a known nlmixr2 convergence failure for this workshop's pk_data. Use whenever writing R code for PK/PD analysis, NCA, population PK modeling, dose simulation, or diagnostics in this project.
---

# PK/PD Analysis Methodology

## Workshop structure

This workshop is split across four Quarto documents (created by `create_workshop_docs.R`), not one: `01_data_setup.qmd`, `02_nca.qmd`, `03_modeling.qmd`, `04_simulation.qmd`. Each has its section headers already in place — write code into the existing sections rather than proposing a new outline.

Infrastructure scripts already exist for these tasks — point to them rather than writing new code for the same job: `check_packages.R` (installs anything missing), `generate_pk_datasets.R` (loads `pk_data` — this is how it gets into the session in the first place), `modeling_helpers.R` (see Document 3 below), `validation/validate_data.R` and `validation/quality_check.R` (data QA, run directly, not built via prompting).

**Run `check_packages.R` proactively, don't wait for a failure to reveal it's needed.** If a package fails to load or `pk_data` doesn't exist, the fix is almost always "run `check_packages.R`," not troubleshooting the symptom in place. (`generate_pk_datasets.R` now runs it automatically as a safety net, but don't rely on that — suggest it explicitly at the start of a session.)

## `pk_data` schema

**mAb dataset** (`dataset_choice <- "mab"`): `SUBJID`, `ARM`, `COHORT`, `DOSE`, `AGE`, `SEX`, `WEIGHT`, `TIME_DAYS`, `CONC` (NA = BLQ), `BLQ_FLAG`, `CRP` (PD endpoint). One row per subject per timepoint, `DOSE` repeated on every row — this is **not** modeling-ready, see below.

**Warfarin dataset** (`dataset_choice <- "warfarin"`, `nlmixr2data::warfarin`): `id`, `time`, `amt`, `dv`, `dvid` (`"cp"` = concentration, `"pca"` = anticoagulant effect — **mixed in the same long dataset**), `evid`, `wt`, `age`, `sex`. No `cmt` column — add one (1 for all rows is fine for a 1-compartment oral model) if a package requires it. Already has `amt`/`evid` dosing-event structure, unlike the mAb dataset.

## Document 3 (Population PK Modeling & Diagnostics): use the provided helpers, don't regenerate reshaping

`source("modeling_helpers.R")` provides three tested functions — **use them instead of writing the reshaping logic from scratch**, since that's exactly where a silent dose-unit bug was found during testing:

```r
source("modeling_helpers.R")
nlmix_data <- reshape_for_nlmixr(pk_data)                       # handles mAb reshape + dose scaling, or Warfarin PK filter
fit <- fit_pk_model(nlmix_data, drug_class = "mab")             # or "warfarin" -- uses correct starting values
check_fit(fit)                                                  # prints estimates + sanity checks
```

`reshape_for_nlmixr()` builds one dosing row (`EVID=1`, `AMT`=dose, `CMT`=1/depot) and one row per observation (`EVID=0`, `DV`=concentration, `CMT`=2/central) per subject for the mAb dataset (auto-scaling `AMT` by 1000 to match how `CONC` was simulated), or filters to `dvid == "cp"` for Warfarin (dropping the `"pca"` PD rows, which would otherwise corrupt a PK fit). `fit_pk_model()` uses the correct drug-class-scaled starting values (see below) automatically. Still write your own `ini()`/`model()` block with Assistant if you want to see how nlmixr2 model syntax works — the helpers exist to remove the risky plumbing, not the learning.

## `nlmixr2` starting values — use drug-class-appropriate values

Generic oral starting values (`CL~10`, `V~80`, `Ka~0.8`) are scaled for small molecules, not biologics. If a mAb model fails to converge, mismatched starting values are the first thing to check — but SAEM often converges fine even from imperfect starting points, so don't assume failure. Use mAb-scale values as the default for that dataset regardless:

```r
# mAb (use nlmix_data from reshape_for_nlmixr() above, or write your own model to compare)
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

## Fitted `nlmixr2` objects already have DV/PRED/IPRED — don't reach for NONMEM-style helpers

A fitted `nlmixr2` object is itself a data frame with `DV`, `PRED`, `IPRED`, `CWRES`, etc. as columns (`names(fit)` shows them directly; `ggplot(fit, aes(PRED, DV))` just works). Don't default to NONMEM-style helper functions like `augPred()` — that's a different package/API shape and will error (e.g. `EVID not found`) on an `nlmixr2` fit object. Check `names(fit)` first if unsure what's available.

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

## `nlmixr2` install failure (tested, confirmed in this environment)

`check_packages.R` may fail to install `nlmixr2` with a compilation error in its
`Deriv` dependency (`error: 'R_ClosureFormals' was not declared in this scope`,
`Rf_allocLang` also missing). This is an R-internals API mismatch: the current
CRAN release of `Deriv` (4.3.0+) requires R-internals functions not present in
this environment's R 4.4.0 headers, and Posit Package Manager has no compiled
binary for `Deriv` on this platform/R-version combination to fall back on
(confirmed via `available.packages(..., type = "binary")` — it still resolves
to `src/contrib`, so this is not a missing-`HTTPUserAgent` issue — that was
checked and ruled out, don't re-test it).

Fix: install an older, pure-R (pre-C++) `Deriv` from the CRAN archive first,
then let the rest of the `nlmixr2` chain install normally:
```r
download.file("https://cran.r-project.org/src/contrib/Archive/Deriv/Deriv_4.1.2.tar.gz",
              "/tmp/Deriv_4.1.2.tar.gz")
install.packages("/tmp/Deriv_4.1.2.tar.gz", repos = NULL, type = "source")
install.packages(c("xgxr", "nlmixr2plot", "nlmixr2"), dependencies = FALSE)
```
Run `source("check_packages.R")` again afterward to confirm `nlmixr2` now
shows `[ok]`. `check_packages.R` already attempts a simpler version of this fix
automatically (install Deriv 4.1.2, no `dependencies = FALSE` step) — if that's
not enough on a given environment, use the fuller sequence above. This is a
stopgap recorded here, not a permanent fix baked into the base image — it will
recur on a genuinely fresh environment until someone updates the install
script itself to pin a compatible `Deriv` version by default.

## Deploying to Posit Connect — two non-interoperable credential paths

If a document gets deployed (e.g. after the Dashboard step), there are two separate, non-interoperable ways to do it: Positron's Publisher UI ("Publish" button) and the R `rsconnect` package driven from the console. **Ask which one the user wants before starting** — don't assume. If someone already clicked "Publish" and it's incomplete, check for a leftover config file before starting a fresh `rsconnect`-driven deployment from scratch; re-derive the server URL from it rather than asking the user to re-enter it.

## Package installs

**See "`nlmixr2` install failure" above** if `Deriv`/`nlmixr2` fails to compile — that section has the full diagnosis and fix.

`check_packages.R` pins the repo to a **dated** Posit Package Manager snapshot (not `/latest` — see the gotcha above for why) at the top of the script. Don't second-guess or change that pin mid-session, and don't pass `repos=` to `install.packages()` or call `options(repos=...)` elsewhere — that would override it and risks pulling from a public CRAN mirror instead. If the pinned date needs updating (e.g. a package genuinely needs a newer version), that's a deliberate edit to `check_packages.R` itself, not a one-off override in analysis code.
