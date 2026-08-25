# Helper functions for Document 3 (Population PK Modeling & Diagnostics).
# Source this once: source("modeling_helpers.R")
#
# These handle the fiddly, error-prone data reshaping so you can focus on
# the actual modeling. The reshaping logic here is deterministic and
# tested -- don't ask Posit Assistant to regenerate it, just call it.

library(tidyverse)

#' Reshape pk_data into nlmixr2's event-based format.
#' Auto-detects which dataset is loaded (mAb vs Warfarin) and applies the
#' correct dose scaling / column names for each.
reshape_for_nlmixr <- function(pk_data) {
  if ("SUBJID" %in% names(pk_data)) {
    # mAb dataset -- one row per timepoint, needs dosing rows constructed,
    # and AMT scaled by 1000 to match how CONC was simulated.
    dosing_rows <- pk_data |>
      distinct(SUBJID, DOSE, WEIGHT) |>
      mutate(TIME_DAYS = 0, AMT = DOSE * 1000, EVID = 1, DV = NA_real_, CMT = 1)

    obs_rows <- pk_data |>
      mutate(AMT = 0, EVID = 0, DV = CONC, CMT = 2) |>
      filter(!is.na(DV))

    bind_rows(dosing_rows, obs_rows) |>
      arrange(SUBJID, TIME_DAYS) |>
      rename(ID = SUBJID, TIME = TIME_DAYS, WT = WEIGHT)

  } else {
    # Warfarin dataset -- already event-structured, just drop the PD rows.
    pk_data |>
      filter(dvid == "cp") |>
      rename(ID = id, TIME = time, DV = dv, AMT = amt, EVID = evid, WT = wt)
  }
}

#' Fit a 1-compartment model with drug-class-appropriate starting values.
#' drug_class: "mab" or "warfarin" (matches dataset_choice in generate_pk_datasets.R).
#' NOTE: nlmixr2's ini() block uses non-standard evaluation -- it cannot be
#' parametrized with bquote()/substitute(). Two literal model functions,
#' not one dynamically-built one.
fit_pk_model <- function(nlmix_data, drug_class = c("mab", "warfarin")) {
  drug_class <- match.arg(drug_class)
  library(nlmixr2)

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

  model <- if (drug_class == "mab") model_mab else model_warfarin
  nlmixr(model, nlmix_data, est = "saem", control = saemControl(print = 0))
}

#' Quick sanity check on a fitted model -- is this worth trusting?
check_fit <- function(fit) {
  p <- exp(fixef(fit))
  cat("Parameter estimates (natural scale):\n")
  print(p)

  cov_ok <- !inherits(try(fit$cov, silent = TRUE), "try-error") && !is.null(fit$cov)
  cat("\nCovariance calculated:", cov_ok,
      if (!cov_ok) "(common with sparse data -- point estimates above are still valid)\n" else "\n")
  cat("(You may see an on-screen 'Error calculating covariance via linearization' message\n")
  cat(" even when this says TRUE -- nlmixr2 tried that method, it failed, then a fallback\n")
  cat(" method succeeded. Trust fit$cov, not the printed message.)\n")

  cat("\nSanity checks:\n")
  cat("  Ka in [0.1, 2]:  ", p["tKa"] > 0.1 && p["tKa"] < 2, "\n")
  cat("  CL, V both > 0:  ", p["tCL"] > 0 && p["tV"] > 0, "\n")
  cat("\nIf CL/V look wildly implausible (near-zero or huge) but Ka looks fine,\n")
  cat("suspect a dose-unit mismatch in reshape_for_nlmixr(), not a bad fit.\n")
}
