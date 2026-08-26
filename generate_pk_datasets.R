# Run this once, before creating any Quarto document, to get your PK/PD
# dataset loaded into your R session. Edit `dataset_choice` below, then
# source this whole script (Ctrl+Shift+S / Cmd+Shift+S, or run line by line).

dataset_choice <- "mab"  # "mab" or "warfarin"

# Safety net: if check_packages.R was skipped, catch it here rather than
# failing deeper into a later document under more time pressure.
if (file.exists("check_packages.R")) source("check_packages.R")

library(tidyverse)

if (dataset_choice == "mab") {

  # ABC-MAB-101: Phase 1 SAD study of a subcutaneous anti-IL6 monoclonal
  # antibody. 3 dose cohorts, 24 subjects, PK + CRP-suppression PD.
  # NCA is the natural primary analysis here (long half-life biologic);
  # population PK modeling is optional and may not converge well.

  set.seed(20260825)

  n_per_cohort <- 8
  doses <- c(100, 300, 600)  # mg, SC
  cohort_names <- c("Cohort 1: 100mg SC", "Cohort 2: 300mg SC", "Cohort 3: 600mg SC")
  sampling_days <- c(0, 1, 3, 7, 14, 21, 28, 42, 56)

  subjects <- tibble(
    SUBJID = sprintf("MAB101-%03d", 1:(length(doses) * n_per_cohort)),
    COHORT = rep(1:length(doses), each = n_per_cohort),
    DOSE   = rep(doses, each = n_per_cohort),
    ARM    = rep(cohort_names, each = n_per_cohort)
  ) |>
    mutate(
      AGE    = pmin(pmax(round(rnorm(n(), 42, 12)), 21), 68),
      SEX    = sample(c("M", "F"), n(), replace = TRUE, prob = c(0.6, 0.4)),
      WEIGHT = round(pmin(pmax(ifelse(SEX == "M", rnorm(n(), 82, 12), rnorm(n(), 68, 10)), 50), 120), 1),
      eta_CL = rnorm(n(), 0, 0.25),
      eta_V  = rnorm(n(), 0, 0.20),
      eta_Ka = rnorm(n(), 0, 0.30)
    )

  CL_pop <- 0.3; V_pop <- 5.0; Ka_pop <- 0.35; F_sc <- 0.70  # typical mAb SC parameters

  pk_data <- subjects |>
    crossing(TIME_DAYS = sampling_days) |>
    mutate(
      CL_i = CL_pop * exp(eta_CL) * (WEIGHT / 70)^0.75,
      V_i  = V_pop  * exp(eta_V)  * (WEIGHT / 70),
      Ka_i = Ka_pop * exp(eta_Ka),
      ke   = CL_i / V_i,
      CONC_PRED = if_else(
        TIME_DAYS == 0, 0,
        (F_sc * DOSE * 1000 * Ka_i / V_i / (Ka_i - ke)) * (exp(-ke * TIME_DAYS) - exp(-Ka_i * TIME_DAYS))
      ),
      CONC_OBS = pmax(CONC_PRED * (1 + rnorm(n(), 0, 0.15)) + rnorm(n(), 0, 0.05), 0),
      LLOQ = 0.1,
      BLQ_FLAG = if_else(CONC_OBS < LLOQ, "BLQ", "Quantifiable"),
      CONC = if_else(BLQ_FLAG == "BLQ", NA_real_, CONC_OBS)
    ) |>
    group_by(SUBJID) |>
    mutate(
      CRP_baseline = pmax(first(rnorm(1, 5.5, 1.5)), 2.0),
      Imax = 0.90, IC50 = 5.0,
      CRP_effect = Imax * CONC_OBS / (IC50 + CONC_OBS),
      CRP = round(CRP_baseline * (1 - CRP_effect) * exp(rnorm(n(), 0, 0.20)), 2)
    ) |>
    ungroup() |>
    select(SUBJID, ARM, COHORT, DOSE, AGE, SEX, WEIGHT, TIME_DAYS, CONC, BLQ_FLAG, CRP)

  dir.create("data", showWarnings = FALSE)
  write_csv(pk_data, "data/mab_pkpd.csv")
  cat("Loaded ABC-MAB-101 into `pk_data` and wrote data/mab_pkpd.csv (", nrow(pk_data), "rows )\n")
  cat("Primary analysis: NCA. PopPK (nlmixr2) is optional and may not converge well for this drug.\n")

} else if (dataset_choice == "warfarin") {

  # nlmixr2data::warfarin - a well-behaved, richly-sampled oral small
  # molecule. Good for a full population PK model that actually converges.
  # Unlike "mab" (fully self-generated), this dataset ships inside the
  # nlmixr2data package, so we install it here if it isn't already available.

  if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
    install.packages("nlmixr2data")
  }
  library(nlmixr2data)
  data(warfarin)
  pk_data <- as_tibble(warfarin)

  dir.create("data", showWarnings = FALSE)
  write_csv(pk_data, "data/warfarin.csv")
  cat("Loaded nlmixr2data::warfarin into `pk_data` and wrote data/warfarin.csv (", nrow(pk_data), "rows )\n")
  cat("Primary analysis: population PK modeling (nlmixr2) - this dataset is well suited to it.\n")

} else {
  stop("Set dataset_choice to 'mab' or 'warfarin' at the top of this script.")
}

glimpse(pk_data)
