# Quick, standalone data validation for pk_data. Run any time -- takes a
# few seconds, no compiling. Not a prompting exercise -- just run it.
# Works for either dataset (mAb or Warfarin) -- detects which is loaded.

library(tidyverse)
library(pointblank)

if (!exists("pk_data")) {
  stop("pk_data not found in this session -- run generate_pk_datasets.R first.")
}

is_mab_schema <- "SUBJID" %in% names(pk_data)

if (is_mab_schema) {
  agent <- create_agent(pk_data, label = "mAb dataset validation") |>
    col_vals_gte(columns = "CONC", value = 0, na_pass = TRUE) |>
    col_vals_not_null(columns = c("SUBJID", "DOSE", "TIME_DAYS")) |>
    col_vals_in_set(columns = "BLQ_FLAG", set = c("BLQ", "Quantifiable")) |>
    col_vals_between(columns = "AGE", left = 18, right = 100) |>
    interrogate()
} else {
  agent <- create_agent(pk_data, label = "Warfarin dataset validation") |>
    col_vals_gte(columns = "dv", value = 0, na_pass = TRUE) |>
    col_vals_not_null(columns = c("id", "time", "dvid")) |>
    col_vals_in_set(columns = "dvid", set = c("cp", "pca")) |>
    interrogate()
}

print(agent)

if (all_passed(agent)) {
  cat("\n[ok] All validation checks passed.\n")
} else {
  cat("\n[warning] Some validation checks failed -- see the report above.\n")
}
