# Quick sanity-check summary of pk_data -- a final gut-check before you
# present results. Run any time; no package beyond tidyverse needed.
# Works for either dataset (mAb or Warfarin) -- detects which is loaded.

library(tidyverse)

if (!exists("pk_data")) {
  stop("pk_data not found in this session -- run generate_pk_datasets.R first.")
}

is_mab_schema <- "SUBJID" %in% names(pk_data)
conc_col <- if (is_mab_schema) "CONC" else "dv"
id_col   <- if (is_mab_schema) "SUBJID" else "id"

cat("=== Quick Quality Check ===\n\n")

neg_conc <- sum(pk_data[[conc_col]] < 0, na.rm = TRUE)
cat(ifelse(neg_conc == 0, "[ok]     ", "[warning]"),
    "No negative concentrations (", neg_conc, "found)\n")

na_count <- sum(is.na(pk_data[[conc_col]]))
cat("[info]    ", na_count, "NA/BLQ concentration values (",
    round(100 * na_count / nrow(pk_data), 1), "%)\n")

n_subj <- n_distinct(pk_data[[id_col]])
cat("[info]    ", n_subj, "unique subjects\n")

cat("\nIf anything above looks wrong, dig into it before presenting your results.\n")
