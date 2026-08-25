# Run this once, early in the workshop, to scaffold your four working
# Quarto documents with section headers already in place. This is
# deterministic setup, not a prompting exercise -- just run it.
# Safe to re-run: it will not overwrite a file that already exists.

docs <- list(
  list(
    file = "01_data_setup.qmd",
    title = "PK/PD Workshop - Part 1: Data Setup & Exploration",
    sections = c("Packages", "Load Data", "Explore & Visualize (EDA)")
  ),
  list(
    file = "02_nca.qmd",
    title = "PK/PD Workshop - Part 2: Non-Compartmental Analysis",
    sections = c("Packages", "NCA Summary")
  ),
  list(
    file = "03_modeling.qmd",
    title = "PK/PD Workshop - Part 3: Population PK Modeling & Diagnostics",
    sections = c("Packages", "Population PK Modeling", "Diagnostics")
  ),
  list(
    file = "04_simulation.qmd",
    title = "PK/PD Workshop - Part 4: Dose Simulation & Decision",
    sections = c("Packages", "Dose Simulation", "Dose Decision Reporting")
  )
)

for (doc in docs) {
  if (file.exists(doc$file)) {
    cat("[skip]   ", doc$file, "already exists -- not overwriting\n")
    next
  }

  header <- c(
    "---",
    paste0('title: "', doc$title, '"'),
    "format:",
    "  html:",
    "    toc: true",
    "---",
    "",
    "<!-- pk_data should already be in your R session from generate_pk_datasets.R -->",
    "<!-- Ask Posit Assistant to write the code for each section below, then run the chunk it inserts. -->",
    ""
  )

  body <- unlist(lapply(doc$sections, function(s) c(paste0("## ", s), "")))

  writeLines(c(header, body), doc$file)
  cat("[created]", doc$file, "\n")
}

cat("\nDone. Open each file in Positron and work through it with Posit Assistant.\n")
