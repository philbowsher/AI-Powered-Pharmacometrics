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
    "<!-- Ask Posit Assistant to write the code for each section below, then run the chunk it inserts. -->",
    "<!-- The chunk below makes this document self-contained: it loads pk_data even if you open -->",
    "<!-- this file fresh (e.g. to render it or convert it to a Dashboard), not just interactively. -->",
    "",
    "```{r}",
    "#| echo: false",
    "if (!exists(\"pk_data\")) {",
    "  if (file.exists(\"data/mab_pkpd.csv\")) {",
    "    pk_data <- readr::read_csv(\"data/mab_pkpd.csv\", show_col_types = FALSE)",
    "  } else if (file.exists(\"data/warfarin.csv\")) {",
    "    pk_data <- readr::read_csv(\"data/warfarin.csv\", show_col_types = FALSE)",
    "  } else {",
    "    stop(\"No dataset found in data/ -- run generate_pk_datasets.R first.\")",
    "  }",
    "}",
    "```",
    ""
  )

  body <- unlist(lapply(doc$sections, function(s) c(paste0("## ", s), "")))

  writeLines(c(header, body), doc$file)
  cat("[created]", doc$file, "\n")
}

cat("\nDone. Open each file in Positron and work through it with Posit Assistant.\n")
