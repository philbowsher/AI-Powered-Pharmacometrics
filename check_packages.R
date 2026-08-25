# Run this FIRST, before anything else in the workshop.
# It checks whether the packages used later are installed, and installs
# any that are missing. Takes a few seconds if everything's already there;
# a minute or two if it needs to install anything.
#
# Environment this was built/tested against (to replicate in a new session):
#   R 4.4.0, x86_64-pc-linux-gnu (Ubuntu 22.04 "jammy")
#   CRAN repo:        https://packagemanager.posit.co/cran/__linux__/jammy/latest
#   Bioconductor repo: https://p3m.dev/bioconductor
#   (already the default here -- don't override with cloud.r-project.org)
#
# Python: NOT required. nlmixr2/rxode2 use Rcpp + the symengine R package
# (a C++ library binding), not Python -- there is no reticulate/Python
# dependency anywhere in this workshop's R toolchain, despite what the
# console messages ("loading into symengine environment...") might suggest.
# A system Python (3.12.3, at /opt/python/3.12.3) is present in this
# environment for unrelated purposes (e.g. Quarto's Python engine) and is
# not needed for anything in this workshop.

required_packages <- c(
  "tidyverse", "patchwork",        # core data wrangling + plotting
  "gt", "gtsummary",                # tables
  "PKNCA",                          # NCA
  "nlmixr2",                        # population PK modeling
  "vpc",                            # visual predictive check
  "mrgsolve",                       # dose simulation
  "pointblank",                     # data validation
  "nlmixr2data",                    # only needed if you pick the Warfarin dataset
  "bslib", "shiny", "shinydashboard", "plotly", "DT"  # dashboard/Shiny step
)

cat("Checking", length(required_packages), "packages...\n\n")

missing <- character(0)

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("[ok]     ", pkg, "\n")
  } else {
    cat("[missing]", pkg, "\n")
    missing <- c(missing, pkg)
  }
}

cat("\n")

if (length(missing) > 0) {
  cat(length(missing), "package(s) missing. Installing now...\n\n")
  install.packages(missing)
  cat("\nDone. Re-run this script to confirm everything now shows [ok].\n")
} else {
  cat("All packages available. You're ready for the workshop.\n")
}
