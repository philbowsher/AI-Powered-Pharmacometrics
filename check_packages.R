# Run this FIRST, before anything else in the workshop.
# It checks whether the packages used later are installed, and installs
# any that are missing. Takes a few seconds if everything's already there;
# a minute or two if it needs to install anything.
#
# Environment this was built/tested against (to replicate in a new session):
#   R 4.4.0, x86_64-pc-linux-gnu (Ubuntu 22.04 "jammy")
#   CRAN repo:        https://packagemanager.posit.co/cran/__linux__/jammy/<dated snapshot, see below>
#   Bioconductor repo: https://p3m.dev/bioconductor
#   (pinned below -- don't override with cloud.r-project.org)
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

# Pin the package repo to a DATED snapshot, not "/latest". /latest is a
# moving target: this environment was validated against R 4.4.0, but a
# live workshop run hit a case where upstream CRAN had drifted a
# dependency (Deriv) past compatibility with that R build, with no
# prebuilt binary available to fall back on. A dated snapshot makes this
# reproducible over time. Update the date below if you re-validate later
# -- this is a deliberate pinned override, not a switch to a public
# CRAN mirror (that's still forbidden).
options(repos = c(
  CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/2026-08-25",
  Bioc = "https://p3m.dev/bioconductor"
))

# Known gotcha (hit during a live workshop run): CRAN's current Deriv
# (>= 4.3.0, a dependency of nlmixr2) ships C++ code calling R-internals
# symbols (R_ClosureFormals, Rf_allocLang) not declared in this R 4.4.0
# installation's headers, so it fails to compile from source -- and no
# prebuilt binary exists for this version on this R/platform combination
# via the package manager (confirmed: available.packages(type="binary")
# still resolves to src/contrib even with a correct platform user-agent,
# so this is not a binary-flag/HTTPUserAgent issue -- that was checked and
# ruled out). Fix: install the last pure-R release (4.1.2, predates the
# C++ rewrite) from the CRAN archive before nlmixr2 pulls in the broken one.
if (!requireNamespace("Deriv", quietly = TRUE)) {
  tryCatch({
    download.file(
      "https://cran.r-project.org/src/contrib/Archive/Deriv/Deriv_4.1.2.tar.gz",
      destfile = file.path(tempdir(), "Deriv_4.1.2.tar.gz"), quiet = TRUE
    )
    install.packages(file.path(tempdir(), "Deriv_4.1.2.tar.gz"), repos = NULL, type = "source")
    # dependencies = FALSE avoids re-triggering the same failure by pulling
    # a fresh Deriv back in while resolving these two -- tested live, needed.
    install.packages(c("xgxr", "nlmixr2plot"), dependencies = FALSE)
    cat("[fixed]   Pre-installed Deriv 4.1.2 + xgxr/nlmixr2plot (avoids a known nlmixr2 install failure)\n\n")
  }, error = function(e) {
    message("Could not pre-install Deriv 4.1.2 -- if nlmixr2 fails to install with ",
            "R_ClosureFormals/Rf_allocLang compile errors, see the note above this line.")
  })
}

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
  cat(length(missing), "package(s) missing. Installing one at a time so failures", 
      "are diagnosed individually instead of as one wall of compiler noise...\n\n")

  failed <- character(0)
  for (pkg in missing) {
    result <- tryCatch({
      install.packages(pkg)
      "ok"
    }, error = function(e) e, warning = function(w) w)

    if (identical(result, "ok") && requireNamespace(pkg, quietly = TRUE)) {
      cat("[installed]", pkg, "\n")
    } else {
      failed <- c(failed, pkg)
      cat("[FAILED]   ", pkg, "\n")
      cat("           -> If this is a compile error (undefined symbols, missing headers),\n")
      cat("              it's likely a source/R-version mismatch with no binary available.\n")
      cat("              Try an older version from the CRAN archive, e.g.:\n")
      cat("              https://cran.r-project.org/src/contrib/Archive/", pkg, "/\n", sep = "")
      cat("              (see check_packages.R's Deriv fix above for a worked example)\n")
    }
  }

  cat("\n")
  if (length(failed) > 0) {
    cat("Still missing:", paste(failed, collapse = ", "), "-- see hints above.\n")
  } else {
    cat("All packages now installed.\n")
  }
} else {
  cat("All packages available. You're ready for the workshop.\n")
}
