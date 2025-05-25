# ---- 0. Point to your R 4.x library first (optional) ----
# .libPaths("/home/user/R/x86_64-pc-linux-gnu-library/4.5")

# ---- 1. CRAN packages to install ----
cran_pkgs <- c(
  "dplyr", "R.matlab", "ggplot2", "stringr", "purrr", "jsonlite",
  "tidyr", "patchwork", "scales", "fdatest", "xml2", "yaml",
  "rlang", "parallel", "foreach", "doParallel", "glue", "readr"
)

# ---- 2. Helper: install if missing ----
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg,
                     repos        = "https://cloud.r-project.org",
                     dependencies = TRUE)
  }
}

# ---- 3. Reinstall pkgload & usethis under R4.x ----
deps <- c("pkgload", "usethis")
for (p in deps) {
  if (p %in% rownames(installed.packages())) {
    message("Removing old ", p, " …")
    remove.packages(p)
  }
  message("Installing fresh ", p, " …")
  install.packages(p,
                   repos        = "https://cloud.r-project.org",
                   dependencies = TRUE)
}

# ---- 4. Install all other CRAN packages ----
invisible(lapply(cran_pkgs, install_if_missing))

# ---- 5. Install (or upgrade) devtools ----
if (!requireNamespace("devtools", quietly = TRUE)) {
  devtools::install.packages("devtools",
                             repos        = "https://cloud.r-project.org",
                             dependencies = TRUE,
                             force = TRUE)
} else {
  message("Updating devtools …")
  # Or the development version from GitHub:
  # install.packages("devtools")
  devtools::install_github("r-lib/devtools")
}

# ---- 6. Install epimod from GitHub if needed ----
if (!requireNamespace("epimod", quietly = TRUE)) {
  # remove.packages("epimod")
  devtools::install_github("https://github.com/qBioTurin/epimod", ref="epimod_pFBA")
}

# ---- 7. Load all libraries ----
all_pkgs <- c(cran_pkgs, deps, "devtools", "epimod")
invisible(lapply(all_pkgs, library, character.only = TRUE))

# ---- 5. Download epimod containers ----
downloadContainers()
