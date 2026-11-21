### ============================================================
### Step 0: install the R packages that every script needs
###
### Installs only what is missing (the packages below and everything they
### depend on) from CRAN into your personal R library. The personal library
### is always writable, so this also works on a lab computer where the main
### R folder is read-only. On Windows and macOS it installs ready-made
### (binary) packages, so nothing is compiled and no questions are asked.
### tidyverse is only there because the authors' script loads it.
###
### Run from the repo root:  Rscript code/novel/00_install_packages.R
### ============================================================

needed <- c("readxl", "dplyr", "tidyr", "metafor", "clubSandwich", "ggplot2", "ranger",
            "tidyverse")

user_lib <- Sys.getenv("R_LIBS_USER")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

repos <- getOption("repos")
if (is.null(repos) || !"CRAN" %in% names(repos) || repos[["CRAN"]] == "@CRAN@")
  options(repos = c(CRAN = "https://cloud.r-project.org"))
type <- if (.Platform$pkgType == "source") "source" else "binary"

### Needed packages that can't be loaded yet
not_loadable <- function() {
  needed[!vapply(needed, function(p) suppressMessages(requireNamespace(p, quietly = TRUE)), logical(1))]
}
### The needed packages plus everything they depend on, minus what is installed
still_missing <- function(db) {
  deps <- unlist(tools::package_dependencies(needed, db = db, recursive = TRUE), use.names = FALSE)
  setdiff(unique(c(needed, deps)), rownames(installed.packages()))
}

for (attempt in 1:2) {
  if (length(not_loadable()) == 0) break       # everything already works: nothing to download
  if (attempt == 2) {
    ### Some Windows networks block R's default download method (an SSL
    ### error); try once more with Windows' own method
    if (.Platform$OS.type != "windows") break
    options(download.file.method = "wininet")
  }
  db <- suppressWarnings(available.packages(type = type))
  if (nrow(db) == 0) next                      # couldn't reach CRAN with this method
  todo <- still_missing(db)
  if (length(todo) == 0) break
  cat("Installing:", paste(todo, collapse = ", "), "\n")
  install.packages(todo, lib = user_lib, type = type, dependencies = FALSE)
}

missing <- not_loadable()
if (length(missing) > 0)
  stop("Could not install: ", paste(missing, collapse = ", "),
       ". Check the internet connection, then run this file again.", call. = FALSE)
cat("All packages ready:", paste(needed, collapse = ", "), "\n")
