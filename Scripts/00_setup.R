# ---- 00. Install & Load libraries ----
# cran
cran_pkgs <- c("here", "dplyr", "tidyr", "pander", "devtools","ggplot2")

installed <- rownames(installed.packages())
to_install_cran <- setdiff(cran_pkgs, installed)
if (length(to_install_cran) > 0) install.packages(to_install_cran)

# load cran pkg
lapply(cran_pkgs, library, character.only = TRUE)

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("BiocParallel", "mzR", "CompoundDb", "MsBackendMgf", 
               "AnnotationHub", "Spectra","MsBackendMassbank","MetaboAnnotation")

for (pkg in bioc_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, force = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Load AnnotationHub and set global object
ah <- AnnotationHub()

# GitHub packages (if not already installed)
# if (!require("MetaboAnnotation")) {
#   devtools::install_github("rformassspectrometry/MetaboAnnotation")
#   library(MetaboAnnotation)
# }


# ---- Set working directory ----
setwd(here())
message("Working directory set to: ", getwd())

dir.create(here::here("data","raw"),recursive= TRUE,showWarning = FALSE)
dir.create(here::here("data","xcms"),recursive= TRUE,showWarning = FALSE)
dir.create(here::here("output","ms1match"),recursive= TRUE,showWarning = FALSE)
# Save session info
writeLines(capture.output(sessionInfo()), "session_info.txt")
