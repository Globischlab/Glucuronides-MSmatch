###############################################
### Load or install packages ##################
###############################################

if (!require("here")) install.packages("here")
library(here)

if (!require("pander")) install.packages("pander")
library(pander)

if (!require("devtools")) install.packages("devtools")

if (!require("Spectra")) devtools::install_github("rformassspectrometry/Spectra")
library(Spectra)

if (!require("MetaboAnnotation")) devtools::install_github("rformassspectrometry/MetaboAnnotation")
library(MetaboAnnotation)

if (!require("MsBackendMassbank")) devtools::install_github("rformassspectrometry/MsBackendMassbank")
library(MsBackendMassbank)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("mzR")) BiocManager::install("mzR")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("AnnotationHub")) BiocManager::install("AnnotationHub")
library(AnnotationHub)
ah <- AnnotationHub()
library(MsBackendMgf)


# Set working directory to repo root
setwd(here())
# Verify
print(paste("Working directory set to:", getwd()))
if (!dir.exists(here("output"))) {
  dir.create(here("output"))
}
# Save session info
writeLines(capture.output(sessionInfo()), "session_info.txt")
