###############################################
### Load or install packages ##################
###############################################
if (!require("RMariaDB")) install.packages("RMariaDB")
library(RMariaDB)

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
