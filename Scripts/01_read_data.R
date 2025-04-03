#################### Read mzml data fuction #################### 
read_mzml <- function(file_path){
    list_of_files <- list.files(path = file_path,
                            recursive = TRUE,
                            pattern = "\\.mzML$",
                            full.names = TRUE)
    sps <- Spectra(list_of_files)
    sps <- filterRt(sps, c(60,1140))
    return(sps)
}

# Read control mzml files & filter precursor mz range
NL <- 176.0321 # Neutral loss of glucuronid acid
proton <- 1.00728
# positive mode
# Read control samples
sps_ctrl <- read_mzml(here("data","positive_ion_mode","ctrl"))
# filter mz range of glucuronide in the control samples
sps_ctrl <- filterPrecursorMz(sps_crt, c((NL+15),1200))

# read enzymatic treatment mzml files
sps_ea <- read_mzml(here("data","positive_ion_mode","ea"))

