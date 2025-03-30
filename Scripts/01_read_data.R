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
sps_ctrl <- read_mzml(here("data","ctrl"))
sps_ctrl <- filterPrecursorMz(sps_crt, c((NL+15),1200))

# ea mzml files
sps_ea <- read_mzml(here("data","ea"))

