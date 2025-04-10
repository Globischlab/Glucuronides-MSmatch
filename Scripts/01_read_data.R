# ---- 01.read spectra data from mzml files  ----
# function to parse all the mzml files in data/raw/
read_mzml <- function(file_path){
    list_of_files <- list.files(path = file_path,
                            recursive = TRUE,
                            pattern = "\\.mzML$",
                            full.names = TRUE)
    sps <- Spectra(list_of_files)
    sps <- filterMsLevel(sps, msLevel = 2L)
    sps <- filterRt(sps, c(60,1140))
    return(sps)
}


NL <- 176.0321 # Neutral loss of glucuronid moeity

# ---- Load data ----

# load spectra from control samples
sps_ctrl <- read_mzml(here("data","raw","positive","Ctrl")) # positive mode
# filter mz range of glucuronide in the control samples
sps_ctrl <- filterPrecursorMzRange(sps_ctrl, c((NL+15),1200))

# load spectraenzymatic treatment mzml files
sps_ea <- read_mzml(here("data","raw","positive","EA"))

