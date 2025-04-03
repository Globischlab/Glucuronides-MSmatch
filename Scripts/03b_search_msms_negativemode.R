######################## 3b. ms/ms screening negative mode only ###################################
##################################################################################################
####################################### Function to screen ms/ms #################################
##################################################################################################

#Fuction: Normalize intensities
norm_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}

#Fuction: Normalize MS/MS, 
low_int <- function(x, ...) {
    x > max(x, na.rm = TRUE) * (5 / 100 ) #remove the Int < 5%
}

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
sps_ctrl <- read_mzml(here("data","negative_ionization","ctrl"))
sps_ctrl <- filterPrecursorMz(sps_crt, c((NL+15),1200))

# negative mode only -- search for glucuronides FPs
# normalize sps
sps_crt_normalized <- addProcessing(sps_crt, norm_int)
sps_crt_normalized <- filterIntensity(sps_crt_normalized, intensity = low_int)

# search for glucuronid acid characteristic fragmentations
# 85.0295,113.0244,175.0248
has_gluc_fp_all <- containsMz(sps_crt_normalized,
    mz = c(85.0295,113.0244,175.0248), # ms/ms fragmentation patterns
    tolerance = 0.005, ppm = 5, 
    which = "all") # which = "all": all fragments should match; which = "any": one of the fragments should be matched 
# subset sps
sps_crt_sub_all <- sps_crt_normalized[has_gluc_fp_all]
spec_df_sub_all <- spectraData(sps_crt_sub_all, c("msLevel","rtime","dataOrigin","precursorMz","scanIndex"))

# clean up output
spec_df_sub_all$precursorMz <-round(spec_df_sub_all$precursorMz,3)
spec_df_sub_all<-spec_df_sub_all[!duplicated(spec_df_sub_all$precursorMz),]
spec_df_sub_all$ID <- seq.int(nrow(spec_df_sub_all))
write.csv(spec_df_sub_all, 
           here("output",crt_spec_hasFPall.csv))