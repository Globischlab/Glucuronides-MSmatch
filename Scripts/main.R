###############################################
############### Path setting ##################
###############################################

# control mzml files directory
crt_path <-
# EA mzml files directory 
ea_path <- 
# save the outputs
output_dic <-



###############################################
######### Parameter settings ##################
###############################################
int_tresh <- 5 # remove MS2 peaks Int < 5%
NL <- 176.0321 # Neutral loss of glucuronid acid
proton <- 1.00728
adduct <- "[M-H]-"


###############################################
######### Parameter settings for MS1 ##########
###############################################
mz_ppm <- 5 # mz search tolerance
parm <- Mass2MzParam(adducts = adduct,
                        tolerance = 0.005, ppm = mz_ppm)


###############################################
######### Parameter settings for MS2 ##########
###############################################
dp_tresh <- 0.5 # threashold for MS2 matching scores
mz_ppm <- 5 # mz search tolerance
parm_ms2 <- MatchForwardReverseParam(ppm = mz_ppm, requirePrecursor = FALSE,
                           THRESHFUN = function(x) which(x >= dp_tresh)
                          #THRESHFUN = select_top_match
                           )

###############################################
############# Call functions ##################
###############################################

# load mzml files
sps_crt <- load_mzml(crt_path)
sps_crt <- filterPrecursorMz(sps_crt, c((NL+15),1200)) # filter precursor mz range
sps_ea <- load_mzml(ea_path)

# search aglycon
mtches_aglycon <- search_aglycon(sps_crt,sps_ea,NL ,polarity=0, adduct,parm)

# negative mode only -- search for glucuronides FP

sps_crt_normalized <- addProcessing(sps_crt, norm_int)
sps_crt_normalized <- filterIntensity(sps_crt_normalized, intensity = low_int)
has_gluc_fp_all <- containsMz(sps_crt_normalized,mz = c(85.0295,113.0244,175.0248), tolerance = 0.005, ppm = 5, which = "all") # which = "any"
sps_crt_sub_all <- sps_crt_normalized[has_gluc_fp_all]
spec_df_sub_all <- spectraData(sps_crt_sub_all, c("msLevel","rtime","dataOrigin","precursorMz","scanIndex"))
#spec_df_sub_all$precursorMz <-round(spec_df_sub_all$precursorMz,4)
spec_df_sub_all$precursorMz <-round(spec_df_sub_all$precursorMz,3)
spec_df_sub_all<-spec_df_sub_all[!duplicated(spec_df_sub_all$precursorMz),]
spec_df_sub_all$ID <- seq.int(nrow(spec_df_sub_all))
write.csv(spec_df_sub_all ,paste(output_dic,"/crt_spec_hasFPall_table.csv", sep = ""))

# load library
query(ah, "MassBank")
mbank <- ah[["AH116166"]] |>
  Spectra()
mbank <- addProcessing(mbank, norm_int)
mbank_neg <-filterPolarity(mbank,polarity = 0)
mbank_neg

gnps <- Spectra("GNPS-LIBRARY.mgf",source = MsBackendMgf())
gnps <- addProcessing(gnps, norm_int)



# negative mode only -- search for glucuronides MS2

msms_match_glucuronides <- msms_match(sps_crt_sub_all,spec_df_sub_all,mbank_neg,parm_ms2)

msms_match_glucuronides <- msms_match(sps_crt_sub_all,spec_df_sub_all,gnps,parm_ms2)

# search for aglycon

# clean data
mtches_aglycon$precursorMz <-round(mtches_aglycon$precursorMz,3)
mtches_aglycon<-mtches_aglycon[!duplicated(mtches_aglycon$precursorMz),]
mtches_aglycon

msms_match_aglycon <- msms_match(sps_ea,mtches_aglycon,mbank_neg,parm_ms2)
msms_match_aglycon <- msms_match(sps_ea,mtches_aglycon,gnps,parm_ms2)
