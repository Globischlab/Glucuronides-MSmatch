# ---- 02. Search for aglycon  ----
# Function to search aglycon between control samples and ea samples
search_aglycon <- function(sps_ctrl,sps_ea,NL,polarity, adduct,parm) {
    # get precursor mz list from crt files
    spec_df_ea <- spectraData(sps_ea, c("rtime","precursorMz"))
    write.csv(spec_df_ea, here("output","sps_ea_table.csv"))
    # retrive crt mass spectra as dataframe
    spec_df <- spectraData(sps_ctrl, c("rtime","precursorMz"))
    spec_df$aglycon_theo <- (spec_df$precursorMz-NL)
    spec_df$exactmass <- if(polarity == 1) {
        round(spec_df$aglycon_theo - proton, 4)
    } else if(polarity == 0) {
        round(spec_df$aglycon_theo + proton, 4)
    } else {
        message("No polarity")
        NA  # Return NA or handle missing case
    }
    write.csv(spec_df,here("output", "sps_ctrl_table.csv"))

    # match the aglycons
    mtches <- matchValues(spec_df_ea, spec_df, parm,mzColname = "precursorMz")
    mtches_aglycon <- matchedData(mtches)
    mtches_aglycon <- mtches_aglycon[!is.na(mtches_aglycon$ppm_error),]
    mtches_aglycon <- mtches_aglycon[order(mtches_aglycon$ppm_error,decreasing = FALSE),]
    mtches_aglycon<- mtches_aglycon[!duplicated(mtches_aglycon$target_exactmass),]
    mtches_aglycon$ID <- seq.int(nrow(mtches_aglycon))
    write.csv(mtches_aglycon , here("output","aglycon_match_table.csv"))
    return(mtches_aglycon)
}

# function for ms1 match with database
MS1_match <- function(df_mtches, df_db,db_param,mz){
    MS1_match <- matchValues(df_mtches, 
                         df_db, 
                         db_param, mzColname= mz)
    MS1_matched <- matchedData(MS1_match)
    MS1_matched <- MS1_matched[!is.na(MS1_matched$ppm_error),]
    MS1_matched <- MS1_matched[order(MS1_matched$ppm_error,decreasing = FALSE),]
    MS1_matched <- MS1_matched[!duplicated(MS1_matched$ID),]
    MS1_matched$ID <- seq.int(nrow(MS1_matched))
    return(MS1_matched)
}


# set parameters and apply to data
proton <- 1.00728
NL <- 176.0321 # Neutral loss of glucuronid acid
adduct <- "[M+H]+" #positive mode
# Negative mode
#adduct <- "[M-H]-" 

# ---- map mass between theoretical aglycon in control samples to mass in EA samples ----
parm <- Mass2MzParam(adducts = adduct,
                        tolerance = 0.005, ppm =10)

# Compare the mz of estimated aglycons to mz of aglycons after EA
mtches_aglycon <- search_aglycon(sps_ctrl,sps_ea,NL ,polarity=1, adduct,parm)

# clean up output
df_mtches <- as.data.frame(mtches_aglycon)
df_mtches <- df_mtches %>%
  select(
    rtime, precursorMz,
    target_rtime, target_precursorMz,
    target_aglycon_theo, target_exactmass,
    adduct, ppm_error, ID
  ) %>%
  dplyr::rename(
    EA_rt = rtime,
    EA_mz = precursorMz,
    rt = target_rtime,
    mz = target_precursorMz,
    ppm_aglycon = ppm_error
  )

write.csv(df_mtches, here("output","aglycon_match_cleanup.csv"))

# ---- precursor mz match with databases ----
# Search mz in database -- MS1 match
parm2 <- Mass2MzParam(adducts = adduct,
                     tolerance = 0.005, ppm = 10)
                     
# with hmdb
df_hmdb <- read.csv(here("data","hmdb_cleanup_v02062023.csv"),header = TRUE, sep = ",")
# with pubchem
df_pubchem <- read.csv(here("data","PubChem_compound_text_glucuronide.csv"),header = TRUE, sep=",")



# EA_mz for alycon match in hmdb
df_ms1match_hmdb <- MS1_match(df_mtches,df_hmdb,parm2,"EA_mz")
df_ms1match_hmdb <- df_ms1match_hmdb[,c("EA_rt","EA_mz","rt","mz","target_aglycon_theo","adduct","ppm_aglycon","ID",
                                                   "target_accession", "target_exactmass","target_name","target_chemical_formula","ppm_error")]
write.csv(df_ms1match_hmdb, 
           here("output",  "mtches_aglycon_hmdb.csv"))


# mz for glucuronides matches in pubchem and hmdb
df_ms1match_pubchem <- MS1_match(df_mtches,df_pubchem ,parm2, "mz")
write.csv(df_ms1match_pubchem, 
           here("output",  "mtches_gluco_pubchem.csv"))

df_ms1match_gluco_hmdb <- MS1_match(df_mtches,df_hmdb ,parm2, "mz")

df_ms1match_gluco_hmdb <- df_ms1match_gluco_hmdb [,c("EA_rt","EA_mz","rt","mz","target_aglycon_theo","adduct","ppm_aglycon","ID",
                                                   "target_accession", "target_exactmass","target_name","target_chemical_formula","ppm_error")]

write.csv(df_ms1match_gluco_hmdb, 
           here("output",  "mtches_gluco_hmdb.csv"))