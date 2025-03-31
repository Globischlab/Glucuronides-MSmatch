#################### 2. Search for aglycon ####################

##################################################################################################
############# Function to search aglycon between control samples and ea samples ##################
##################################################################################################

# Function to search aglycon between control samples and ea samples
search_aglycon <- function(sps_ctrl,sps_ea,NL,polarity, adduct,parm) {
    # get precursor mz list from crt files
    spec_df_ea <- spectraData(sps_ea, c("msLevel","rtime","dataOrigin","precursorMz","scanIndex"))
    write.csv(spec_df_ea, here("output","sps_ea_table.csv"))
    # get precursor mz list from crt file
    sps_ctrl <- filterPrecursorMz(sps_ctrl, c((NL+15),1200))
    # retrive crt mass spectra as dataframe
    spec_df <- spectraData(sps_ctrl, c("msLevel","rtime","dataOrigin","precursorMz","scanIndex"))
    spec_df$aglycon_theo <- (spec_df$precursorMz-NL)
    if(polarity == 1){# positive mode
        spec_df$exactmass <-round(spec_df$aglycon_theo-proton,4)

    } 
    if(polarity == 0){# negative mode
        spec_df$exactmass <-round(spec_df$aglycon_theo+proton,4)
    }
    else{message("No polarity")}
    
    write.csv(spec_df,here("output", "sps_ctrl_table.csv"))

    # match the aglycons
    mtches <- matchValues(spec_df_ea, spec_df, parm,mzColname = "precursorMz")
    mtches_aglycon <- matchedData(mtches)[whichQuery(mtches),]
    mtches_aglycon <- mtches_aglycon[!is.na(mtches_aglycon$score),]
    mtches_aglycon <- mtches_aglycon[order(mtches_aglycon$ppm_error,decreasing = FALSE),]
    mtches_aglycon<- mtches_aglycon[!duplicated(mtches_aglycon$target_exactmass),]
    mtches_aglycon$ID <- seq.int(nrow(mtches_aglycon))
    write.csv(mtches_aglycon , here("output","aglycon_match_table.csv"))
    return(mtches_aglycon)
}

# set parameters and apply to data
proton <- 1.00728
NL <- 176.0321 # Neutral loss of glucuronid acid
adduct <- "[M+H]+"
# Negative mode
#adduct <- "[M-H]-" 
parm <- Mass2MzParam(adducts = adduct,
                        tolerance = 0.005, ppm = 5)

# Compare the mz of estimated aglycons to mz of aglycons after EA
mtches_aglycon <- search_aglycon(sps_ctrl,sps_ea,NL ,polarity=1, adduct,parm)

# Search mz in database -- MS1 match
# with hmdb
df_hmdb <- read.csv(here("data","databases","hmdb_cleanup_v02062023.csv"),header = TRUE, sep = ",")
# with pubchem
df_pubchem <- read.csv(here("data","databases","pubchem_glucuronides.csv"),header = TRUE，sep=",")

parm2 <- Mass2MzParam(adducts = adduct,
                     tolerance = 0.001, ppm = 5)

MS1_match <- function(df_match, df_db,parm2,mz_name){
    df_match <- df_match[,c("mz","rt","target_mz","target_rt")]#,"ID")]
    MS1_match <- matchValues(df_match, 
                         df_db, 
                         parm2, mzColname= mz_name)
    MS1_matched <- matchedData(MS1_match)[whichQuery(MS1_match),]
    MS1_matched <- MS1_matched[!is.na(MS1_matched$score),]
    MS1_matched <- MS1_matched[order(MS1_matched$ppm_error,decreasing = FALSE),]
    MS1_matched <- MS1_matched[!duplicated(MS1_matched$mz),]
    MS1_matched$ID <- seq.int(nrow(MS1_matched))
    MS1_matched <- MS1_matched[,c("mz","rt","target_mz","target_rt","ID","target_accession","target_exactmass","target_name","target_chemical_formula","ppm_error")]
    write.csv2(MS1_matched, 
           here(output, paste0("ms1mtch_", substitute(df_db),"_",substitute(df_match), ".csv")))
           
    return(MS1_matched)

}

df_ms1match_hmdb <- MS1_match(mtches_aglycon,df_hmdb,parm2, mz_name = "mz")
df_ms1match_pubchem <- MS1_match(mtches_aglycon,df_pubchem ,parm2, mz_name = "target_mz")
df_ms1match_gluco_hmdb <- MS1_match(mtches_aglycon,df_hmdb ,parm2, mz_name = "target_mz")
