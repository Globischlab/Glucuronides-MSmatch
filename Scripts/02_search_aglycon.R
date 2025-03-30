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