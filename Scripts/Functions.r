# Parse the directories and Load mzml files
load_mzml <- function(file_path){
    list_of_files <- list.files(path = file_path,
                            recursive = TRUE,
                            #pattern = "\\.mzML$",
                            full.names = TRUE)
    sps <- Spectra(list_of_files)
    sps <- filterRt(sps, c(60,1140))
    return(sps)
}


#Normalize MS/MS, remove the Int < int_tresh
low_int <- function(x, ...) {
    x > max(x, na.rm = TRUE) * (int_tresh / 100 )
}


# normalize intensities
norm_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}


# Compare the mz of estimated aglycons to mz of aglycons after EA
search_aglycon <- function(sps_crt,sps_ea,NL,polarity, adduct,parm) {
    # get precursor mz list from crt file
    sps_crt <- filterPrecursorMz(sps_crt, c((NL+15),1200))
    # retrive crt mass spectra as dataframe
    spec_df <- spectraData(sps_crt, c("msLevel","rtime","dataOrigin","precursorMz","scanIndex"))
    spec_df$aglycon_theo <- (spec_df$precursorMz-NL)
    if(polarity == 1){# positive mode
        spec_df$exactmass <-round(spec_df$aglycon_theo-proton,4)
    } 
    if(polarity == 0){# negative mode
        spec_df$exactmass <-round(spec_df$aglycon_theo+proton,4)
    }
    else{message("No polarity")}
    write.csv(spec_df,paste(output_dic,"/crt_spec_table.csv", sep = ""))
    # get precursor mz list from crt files
    spec_df_ea <- spectraData(sps_ea, c("msLevel","rtime","dataOrigin","precursorMz","scanIndex"))
    write.csv(spec_df_ea,paste(output_dic,"/ea_spec_table.csv", sep = ""))
    # match the aglycons
    
    mtches <- matchValues(spec_df_ea, spec_df, parm,mzColname = "precursorMz")
    print(mtches)
    mtches_aglycon <- matchedData(mtches)[whichQuery(mtches),]
    mtches_aglycon <- mtches_aglycon[!is.na(mtches_aglycon$score),]
    mtches_aglycon <-mtches_aglycon[order(mtches_aglycon$ppm_error,decreasing = FALSE),]
    mtches_aglycon<-mtches_aglycon[!duplicated(mtches_aglycon$target_exactmass),]
    mtches_aglycon$ID <- seq.int(nrow(mtches_aglycon))
    print(mtches_aglycon)
    write.csv(mtches_aglycon , paste(output_dic,"/aglycon_match_table.csv", sep = ""))
    return(mtches_aglycon)
}


#MS2 match with experimental MS2 databases
msms_match <- function(sps,mtches_aglycon,db,param){
    # setting up the sub directory

sub_dir <- paste0((substitute(db)))

# check if sub directory exists 
if (file.exists(sub_dir)){    
} 
else {		
		dir.create(file.path(output_dic, sub_dir))		
}        
    # Normalize sps
    sps_normalized <- addProcessing(sps, norm_int)
    sps_normalized <- filterIntensity(sps_normalized, intensity = low_int)
    db_normalized <- addProcessing(db,norm_int)
    for(i in seq_len(nrow(mtches_aglycon))){
        mz <- mtches_aglycon$precursorMz[i]
        id <- mtches_aglycon$ID[i]
        sps_ms <- filterPrecursorMzValues(sps_normalized,mz,ppm = 10, tolerance = 0.005)
        if (length(sps_ms)>1){
            mtch <- matchSpectra(sps_ms, db_normalized, param)
            mtch_sub <- mtch[whichQuery(mtch)]
            df_mtch_sub <- apply(spectraData(mtch_sub),2,as.character)
            if(length(df_mtch_sub) == 0){
                message("No hit with mz = ", mz)
            }
            else{write.csv2(df_mtch_sub,paste(file.path(output_dic, sub_dir,paste0("ms2mtch_",id, ".csv"))))}
        }
        
    }
}


