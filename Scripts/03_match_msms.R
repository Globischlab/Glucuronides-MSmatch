#################### 3. ms/ms screening ####################
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

#Fuction: Match ms with selected m/z 
msms_match <- function(sps,df_match,db,param){
    # set up the sub directory

sub_dir <- paste0((substitute(db)))

# check if sub directory exists 
if (!dir.exists(here("output",sub_dir))){    
} 
else {		
		dir.create(here("output", sub_dir))		
}        
    # Normalize sps
    sps_normalized <- addProcessing(sps, norm_int)
    sps_normalized <- filterIntensity(sps_normalized, intensity = low_int)
    db_normalized <- addProcessing(db,norm_int)
    # filter sps with mz
    for(i in seq_len(nrow(df_match))){
        mz <- df_match$precursorMz[i]
        id <- df_match$ID[i]
        rt <- df_match$rtime[i]
        sps_ms <- filterPrecursorMzValues(sps_normalized,mz,ppm = 10, tolerance = 0.005)
        # screen ms/ms
        if (length(sps_ms)>2){
            mtch <- matchSpectra(sps_ms, db_normalized, param)
            mtch_sub <- mtch[whichQuery(mtch)]
            df_mtch_sub <- apply(spectraData(mtch_sub),2,as.character)
            if(length(df_mtch_sub) == 0){
                message("No hit with mz = ", mz)
            }
            else{write.csv2(df_mtch_sub, 
           here(output, sub_dir, paste0("ms2mtch_", id, ".csv")))}
        }
        
    }
}

parm_ms2 <- MatchForwardReverseParam(ppm = 5, requirePrecursor = FALSE #TRUE,
                           THRESHFUN = function(x) which(x >= 0.7)
                          #THRESHFUN = select_top_match
                           )

# load library
query(ah, "MassBank")
mbank <- ah[["AH116166"]] |>
  Spectra()
mbank <- addProcessing(mbank, norm_int)
mbank_sub <-filterPolarity(mbank,polarity = 1)

# clean data (optional)
mtches_aglycon$precursorMz <-round(mtches_aglycon$precursorMz,3)
mtches_aglycon<-mtches_aglycon[!duplicated(mtches_aglycon$precursorMz),]

# search ms/ms
msms_match_aglycon <- msms_match(sps_ea,mtches_aglycon,mbank_sub,parm_ms2)