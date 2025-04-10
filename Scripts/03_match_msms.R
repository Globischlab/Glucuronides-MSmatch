# ---- 03. ms/ms screening   ----
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

maxTic <- function(x, ...){
    tic <- vapply(x,function(z) sum(z[,"intensity"], na.rm = TRUE),
    numeric(1))
    x[[which.max(tic)]]

}

    # set up the sub directory
make_dir <- function(db) {
    # get the name of the database
    sub_dir <- deparse((substitute(db)))
    
    # create the subfolder path
    dir_path <- here::here("output",sub_dir)

    # check if directory exists, if not creat one
    if (!dir.exists(dir_path)){
        dir.create(dir_path, recursive = TRUE)
        message("Create directory: ", dir_path)
    } else {
       message("Directory already exists", dir_path)
    }

    # Rturn the path
    invisible(dir_path)
}

#Fuction: Match ms with selected m/z 
msms_match <- function(sps,df_match,db,param, polarity){
     # Normalize sps for control samples
    sub_dir <- deparse(substitute(db))
    sps_normalized <- addProcessing(sps, norm_int)
    sps_normalized <- filterIntensity(sps_normalized, intensity = low_int)

    # Normalize sps for libraries
    db <- filterPolarity(db, polarity = polarity)
    db <- addProcessing(db,norm_int)

    # filter sps with mz
    for(i in seq_len(nrow(df_match))){
        mz <- df_match$mz[i]
        id <- df_match$ID[i]
        rtime <- df_match$rt[i]
        sps_ms <- filterValues(sps_normalized,
                        spectraVariables = c("precursorMz", "rtime"),
                        values = c(mz,rt),
                        tolerance = c(0.005, 20),
                        ppm = c(10,0), match = "all")
        # screen ms/ms
        sps_agg <- combineSpectra(sps_ms, FUN = maxTic, minProp =5)
        if (length(sps_ms)>1){
            mtch <- matchSpectra(sps_ms, db, param)
            mtch_sub <- mtch[whichQuery(mtch)]
            df_mtch_sub <- apply(spectraData(mtch_sub),2,as.character)
            if(length(df_mtch_sub) == 0){
                message("No hit with mz = ", mz)
            }
            else{write.csv(df_mtch_sub, 
                here("output", sub_dir, paste0("ms2mtch_gluc_", id, ".csv")))
                }
        }
    }
}


# parameters for ms/ms matching
parm_ms2 <- MatchForwardReverseParam(ppm = 10, requirePrecursor = TRUE,
                           THRESHFUN = function(x) which(x >= 0.6)
                          #THRESHFUN = select_top_match
                           )

# load library
query(ah, "MassBank")
mbank <- ah[["AH119519"]] |>
  Spectra()

# create a subfolder
path_subfolder <- make_dir(mbank)

# search ms/ms
## ms1 matched in hmdb
df_ms1match_gluco_hmdb <- read.csv(here("output","mtches_gluco_hmdb.csv"))

msms_match_gluco <- msms_match(sps_ctrl,df_ms1match_gluco_hmdb,mbank,parm_ms2,polarity = 1)

# ms1 matched in pubchem
df_ms1match_gluco_pubchem <- read.csv((here("output","mtches_gluco_pubchem.csv")))

msms_match_gluco_pubchem <- msms_match(sps_ctrl,df_ms1match_pubchem,mbank,parm_ms2,polarity = 1)
