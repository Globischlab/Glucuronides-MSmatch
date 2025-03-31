library(MetaboAnnotation)
#mz_ppm <- 5 # mz search tolerance
NL <- 176.0321
proton <- 1.00728
parm <- Mass2MzParam(adducts = "[M-H]-",
                     tolerance = 0.001, ppm = 5)

ctr_xcms_path <- "F:/Ioanna/Glucuronides/From Ioanna with LOVE/Feces_NEG/final_peaks_Feces_NEG_Control_filter.csv"
ea__xcms_path <-  "F:/Ioanna/Glucuronides/From Ioanna with LOVE/Feces_NEG/final_peaks_Feces_NEG_EA_filter.csv"
output_dic <-"F:/Ioanna/Glucuronides/From Ioanna with LOVE/Feces_NEG/Output"
ctr_df <- read.csv(ctr_xcms_path,header = TRUE, sep =",")
ea_df <- read.csv(ea__xcms_path,header = TRUE, sep =",")
ctr_df$aglycon_est <- (ctr_df$mz-NL)
ctr_df$exactmass <-round(ctr_df$aglycon_est+proton,4)
mtches <- matchValues(ea_df, ctr_df, parm,mzColname = "mz")
mtches_aglycon <- matchedData(mtches)[whichQuery(mtches),]
mtches_aglycon <- mtches_aglycon[!is.na(mtches_aglycon$score),]
mtches_aglycon <-mtches_aglycon[order(mtches_aglycon$ppm_error,decreasing = FALSE),]
mtches_aglycon$target_exactmass <- round(mtches_aglycon$target_exactmass,3)
mtches_aglycon<-mtches_aglycon[!duplicated(mtches_aglycon$target_exactmass),]
mtches_aglycon$ID <- seq.int(nrow(mtches_aglycon))
mtches_aglycon <- mtches_aglycon[,c("mz","rt","target_mz","target_rt","target_aglycon_est","target_exactmass","adduct","ppm_error")]#,"ID")]
write.csv(mtches_aglycon , paste(output_dic,"/aglycon_match_table_xcms2.csv", sep = ""))
df_hmdb = read.csv("F:/Fan/hmdb_cleanup_v02062023.csv",header = TRUE, sep = ",")

mtches_aglycon <- mtches_aglycon[,c("mz","rt","target_mz","target_rt")]#,"ID")]
MS1_match <- matchValues(mtches_aglycon, 
                         df_hmdb, 
                         parm, mzColname="mz")
MS1_matched <- matchedData(MS1_match)[whichQuery(MS1_match),]
MS1_matched <- MS1_matched[!is.na(MS1_matched$score),]
MS1_matched <- MS1_matched[order(MS1_matched$ppm_error,decreasing = FALSE),]
MS1_matched <- MS1_matched[!duplicated(MS1_matched$mz),]
MS1_matched$ID <- seq.int(nrow(MS1_matched))
MS1_matched <- MS1_matched[,c("mz","rt","target_mz","target_rt","ID","target_accession","target_exactmass","target_name","target_chemical_formula","ppm_error")]
write.csv(MS1_matched , paste(output_dic,"/aglycon_match_hmdb_filtered.csv", sep = ""))


MS1_match_gluc <- matchValues(mtches_aglycon, 
                         df_hmdb, 
                         parm, mzColname="target_mz")
MS1_matched_gluc <- matchedData(MS1_match_gluc)[whichQuery(MS1_match_gluc),]
MS1_matched_gluc <- MS1_matched_gluc[!is.na(MS1_matched_gluc$score),]
MS1_matched_gluc <- MS1_matched_gluc[order(MS1_matched_gluc$ppm_error,decreasing = FALSE),]
MS1_matched_gluc <-MS1_matched_gluc[!duplicated(MS1_matched_gluc$target_mz),]
MS1_matched_gluc <- MS1_matched_gluc[,c("mz","rt","target_mz","target_rt","target_aglycon_est","ID","target_accession","target_exactmass","target_name","target_chemical_formula","ppm_error")]
write.csv(MS1_matched_gluc , paste(output_dic,"/gluco_match_hmdb_filtered.csv", sep = ""))

df_pubchem = read.csv("F:/Fan/PubChem_compound_text_glucuronide.csv",header = TRUE, sep = ",")

MS1_match_gluc <- matchValues(mtches_aglycon, 
                              df_pubchem, 
                              parm, mzColname="target_mz")
MS1_matched_gluc <- matchedData(MS1_match_gluc)[whichQuery(MS1_match_gluc),]
MS1_matched_gluc <- MS1_matched_gluc[!is.na(MS1_matched_gluc$score),]
MS1_matched_gluc <- MS1_matched_gluc[order(MS1_matched_gluc$ppm_error,decreasing = FALSE),]
#MS1_matched_gluc <-MS1_matched_gluc[!duplicated(MS1_matched_gluc$target_mz),]
#MS1_matched_gluc <- MS1_matched_gluc[,c("mz","rt","target_mz","target_rt","target_aglycon_est","ID",,"target_exactmass","target_cmpdname","target_cid","ppm_error")]
write.csv(MS1_matched_gluc , paste(output_dic,"/gluco_match_pubchem_filtered.csv", sep = ""))


