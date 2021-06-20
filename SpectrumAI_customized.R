# Original script was downloaded from:
# https://github.com/yafeng/SpectrumAI


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

## mzR will be installed automatically as MSnbase depends on it
#BiocManager::install(c("MSnbase", "protViz"))

setwd("E:/Lazaro/Segundo TMT 11/Segundo with isotope correction/20191128_MM_TMT11_IsotopeCorrection/Spectrum AI/SpectrumAI scripts and tables")  #set your working directory
mzml_path = "E:/Lazaro/Segundo TMT 11/TMT11 raw data" # set file path to which raw files are located
infile_name ="2_MMSegundo_SAAV_SpectrumAI_results.xlsx"  # PSM table file name

Frag.ions.tolerance= 0.02 # 0.02 Da tolerance for MS2 fragment ions mass accuracy.
relative=FALSE # set TRUE if ppm value is used for Frag.ions.tolerance

# or you can use ppm threshold
# Frag.ions.tolerance= 10 # 10 ppm tolerance for MS2 fragment ions mass accuracy.
# relative=TRUE

library(openxlsx)
library(readxl)
library(protViz)
library(MSnbase)
library(stringr)

DF <- psm.table[i,]
SpectraFile <- spectra_file

InspectSpectrum <- function (DF, SpectraFile){
  
  start.time <- Sys.time()
  ScanNum=as.integer(DF$ScanNum)
  peptide=as.character(DF$Peptide)
  sub_pos=as.integer(DF$sub_pos)
  seq=DF$Sequence
  
  exp_peaks <- peaks(SpectraFile, scan=ScanNum)
  predicted_peaks = predict_MS2_spectrum(Peptide =  as.character(DF$Peptide))
  match_ions = match_exp2predicted(exp_peaks, predicted_peaks, tolerance =Frag.ions.tolerance, relative = relative )
  
  maxintensity=max(exp_peaks[,2])
  average_intensity=mean(exp_peaks[,2])
  median_intensity=median(exp_peaks[,2])
  
  DF$sum.fragmentions.intensity=sum(exp_peaks[,2])
  DF$maxintensity=maxintensity
  DF$average_intensity=average_intensity
  DF$median_intensity=median_intensity
  
  DF$status="checked"
  
  if (nrow(match_ions)!=0){
    DF$matched_ions=paste(unique(match_ions$ion),collapse = ",")
    DF$sum.matchedions.intensity=sum(match_ions$intensity)
    
    supportions_intensity=0
    ions_support="NO"
    supportions=""
    
    for (j in 1:nrow(match_ions)){
      type=match_ions[j,]$type
      pos=match_ions[j,]$pos
      ion=match_ions[j,]$ion
      if (type=="b" & pos>=sub_pos){
        ions_support="YES"
        supportions_intensity=supportions_intensity+match_ions[j,]$intensity
        supportions=paste0(supportions,ion,",")
      }else if (type=="y" & pos>nchar(seq)-sub_pos){
        ions_support="YES"
        supportions_intensity=supportions_intensity+match_ions[j,]$intensity
        supportions=paste0(supportions,ion,",")
      }
      
      DF$ions_support=ions_support
      DF$support_ions=supportions
      DF$sum.supportions.intensity=supportions_intensity
      
      #check if it is a noise peak or isotope peak supporting mutant ions
      if (DF$sum.supportions.intensity < DF$median_intensity){DF$ions_support <- "NO"}
      
      flanking_ions_left=c()
      flanking_ions_right=c()
      flanking_ions=c()
      flanking_ions_support="NO"
      n1=DF$peptide_length
      n2=sub_pos
      if (n2 ==1){
        flanking_ions=c("b1",paste0("y",as.character(n1-1)))
        flanking_ions=intersect(flanking_ions,match_ions$ion)
        if (length(flanking_ions)>0){
          flanking_ions_support="YES"
        }
      }else if (n2 == n1){
        flanking_ions=c("y1",paste0("b",as.character(n1-1)))
        flanking_ions=intersect(flanking_ions,match_ions$ion)
        if (length(flanking_ions)>0){
          flanking_ions_support="YES"
        }
      }else {
        flanking_ions_left=c(paste0("b",as.character(n2-1)))
        flanking_ions_left=c(flanking_ions_left,paste0("y",as.character(n1-n2+1)))
        
        flanking_ions_right=c(paste0("b",as.character(n2)))
        flanking_ions_right=c(flanking_ions_right,paste0("y",as.character(n1-n2)))
        
        
        flanking_ions_left=intersect(flanking_ions_left,match_ions$ion)
        flanking_ions_right=intersect(flanking_ions_right,match_ions$ion)
        
        flanking_ions=union(flanking_ions_left,flanking_ions_right)
        if ( length(flanking_ions_left)>0 & length(flanking_ions_right)>0){
          flanking_ions_support="YES"
        }
      }
      
      DF$flanking_ions_support=flanking_ions_support
      DF$flanking_ions=paste(flanking_ions,collapse = ",")
      DF$sum.flanking.ions.intensity=sum(match_ions[match_ions$ion %in% flanking_ions,]$intensity)
      
      if (DF$sum.flanking.ions.intensity < DF$median_intensity){DF$flanking_ions_support <- "NO"}
      
      #fragmentation is not preferable at Cterm side of proline, so only require supporting ions
      if (grepl("P",substr(seq, sub_pos-1, sub_pos))){
        DF$flanking_ions_support=DF$ions_support
      }
    }
  }
  end.time <- Sys.time()
  DF$Time.taken <- paste(end.time - start.time)
  return(DF)
}

start.time <- Sys.time()

source('SpectrumAI_SpectraFunctions.R')
options(stringsAsFactors = FALSE)

df.psm=readxl::read_xlsx(infile_name, sheet=1)
#sapply(df.psm, class)
#Before running the next command, double check the header names in the input PSM table
#The df.psm dataframe should have at least the following columns with exactly same names (the order can be different): 
# "SpectraFile", "ScanNum", "Peptide",  "sub_pos" 

#df.psm <- df.psm[which(df.psm$Comment != "" & df.psm$Comment != "mzML file not found" & df.psm$Comment != "mzML file not found"),]

# replace file names wherever necessary
df.psm[grep("20180801_MM_TMT_B12_F6_R1", df.psm$SpectraFile),"SpectraFile"] <-"20180801_MM_TMT_B12_F6_R1_20180801115227" 
df.psm[grep("20180509_MM_TMT_B2_61min_R", df.psm$SpectraFile),"SpectraFile"] <-"20180509_MM_TMT_B2_61min_R2" 
df.psm[grep("20180719_MM_TMT_B9_F15_R1", df.psm$SpectraFile),"SpectraFile"] <-"20180719_MM_TMT_B8_F15_R1" 
df.psm[grep("20180509_MM_TMT_B2_F19_R1", df.psm$SpectraFile),"SpectraFile"] <-"20180509_MM_TMT_B1_F19_R1" 


df.psm.finished <- df.psm[which(df.psm$status =="checked" & df.psm$Comment !="mzML file not found"),]
df.psm.acetylated <- df.psm[which(df.psm$Comment =="acetylated peptide"),]
df.psm <- df.psm[which(df.psm$status !="checked" & df.psm$Comment !="acetylated peptide"),]

raw.files <- unique(df.psm$SpectraFile)
length(raw.files)
folders <- unique(unlist(lapply(raw.files, function(x) {strsplit(x, split = "_")[[1]][1]})))
length(folders)

df.psm$Comment <- "none"

df.psm.all.result <- matrix(nrow=0, ncol=ncol(df.psm))
all.no.peptides <- nrow(df.psm)
for (f in 1:length(folders)) {
  raw.files.inside <- raw.files[grepl(folders[f], raw.files)]
  mzml_path.folder = paste(mzml_path, folders[f], sep="/")
  df.psm.folder <- df.psm[(df.psm$SpectraFile %in% raw.files.inside),]
  no.peptides <- nrow(df.psm.folder)
  df.psm.folder.result <- matrix(nrow=0, ncol=ncol(df.psm.folder))
  colnames(df.psm.folder.result) <- colnames(df.psm.folder)
  outfile_name <- paste(c(folders[f],"SpectrumAI", "Results",no.peptides,"Peptides"), collapse="_")
  if (!dir.exists(mzml_path.folder)) {
    folderlist <- list.dirs(mzml_path)
    if (length(grep(folders[f], folderlist)) ==0) {
      df.psm.folder[,"Comment"] <- "Directory does not exist"
    } else {
      mzml_path.folder <- folderlist[grep(folders[f], folderlist)]
    }
  } 
  if (!dir.exists(mzml_path.folder)) {
    next
  } else {
    for (r in 1:length(raw.files.inside)) {
      raw.f <- raw.files.inside[r]
      psm.table <- df.psm.folder[df.psm.folder$SpectraFile == raw.f,]
      mzml_file=file.path(mzml_path.folder,raw.f)
      mzml_file <- paste(mzml_file, "mzML", sep=".")
      if (!file.exists(mzml_file)) {
        mzml_file.2 <- paste(c("E:/Lazaro/Segundo TMT 11/TMT11 raw data/20180509/", raw.f,".mzML"), collapse="")
        if ((!file.exists(mzml_file.2))) {
          psm.table$Comment <- "mzML file not found"
          df.psm.folder.result <- rbind(df.psm.folder.result,psm.table)
          next
        } else if (file.exists(mzml_file.2)) {
          spectra_file = openMSfile(mzml_file.2, verbose=T)
          
          for (i in 1:nrow(psm.table)) {
            if (is.na(psm.table[i,"sub_pos"])) {
              psm.table[i,"Comment"] <- "sub pos NA"
              psm.table[i,"sub_pos"] <- 0
              #next
            } else if (psm.table[i,"sub_pos"] ==0) {
              psm.table[i,"Comment"] <- "sub pos 0"
              #next
            } else if (grepl("Acetyl", psm.table[i,"Peptide"])) {
              psm.table[i,"Comment"] <- "acetylated peptide"
              next
            } else if (psm.table[i,"peptide_length"] < psm.table[i,"sub_pos"]) {
              psm.table[i,"Comment"] <- paste("sub pos too long", psm.table[i,"sub_pos"], sep="-")
              psm.table[i,"sub_pos"] <- 0
              #next
            }
            psm.table[i,] <- InspectSpectrum(psm.table[i,], spectra_file)
          }
        }
      } else {
        spectra_file = openMSfile(mzml_file, verbose=T)
        
        for (i in 1:nrow(psm.table)) {
          if (is.na(psm.table[i,"sub_pos"])) {
            psm.table[i,"Comment"] <- "sub pos NA"
            psm.table[i,"sub_pos"] <- 0
            #next
          } else if (psm.table[i,"sub_pos"] ==0) {
            psm.table[i,"Comment"] <- "sub pos 0"
            #next
          } else if (grepl("Acetyl", psm.table[i,"Peptide"])) {
            psm.table[i,"Comment"] <- "acetylated peptide"
            next
          } else if (psm.table[i,"peptide_length"] < psm.table[i,"sub_pos"]) {
            psm.table[i,"Comment"] <- paste("sub pos too long", psm.table[i,"sub_pos"], sep="-")
            psm.table[i,"sub_pos"] <- 0
            #
          } 
          psm.table[i,] <- InspectSpectrum(psm.table[i,], spectra_file)
        }
      }
      df.psm.folder.result <- rbind(df.psm.folder.result,psm.table)
    }
  }
  df.psm.all.result <- rbind(df.psm.all.result, df.psm.folder.result)
  #write.table(df.psm.folder.result,paste(outfile_name, "txt", sep="."),sep="\t",quote=F,row.names=F)
}

df.psm.all.result2 <- df.psm.all.result
#sapply(df.psm.all.result2, class)
for (i in 1:nrow(df.psm.all.result2)) {
  if (grepl("none", df.psm.all.result2[i,"Comment"])) { next}
  if (grepl("too long", df.psm.all.result2[i,"Comment"])) {
    no <- strsplit(as.character(df.psm.all.result2[i,"Comment"]), split = "long-")[[1]][2]
    df.psm.all.result2[i,"sub_pos"] <- as.numeric(as.character(no))
  }
  else if (grepl("sub pos NA", df.psm.all.result2[i,"Comment"])) {
    df.psm.all.result2[i,"sub_pos"] <- NA
  }
}

checked <- df.psm.all.result2[df.psm.all.result2$status =="checked",]
skipped <- df.psm.all.result2[df.psm.all.result2$status !="checked",]
not.detected <- checked[checked$ions_support !="YES",]


every.res <- rbind(df.psm.finished, df.psm.acetylated, df.psm.all.result2)

write.table(every.res, paste(resultfile_name, "txt", sep="."),sep="\t",quote=F,row.names=F)

wb <- loadWorkbook(infile_name, xlsxFile = NULL, isUnzipped = FALSE)
addWorksheet(wb, "Result2")

writeData(
  wb,
  sheet = "Result2",
  every.res)
saveWorkbook(wb, infile_name, overwrite = TRUE)

end.time <- Sys.time()
