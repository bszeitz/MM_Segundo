# Original script was downloaded from:
# https://github.com/yafeng/SpectrumAI
# Published in Zhu et al., 2018, PMID: 29500430

# Script accessed from: https://github.com/yafeng/SpectrumAI/blob/master/Spectra_functions.R

require(mzR)
require(protViz)
require(stringr)

predict_MS2_spectrum <- function (Peptide, product_ion_charge = 1){
  pepMSGF=gsub("[^A-Z]","",Peptide)
  size = nchar(pepMSGF)
  pepMSGFMods=Peptide
  pepMSGFMods=str_replace_all(pepMSGFMods,pattern = "([\\+,0-9,:,\\.]+)" ,replacement = "\\1:" )
  
  l='none'
  while (l != pepMSGFMods){
    l=pepMSGFMods
    pepMSGFMods=sub("([^\\+,\\-,0-9,:])([^\\+,\\-,0-9,:])", "\\1:\\2", pepMSGFMods, perl=TRUE)
  }
  pepMSGFMods=t(str_split(pepMSGFMods,pattern = ":",simplify = T))
  pepMSGFMods=data.frame(X1=str_extract(pepMSGFMods,"[A-Z]+"),X2=str_extract(pepMSGFMods,"[^A-Z]+"), stringsAsFactors = F)
  pepMSGFMods[,2][is.na(pepMSGFMods[,2])] <- 0
  pepMSGFMods[,2]=  as.double(pepMSGFMods[,2])
  if (is.na(pepMSGFMods[1,1])){ # add mass of N-term modification to the first amino acid if there is one
    pepMSGFMods[2,2] = pepMSGFMods[1,2] + pepMSGFMods[2,2]
  }
  
  pepMSGFMods=na.omit(pepMSGFMods) # remove the line of N-term Mass  
  pepMSGFWeights = protViz::aa2mass(pepMSGF)[[1]]
  pepMSGFWeights =  pepMSGFWeights + t(as.double(pepMSGFMods[,2]))
  
  ions=fragmentIon(pepMSGFWeights)
  #convert it data frame
  ions <- data.frame(mz=c(ions[[1]]$b,ions[[1]]$y),
                     ion=c(paste0(rep("b",size),1:size),paste0(rep("y",size),1:size)),
                     type=c(rep("b",size),rep("y",size)),
                     pos =rep(1:size,2),
                     z=rep(1,size*2))
  
  proton_mono_mass = 1.007276
  if (product_ion_charge>1){
    ions2 = ions
    ions2$mz = (ions$mz + proton_mono_mass)/2
    ions2$z = 2
    ions = rbind(ions,ions2)
  }
  return (ions)
}

match_exp2predicted <- function (exp_peak,pred_peak,tolerance = 0.02, relative=FALSE){
  pred_peak$error = apply(pred_peak,1,function(x) min(abs(as.numeric(x[1])-exp_peak[,1])))
  pred_peak$intensity = apply(pred_peak,1,function(x) exp_peak[which.min(abs(as.numeric(x[1]) - exp_peak[,1])),2])
  pred_peak$ppm = round(pred_peak$error/pred_peak$mz*1000000,2)
  
  if (relative){
    match_ions = pred_peak[pred_peak$ppm <tolerance,]
  }else{match_ions = pred_peak[pred_peak$error <tolerance,]
  }
  
  return (match_ions)
}