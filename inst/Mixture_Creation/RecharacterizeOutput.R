#-----------------------------------------------------------------------------------------------------#
# READ IN COMMAND LINE ARGUMENTS                                                                      #
#-----------------------------------------------------------------------------------------------------#

args<-commandArgs(TRUE)

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  file = "g10h90"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#-------------------------------------------------------------------------#
# Setting the Input Directory                                             #
#-------------------------------------------------------------------------#

library(gtools)

setwd("/netscr/drwilson/2018-04-05 Paper 1/Programs/")
source("geneModel_multcell_edit.R")
source("loadData_edit.R")
source("isoDeconv_geneModel_revised.R")
source("pdist_gen.R")
source("EffectiveLength_Total.R")

setwd("/netscr/drwilson/2018-04-05 Paper 1/MappedData/Real_Data_Opt/GeneModels/")

#-------------------------------------------------------------------------#
# CALL the geneMod Functions                                              #
#-------------------------------------------------------------------------#
output = sprintf("sigred_%s_geneMod.RData",file)

load(output)

#-------------------------------------------------------------------#
# EDIT TO GROUP CELL TYPES                                          #
#-------------------------------------------------------------------#
info_mat = sig_geneMod[["Sample_Info"]]
cellTypes = unique(info_mat$Cell_Type)

ctList = list()

for(j in 1:length(cellTypes)){
  idx = which(info_mat$Cell_Type==cellTypes[j])
  ctList[[cellTypes[j]]] = list(samps = info_mat$Label[idx], tots = info_mat$Total[idx])
}

idx2consider = which(names(sig_geneMod)!="Sample_Info")
for(k in idx2consider){
  for(l in 1:length(cellTypes)){
    samps2use = ctList[[l]]$samps
    tots      = ctList[[l]]$tots

    y_vecs  = paste("sig_geneMod[[k]]$y",samps2use,sep = "_")
    y_vecsc = paste(y_vecs,collapse = ",")
    nExon = eval(parse(text=sprintf("length(%s)",y_vecs[1])))
    textcmd = sprintf("matrix(c(%s),nrow=nExon,ncol=length(samps2use))",y_vecsc)
    expMat  = eval(parse(text=textcmd))

    totmg   = tots-colSums(expMat)
    expMat2 = rbind(totmg,expMat)

    if(cellTypes[l]!="mix"){
      sig_geneMod[[k]][[cellTypes[l]]] = list(cellType=cellTypes[l],rds_exons=expMat2)
    } else {
      sig_geneMod[[k]][[cellTypes[l]]] = list(cellType=cellTypes[l],rds_exons_t=expMat2)
    }
    
  }
}

save(sig_geneMod,file=sprintf("rechar_%s_geneMod.RData",file))