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
# Place holders to establish what needs to be input into the functions:
# countData= c("mf_sbn_g60h40_counts.txt",
#              "refsamp_sbn_gm12878_counts.txt",
#              "refsamp_sbn_hmec_counts.txt")
# fragSizeFile = "full_g60h40_fraglens.txt"
# output = "unred_g60h40_geneMod.RData"

sys_statement1 = sprintf("countData=c(\"mf_%s_counts.txt\",
                         \"CD4_SRR1550995_SBC_counts.txt\",
                         \"CD4_SRR1551016_SBC_counts.txt\",
                         \"CD4_SRR1551098_SBC_counts.txt\",
                         \"CD8_SRR1551024_SBC_counts.txt\",
                         \"CD8_SRR1551051_SBC_counts.txt\",
                         \"CD8_SRR1551065_SBC_counts.txt\")",file)
eval(parse(text=sys_statement1))
labels = c("mix","CD4_ref1","CD4_ref2","CD4_ref3",
           "CD8_ref1","CD8_ref2","CD8_ref3")
cellTypes = c("mix","CD4","CD4","CD4",
              "CD8","CD8","CD8")
fragSizeFile = sprintf("/netscr/drwilson/2018-04-05 Paper 1/Frag_Lengths/%s_fraglens.txt",file)
bedFile = "/netscr/drwilson/Reference_Annotations/Homo_sapiens/Homo_sapiens.GRCh37.66.nonoverlap.exon.bed"
knownIsoforms = "/netscr/drwilson/Reference_Annotations/Homo_sapiens/Homo_sapiens.GRCh37.66.nonoverlap.exon.knownIsoforms.RData"
output = sprintf("unred_%s_geneMod.RData",file)
readLen=50
lmax=600
eLenMin=1

# variables: countData,labels,cellTypes,total_cts,bedFile,knownIsoforms,fragSizeFile,output,readLen,lmax,eLenMin

dev_compiled_geneMod(countData=countData,labels = labels,total_cts = total_cts, cellTypes=cellTypes, bedFile=bedFile,knownIsoforms=knownIsoforms,
                     fragSizeFile=fragSizeFile,output=output,readLen=readLen,lmax=lmax,eLenMin=eLenMin)