args<-commandArgs(TRUE)

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  inputDir = "/netscr/drwilson/2018-04-05 Paper 1/MappedData/CD4"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#--------------------------------------------------------------------------------------
# ISOFORM Software Library Access
#--------------------------------------------------------------------------------------
library(isoform, lib.loc="/nas02/home/d/r/drwilson/R/Rlibraries/") 

#------------------------------------------------------------------------------------------------
# Set bedFile (once identified and downloaded)               
#------------------------------------------------------------------------------------------------
bedFile = "/netscr/drwilson/Reference_Annotations/Homo_sapiens/Homo_sapiens.GRCh37.66.nonoverlap.exon.bed"

#------------------------------------------------------------------------------------------------
# Set Working Directory
#------------------------------------------------------------------------------------------------
setwd(inputDir)

#------------------------------------------------------------------------------------------------
# Loop across all replicates and files:
#        For loop assumes that all files within a class (ECC1, HMEC, GM12878) have same bed File.
#------------------------------------------------------------------------------------------------

cmd  = "ls *_sorted_by_name_uniq_filtered.bam"
ffs  = system(cmd, intern=TRUE)
length(ffs)
head(ffs)
sams = gsub("_sorted_by_name_uniq_filtered.bam", "", ffs)

for(i in 1:length(ffs)){

  sam1 = sams[i]
  cat(i, sam1, date(), "\n")
  
  bamFile = ffs[i]
  outFile = sprintf("%s_counts.txt", sam1)
  
  countReads(bamFile, bedFile, outFile)
}