#--------------------------------------------------------#
# RUN on a single Cluster                                #
#--------------------------------------------------------#
# Code necessary for running a single cluster.

#--------------------------------------------------------#
# Read in Command Line Arguments                         #
#--------------------------------------------------------#
args<-commandArgs(TRUE)

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  stop("Arguments Must Be Provided!\n")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#--------------------------------------------------------#
# LOAD Necessary Code                                    #
#--------------------------------------------------------#
#--- Load Necessary Functions ---#
setwd("/netscr/drwilson/2018-04-05 Paper 1/Programs/")
source("Production Functions (Err_Control_MI_v2).R")
library(alabama)

#--- Load Data ---#
curr.dir = sprintf("/netscr/drwilson/2018-04-05 Paper 1/MappedData/Real_Data_Opt/Compiled_Output/")
setwd(curr.dir)
load(file_name)

#--------------------------------------------------------#
# Establish input break ups                              #
#--------------------------------------------------------#
# Save Directory:
setwd(sprintf("./%s/",file))

# Cell-Types:
cellTypes = c("cd4","cd8")

curr.clust.opt = tmp.data[c(start.pt:end.pt)]
curr.clust.out = STG.Update_Cluster.All(all_data=curr.clust.opt,cellTypes = cellTypes,optimType="nlminb",simple.Init=FALSE,initPts=c(0.5))

outFile = sprintf("./%s/%s_c%s_%s.RData",file,file,start.pt,end.pt)

save(curr.clust.out,file=outFile)
