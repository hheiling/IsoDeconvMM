#-------------------------------------------------------------------#
# READ in COMMAND LINE ARGUMENTS                                    #
#-------------------------------------------------------------------#
args<-commandArgs(TRUE)

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  file = "CD4_0_CD8_100"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#-------------------------------------------------------------------#
# REVISED SIMULATION CODE                                           #
#-------------------------------------------------------------------#
library(gtools)
library(alabama)

setwd("/netscr/drwilson/2018-04-05 Paper 1/Programs/")
  source("PS Production Functions (Err Control).R")

setwd("/netscr/drwilson/2018-04-05 Paper 1/MappedData/GeneModels/")

#-----------------------------------------------------------#
# CALL Pure Sample                                          #
#-----------------------------------------------------------#
cellTypes = c("cd4","cd8")

f.curr = sprintf("./rechar_%s_geneMod.RData",file)
load(f.curr)

sim.out = sig_geneMod[which(names(sig_geneMod)!="Sample_Info")]

# Clusters with single isoforms:
#  EXCLUDE THEM FOR THE MOMENT!.
dim_mat = matrix(0,nrow=length(sim.out),ncol=2)
excl_clust = c()
excl_clust2 = c()
for(i in 1:length(sim.out)){
  dim_mat[i,] = dim(sim.out[[i]][["X"]])
  if(all(dim_mat[i,]==c(1,1))){
    excl_clust = c(excl_clust,i)
  }
  if(dim_mat[i,2] == 1){
    excl_clust2 = c(excl_clust2,i)
  }
}

excl_clust_union = union(excl_clust,excl_clust2)
if(length(excl_clust_union)>0){
  sim.new = sim.out[-c(excl_clust_union)]
} else {
  sim.new = sim.out
}


# Optimize the Pure Sample Functions:
tmp.data = Pure.apply.fun(data.list = sim.new,cellTypes = cellTypes,corr_co = 1)
save_cmd = sprintf("save(tmp.data,file=\"../Compiled_Output/%s_pure_est.RData\")",file)
eval(parse(text=save_cmd))

message("Sample ",file,"Completed")
