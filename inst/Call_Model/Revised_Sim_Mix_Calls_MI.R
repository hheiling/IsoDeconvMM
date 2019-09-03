#--------------------------------------------------------#
# RUN on a single Cluster                                #
#--------------------------------------------------------#
# Code necessary for running a single cluster.

#--------------------------------------------------------#
# LOAD Necessary Code                                    #
#--------------------------------------------------------#
#--- Load Necessary Functions ---#
setwd("/netscr/drwilson/2018-04-05 Paper 1/Programs/")
source("Production Functions (Err_Control_MI_v2).R")

#--- Set Directory ---#
curr.dir = sprintf("/netscr/drwilson/2018-04-05 Paper 1/MappedData/Real_Data_Opt/Compiled_Output/")
setwd(curr.dir)

files = list.files(pattern = "pure_est.RData")
samps = gsub(pattern = "_pure_est.RData",replacement = "",x = files)

for(i in 1:length(files)){
  #--------------------------------------------------------#
  # LOAD CURRENT FILE                                      #
  #--------------------------------------------------------#
  #--- Set Directory ---#
  # Combat any directory switches 
  curr.file = sprintf("./%s_pure_est.RData",samps[i])
  load(curr.file) 
  
  #--------------------------------------------------------#
  # Establish input break ups                              #
  #--------------------------------------------------------#
  # Cell-Types:
  cellTypes = c("cd4","cd8")
  
  # Data Set Necessities:
  clust.start = 1
  clust.end = length(tmp.data)
  by.value = 15
  
  start.pts = seq(from = 1,to = clust.end,by = by.value)
  end.pts = c((start.pts[-1]-1),clust.end)
  
  setwd("/netscr/drwilson/2018-04-05 Paper 1/Programs/")
  for(m in 1:length(start.pts)){
    fname = sprintf("%s_pure_est.RData",samps[i])
    cloud.job.cmd = sprintf("bsub -o ./Run_Output/job_sub_%s.txt R CMD BATCH --no-save --no-restore '--args start.pt=%s end.pt=%s file_name=\"%s\" file=\"%s\"' Revised_Sim_MixCode_SI.R ./Run_Output/rsub_mix_%s.txt",
                            samps[i],start.pts[m],end.pts[m],fname,samps[i],samps[i])
    system(cloud.job.cmd)
  }
  curr.dir = sprintf("/netscr/drwilson/2018-04-05 Paper 1/MappedData/Real_Data_Opt/Compiled_Output/")
  setwd(curr.dir)
}
