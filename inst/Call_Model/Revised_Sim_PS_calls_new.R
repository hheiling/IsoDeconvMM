#-------------------------------------------------------------------#
# REVISED SIMULATION CODE                                           #
#-------------------------------------------------------------------#
library(gtools)

setwd("/netscr/drwilson/2017-02-17 Paper 1/Programs/")

#-------------------------------------------------------------------#
files = c("CD4_0_CD8_100",
          "CD4_10_CD8_90",
          "CD4_20_CD8_80",
          "CD4_30_CD8_70",
          "CD4_40_CD8_60",
          "CD4_50_CD8_50",
          "CD4_60_CD8_40",
          "CD4_70_CD8_30",
          "CD4_80_CD8_20",
          "CD4_90_CD8_10",
          "CD4_100_CD8_0")

for(i in 1:length(files)){
  cmd1 = sprintf("cd /netscr/drwilson/2018-04-05\ Paper\ 1/Programs/\n\nbsub -o ./Run_Output/PS_calls_realdata.txt R CMD BATCH --no-save --no-restore '--args file=\"%s\"' Revised_Sim_PS_codeNew.R ./Run_Output/%s_ps.Rout",files[i],files[i])
  system(cmd1)
}