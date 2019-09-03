#-----------------------------------------------------------------------#
# Compile Files                                                         #
#-----------------------------------------------------------------------#
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

for(j in 1:length(files)){
  dir_name = sprintf("/netscr/drwilson/2018-04-05 Paper 1/MappedData/Real_Data_Opt/Compiled_Output/%s/",files[j])
  setwd(dir_name)
  
  comp.out= NULL
  curr.clust.out = NULL
  
  #---- Set up new pattern ----#
  est.chunks = list.files(path="./",pattern = sprintf("%s_c%s_",files[j],files[j]))
  
  message("File ",files[j])
  
  #---- Populate Output Dataset ----#
  comp.out = list()
  r = 1
  
  for(i in 1:length(est.chunks)){
    load(est.chunks[i])
    nl = length(curr.clust.out)
    for(m in 1:nl){
      comp.out[[r]]=curr.clust.out[[m]]
      r = r+1
    }
  }
  
  save_file = sprintf("save(comp.out,file = \"./FinCompiled_%s.RData\")",files[j])
  eval(parse(text=save_file))
}
