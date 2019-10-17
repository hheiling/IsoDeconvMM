#-------------------------------------------------------------#
# Fragment Length Files: Generation                           #
#-------------------------------------------------------------#

# Specify Directory where read files are kept:
setwd("/netscr/drwilson/2018-04-05 Paper 1/MappedData/")

# BAM Files where sequenced reads are located:
inputFiles = c("./Mixtures/merged/mf_CD4_0_CD8_100.bam",
               "./Mixtures/merged/mf_CD4_10_CD8_90.bam",
               "./Mixtures/merged/mf_CD4_20_CD8_80.bam",
               "./Mixtures/merged/mf_CD4_30_CD8_70.bam",
               "./Mixtures/merged/mf_CD4_40_CD8_60.bam",
               "./Mixtures/merged/mf_CD4_50_CD8_50.bam",
               "./Mixtures/merged/mf_CD4_60_CD8_40.bam",
               "./Mixtures/merged/mf_CD4_70_CD8_30.bam",
               "./Mixtures/merged/mf_CD4_80_CD8_20.bam",
               "./Mixtures/merged/mf_CD4_90_CD8_10.bam",
               "./Mixtures/merged/mf_CD4_100_CD8_0.bam",
               "./CD4/_count/SRR1550995_SBC_sorted_by_name_uniq_filtered.bam",
               "./CD4/_count/SRR1551016_SBC_sorted_by_name_uniq_filtered.bam",
               "./CD4/_count/SRR1551098_SBC_sorted_by_name_uniq_filtered.bam",
               "./CD8/_count/SRR1551024_SBC_sorted_by_name_uniq_filtered.bam",
               "./CD8/_count/SRR1551051_SBC_sorted_by_name_uniq_filtered.bam",
               "./CD8/_count/SRR1551065_SBC_sorted_by_name_uniq_filtered.bam")

# Associated Labels for output files:
outputLabels = c("cd4_0_cd8_100",
                 "cd4_10_cd8_90",
                 "cd4_20_cd8_80",
                 "cd4_30_cd8_70",
                 "cd4_40_cd8_60",
                 "cd4_50_cd8_50",
                 "cd4_60_cd8_40",
                 "cd4_70_cd8_30",
                 "cd4_80_cd8_20",
                 "cd4_90_cd8_10",
                 "cd4_100_cd8_0",
                 "cd4_r1",
                 "cd4_r2",
                 "cd4_r3",
                 "cd8_r1",
                 "cd8_r2",
                 "cd8_r3")

# If any files are to be combined, list them in separate units:
comboList = list()
comboList[[1]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_0_cd8_100_lengths.txt")
comboList[[2]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_10_cd8_90_lengths.txt")
comboList[[3]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_20_cd8_80_lengths.txt")
comboList[[4]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_30_cd8_70_lengths.txt")
comboList[[5]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_40_cd8_60_lengths.txt")
comboList[[6]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_50_cd8_50_lengths.txt")
comboList[[7]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_60_cd8_40_lengths.txt")
comboList[[8]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_70_cd8_30_lengths.txt")
comboList[[9]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_80_cd8_20_lengths.txt")
comboList[[10]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_90_cd8_10_lengths.txt")
comboList[[11]] = c("cd4_r1_lengths.txt","cd4_r2_lengths.txt","cd4_r3_lengths.txt",
                   "cd8_r1_lengths.txt","cd8_r2_lengths.txt","cd8_r3_lengths.txt",
                   "cd4_100_cd8_0_lengths.txt")

# Combo Output Labels:
comboLabels = c("cd4_0_cd8_100",
                 "cd4_10_cd8_90",
                 "cd4_20_cd8_80",
                 "cd4_30_cd8_70",
                 "cd4_40_cd8_60",
                 "cd4_50_cd8_50",
                 "cd4_60_cd8_40",
                 "cd4_70_cd8_30",
                 "cd4_80_cd8_20",
                 "cd4_90_cd8_10",
                 "cd4_100_cd8_0")


fragLengths<-function(Input_Files,outputLabels,comboList,comboLabels,useCombo){
  outLabels = paste(outputLabels,"_lengths.txt",sep="")
  for(i in 1: length(Input_Files)){
    cmd1 = sprintf("samtools view -f 65 %s | awk '{print ($8>=$4) ? $8-$4+51 : $4-$8+51}' > %s",inputFiles[i],outLabels[i])
    system(cmd1)
  }
  
  if(missing(comboList) && useCombo==0){
    for(j in 1:length(outLabels)){
      cmd2_a = sprintf("cat %s | sort -n | uniq -c > %s_fraglens.txt",outLabels[j],outputLabels[j])
      system(cmd2_a)
    }
  } else if(useCombo==1 && missing(comboList)){
    stop("If you wish to combine files, you must list which files are to be combined!")
    
  } else { # comboList specified, and useCombo = 1
    
    ftc = unlist(lapply(X = comboList,FUN = function(x) {return(paste(x,collapse=" "))}))
    
    for(k in 1:length(comboList)){
      cmd2_b = sprintf("cat %s | sort -n | uniq -c > %s_fraglens.txt",ftc[k],comboLabels[k])
      system(cmd2_b)
    }
    
  }
}

fragLengths(Input_Files = inputFiles,outputLabels=outputLabels,comboList = comboList,comboLabels=comboLabels,useCombo=1)