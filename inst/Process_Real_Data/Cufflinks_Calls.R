#---------------------------------------------------------#
# TOPHAT CALLS:
#---------------------------------------------------------#
# Automated TOPHAT calls for sequence alignment.

CL_calls<-function(input_direc,j){
  cmd1 = sprintf("#! /bin/bash\n\ncd /netscr/drwilson/SRA_Data_dwnld/files/%s/\n\nbsub -o clinks_%s_wn.out -n 12 -R \"span[hosts=1]\" -x -q whole_node cufflinks -p 12 -o ./CLOUT2/ -g /netscr/drwilson/Reference_Annotations/Homo_sapiens/Homo_sapiens.GRCh37.66.updated.gtf %s_SBCv2.bam",
                 input_direc,input_direc,input_direc)
  write(cmd1,file = "tmp_out.txt")
  
  cmd2 = sprintf("bsub -o ./Run_Output/CLINKS_call_%s.txt < tmp_out.txt",j)
  system(cmd2)
  
  system("rm tmp_out.txt")
}

#---------------------------------------------------------#
# CALLING THE FUNCTIONS                                   #
#---------------------------------------------------------#
setwd("/netscr/drwilson/2018-04-05 Paper 1/Programs/")

input_direc = c("SRR1550995",
                "SRR1551016",
                "SRR1551024",
                "SRR1551051",
                "SRR1551071",
                "SRR1551079",
                "SRR1551065",
                "SRR1551098")

for(j in 1:length(input_direc)){
  CL_calls(input_direc = input_direc[j],j)
}