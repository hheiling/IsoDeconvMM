#-------------------------------------------------------------------#
# Generate Downsampled Files:                                       #
#-------------------------------------------------------------------#

downsample_BAM<-function(infile,outlabel,directory,props.out,desct){
  cmd1 = sprintf("samtools view %s | wc -l > tmp_count.txt",infile)
  system(cmd1)
  total_ct = scan("tmp_count.txt")
  
  for(j in props.out){
    prob_val = (j/100)*(desct)/total_ct
    cmd2 = sprintf("java -jar /nas02/apps/picard-1.88/picard-tools-1.88/DownsampleSam.jar INPUT=%s OUTPUT=%s%s_%s.bam PROBABILITY=%s VALIDATION_STRINGENCY=LENIENT RANDOM_SEED=null",infile,directory,outlabel,j,prob_val)
    cmd2
    system(cmd2)
  }
  
  system("rm tmp_count.txt")
}

#---------------- APPLY DOWNSAMPLE_BAM -----------------------------#
# Set the directory where input files are located:
setwd("/netscr/drwilson/2018-04-05 Paper 1/MappedData/")

# Insert Input Files:
inFiles = c("/CD4/_count/SRR1551071_SBC_sorted_by_name_uniq_filtered.bam",
            "/CD8/_count/SRR1551079_SBC_sorted_by_name_uniq_filtered.bam")

# Insert Output Labels:
outlabels = c("dsCD4pos",
              "dsCD8pos")

# Set Output Directory:
outDirec = "/netscr/drwilson/2018-04-05 Paper 1/MappedData/Mixtures/"

# Proportion Lists:
prop.vecs = list()
prop.vecs[[1]] = seq(from = 10,to=100,by=10)
prop.vecs[[2]] = seq(from = 10,to=100,by=10)

for(i in 1:length(inFiles)){
  downsample_BAM(infile = inFiles[i],outlabel = outlabels[i],
                 directory=outDirec,props.out=prop.vecs[[i]])
}
