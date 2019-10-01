#----------------------------------------------------------#
# MERGING DOWNSAMPLED DATA                                 #
#----------------------------------------------------------#

merge_mat = matrix(c("dsCD4pos_10.bam","dsCD8pos_90.bam","mf_CD4_10_CD8_90.bam",
                     "dsCD4pos_20.bam","dsCD8pos_80.bam","mf_CD4_20_CD8_80.bam",
                     "dsCD4pos_30.bam","dsCD8pos_70.bam","mf_CD4_30_CD8_70.bam",
                     "dsCD4pos_40.bam","dsCD8pos_60.bam","mf_CD4_40_CD8_60.bam",
                     "dsCD4pos_50.bam","dsCD8pos_50.bam","mf_CD4_50_CD8_50.bam",
                     "dsCD4pos_60.bam","dsCD8pos_40.bam","mf_CD4_60_CD8_40.bam",
                     "dsCD4pos_70.bam","dsCD8pos_30.bam","mf_CD4_70_CD8_30.bam",
                     "dsCD4pos_80.bam","dsCD8pos_20.bam","mf_CD4_80_CD8_20.bam",
                     "dsCD4pos_90.bam","dsCD8pos_10.bam","mf_CD4_90_CD8_10.bam"),byrow=TRUE,ncol=3)

merge_func<-function(merge_mat,directory){
  for(i in 1:nrow(merge_mat)){
    cmd1 = sprintf("java -jar /nas02/apps/picard-2.10.3/picard-2.10.3/picard.jar MergeSamFiles I=%s I=%s O=%s/%s SORT_ORDER=queryname",merge_mat[i,1],merge_mat[i,2],directory,merge_mat[i,3])
    system(cmd1)
  }
}

#----------------------------------------------------------#
# APPLY THE FUNCTION                                       #
#----------------------------------------------------------#
setwd("/netscr/drwilson/2018-04-05 Paper 1/MappedData/Mixtures/")
directory = "merged"


merge_func(merge_mat = merge_mat,directory = directory)