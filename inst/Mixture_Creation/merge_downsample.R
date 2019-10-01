#----------------------------------------------------------#
# MERGING DOWNSAMPLED DATA                                 #
#----------------------------------------------------------#

merge_mat = matrix(c("dsM11_10.bam","dsM21_90.bam","mf_m11_10_m21_90.bam",
                     "dsM11_20.bam","dsM21_80.bam","mf_m11_20_m21_80.bam",
                     "dsM11_30.bam","dsM21_70.bam","mf_m11_30_m21_70.bam",
                     "dsM11_40.bam","dsM21_60.bam","mf_m11_40_m21_60.bam",
                     "dsM11_50.bam","dsM21_50.bam","mf_m11_50_m21_50.bam",
                     "dsM11_60.bam","dsM21_40.bam","mf_m11_60_m21_40.bam",
                     "dsM11_70.bam","dsM21_30.bam","mf_m11_70_m21_30.bam",
                     "dsM11_80.bam","dsM21_20.bam","mf_m11_80_m21_20.bam",
                     "dsM11_90.bam","dsM21_10.bam","mf_m11_90_m21_10.bam"),byrow=TRUE,ncol=3)

merge_func<-function(merge_mat,directory){
  for(i in 1:nrow(merge_mat)){
    cmd1 = sprintf("samtools merge %s/%s %s %s",directory,merge_mat[i,3],merge_mat[i,1],merge_mat[i,2])
    system(cmd1)
  }
}

#----------------------------------------------------------#
# APPLY THE FUNCTION                                       #
#----------------------------------------------------------#
setwd("/netscr/drwilson/2015-08-10_Macrophages/Downsampled/")
directory = "/netscr/drwilson/2015-08-10_Macrophages/Downsampled"


merge_func(merge_mat = merge_mat,directory = directory)