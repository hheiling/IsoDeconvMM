#---------------------------------------------------------#
# TOPHAT CALLS:
#---------------------------------------------------------#
# Automated TOPHAT calls for sequence alignment.

TH_calls<-function(PE_reads_mat,input_direc,j){
  cmd1 = sprintf("#! /bin/bash\n\ncd /netscr/drwilson/SRA_Data_dwnld/files/%s/\n\nbsub -n 12 -R \"span[hosts=1]\" -x -q whole_node STAR --runThreadN 12 --genomeDir /proj/seq/data/STAR_genomes/hg19/ --sjdbGTFfile /netscr/drwilson/Reference_Annotations/Homo_sapiens/Homo_sapiens.GRCh37.66.updated.gtf --sjdbOverhang 49 --readFilesIn %s %s --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --alignEndsType EndToEnd",
                 input_direc,PE_reads_mat[1],PE_reads_mat[2])
  write(cmd1,file = "tmp_out.txt")
  
  cmd2 = sprintf("bsub -o ./Run_Output/STAR_call_%s.txt < tmp_out.txt",j)
  system(cmd2)
  
  system("rm tmp_out.txt")
}

#---------------------------------------------------------#
# CALLING THE FUNCTIONS                                   #
#---------------------------------------------------------#
setwd("/netscr/drwilson/2018-04-05 Paper 1/Programs/")
PE_reads_mat = matrix(c("SRR1550989_1.fastq","SRR1550989_2.fastq",
                        "SRR1551000_1.fastq","SRR1551000_2.fastq",
                        "SRR1551011_1.fastq","SRR1551011_2.fastq",
                        "SRR1551023_1.fastq","SRR1551023_2.fastq",
                        "SRR1551043_1.fastq","SRR1551043_2.fastq",
                        "SRR1551050_1.fastq","SRR1551050_2.fastq",
                        "SRR1551111_1.fastq","SRR1551111_2.fastq",
                        "SRR1551006_1.fastq","SRR1551006_2.fastq",
                        "SRR1551012_1.fastq","SRR1551012_2.fastq",
                        "SRR1551017_1.fastq","SRR1551017_2.fastq",
                        "SRR1551037_1.fastq","SRR1551037_2.fastq",
                        "SRR1551086_1.fastq","SRR1551086_2.fastq",
                        "SRR1551099_1.fastq","SRR1551099_2.fastq",
                        "SRR1551105_1.fastq","SRR1551105_2.fastq"),ncol=2,byrow = TRUE)

# PE_reads_mat = matrix(c("SRR1550995_1.fastq","SRR1550995_2.fastq",
#                         "SRR1551016_1.fastq","SRR1551016_2.fastq",
#                         "SRR1551024_1.fastq","SRR1551024_2.fastq",
#                         "SRR1551051_1.fastq","SRR1551051_2.fastq",
#                         "SRR1551071_1.fastq","SRR1551071_2.fastq",
#                         "SRR1551079_1.fastq","SRR1551079_2.fastq",
#                         "SRR1551065_1.fastq","SRR1551065_2.fastq",
#                         "SRR1551098_1.fastq","SRR1551098_2.fastq"),ncol=2,byrow = TRUE)

insert_len = 0

input_direc = c("SRR1550989",
                "SRR1551000",
                "SRR1551011",
                "SRR1551023",
                "SRR1551043",
                "SRR1551050",
                "SRR1551111",
                "SRR1551006",
                "SRR1551012",
                "SRR1551017",
                "SRR1551037",
                "SRR1551086",
                "SRR1551099",
                "SRR1551105")

for(j in 1:length(input_direc)){
  TH_calls(PE_reads_mat = PE_reads_mat[j,],input_direc = input_direc[j],j)
}