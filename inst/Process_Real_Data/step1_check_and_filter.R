library(asSeq, lib="/nas02/home/d/r/drwilson/R/Rlibraries/")

args<-commandArgs(TRUE)

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  inputDir = "/netscr/drwilson/2018-04-05 Paper 1/MappedData/CD4"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


#------------------------------------------------------------------
# Generating BAM file List            
#------------------------------------------------------------------

#Set working directory to folder where BAM files are located:
setwd(inputDir)

#Generate initial list of files: (NOTE- .bam.bai files included)
init_list = list.files(pattern=".bam")

#Generate list of .bai files
int_list = list.files(pattern=".bai")

#Final List (excluding .bai files):
BAM_list = setdiff(init_list,int_list)

# -----------------------------------------------------------------
# check the bam files
# -----------------------------------------------------------------

#Checks length of BAM files to ensure all has run properly:
length(BAM_list)

#Displays BAM list as another check for errors:
BAM_list

bam2use = BAM_list

#Loop across all BAM files in the folder:
for(i in 1:length(BAM_list)){

  bami = bam2use[i]
  
  #Generates a name that can be used in the files by stripping 
  #off the .BAM extension:
  sami = substr(bam2use[i],start=1,stop=nchar(bam2use[i])-4)
  
  # ----------------------------------------------------------
  # counting
  # ----------------------------------------------------------
  ctF  = sprintf("_count/count_%s.txt", sami)
  cmd1 = sprintf("samtools view %s | wc -l >> %s\n", bam2use[i], ctF)
  system(cmd1)
  
  # ----------------------------------------------------------
  # sorting
  # ----------------------------------------------------------
  cmd2 = sprintf("samtools sort -n %s _count/%s_sorted_by_name", bam2use[i], sami)
  system(cmd2)
  bamF = sprintf("_count/%s_sorted_by_name.bam", sami)
  
  # ----------------------------------------------------------
  # getUnique and filtering
  # ----------------------------------------------------------
  prepareBAM(bamF, sprintf("_count/%s_sorted_by_name", sami), sortIt=FALSE)
  
  system(sprintf("rm %s", bamF))
  
  # ----------------------------------------------------------
  # counting again
  # ----------------------------------------------------------
  cmd3   = sprintf("samtools view _count/%s_sorted_by_name_uniq_filtered.bam | wc -l >> %s\n", sami, ctF)
  system(cmd3)
}


