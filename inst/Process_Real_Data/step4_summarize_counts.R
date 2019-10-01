#---------------------------------------------------------------------------------
# Set Working Directory [EDIT]
#---------------------------------------------------------------------------------

setwd("/netscr/drwilson/2015-06-19 Refined Half Simulate/Mixture Files/_counts/")

#---------------------------------------------------------------------------------
# Load/Organize Gene List into Matrix Form
#---------------------------------------------------------------------------------

load("/netscr/drwilson/Homo_sapiens.GRCh37.66.nTE.RData")
dim(nTE)
nTE[1:2,]
length(unique(nTE$geneId))

# Lists all count files in the current directory. Assigns number of such files
# to the idx2 random variable. To be used in looping.

ffs  = list.files(pattern="_counts.txt")
nn   = length(ffs)
ffs[1:2]
nn
idx2<-nn

# Creates matrix of following form:
#  -> One geneID per row
#  -> One Cell Line Sample per column
#  -> Count for gene i in sample j is cell value

# Labels rows by GeneID, labels columns (samples) by file name.

sams  = gsub("_counts.txt", "", ffs)
couts = matrix(0, nrow=nrow(nTE), ncol=(idx2-1+1))
colnames(couts) = sams[1:idx2]
rownames(couts) = nTE$geneId

#---------------------------------------------------------------------------------
# READ IN SAMPLE DATA:
#
# Reads in output from step 3 (countReads) and organizes it into a matrix form.
# Column 1: Read Count
# Column 2: Reference site (Transcript Cluster, Gene ID, Exons involved) 
#---------------------------------------------------------------------------------

for(idx in 1:idx2){
  
  cat(idx, date(), "\n")
  
  f1   = ffs[idx]
  dat  = scan(f1, what=character(0))  
  dat  = matrix(dat, ncol=2, byrow=TRUE);
 
  colNames = c("count", "exons")
  cN = sprintf("%s and %s", colNames[1], colNames[2])
 
# Error Checking: if a row does not have two columns, the procedure
# halts and warns the user.

  if(ncol(dat) != 2){
    stop(countFile, " should have 2 columns: ", cN, "\n")
  }
  
  colnames(dat) = colNames
  dim(dat)
  dat[1:2,]
  
#---------------------------------------------------------------------------------
# OBTAIN TRANSCRIPT CLUSTERS AND GENE IDS
#---------------------------------------------------------------------------------
  
# Splits each string from column 2 (reference sites) of the dat matrix into 3
# components: (i) Transcript cluster, (ii) Ensembl Gene ID, (iii) Exon Cluster.
# Output is a list of exon sets, one list for each line of original dat matrix.

  groupIDs = strsplit(dat[,"exons"], split=";|\\|", perl=TRUE)
  
# Creates function that will be used to extract the UNIQUE gene IDs from within
# each exon set.

  splitFun <- function(vx){
    unique(vx[grep("ENSG", vx)])
  }
  
#Extract unique geneIDs for each exon set read count.

  date()
  geneIDs = lapply(groupIDs, splitFun)
  date()
  geneIDs[1:2]
  
#---------------------------------------------------------------------------------
# DATA CHECKING:
#     This section assesses the number of unique gene IDs present for each exon 
# set/read count. If more than one unique geneID is present for an exon set,
# its reads are excluded. If only one geneID is present, these reads are assigned
# to that gene.
#---------------------------------------------------------------------------------
 
  ngIDs   = sapply(geneIDs, length)
  table(ngIDs)
  
  w2check = which(ngIDs > 1)
  
  chkGIDs <- function(g1){
    gkp    = ""
    ncombo = sum(!grepl(":", g1, fixed=TRUE))
    
    if(ncombo <= 1){
      g2s = strsplit(g1, split=":", fixed=TRUE)
      gus = unique(unlist(g2s))
      foundONE = FALSE
      
      for(gu1 in gus){
        if (all(grepl(gu1, g1))){
          if(foundONE){  
            gkp = ""
            break
          }else{
            foundONE = TRUE
            gkp = gu1
          }
        }
      }
    }
    
    gkp
  }

  gIDchk  = sapply(geneIDs[w2check], chkGIDs)
  length(gIDchk)
  gIDchk[1:4]
  
  geneIDs[w2check] = gIDchk
  n1      = length(geneIDs)
  geneIDs = unlist(geneIDs)
  if(n1 != length(geneIDs)){ stop("non-unique geneIDs\n") }
  
  gID2rm = w2check[which(gIDchk=="")]
  str1   = "combinations are skipped because they belong to different genes"
  message(length(gID2rm), " exon ", str1)
  
  dim(dat)
  if(length(gID2rm) > 0){
    dat     = dat[-gID2rm,]
    geneIDs = geneIDs[-gID2rm]
  }
  
  dim(dat)
  length(unique(geneIDs))
  
#---------------------------------------------------------------------------------
# RECORDING COUNTS:
#---------------------------------------------------------------------------------
 
# Converts counts from data matrix from character to numeric. Sums across all
# exon sets with the same, unique (after previous step), geneID.

  cts     = as.numeric(dat[,"count"])
  geneCts = tapply(cts, geneIDs, sum)
 
# Matches geneIDs from the geneCts matrix above with the nTE data set geneIDs
# for proper almagamation.
 
  mat1    = match(names(geneCts), nTE$geneId)
  wNotNA  = which(! is.na(mat1))
 
# Provides information on the number of geneIDs skipped during the process
# due to: (i) having two geneIDs present in read. 

  pp1 = round(sum(geneCts[-wNotNA])/sum(geneCts),4)
  nn1 = length(geneCts) - length(wNotNA)
  message(100*pp1, "% of reads @ ", nn1, " exon combinations are skipped\n")

# Updates original countmatrix with counts from the current sample passed through
# the for loop.

  couts[mat1[wNotNA],idx-1+1] = geneCts[wNotNA]
}

# Outputs a dataset with the gene counts matrix.

outF = sprintf("gene_counts_%d_%d.txt", 1, idx2)

write.table(couts, file = outF, append = FALSE, quote = FALSE, sep = "\t",
row.names = TRUE, col.names = TRUE)


