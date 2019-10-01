#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 1 - COPY NUMBER CODE                       #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
CN_Calc <- function(z,subjFile,outVal){
  wChr   = as.character(z[2])
  chrObj = subset(subjFile$CN_segs,subset = (Chr==wChr))
  
  gs = matrix(as.numeric(z[3]),ncol=1,nrow=nrow(chrObj))
  ge = matrix(as.numeric(z[4]),ncol=1,nrow=nrow(chrObj))
  
  inInt   = !(gs>chrObj$End_bp_hg19|ge<chrObj$Start_bp_hg19)
  overlap = rowMins(cbind(ge,chrObj$End_bp_hg19))-
            rowMaxs(cbind(gs,chrObj$Start_bp_hg19))
  
  if(all(inInt==FALSE)){
    return(NA)
  } else {
    overlap[which(inInt==FALSE)] = -10000
    
    R2use = which.max(overlap)
    
    if(outVal=="CNT"){
      return(chrObj$CN_1[R2use]+chrObj$CN_2[R2use])
    } else if(outVal=="CN1"){
      return(chrObj$CN_1[R2use])
    } else if(outVal=="CN2"){
      return(chrObj$CN_2[R2use])
    } else {
      stop("outVal specification is invalid! Use CNT, CN1, or CN2 only.")
    }
    
  }
}

rowMins <- function(matrix){
  return(as.matrix(apply(X = matrix,MARGIN = 1,FUN=min)))
}

rowMaxs<- function(matrix){
  return(as.matrix(apply(X = matrix,MARGIN = 1,FUN=max)))
}

#-------------------------------------------------------------------#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#                SECTION 2 - CYCLE CODE                             #
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#-------------------------------------------------------------------#
#GeneAnno = load("/netscr/drwilson/Reference_Annotations/GeneLists.RData")
GeneAnno = load("/nas/longleaf/home/drwilson/Paper 2/Annotations/GeneLists.RData")

sigGene = FALSE
CNvers  = "CN1"

if(sigGene){
  GenesFile = redGenelist
} else {
  GenesFile = fullGenelist
}

### Load Paul's Files
setwd("/pine/scr/p/l/pllittle/BRCA/ASCAT/")

subjList = list.files()

subjList  = list.files()
subjName  = rep("",length(subjList))
subjPur   = rep(0,length(subjList))

CN_Dat = matrix(NA,nrow=length(subjList),ncol=nrow(GenesFile))

for(i in 1:length(subjList)){
  direc = subjList[i]
  setwd(direc)
  
  if(("pur_plo.rds"%in%list.files())){
    csubjFile = readRDS("pur_plo.rds")
    subjName[i] = substr(direc,start = 1,stop = 16)
    subjPur[i]  = csubjFile$purity
    
    if(csubjFile$purity<1.0){
      csubjOut  = apply(X = GenesFile[,-c(2)],MARGIN = 1,
                        FUN = CN_Calc,subjFile=csubjFile,outVal=CNvers)
      
      CN_Dat[i,] = csubjOut
    }
  }
  
  if(i%%100==0){message("Subject ",i," Complete!")}
  setwd("../")
}

colnames(CN_Dat) = GenesFile$Gene
names(subjPur) = subjName

CN_Dat = data.frame(FileName = subjList,barcode=subjName,CN_Dat)

save(CN_Dat,
     file=sprintf("/nas/longleaf/home/drwilson/Paper 2/All_%s_Data.RData",CNvers))