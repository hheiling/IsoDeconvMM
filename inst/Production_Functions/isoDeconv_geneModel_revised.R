#------------------------------------------------------------------------------------------------------------------------------------#
# ISODECONV                                                                                                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# NAME:                                                                                                                              #
#    isoDeconv_d2.R                                                                                                                  #
# DATE:                                                                                                                              #
#    7/23/2003                                                                                                                       #
# PROGRAMMER:                                                                                                                        #
#    Wei Sun, PhD                                                                                                                    #
# ALTERED BY:                                                                                                                        #
#    Douglas Roy Wilson, Jr.                                                                                                         #
# DESCRIPTION:                                                                                                                       # 
#    Edits Dr. Wei Sun's geneModel creation program to acommodate multiple cell types. Both the wrapper for the geneModel function   #
#    and the geneModel function itself are in need of several edits. All comments provided were written by Doug Wilson in an attempt #
#    to understand and thoroughly document the programs.                                                                             #
#------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------------------------------------------------------#
# SETTING VARIABLES                                                                                                                  #
#------------------------------------------------------------------------------------------------------------------------------------#
#  - countData:
#       Vector of file names of the summarized count files (.txt files). Read counts summarized by exon set.
#  - labels:
#       Vector of labels which will be used to name the read counts within clusters. 
#  - cellTypes:
#       Factor vector which details both the number and cell type for each read file. Must follow the orders
#       of the countData and labels vectors.
#  - bedFile:
#       BED File which contains the genome information (gene labels, etc.)
#  - knownIsoforms:
#       File which contains for each gene cluster information on the isoforms used by the cluster.
#  - fragSizeFile:
#       A single file which computes the distribution of reads of various lengths across all of the read experiments.
#  - output:
#       A label for the output file.
#  - readLen:
#       The read length for each of the experiments. Assumes homogeneous read length.
#  - lmax: 
#       Maximum read length known from prior knowledge.
#  - elenmin:
#       Minimum value of effective length which allows for errors in the sequencing/mapping process
#  - mix_sams:
#       Listing of names for use in the generation of exon set counts for the mixture files only. Restricts output to only exon sets
#       defined by the mixture samples to reduce modeling complexity.

dev_compiled_geneMod <- function(countData,labels,cellTypes,total_cts,bedFile,knownIsoforms,fragSizeFile,output,readLen,lmax,eLenMin){
  
#------------------------------------------------------------------------------------------------------------------------------------#
# LOADING THE DATA                                                                                                                   #
#------------------------------------------------------------------------------------------------------------------------------------#

# Calls to Dr. Sun's loadData function and generates one GeneModel for each sample
# present in the countData vector.

for(i in 1:length(countData)){
  ct_datai = countData[i]
  outFile = labels[i]
  cmdi = sprintf("%s = loadData_djEdit(ct_datai,bedFile,readLen,lmax)",outFile)
  eval(parse(text=cmdi))
}

#------------------------------------------------------------------------------------------------------------------------------------#
# LIST APPEND: ALL SAMPLES                                                                                                           #
#------------------------------------------------------------------------------------------------------------------------------------#

sam_names = labels

# list_append function: concatenates all geneModels into a large list, labeled by sample names.

list_append <- function(list_names){
  nl = length(list_names)
  concat_list = list()
  for(i in 1:nl){
    cmdi = sprintf("concat_list[[i]] = %s",list_names[i])
    eval(parse(text=cmdi))
  }
  names(concat_list) = list_names
  return(concat_list)
}

concat_list = list_append(sam_names)

#------------------------------------------------------------------------------------------------------------------------------------#
# OBTAINING LIST OF TRANSCRIPT CLUSTER NAMES                                                                                         #
#------------------------------------------------------------------------------------------------------------------------------------#

tnames = NULL

for(i in 1:length(sam_names)){
  call1 = sprintf("tclust_i = names(%s)",sam_names[i])
  eval(parse(text=call1))
  call2 = sprintf("tclust_iand1 = names(%s)",sam_names[i+1])
  eval(parse(text=call2))
  if(i==1){tnames = union(tclust_i,tclust_iand1)}
  else{tnames = union(tnames,tclust_i)}
}

tnames = tnames[order(tnames)]

# Create a list of length length(tnames) wherein the i-th element of the list is a list composed of #samples+1 datasets and a row vector.
# The vector is info_status and confirms whether or not the info dataset is complete and accurate.
# The first dataset is $info, which contains information on all exons in the cluster. The subsequent datasets are the counts at each
# exon set in each sample.

concat_geneMod = list()

for(i in 1:length(tnames)){
  concat_geneMod[[i]] = list(info_status="Not Checked",info=NULL)
  incurr_clust = rep(0,length(sam_names))
  z = rep(TRUE,length(sam_names)-1)
  for(j in 1:length(sam_names)){
    tclust_currsamp = names(concat_list[[j]])
    if(tnames[i] %in% tclust_currsamp){
      cmdij = sprintf("concat_geneMod[[i]]$%s = concat_list[[j]][[tnames[i]]]$count",sam_names[j])
      eval(parse(text=cmdij))
      incurr_clust[j] = 1
    } else {
      cmdij = sprintf("concat_geneMod[[i]]$%s = data.frame(count = as.numeric(NULL),exons=as.character(NULL))",sam_names[j])
      eval(parse(text=cmdij))
      incurr_clust[j] = 0
    }
  }
  sam_clust = which(incurr_clust==1)
  fsamp = sam_clust[1]
  concat_geneMod[[i]]$info = concat_list[[fsamp]][[tnames[i]]]$info
  for(k in 1:length(sam_clust)){
    fsamp_k = sam_clust[k]
    dim_orig = dim(concat_geneMod[[i]]$info)
    dim_final = dim(concat_list[[fsamp_k]][[tnames[i]]]$info)
    if(all(dim_orig == dim_final)){
      z[k-1] = all(concat_geneMod[[i]]$info==concat_list[[fsamp_k]][[tnames[i]]]$info)
    } else {
      z[k-1] = FALSE
    }
  }
  if (all(z)){
    concat_geneMod[[i]]$info_status = "PASS"
  } else {
    concat_geneMod[[i]]$info_status = "FAIL"
  }
}

names(concat_geneMod) = tnames

#-------------------------------------------------------------------------------------------------------------------------------------------#
# CONCAT_GENEMOD:                                                                                                                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FORMAT:                                                                                                                                   #
#     List                                                                                                                                  #
# I-TH ELEMENT:                                                                                                                             #
#     - Label     : "Transcript Cluster I"                                                                                                  #                                                                                                   #
#     - Element 1 : $info_status (vector)                                                                                                   #
#         - "Not Checked" -- error in the program and $info data set was not checked for correctness.                                       #
#         - "PASS" -- $info dataset passes checks for accuracy                                                                              #
#         - "FAIL" -- $info dataset failed checks for accuracy                                                                              # 
#     - Element 2 : $info dataset                                                                                                           #
#         - chr   : contains the chromosome location of the exon being considered                                                           #
#         - start : contains the starting location of the exon being considered                                                             #
#         - end   : contains the ending location of the exon being considered                                                               #
#         - exon  : contains the unique exon id for the exon being considered                                                               #
#     - Element 3 : $sample_name dataset                                                                                                    #   
#         - exons : Contains the IDs of the exons comprising the exon set                                                                   #
#         - count : contains the number of reads at the exon set in question for "sample_name"                                              #
#-------------------------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------------------------------#
# GENERATING PDDIST:                                                                                                                        #
#-------------------------------------------------------------------------------------------------------------------------------------------#
# Generates an estimate of the distribution of fragment lengths necessary for the computation of effective length.

pdDist_gen <- function(fragSizeFile,lmax){
  md = read.table(fragSizeFile)
  if (ncol(md) != 2) {
    stop(fragSizeFile, " should have 2 columns for Freq and Len\n")
  }
  names(md) = c("Freq", "Len")
  pd = rep(0, lmax)
  w2 = which(md$Len <= lmax)
  pd[md$Len[w2]] = md$Freq[w2]
  pdDist = pd/sum(pd)
  return(pdDist)
}

pdDist = pdDist_gen(fragSizeFile,lmax)

#-------------------------------------------------------------------------------------------------------------------------------------------#
# CHECKING Presence of Isoforms File:                                                                                                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# IsoDeconv requires a list of known isoforms in order to model intra-sample heterogeneity. Stops program if file not present.              #
#-------------------------------------------------------------------------------------------------------------------------------------------#

if(!is.null(knownIsoforms)){
  load(knownIsoforms)
}
else {
  stop("Isoforms file Must be present!")
}

#-------------------------------------------------------------------------------------------------------------------------------------------#
# CALL TO GENEMODEL:                                                                                                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Calls to the geneModel function for each transcript cluster.                                                                              #
#-------------------------------------------------------------------------------------------------------------------------------------------#

fin_geneMod = list()
nms = names(concat_geneMod)
sep = "------"
mix_sams = sam_names[which(tolower(cellTypes)=="mix")]

for (i in 1:length(nms)) {
  if (i%%100 == 0) {
      message(sep, i, "  ", date(), sep)
  }
  nm1 = nms[i]
  ge1 = concat_geneMod[[nm1]]
  if (!is.null(knownIsoforms)) {
    isoforms = isoAll[[nm1]]
  } 
  # Call to geneModel function: ge1 is a list with components $count and $info for a single transcript
  # cluster. isoforms is a matrix which details all of the isoforms used in a transcript cluster with
  # indicators for whether a particular exon is used by an isoform.
  
  gm1 = geneModel_multcell_Edit(ge1, d = readLen, pdDist, isoforms, lmax, 
                  eLenMin, verbose=1,sam_names=sam_names,mix_sams=mix_sams)
  fin_geneMod[[nm1]] = gm1
}

cellType_count = sum(unique(tolower(cellTypes))!="mix")
info_mat = data.frame(Label = labels, Cell_Type = tolower(cellTypes), Total = total_cts,stringsAsFactors = FALSE)
fin_geneMod["Sample_Info"] = list(info = info_mat,tclust_tot=length(fin_geneMod),
                                      cellType_count=cellType_count)

# #-------------------------------------------------------------------#
# # EDIT TO GROUP CELL TYPES                                          #
# #-------------------------------------------------------------------#
# info_mat = fin_geneMod[["Sample_Info"]]$info
# cellTypes = unique(info_mat$Cell_Type)
# 
# ctList = list()
# 
# for(j in 1:length(cellTypes)){
#   idx = which(info_mat$Cell_Type==cellTypes[j])
#   ctList[[cellTypes[j]]] = list(samps = info_mat$Label[idx], tots = info_mat$Total[idx])
# }
# 
# idx2consider = which(names(fin_geneMod)!="Sample_Info")
# for(k in idx2consider){
#   for(l in 1:length(cellTypes)){
#     samps2use = ctList[[l]]$samps
#     tots      = ctList[[l]]$samps
#     
#     y_vecs  = paste("fin_geneMod[[k]]$y",samps2use,sep = "_")
#     y_vecsc = paste(y_vecs,collapse = ",")
#     nExon = eval(parse(text=sprintf("length(%s)",y_vecs[1])))
#     textcmd = sprintf("matrix(c(%s),nrow=nExon,ncol=length(samps2use))",y_vecsc)
#     expMat  = eval(parse(text=textcmd))
#     
#     totmg   = tots-colSums(expMat)
#     expMat2 = rbind(totmg,expMat)
#     
#     fin_geneMod[[k]][[cellTypes[l]]] = expMat2
#   }
# }

save(fin_geneMod, file = output)
}

#-------------------------------------------------------------------#
# EDIT TO GROUP CELL TYPES                                          #
#-------------------------------------------------------------------#
# info_mat = fin_geneMod[["Sample_Info"]]$info
# cellTypes = unique(info_mat$Cell_Type)
# 
# ctList = list()
# 
# for(j in 1:length(cellTypes)){
#   idx = which(info_mat$Cell_Type==cellTypes[j])
#   ctList[[cellTypes[j]]] = list(samps = info_mat$Label[idx], tots = info_mat$Total[idx])
# }
# 
# idx2consider = which(names(fin_geneMod)!="Sample_Info")
# for(k in idx2consider){
#   for(l in 1:length(cellTypes)){
#     samps2use = ctList[[l]]$samps
#     tots      = ctList[[l]]$samps
#     
#     y_vecs  = paste("y",samps2use,sep = "_")
#     y_vecsc = paste(y_vecs,collapse = ",")
#     nExon = eval(parse(text=sprintf("length(%s)",y_vecs[1])))
#     textcmd = sprintf("matrix(c(%s),nrow=nExon,ncol=length(samps2use))",y_vecsc)
#     expMat  = eval(parse(text=textcmd))
#     
#     totmg   = tots-colSums(expMat)
#     expMat2 = rbind(totmg,expMat)
#     
#     fin_geneMod[[k]][[cellTypes[l]]] = expMat2
#   }
# }