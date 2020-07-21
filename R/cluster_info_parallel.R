# parallelization of cluster_info

#' @importFrom stringr str_split
#' @import parallel
#' @import doParallel
#' @importFrom pbapply pboptions pblapply
cluster_info_parallel = function(countData, labels, cellTypes, bedFile, discrim_genes, 
                                 discrim_clusts, readLen, lmax, max_cores = 1, cluster_type){
  
  #------------------------------------------------------------------------------------------------------------------------------------#
  # LOADING THE DATA                                                                                                                   #
  #------------------------------------------------------------------------------------------------------------------------------------#
  
  # Calls to Dr. Sun's loadData function and generates one GeneModel for each sample
  # present in the countData vector.
  
  sam_names = labels
  names(countData) = sam_names
  
  #------------------------------------------------------------------------------------------------------------------------------------#
  # LOAD DATA                                                                                                          #
  #------------------------------------------------------------------------------------------------------------------------------------#
  
  
  data_extract = function(idx, labels, countData_all, bedFile, readLen, lmax){
    
    sam_names = labels[idx]
    # Select elements of count data list object to use in single cluster
    countData = countData_all[idx]
    
    concat_list = list()
    
    # See "R/loadData_edit.R" file for loadData_djEdit() function code
    
    for(i in 1:length(countData)){
      ct_datai = countData[[i]]
      outFile = sam_names[i]
      cat("\n Loading data for ", outFile, "\n")
      concat_list[[i]] = loadData_djEdit(ct_datai,bedFile,readLen,lmax)
      
    }
    
    names(concat_list) = sam_names
    
    return(concat_list)
    
  } # End data_extract function
  
  
  # If max_cores > 1, incorporate parallelization options
  if(max_cores > 1){
    # If max_cores > length(countData) then only use length(countData) number of cores
    if(max_cores > length(countData)){
      num_cores = length(countData)
    }else{
      num_cores = max_cores
    }
    
    size = ceiling(length(countData)/num_cores)
    full_idx = 1:length(countData)
    chunks = split(full_idx, ceiling(seq_along(full_idx)/size))
    
    if(cluster_type == "Socket"){
      cl = parallel::makeCluster(num_cores)
      registerDoParallel(cl)
      vars = c("countData","labels","bedFile","readLen","lmax","loadData_djEdit")
      # vars = list(countData, labels, bedFile, readLen, lmax)
      clusterExport(cl, varlist=vars, envir = environment())
      # clusterExport(cl, varlist = vars, envir = environment())
    }else if(cluster_type == "Fork"){
      cl = parallel::makeForkCluster(num_cores, outFile = "")
      registerDoParalle(cl)
    }
    
    pboptions(type="timer")
    info_output = pblapply(chunks, FUN = data_extract, countData_all = countData, labels = labels,
                           bedFile = bedFile, readLen = readLen, lmax = lmax, cl = cl)
    parallel::stopCluster(cl)
    
    concat_list = info_output[[1]]
    for(i in 2:length(info_output)){
      concat_list = c(concat_list, info_output[[i]])
    }
    
  }else{ # max_cores = 1 (no parallelization allowed)
    
    concat_list = data_extract(idx = 1:length(countData), labels = labels, 
                               countData_all = countData,bedFile = bedFile, 
                               readLen = readLen, lmax = lmax)
    
  } # End max_cores if-else
  
  #------------------------------------------------------------------------------------------------------------------------------------#
  # OBTAINING LIST OF TRANSCRIPT CLUSTER NAMES                                                                                         #
  #------------------------------------------------------------------------------------------------------------------------------------#
  
  tnames = NULL
  tnames = names(concat_list[[1]])
  for(i in 2:length(concat_list)){
    tclust_i = names(concat_list[[i]])
    tnames = union(tnames, tclust_i)
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
  
  #-----------------------------------------------------------------------------#
  # ESTABLISHING CLUSTERS WITH HIGHEST LEVELS OF DISCRIMINATORY CAPABILITY      #
  #-----------------------------------------------------------------------------#
  # Limit clusters examined to those with discriminatory genes 
  # (discrim_genes or discrim_clusts)                                                                                                                       #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Gene names for cluster given in "info" matrix                               #
  #-----------------------------------------------------------------------------#
  
  #------------------ Identify Highly Discriminatory Clusters -------------------#
  
  # User inputs discrim_genes or discrim_clusts information
  # Find names of clusters that contain these discriminatory genes
  
  if(!is.null(discrim_clusts)){
    
    discrim_clusters = discrim_clusts
    
  }else if(!is.null(discrim_genes)){
    
    all_clusters = names(concat_geneMod)
    idx_clust_tmp = numeric(length(all_clusters))
    
    
    for(clust in all_clusters){
      clust_genes = unique(unlist(str_split(concat_geneMod[[clust]]$info$gene, pattern = ":")))
      if(any(clust_genes %in% discrim_genes)){
        idx_clust_tmp[which(all_clusters == clust)] = 1
      }
    }
    
    idx_clust = which(idx_clust_tmp==1)
    discrim_clusters = unique(all_clusters[idx_clust])
    
  }
  
  cat("Length discrim_clusters: ", length(discrim_clusters), "\n")
  
  concat_geneMod = concat_geneMod[discrim_clusters]
  
  return(concat_geneMod)
  
}

