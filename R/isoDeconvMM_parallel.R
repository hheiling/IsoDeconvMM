# parallelization of IsoDeconvMM

IsoDeconvMM_parallel = function(directory = NULL, mix_files, pure_ref_files,
                       fraglens_files, bedFile, knownIsoforms, 
                       discrim_genes = NULL, discrim_clusts = NULL,
                       readLen, lmax = 600, eLenMin = 1, mix_names = NULL,
                       initPts = NULL, max_cores = 1, cluster_type = c("Fork","Socket"),
                       optim_options = optimControl(), ClustInfoOnly = FALSE){
  
  # time0 = proc.time()
  
  if(is.null(discrim_genes) & is.null(discrim_clusts)){
    stop("Either discrim_genes or discrim_clusts must be specified \n")
  }else if(!is.null(discrim_genes) & !is.null(discrim_clusts)){
    warning("Both discrim_genes and discrim_clusts specified, ",
            "will only use discrim_clusts \n", immediate. = T)
  }
  
  if(length(cluster_type) > 1) cluster_type = cluster_type[1]
  if(!(cluster_type %in% c("Fork","Socket"))){
    stop("cluster_type must by either 'Fork' or 'Socket'")
  }
  
  #-----------------------------------------------------------------------------#
  # Download all needed files, convert input items to useful formats            #
  #-----------------------------------------------------------------------------#
  
  if(!is.null(directory)){
    setwd(directory)
  }
  
  countData_mix = mix_files
  countData_pure = pure_ref_files[,1] 
  
  cellTypes_pure = as.character(pure_ref_files[,2])
  # ctpure_names = pure cellTypes names
  ctpure_names = levels(as.factor(cellTypes_pure))
  
  labels_pure = character(length(countData_pure))
  
  for(type in ctpure_names){
    locations = which(cellTypes_pure == type)
    labels_type = cellTypes_pure[locations]
    num_ref = length(labels_type)
    labels_ref = str_c(labels_type, "_ref", 1:num_ref)
    labels_pure[locations] = labels_ref
  }
  
  if(!is.null(mix_names)){
    labels_mix = mix_names
  }else{
    labels_mix = str_c("mix",1:length(mix_files))
  }
  
  # Download all count text files, compute total counts for each file
  # Output: list with elements total_cts, counts_list
  ## See comp_total_cts() function under "Internal isoDeconvMM Functions" heading later in this document
  pure_input = comp_total_cts(directory = directory, countData = countData_pure)
  pure_counts = pure_input$counts_list
  
  mix_input = comp_total_cts(directory = directory, countData = countData_mix)
  mix_counts = mix_input$counts_list
  
  countData = pure_counts
  countData[(length(pure_counts)+1):(length(pure_counts)+length(mix_counts))] = mix_counts
  
  if(!(length(fraglens_files) %in% c(1, length(mix_files)))){
    stop("Length of fraglens_files must be equal to 1 or length of mix_files")
  }
  
  fraglens_list = list()
  for(i in 1:length(fraglens_files)){
    fraglens = read.table(fraglens_files[i], as.is = T)
    if (ncol(fraglens) != 2) {
      stop(fraglens_files[i], " should have 2 columns for Freq and Len\n")
    }
    names(fraglens) = c("Freq", "Len")
    
    fraglens_list[[i]] = fraglens
  }
  
  
  # CHECKING Presence of Isoforms File:
  # IsoDeconv requires a list of known isoforms in order to model intra-sample heterogeneity. Stops program if file not present.
  # Load knownIsoforms .RData object:
  
  if(!is.null(knownIsoforms)){
    assign("isoAll", get(load(knownIsoforms)))
  }else{
    stop("knownIsoforms list object must be present!")
  }
  
  # Download .bed file
  bedFile_info = read.table(sprintf("%s", bedFile), sep = "\t", as.is = TRUE)
  
  bf_colNames = c("chr", "start", "end", "exon", "score", "strand")
  
  if (ncol(bedFile_info) != 6) {
    cN = paste(bf_colNames, collapse = ", ")
    stop(bedFile, " should have 6 columns: ", cN, "\n")
  }
  
  names(bedFile_info) = bf_colNames
  
  print("Finished loading supporting files")
  # time1 = proc.time()
  
  # Check initPts is specified correctly OR create initPts matrix if not specified
  
  if(is.null(initPts)){
    initPts = matrix(1/length(ctpure_names), nrow = 1, ncol = length(ctpure_names))
    colnames(initPts) = ctpure_names
  }else if(class(initPts) == "matrix"){
    if(ncol(initPts) != length(ctpure_names)){
      stop("number of columns of initPts must be equal to the number of pure reference cell types")
    }
    if(!all(colnames(initPts) %in% ctpure_names)){
      stop("colnames of initPts must match pure cell type names in pure_ref_files; see documentation for more details")
    }
    if(sum(initPts > 1) > 0 | sum(initPts < 0) > 0 | (nrow(initPts) > 1 & sum(rowSums(initPts) > 1) > 0)){
      stop("initPts specifies probabilities, which cannot be below 0 or above 1 or sum to more than 1")
    }
    
    # Re-arrange column order of initPts matrix to match order of ctpure_names
    initPts = initPts[,ctpure_names]
    if(!is.matrix(initPts)){
      # If only one initial point, convert above result from vector to matrix with one row
      initPts = matrix(initPts, nrow = 1)
    }
    colnames(initPts) = ctpure_names
    
  }else{
    stop("initPts must be a matrix, see documentation for details")
  }
  
  #--------------------------------------------------------------------------------#
  # Step 1
  # concat_geneMod() extracts gene info, counts, and exon information
  # geneModel_creatoin() calls geneModel_multcell_Edit(), which is an edit of Dr. Wei Sun's
  # geneModel creation problem (after edits, now accommodates multiple cell types)
  # geneModel() altered by:                                                                                                                        
  #    Douglas Roy Wilson, Jr. 
  #--------------------------------------------------------------------------------#
  
  labels = c(labels_pure, labels_mix)
  cellTypes = c(cellTypes_pure, rep("mix", times = length(mix_files)))
  
  concat_geneMod = cluster_info_parallel(countData = countData, labels = labels, cellTypes = cellTypes, 
                                bedFile = bedFile_info, discrim_genes = discrim_genes,
                                discrim_clusts = discrim_clusts, readLen = readLen, lmax = lmax,
                                max_cores = max_cores, cluster_type = cluster_type)
  
  final_geneMod = list()
  
  for(j in 1:length(mix_files)){
    cellTypes_sub = c(cellTypes_pure, "mix")
    labels_sub = c(labels_pure, labels_mix[j])
    total_cts = c(pure_input$total_cts, mix_input$total_cts[j])
    if(length(fraglens_list) == 1){
      fragSizeFile = fraglens_list[[1]]
    }else if(length(fraglens_list) == length(mix_files)){
      fragSizeFile = fraglens_list[[j]]
    }
    
    # Need to select specific components of concat_geneMod?
    
    fin_geneMod = geneModel_creation(concat_geneMod = concat_geneMod, fragSizeFile = fragSizeFile, 
                                     cellTypes = cellTypes_sub, labels = labels_sub,
                                     knownIsoforms = isoAll, readLen = readLen, lmax = lmax, 
                                     eLenMin = eLenMin, total_cts = total_cts)
    
    # Perform some checks on the geneMod output
    ## See R/rem_clust.R for rem_clust() code
    sig_geneMod = rem_clust(geneMod = fin_geneMod,co = 5,min_ind = 0)
    
    final_geneMod[[j]] = sig_geneMod
    
  }
  
  print("Finished creation of gene model")
  # time2 = proc.time()
  
  #--------------------------------------------------------------------------------#
  # Step 2
  # Add rds_exons matrix for each cell type to each cluster element
  # rds_exons matrix: columns = samples associated with each cell type,
  # rows (except first) = read count for each exon set in the given gene/cluster 
  # for sample j of cell type k,
  # first row = total read count outside gene/cluster of interest in sample j of 
  # cell type k
  # EDIT TO GROUP CELL TYPES
  #--------------------------------------------------------------------------------#
  
  ## See mod_sig_gM() function under "Internal isoDeconvMM Functions" heading later in this document
  modified_sig_geneMod = mod_sig_gM(significant_geneMod = final_geneMod)
  
  #-----------------------------------------------------------#
  # Step 3
  # Pure Cell Type Parameter Estimation                       
  #-----------------------------------------------------------#
  
  ## See pure_estimation() function under "Internal isoDeconvMM Functions" heading later in this document
  pure_est = pure_estimation(modified_sig_geneMod = modified_sig_geneMod, cellTypes = ctpure_names)
  
  print("Finished pure cell type parameter estimation")
  # time3 = proc.time()
  
  if(is.null(mix_names)){
    names(pure_est) = str_remove(mix_files, ".txt")
  }else{
    names(pure_est) = mix_names
  }
  
  cluster_params = structure(list(Cluster_Details = pure_est, cellTypes = ctpure_names),
                             class = c("clusterInfo"))
  
  if(ClustInfoOnly){
    return(cluster_params)
  } 
  #-------------------------------------------------------------------#
  # Testing purposes only: Simulation    
  # If only have one reference sample for each pure cell type,
  # simulate additional sample counts 
  # Use pure cell type parameter estimates from Step 3 to simulate
  # additional counts
  #-------------------------------------------------------------------#
  
  # Step 4
  #------------------------------------------------------------------------------#
  # Step 4
  # Mixture Cell Type Parameter Estimation
  # Calls STG.Update_Cluster.All() for this estimation procedure
  #------------------------------------------------------------------------------#
  
  IsoDeconv_Output = mix_fit(cluster_params = cluster_params, initPts = initPts, 
                             optim_options = optim_options, max_cores = max_cores)
  
  print("Finished mixture paramter estimation for all samples")
  
  return(IsoDeconv_Output)
  
  # for(i in 1:length(pure_est)){
  #   
  #   tmp.data = pure_est[[i]]
  #   
  #   # Establish input break ups
  #   
  #   # Data Set Necessities:
  #   clust.start = 1
  #   clust.end = length(tmp.data)
  #   by.value = 15
  #   
  #   start.pts = seq(from = 1,to = clust.end,by = by.value)
  #   end.pts = c((start.pts[-1]-1),clust.end)
  #   
  #   cluster_output = list()
  #   for(m in 1:length(start.pts)){
  #     start.pt = start.pts[m]
  #     end.pt = end.pts[m]
  #     
  #     curr.clust.opt = tmp.data[c(start.pt:end.pt)]
  #     ## See R/Production_Functions_MixedSamp.R for STG.Updat_Cluster.All() code
  #     curr.clust.out = STG.Update_Cluster.All(all_data=curr.clust.opt, cellTypes = ctpure_names,
  #                                             optimType=optim_options$optimType, 
  #                                             simple.Init=optim_options$simple.Init, 
  #                                             initPts = initPts)
  #     
  #     cluster_output[[m]] = curr.clust.out
  #   }
  #   
  #   IsoDeconv_Output[[i]] = cluster_output
  #   
  #   cat("Finished mixture param estimation for mix file ", i, "\n")
  # }
  
  # print("Finished mixture paramter estimation for all samples")
  # time4 = proc.time()
  
  #---------------------------------------------------------------------------------------------#
  # Step 5
  # Re-compile Step 4 output such that all information organized as follows:
  # First layer of Final_Compiled_Output list is associated with a mixture file of 
  # Second layer of list is associated with name of a cluster
  # Third layer of list contains all cluster-specific information
  #---------------------------------------------------------------------------------------------#
  
  # Final_Compiled_Output = list()
  # 
  # for(j in 1:length(IsoDeconv_Output)){
  #   
  #   comp.out= NULL
  #   curr.clust.out = NULL
  #   
  #   #---- Set up new pattern ----#
  #   est.chunks = IsoDeconv_Output[[j]]
  #   
  #   message("File ", j)
  #   
  #   #---- Populate Output Dataset ----#
  #   comp.out = list()
  #   
  #   for(i in 1:length(est.chunks)){
  #     curr.clust.out = est.chunks[[i]]
  #     clust_names = names(curr.clust.out)
  #     nl = length(curr.clust.out)
  #     for(m in 1:nl){
  #       comp.out[[clust_names[m]]]=curr.clust.out[[m]]
  #     }
  #   }
  #   
  #   Final_Compiled_Output[[j]] = comp.out 
  #   
  # }
  # 
  # if(is.null(mix_names)){
  #   names(Final_Compiled_Output) = str_remove(mix_files, ".txt")
  # }else{
  #   names(Final_Compiled_Output) = mix_names
  # }
  # 
  # # Output timings
  # # time_record = rbind(time0[1:3],time1[1:3],time2[1:3],time3[1:3],time4[1:3])
  # # rownames(time_record) = c("Start","End Loading","End Gene Model","End Pure Fit","End Mixture Fit")
  # # 
  # # Final_Compiled_Output[["Time"]] = time_record
  # 
  # return(Final_Compiled_Output)
  
  
} # End isoDeconvMM() function
