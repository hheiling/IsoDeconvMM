
#' @title Mixture Proportion Estimation Algorithm
#'  
#' @description Calculates the proportions of pure cell type components in heterogeneous cell type 
#' samples using the S4 object with class "clusterInfo" created in \code{\link{IsoDeconvMM}} when 
#' the argument \code{ClustInfoOnly} is set to \code{TRUE}.
#' 
#' @param clusterInfo list object of class "clusterInfo" created in  \code{\link{IsoDeconvMM}} when 
#' the argument \code{ClustInfoOnly} is set to \code{TRUE}.
#' @param mix_names character vector of 
#'
#' @import doParallel
#' @importFrom pbapply pboptions pbapply
#' @export
mix_fit = function(clusterInfo, mix_names = NULL, initPts = NULL, optim_options = optimControl(),
                   max_cores = 1, cluster_type = c("Fork","Socket")){
  
  if(length(cluster_type) > 1) cluster_type = cluster_type[1]
  if(!(cluster_type %in% c("Fork","Socket"))){
    stop("cluster_type must be either 'Fork' or 'Socket'")
  }
  
  if(is.null(mix_names)){
    pure_est = clusterInfo$Cluster_Details
  }else{
    tmp_all = clusterInfo$Cluster_Details
    pure_est = tmp_all[which(names(tmp_all) %in% mix_names)]
  }
  
  cellTypes = clusterInfo$cellTypes
  
  # Check initPts is specified correctly OR create initPts matrix if not specified
  if(is.null(initPts)){
    initPts = matrix(1/length(cellTypes), nrow = 1, ncol = length(cellTypes))
    colnames(initPts) = cellTypes
  }else if(class(initPts) == "matrix"){
    if(ncol(initPts) != length(cellTypes)){
      stop("number of columns of initPts must be equal to the number of pure reference cell types")
    }
    if(!all(colnames(initPts) %in% cellTypes)){
      stop("colnames of initPts must match pure cell type names in pure_ref_files; see documentation for more details")
    }
    if(sum(initPts > 1) > 0 | sum(initPts < 0) > 0 | (nrow(initPts) > 1 & sum(rowSums(initPts) > 1) > 0)){
      stop("initPts specifies probabilities, which cannot be below 0 or above 1 or sum to more than 1")
    }
    
    # Re-arrange column order of initPts matrix to match order of cellTypes
    initPts = initPts[,cellTypes]
    if(!is.matrix(initPts)){
      # If only one initial point, convert above result from vector to matrix with one row
      initPts = matrix(initPts, nrow = 1)
    }
    colnames(initPts) = cellTypes
    
  }else{
    stop("initPts must be a matrix, see documentation for details")
  }
  
  # Step 4
  #------------------------------------------------------------------------------#
  # Step 4
  # Mixture Cell Type Parameter Estimation
  # Calls STG.Update_Cluster.All() for this estimation procedure
  #------------------------------------------------------------------------------#
  
  mixFitFun = function(clust_idx, tmp.data, cellTypes, optim_options, initPts){
    
    curr.clust.opt = tmp.data[clust_idx]
    
    curr.clust.out = STG.Update_Cluster.All(all_data=curr.clust.opt, cellTypes = cellTypes,
                                            optimType=optim_options$optimType, 
                                            simple.Init=optim_options$simple.Init, 
                                            initPts=initPts)
    
    return(curr.clust.out)
    
  }
  
  # If specify only some mixture sample names, restrict pure_est object
  if(!is.null(mix_names)){
    pure_est = pure_est[which(names(pure_est) %in% mix_names)]
  }
  
  IsoDeconv_Output = list()
  
  for(i in 1:length(pure_est)){
    
    tmp.data = pure_est[[i]]
    
    if(max_cores > length(tmp.data)){
      num_cores = length(tmp.data)
    }else{
      num_cores = max_cores
    }
    
    clust_names = names(tmp.data)
    
    # Establish input break ups
    size = ceiling(seq_along(length(tmp.data)/num_cores))
    full_idx = 1:length(tmp.data)
    chunks = split(full_idx, ceiling(seq_along(full_idx)/size))
    
    cl = NULL # No parallelization used if cl stays as NULL (num_cores == 1)
    
    if((cluster_type == "Socket") & (num_cores > 1)){
      cl = parallel::makeCluster(num_cores)
      registerDoParallel(cl)
      vars = c("tmp.data","cellTypes","optim_options","initPts","STG.Update_Cluster.All",
               "STG.Update_Cluster.Single","STG.Update_Cluster.SingMI","applyMI",
               "update.iso.STG2","lab_fit_update_v2","compute.pmeans","applyLik",
               "ActLik4test","mix.full.stage.lik.pt")
      clusterExport(cl, varlist = vars, envir = environment())
    }else if((cluster_type == "Fork") & (num_cores > 1)){
      cl = parallel::makeForkCluster(num_cores, outFile = "")
      registerDoParalle(cl)
    }
    
    pboptions(type="timer")
    
    cluster_output = pblapply(chunks, FUN = mixFitFun, tmp.data = tmp.data,
                             cellTypes = cellTypes, optim_options = optim_options,
                             initPts = initPts, cl = cl)
    if(num_cores > 1){
      parallel::stopCluster(cl)
    }
    
    IsoDeconv_Output[[i]] = cluster_output
    
    cat("Finished mixture param estimation for mix file ", i, "\n")
  }
  
  print("Finished mixture paramter estimation for all samples")
  
  return(IsoDeconv_Output)
  
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
  # return(Final_Compiled_Output)
  
  
} # End mix_fit function