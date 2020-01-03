#-----------------------------------------------------------------------#
# SUMMARIZING the Output                                                #
#-----------------------------------------------------------------------#
# Compile the estimates from each cluster, provide a histogram, 
# and summarize the values with a geometric median.

#' @importFrom ICSNP spatial.median
#' @import ggplot2
#' @export
Summarize_Report<-function(final_output, restrict = FALSE, clusters = NULL, bin_num = 20){
  
  mix_names = names(final_output)
  
  if(restrict==TRUE && length(clusters)>0){
    warn("Set restrict to FALSE, but specified clusters for restriction.\n")
    warn("Will restrict to provided clusters.\n")
  }
  
  final_summary = list()
  
  for(j in 1:length(final_output)){
    cdata = final_output[[j]]
    
    #--------------------------------------------------------------#
    # CODING for restrictions                                      #
    #--------------------------------------------------------------#
    if(length(clusters)>0){
      cdata_fin = cdata[clusters]
    } else {
      # Eliminate clusters with "WARN" issues
      cdata_fin = cdata
    }
    
    #--------------------------------------------------------------#
    # EXTRACT cell-types info                                      #
    #--------------------------------------------------------------#
    cto = cdata_fin[[1]][["CellType_Order"]]
    nclust = length(cdata_fin)
    clust_names = names(cdata_fin)
    
    p_mat = matrix(0, nrow = nclust, ncol = length(cto))
    colnames(p_mat) = cto
    rownames(p_mat) = clust_names
    
    for(i in 1:nclust){
      p_mat[i,] = cdata_fin[[i]][["mix"]][["p.est"]]
    }
    
    if(any(p_mat< -1e-5)){message("Warning: Prop < -1e-5")}
    if(any(p_mat>(1+1e-5))){message("Warning: Prop > 1 + 1e-5")}
    
    p_mat[which(p_mat<0)] = 0
    # p_mat[which(p_mat>1)] = 1
    
    #--------------------------------------------------------------#
    # PLOTTING the cell type info                                  #
    #--------------------------------------------------------------#
    
    df = as.data.frame(p_mat)
    
    histograms = lapply(1:2, function(k){
      mtitle = sprintf("Distribution of Estimated Proportions of (%s)",cto[k])
      xlabel= sprintf("Prop. of (%s)",cto[k])
      hg = ggplot(data = df, mapping = aes(df[,k])) + geom_histogram(bins = bin_num) +
        ggtitle(label = mtitle, subtitle = sprintf("Mixture %s", mix_names[j])) + xlab(xlabel)
    })
    
    names(histograms) = cto
    
    #---------------------------------------------------------------#
    # Summarizing the Values                                        #
    #---------------------------------------------------------------#
    fin_est = matrix(0,nrow=1,ncol=length(cto))
    colnames(fin_est) = cto
    
    p.fin = spatial.median(X = p_mat[,-length(cto)])
    
    fin_est[1,] = c(p.fin,1-sum(p.fin))
    
    print(fin_est)
    print(head(p_mat))
    print(class(histograms))
    print(length(histograms))
    print(class(histograms$set1))
    
    final_summary[[mix_names[j]]] = list(p_est = fin_est, p_mat = p_mat, histograms = histograms)
  }
  
  
  #----------------------------------------------------------------#
  # RETURN values                                                  #
  #----------------------------------------------------------------#
  return(final_summary)
}