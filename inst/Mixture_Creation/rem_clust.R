#-------------------------------------------------------------------------------------------#
# REDUCING TRANSCRIPT CLUSTERS                                                              #
#-------------------------------------------------------------------------------------------#
# NAME:                                                                                     #
#   remove_tclust_take2.R                                                                   #
# DATE:                                                                                     #
#   12/30/2013                                                                              #
# VERSION:                                                                                  #
#   R 2.15.3                                                                                #
# PROGRAMMER:                                                                               #
#   Douglas Roy Wilson, Jr.                                                                 #
#-------------------------------------------------------------------------------------------#
# DESCRIPTION:                                                                              #
#   After several failed attempts to make the previous optimization framework function as   #
#   expected, the entire program suite is being rewriten. To this end, the cluster removal  #
#   function is also being rewritten.                                                       #
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
# FUNCTION: rem_clust()                                                                     #
#-------------------------------------------------------------------------------------------#

rem_clust<-function(geneMod,co,min_ind){
  indices = which(names(geneMod)!="Sample_Info")
  geneMod_labels = paste("y_",geneMod[["Sample_Info"]][["info"]]$Label,sep="")
  co_ct = c()
  
  for(i in indices){
    tot_all = rep(0,length(geneMod_labels))
    for(j in 1:length(geneMod_labels)){
      tot_all[j] = sum(geneMod[[i]][[geneMod_labels[j]]])
    }
    if(min_ind==1 && min(tot_all)<co){co_ct=c(co_ct,i)}
    if(min_ind==0 && median(tot_all)<co){co_ct=c(co_ct,i)}
    if(all(dim(as.matrix(geneMod[[i]][["X"]]))==c(1,1))==1){co_ct=c(co_ct,i)}
    if(any(dim(as.matrix(geneMod[[i]][["X"]]))==0)){co_ct = c(co_ct,i)}
  }
  co_ct = unique(co_ct)
  
  if(length(co_ct)!=0){geneMod = geneMod[-co_ct]}
  
  return(geneMod)
}
