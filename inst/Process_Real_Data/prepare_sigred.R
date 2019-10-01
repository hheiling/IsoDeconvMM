rem_clust<-function(geneMod,co,min_ind){
  indices = which(names(geneMod)!="Sample_Info")
  geneMod_labels = paste("y_",geneMod[["Sample_Info"]]$Label,sep="")
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


setwd("/netscr/drwilson/2018-04-05 Paper 1/MappedData/Real_Data_Opt/GeneModels/")

files = list.files(pattern = "unred")
file_labels = gsub(pattern = "unred",replacement = "sigred",x=files)

geneslist = load(file="../../../Programs/EnsemblIds2Use.RData")
analy_genes = finalOut2$Ensembl.ID

for(j in 1:length(files)){
  load(files[j])

  indices2chk = which(names(fin_geneMod)!="Sample_Info")
  indices_tmp = NULL
  indices=NULL
  indices_tmp = rep(0,length(geneMod))
  for(i in indices2chk){
    infodf = fin_geneMod[[i]]$info
    genesi = unique(infodf$gene)
    genesi = unique(unlist(strsplit(x=genesi,split = ":")))
    if(any(genesi %in% analy_genes)){indices_tmp[i]=1}
  }
  indices = which(indices_tmp==1)
  
  sig_geneMod = fin_geneMod[indices]
  sig_geneMod["Sample_Info"] = fin_geneMod["Sample_Info"]
  
  sig_geneMod = rem_clust(geneMod = sig_geneMod,co = 5,min_ind = 0)
  
  message("File: ",file_labels[j]," // nGenes= ",length(sig_geneMod))
  
  save(sig_geneMod,file=file_labels[j])
}
  

  
