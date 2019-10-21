
#' @import cummeRbund 
#' @export
cuffdiff_siggenes<-function(folder, directory){
  
  #-------------------------------------------------------------------#
  # READ-In Data                                                      #
  #-------------------------------------------------------------------#
  setwd(directory)
  cuffdiff_total<- readCufflinks(folder,rebuild = T)
  
  #-------------------------------------------------------------------#
  # Consider Sig Genes                                                #
  #-------------------------------------------------------------------#
  gene_diff<- diffData(genes(cuffdiff_total))
  gene_diff_sub <- subset(gene_diff,(status=="OK"))
  nrow(gene_diff_sub)
  sig_gene_diff <- subset(gene_diff_sub,(significant=='yes'))
  sig_gene_diff <- sig_gene_diff[order(sig_gene_diff$q_value),]
  sig_gene_diff$class = rep("Gene",nrow(sig_gene_diff))
  nrow(sig_gene_diff)
  
  #-------------------------------------------------------------------#
  # Consider Isoform Genes                                            #
  #-------------------------------------------------------------------#
  prom_diff<- distValues(promoters(cuffdiff_total))
  prom_diff_sub <- subset(prom_diff,(status=="OK"))
  nrow(prom_diff_sub)
  prom_diff_sub = prom_diff_sub[order(prom_diff_sub$p_value),]
  prom_diff_sub = prom_diff_sub[which(prom_diff_sub$p_value<0.10),]
  prom_diff_sub$class = rep("Iso",nrow(prom_diff_sub))
  
  splice_diff<- distValues(splicing(cuffdiff_total))
  splice_diff_sub <- subset(splice_diff,(status=="OK"))
  nrow(splice_diff_sub)
  splice_diff_sub = splice_diff_sub[order(splice_diff_sub$p_value),]
  splice_diff_sub = splice_diff_sub[which(splice_diff_sub$p_value<0.10),]
  splice_diff_sub$class = rep("Iso",nrow(splice_diff_sub))
  
  relCDS_diff<- distValues(relCDS(cuffdiff_total))
  relCDS_diff_sub <- subset(relCDS_diff,(status=="OK"))
  nrow(relCDS_diff_sub)
  relCDS_diff_sub = relCDS_diff_sub[order(relCDS_diff_sub$p_value),]
  relCDS_diff_sub = relCDS_diff_sub[which(relCDS_diff_sub$p_value<0.10),]
  relCDS_diff_sub$class = rep("Iso",nrow(relCDS_diff_sub))
  
  #Obtain Gene Names for All significant locations:
  totGeneMat = rbind(sig_gene_diff[1:240,c("gene_id","p_value","class")],
                     prom_diff_sub[1:240,c("gene_id","p_value","class")],
                     splice_diff_sub[1:240,c("gene_id","p_value","class")],
                     relCDS_diff_sub[1:135,c("gene_id","p_value","class")])
  gene_names<-annotation(genes(cuffdiff_total))
  gene_names_sub <- gene_names[,c("gene_id","gene_short_name")]
  testing<-merge(totGeneMat,gene_names_sub,by=c("gene_id"),all.x=TRUE)
  testing<-testing[which(!duplicated(testing$gene_id)),]
  
  return(testing)
}

#' @import cummeRbund
#' @export
EnsemblIds2Use = function(folder, directory){
  
  testing = cuffdiff_siggenes(folder = folder, directory)
  
  #------------------ LOAD AND PROCESS GENE NAMES ---------#
  queryList = testing$gene_short_name
  ensIds    = queryMany(queryList,scopes="symbol",fields=c("ensembl.gene"),
                        species="human",return.as = "records")
  
  ensLink = matrix(NA,ncol=2,nrow=length(ensIds))
  
  for(i in 1:length(ensIds)){
    if(length(ensIds[[i]]$ensembl)==0){
      ensLink[i,1] = ensIds[[i]]$query
    } else {
      if(length(ensIds[[i]]$ensembl)==1){
        ensLink[i,1] = ensIds[[i]]$query
        ensLink[i,2] = ensIds[[i]]$ensembl$gene
      } else {
        ensLink[i,1] = ensIds[[i]]$query
        ensLink[i,2] = ensIds[[i]]$ensembl[[1]]$gene
      }
    }
  }
  
  ensLink = as.data.frame(ensLink)
  colnames(ensLink) = c("gene_short_name","Ensembl.ID")
  
  finalOut = merge(testing,ensLink,by="gene_short_name",all.x=TRUE)
  finalOut2 = finalOut[which(!is.na(finalOut$Ensembl.ID)),]
  
  # Final Output
  finalOut2 = finalOut2[order(finalOut2$class,finalOut2$p_value),]
  
  return(finalOut2)
}
