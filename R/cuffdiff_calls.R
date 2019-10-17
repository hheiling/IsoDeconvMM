
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
