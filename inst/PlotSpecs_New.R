#-------------------------------------------------------------------#
# Summarizing Simulations                                           #
#-------------------------------------------------------------------#
library(ICSNP)
#---------- Programs ------------#
histogram_func<-function(vec){
  mtitle = sprintf("Distribution of Estimated Proportions of (%s)",vec[1])
  xlabel= sprintf("Prop. of (%s)",vec[1])
  hist(x = as.numeric(vec[-1]),main = mtitle,xlab = xlabel)
}

Summarize_Report<-function(cdata,restrict,clusters=NULL){
  #--------------------------------------------------------------#
  # CODING for restrictions                                      #
  #--------------------------------------------------------------#
  if(length(clusters)>0){
    if(restrict==0 && length(clusters)>0){
      warn("Set Restrict to no, but specified clusters for restriction.\n")
      warn("Will restrict to provided clusters.\n")
    }
    cdata_fin = cdata[clusters]
  } else {
    cdata_fin = cdata
  }
  
  #--------------------------------------------------------------#
  # EXTRACT cell-types info                                      #
  #--------------------------------------------------------------#
  cto = cdata_fin[[1]][["CellType_Order"]]
  nclust = length(cdata_fin)
  
  p_mat = matrix(0,nrow = nclust,ncol = length(cto))
  colnames(p_mat) = cto
  
  for(i in 1:nclust){
    p_mat[i,] = cdata_fin[[i]][["mix"]][["p.est"]]
  }
  
  if(any(p_mat< -1e-5)){message("Warning: Prop < -1e-5")}
  if(any(p_mat>(1+1e-5))){message("Warning: Prop > 1 + 1e-5")}
  
  p_mat[which(p_mat<0)] = 0
  
  #--------------------------------------------------------------#
  # PLOTTING the cell type info                                  #
  #--------------------------------------------------------------#
  df = as.data.frame(x = t(p_mat))
  names_vec = paste("clust",1:nclust,"est",sep = ".")
  colnames(df) = names_vec
  df$cellType = cto
  df = df[c("cellType",names_vec)]
  
  apply(X = df,MARGIN = 1,FUN = histogram_func)
  
  #---------------------------------------------------------------#
  # Summarizing the Values                                        #
  #---------------------------------------------------------------#
  fin_est = matrix(0,nrow=1,ncol=length(cto))
  colnames(fin_est) = cto
  
  p.fin = spatial.median(X = p_mat[,-length(cto)])
  
  fin_est[1,] = c(p.fin,1-sum(p.fin))
  
  #----------------------------------------------------------------#
  # RETURN values                                                  #
  #----------------------------------------------------------------#
  return(list(p.est = fin_est,p.mat = p_mat))
}

#----- Load the Genes ------#
test1 = load("D://DougData/Documents/Dissertation/Paper 1 - Cell Type Abundance Estimation/Daily Work/1_30_2018/Mixture_Creation/EnsemblIds2Use.RData")

isoGenes  = unique(finalOut2$Ensembl.ID[which(finalOut2$class=="Iso")])
GeneGenes = unique(finalOut2$Ensembl.ID[which(finalOut2$class=="Gene")]) 

#----- Load it Up with Simulation Data --------#
DirName = sprintf("D:/DougData/Documents/Dissertation/Paper 1 - Cell Type Abundance Estimation/Daily Work/1_30_2018/Results/Real_Data/")
setwd(DirName)

p = seq(from=0,to=100,by=10)

Plot_Dat = data.frame(CD4_True=p)
names(Plot_Dat) = c("True.CD4.Prop")

no.combo = nrow(Plot_Dat)

p.A1   = p.A2 = p.A3  = nClust = bClust = lClust = rep(-1,no.combo)

for(j in p){
  jp = 100-j
  File_Name = sprintf("FinCompiled_CD4_%s_CD8_%s.RData",j,jp)
  load(File_Name)
  
  #------ Consider Missing Data -------#
  message("P = ",j)
  message("No. Clusters: ",length(comp.out))
  
  #------ Consider Warned Clusters -----#
  warn.clust = rep(0,length(comp.out))
  for(z in 1:length(comp.out)){
    warn.clust[z] = comp.out[[z]][["WARN"]]
  }
  table(warn.clust)
  
  Bad_Clust = sum(as.integer((warn.clust==4)|(warn.clust==5)))
  lim_Clust = sum(as.integer((warn.clust==1)))
  
  #------ Summarize the Output ------#
  indices_tmp = NULL
  indicesIso=NULL
  indices_tmp = rep(0,length(comp.out))
  for(i in 1:length(comp.out)){
    infodf = comp.out[[i]]$info
    genesi = unique(infodf$gene)
    genesi = unique(unlist(strsplit(x=genesi,split = ":")))
    if(any(genesi %in% isoGenes)){indices_tmp[i]=1}
  }
  indicesIso = which(indices_tmp==1)
  
  indices_tmp = NULL
  indicesGenes=NULL
  indices_tmp = rep(0,length(comp.out))
  for(i in 1:length(comp.out)){
    infodf = comp.out[[i]]$info
    genesi = unique(infodf$gene)
    genesi = unique(unlist(strsplit(x=genesi,split = ":")))
    if(any(genesi %in% GeneGenes)){indices_tmp[i]=1}
  }
  indicesGenes = which(indices_tmp==1)
  
  clust2use = intersect(which((warn.clust!=4)&(warn.clust!=5)),indicesIso)
  p.A1t = Summarize_Report(cdata=comp.out,restrict = 1,clusters = clust2use)[["p.est"]][1]
  
  clust2use = intersect(which((warn.clust!=4)&(warn.clust!=5)),indicesGenes)
  p.A2t = Summarize_Report(cdata=comp.out,restrict = 1,clusters = clust2use)[["p.est"]][1]

  clust2use = which((warn.clust!=4)&(warn.clust!=5))
  p.A3t = Summarize_Report(cdata=comp.out,restrict = 1,clusters = clust2use)[["p.est"]][1]
  
  #------ Store the Summaries -------#
  idx2use = which((Plot_Dat$True.CD4.Prop==j))
  
  p.A1[idx2use]   = p.A1t
  p.A2[idx2use]   = p.A2t
  p.A3[idx2use]   = p.A3t
  nClust[idx2use] = length(comp.out)
  bClust[idx2use] = Bad_Clust
  lClust[idx2use] = lim_Clust
}

Plot_Dat$Est.Propv1 = p.A1
Plot_Dat$Est.Propv2 = p.A2
Plot_Dat$nClust     = nClust 
Plot_Dat$bClust     = bClust
Plot_Dat$lClust     = lClust

Plot_Dat = Plot_Dat[which(Plot_Dat$Avg_RPG!="Avg. RPG = 400"),]

#---- Generate the Lattice Plots using GGPLOT ----#
library(ggplot2)
AccPlots = ggplot(data=Plot_Dat,aes(x=True.Prop,y=Est.Propv1))+geom_point()+facet_grid(Avg_RPG~Var)+
  geom_abline(slope=1,intercept=0)+ggtitle("Simulation Results",subtitle=sprintf("%s Pure Samples per CT",nkper_i))+
  labs(x="True Prop CTA",y="Est. Prop CTA")
AccPlots