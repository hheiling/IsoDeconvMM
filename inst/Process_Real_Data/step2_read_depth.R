# -----------------------------------------------------------------
# read in sample information
# -----------------------------------------------------------------
args<-commandArgs(TRUE)

if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  inputDir = "/netscr/drwilson/2018-04-05 Paper 1/MappedData/CD4/"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


# -----------------------------------------------------------------
# read in count data
# -----------------------------------------------------------------

setwd("./_count")

files = list.files(pattern="count")
sams  = gsub("count_", "", files, fixed=TRUE)
sams  = gsub(".txt",  "", sams,  fixed=TRUE)
cts   = matrix(NA, nrow=length(sams), ncol=2)

for(i in 1:length(sams)){

  sam1 = sams[i]
  ff1  = files[i]

  ct   = scan(ff1, quiet=TRUE)
  cts[i,] = ct
}

cts = cts/1e6

summary((cts[,1] - cts[,2])/cts[,1] )

# -----------------------------------------------------------------
# summarize
# -----------------------------------------------------------------

setwd("../")

#Create a PDF plot with a certain height and width.
pdf("TReC_vs_QCed_TReC.pdf", width=3, height=3)

#Set margins of the figure (in inches).
par(mar=c(5,4,1,1))

#Creates the plot and places it in the PDF. Plot, here,
# places uncorrected counts as X and corrected counts as Y.

plot(cts[,1], cts[,2], xlab="Total # of Reads (million)", 
  ylab="Total # of Reads after QC", bty="n", cex=0.5, pch=20, col="darkgrey")
points(cts[,1], cts[,2], col="darkgreen", cex=0.5)

b = median(cts[,2])/median(cts[,1])
abline(0, b, col="darkred", lwd=2)
abline(0, 0.70, col="orange", lwd=2)

table = (cts[,2]/cts[,1]<.7)

lg = c(sprintf("y=%.2fx", b), "y=0.70x")
legend("topleft", lg, lty=c(1,1), col=c("darkred", "orange"), bty="n")
dev.off()

db = data.frame(sample=sams, cts, round(cts[,2]/cts[,1],4))

names(db) = c("sample", "TReC", "TReC_after_QC", "ratio")
dim(db)
db[1:2,]

setwd("./")

write.table(db, file = "samples_count.txt", append = FALSE, quote = FALSE, 
  sep = "\t", row.names = FALSE, col.names = TRUE)


