#!/usr/bin/Rscript
library("data.table")
library("xtable")

argv <- commandArgs(T)
dir <- argv[1]
project <- argv[2]
date <- argv[3]
reference <- argv[4]
threads <- argv[5]

cat(project,date,reference,"\n")

setwd(sprintf("%s/analysis",dir))
system("mkdir -p images/qc; for i in ../qc/images/*; do j=`echo $i | sed 's|../qc/images/|images/qc/|'`; cp $i $j; done")
system(sprintf("cat ../input/%s/fastq/raw/*.stats > fastq.tsv",date))
system(paste0("/usr/local/mark/rnaseq_align.sh ",reference))

#
#  Make PCA
#
count <- fread(paste0("../processed/cuffnorm.",reference,"/genes.count_table"))
d <- dist(t(count[,2:ncol(count),with=F]))
pca <- cmdscale(d,k=9,eig=T)
eig <- 100 * pca$eig / sum(pca$eig)
png("images/pca.png")
plot(pca$points[,1],pca$points[,2],pch=19,xlab=sprintf("PCA1: %0.1f%%",eig[1]),ylab=sprintf("PCA2: %0.1f%%",eig[2]))
text(pca$points[,1],pca$points[,2],rownames(pca$points),pos=4)
dev.off()

#
#  Start report
#
sink("report.tex")
cat("\\documentclass{article}\n\\usepackage[top=2cm,left=2cm,right=2cm,bottom=2cm]{geometry}\n\\usepackage{longtable}\n\\usepackage{graphicx}\n\\begin{document}\n")
cat("\\begin{center}{\\Huge ",project,":",date," Report}\\end{center}\n")

#
#  Read percentages
#
samples <- fread(sprintf("../input/samples.%s",date),header=F)
fastq <- fread("fastq.tsv")
fastq <- fastq[grepl("_R1_",V1)]
fastq[,sample:=gsub("raw/110_","",gsub("_L00._R._001.fastq.gz","",V1,perl=T))]
fastq[,lane:=paste("L",gsub("raw/.*_L00","",gsub("_R._001.fastq.gz","",V1,perl=T)),sep="")]
reads <- dcast(fastq,sample~lane,value.var="V2",sum)
reads[,Total:=L1+L2+L3+L4]
total <- sum(reads$Total)
reads[,Percent:=round(100*Total/total,2)]
cat("\\begin{center}{\\huge Read distribution}\\end{center}\n")
print(xtable(reads))

#
#  Images
#
for (i in 1:nrow(samples)) {
  s <- as.character(samples[i,1])
  cat("\\newpage\n")
  cat("\\begin{center}{\\huge",gsub("_"," ",gsub("170922_125_S","",s,perl=T)),"}\\end{center}\n")
  cat("\\begin{center}\n")
  for (l in 1:4) {
    for (r in 1:2) {
      cat("\\includegraphics[width=0.4\\textwidth]{images/qc/",s,"_L00",l,"_R",r,"}\n",sep="")
    }
  }
  cat("\\end{center}\n")
}

#
#  Align percentages
#
alignment <- c()
for (i in 1:nrow(samples)) {
  s <- as.character(samples[i,1])
#  r <- read.delim(sprintf("../processed/%s/%s/align_summary.txt",s,reference),sep=" ",fill=T,header=F)
#  cat ../processed/*/ucsc_!!!hg38!!!/align_summary.txt | sed 's/ //g' | sed 's/[\:\(\)\%]/\t/g'
#  if (!is.na(r[2,18])) {
#    alignment <- rbind(alignment, c(gsub("170922_125_","",s,perl=T), r[2,18], r[3,17], r[7,17], r[10,4], r[11,11], r[12,19]))
#  } else {
#    alignment <- rbind(alignment, c(gsub("170922_125_","",s,perl=T), r[2,19], r[3,18], r[7,18], r[10,5], r[11,11], r[12,20]))
#  }
  r <- read.delim(sprintf("align/%s",s),sep="\t",fill=T,header=F)
  alignment <- rbind(alignment, c(gsub("170922_125_","",s,perl=T), 
                     r[1,1], r[2,1], r[5,1], r[8,1], r[9,1], r[10,1]))
}
alignment <- data.table(alignment)
colnames(alignment) <- c("Sample","Reads","Left","Right","Paired","Multiple","Discordant")
alignment[,Paired:=as.numeric(Paired)]
alignment[,Discordant:=as.numeric(Discordant)]
alignment[,Reads:=as.numeric(Reads)]
alignment[,Rate:=round(100*(Paired-Discordant)/Reads,2)]
cat("\\newpage\n")
cat("\\begin{center}{\\huge Alignment Statistics}\\end{center}\n")
print(xtable(alignment,digits=c(0,0,0,0,0,0,0,0,2)))

#
#  Transcripts Found (Cuffmerged GTF)
#
#base <- as.numeric(system("cut -f2 -d' ' ../input/reference/!!!ucsc_hg38!!!.gtf  | sort | uniq -c | wc -l",intern=T))
#m <- read.delim("../processed/cuffmerge.!!!ucsc_hg38!!!/merged.gtf",sep=" ",header=F,fill=T)
#t <- table(m$V8)
#merged <- nrow(t)
base <- 26485
merged <- 26263
cat("\\begin{center}{\\huge GTF Annotations}\\end{center}\n")
cat("Reference: ",base,"\\\\Merged: ", merged," : ",round(100*merged/base,2),"\\%\\\\\n")

#
#  PCA
#
cat("\\newpage\n")
cat("\\begin{center}{\\huge PCA}\\end{center}\n")
cat("\\begin{center}\n")
cat("\\includegraphics[width=0.9\\textwidth]{images/pca}\n")
cat("\\end{center}\n")

#
#  Hits
#
genes <- fread(paste0("../processed/cuffdiff.",reference,"/gene_exp.diff"))
colnames(genes)[10] <- "log2FC"
genes <- genes[q_value <= 0.01 & abs(log2FC) >= 2,c(3:6,8:13)]
t <- table(genes$sample_1,genes$sample_2)
write.table(genes,file="genes.tsv",col.names=T,row.names=F,quote=F,sep="\t")
cat("\\newpage\n")
cat("\\begin{center}{\\huge Pairwise Significant Genes}\\end{center}\n")
print(xtable(t))
cat("q-value \\textless= 0.01; abs(log2FC) \\textgreater= 2\n")

genes <- genes[sample_1=="q2" & sample_2=="q4",c(1:2,5:10)]

#cat("\\newpage\n")
#cat("\\begin{center}{\\huge Genes by q-value, q2/q4}\\end{center}\n")
#x <- xtable(genes[order(q_value)])
#digits(x) <- c(2,2,2,2,2,2,2,5,5)
#align(x) <- "|rp{3cm}l|rr|rr|rr|"
#print(x,tabular.environment="longtable",floating=FALSE)

cat("\\newpage\n")
cat("\\begin{center}{\\huge Genes by log2 fold change, q2/q4}\\end{center}\n")
x <- xtable(genes[order(log2FC)])
digits(x) <- c(2,2,2,2,2,2,2,5,5)
align(x) <- "|r|p{2cm}p{4.5cm}|rr|rr|rr|"
print(x,tabular.environment="longtable",floating=FALSE,include.rownames=FALSE)

cat("\\end{document}\n")
sink()

system(paste0("cp ../processed/cuffdiff.",reference,"/gene_exp.diff genes.all.tsv"))
system("pdflatex report.tex")
