#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2", version = "3.8")

library("DESeq2")
library("tidyr")
library("stringr")

setwd("C:/Users/tbleazar/OneDrive - MHRA/Documents/R/rnaseq/152deseq")
countfiles = list.files(path=getwd(), pattern="*.htseq.a5.txt")
countfiles
# countfilesfiltered = countfiles[grepl("A549", countfiles)]
# countfilesfiltered = countfilesfiltered[grepl("6hpi", countfilesfiltered)]
# countfilesfiltered = countfilesfiltered[! grepl("WT", countfilesfiltered)]
#countfilesfiltered = countfiles
countfilesfiltered = countfiles[grepl("B_", countfiles)]
countfilesfiltered = countfilesfiltered[! grepl("2B", countfilesfiltered)]
countfilesfiltered = countfilesfiltered[! grepl("4B", countfilesfiltered)]
countfilesfiltered

#sampleNames = countfilesfiltered
#sampleNames = c("A_Mock_S1", "B_Mock_S2", "C_Mock_S3", "D_Mock_S4", "E_In1382_S5", "F_In1382_S6", "G_In1382_S7", "H_In1382_S8")
#sampleNames = c("Day1_1A","Day1_2A","Day2_1A","Day2_2A","Day3_2A","Day4_1A","Day4_2A","Day5_1A","Day5_2A")
sampleNames = str_replace(countfilesfiltered, ".htseq.a5.txt","")
sampleNames = str_replace(sampleNames, "152_","")
sampleNames

sampleFiles = countfilesfiltered

sampleCondition = c("Untreated","Antitoxin_B","Untreated","Untreated","Antitoxin_B","Untreated","Antitoxin_B","Untreated","Antitoxin_B")
treatments = c("Untreated", "Antitoxin_B")

#sampleCondition = c("A549_ni_6hpi","A549_ni_6hpi","A549_ni_6hpi","A549_ni_6hpi","A549_ni_6hpi","A549_Vx_6hpi","A549_Vx_6hpi","A549_Vx_6hpi","A549_Vx_6hpi","A549_Vx_6hpi")
#treatments = c("A549_ni_6hpi","A549_Vx_6hpi")

sampleTable = data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = getwd(), design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition, levels = treatments)
dds <- DESeq(ddsHTSeq)
res <- results(dds)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'

#write.csv(resdata, file="/home/AD/tbleazar/169/results/hsa.A549_ni_6hpi-vs-A549_Vx_6hpi.csv")
write.csv(resdata, file = "152-Untreated-vs-Antitoxin-B.csv")

rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)
library("genefilter")
library("ggplot2")
library("grDevices")
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))
condition <- sampleCondition
scores <- data.frame(pc$x, condition)
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 3)
  + ggtitle("Normalised Expression PCA Plot")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = "right",
    legend.key = element_rect(fill = 'grey'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

#head(scores)
#head(pc)
summary(pc)

#ggsave(pcaplot,file="152-Untreated-vs-Toxin-A.pdf")
#Better to use the plot window then export

#here are some commands now to make a heatmap

#rld <- rlogTransformation(dds, blind=T)
#vsd <- varianceStabilizingTransformation(dds, blind=T)
#hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)


distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), sampleNames)

#summary(mat)
#head(mat)

#heatmap(mat)

#install.packages("gplots")
library(gplots)
pdf("heattest.pdf", width=14, height=14)

heatmap.2(mat, trace = "none", margin = c(22,22))

dev.off()
