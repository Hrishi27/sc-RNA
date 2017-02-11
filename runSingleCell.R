library(scater)
library(data.table)
library(dplyr)
library(readxl)
library(gtools)
library(scran)

###Reading the expression matrix

myReadCount <-fread('Probiotic_VSL3_RNA-DGE_raw_data/SSF-1972_H7TJWBGXY.unq.refseq.umi.dat', header=TRUE, data.table=FALSE)
rownames(myReadCount)=myReadCount[,1]
myReadCount[, 1] = NULL

cnames = unlist(lapply(strsplit(colnames(myReadCount), split = "_"), tail, n = 1))
colnames(myReadCount) = cnames

### Reading the phenoFile

myPheno <-read.table("pheno.txt", header=T,sep="\t")
Wells<-na.omit(myPheno$Well)
remove_column<-setdiff(cnames,Wells)

mycountFrame<-myReadCount[,!(names(myReadCount) %in% remove_column)]

rownames(myPheno) <- myPheno$Well
mysortedPheno <-myPheno[names(mycountFrame),]
pheno_data <-new("AnnotatedDataFrame",mysortedPheno)

########
cellstats = data.frame(cell = colnames(mycountFrame), total_counts = colSums(mycountFrame), 
     genes_detected = colSums(mycountFrame > 0), robust_genes_detected = colSums(mycountFrame > 
         10))

print(cellstats)

#######

umi <- newSCESet(
    countData = mycountFrame,
    phenoData = pheno_data
)

umi <- scater::calculateQCMetrics(umi)
plotQC(umi, type = "exprs-freq-vs-mean")
keep_feature <- rowSums(counts(umi) > 0) > 0
keep_feature <- rowMeans(counts(umi)) >= 3 
umi <-umi[keep_feature,]


###Plots###

hist(umi$total_counts/1e6, xlab="Library sizes (millions)", main="",
     breaks=20, col="grey80", ylab="Number of cells")
hist(umi$total_features, xlab="Number of expressed genes", main="",
     breaks=20, col="grey80", ylab="Number of cells")

plotPCA(umi,colour_by = "Status")

umi <- computeSumFactors(umi, sizes=c(10, 20, 30))
summary(sizeFactors(umi))

umi <- normalize.SCESet(umi)


plot(sizeFactors(umi), colSums(counts(umi))/1e6, log="xy",
   ylab="Library Size (Total Counts in Millions)", xlab="Pooled Size Factor Estimate",
    main="Normalization factor versus library size")


scater::plotQC(umi, type = "highest-expression")

