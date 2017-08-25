library("data.table")
library("ggplot2")
library("DESeq2")
library("matrixStats")

pdf("Rplots.pdf", height=10,width=12)

myReadCount <-fread("SSF-1972_H7TJWBGXY.unq.refseq.umi.dat", header=T,sep="\t",data.table=FALSE)
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

patient_status <- ifelse(grepl("MS",mysortedPheno$SampleName),'MS','HC')


### Individual cell Statistics ####

StatsCell <- data.frame(Well= colnames(mycountFrame), total_counts = colSums(mycountFrame), total_genes = colSums(mycountFrame > 0))

plotStatistics <- function(Frame,string){
    if(string=='count'){
        p1 <- ggplot(Frame, aes(Well, total_counts)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
        return(p1)
    }
    else{
       p1 <- ggplot(Frame, aes(Well, total_genes)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
       return (p1)
    }
 }

 # WellAD <- plotStatistics(StatsCell[1:26,],"count")
 # WellEF <- plotStatistics(StatsCell[27:50,],"count")
 # WellGH <- plotStatistics(StatsCell[51:74,],"count")

 # GeneAD <- plotStatistics(StatsCell[1:26,],"genes")
 # GeneEF <- plotStatistics(StatsCell[27:50,],"genes")
 # GeneGH <- plotStatistics(StatsCell[51:74,],"genes")



# Status <- mysortedPheno$Status
# design = ~Status
# dds = DESeqDataSetFromMatrix(countData = mycountFrame, colData = mysortedPheno, design = design)
# transformedData = varianceStabilizingTransformation(dds)

# #### Library Complexity #####

# libCom <- ggplot(StatsCell,aes(total_counts,total_genes)) + geom_point()

# ####### PCA of top500 genes #####

# rv <- rowVars(assay(transformedData))
# select <- order(rv, decreasing=T)[seq_len(min(1000,length(rv)))]
# pc <- prcomp(t(assay(transformedData)[select,]))

# # set condition
# condition <- Status
# P_type <- patient_status
# scores <- data.frame(pc$x, condition)

# pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition)),shape=factor(P_type))) + geom_point(size = 2) + ggtitle("Principal Components") 

MSPhenoData <- subset(mysortedPheno,grepl("^MS",SampleName))
HCPhenoData <- subset(mysortedPheno,grepl("^HC",SampleName))

MSCountData <- mycountFrame[,as.character(MSPhenoData$Well)]
HCCountData <- mycountFrame[,as.character(HCPhenoData$Well)]


FilteredMS <- MSCountData[rowSums(MSCountData) > 0,]
FilteredHC <- HCCountData[rowSums(HCCountData) > 0,]


plotMA_full = function(res) {
   ymax = max(res$log2FoldChange, na.rm = TRUE)
   ymin = min(res$log2FoldChange, na.rm = TRUE)
     plotMA(res, ylim = c(ymin, ymax))
 }

processPCA <- function(object){
  transformedData = varianceStabilizingTransformation(object)
  rv <- rowVars(assay(transformedData))
  select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
  pc <- prcomp(t(assay(transformedData)[select,]))
  return(pc)
}


### MS Analysis
MSdds = DESeqDataSetFromMatrix(countData = FilteredMS, colData = MSPhenoData,design=~Status)

###PCA analysis

MS_PCA <- processPCA(MSdds)
condition <- MSPhenoData$Status
scores <- data.frame(MS_PCA$x, condition)
pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition)))) + geom_point(size = 2) + ggtitle("Principal Components") 
print(pcaplot)

#### Collapse and DGE
MSddsColl <- collapseReplicates(MSdds, MSdds$SampleName, MSdds$Well)
MSddsColl=estimateSizeFactors(MSddsColl) 
MSddsColl=DESeq(MSddsColl,test = "LRT", reduced = ~ 1)
plotDispEsts(MSddsColl)
MSgroupComp <- results(MSddsColl)
MSgroupComp <- MSgroupComp[order(MSgroupComp$pvalue),]
plotMA_full(MSgroupComp)
print(MSgroupComp)
write.csv(MSgroupComp,'MSgroupComp.csv')


### HC Analysis

HCdds = DESeqDataSetFromMatrix(countData = FilteredHC, colData = HCPhenoData,design=~Status)

###PCA analysis

HC_PCA <- processPCA(HCdds)
condition <- HCPhenoData$Status
scores <- data.frame(HC_PCA$x, condition)
pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition)))) + geom_point(size = 2) + ggtitle("Principal Components") 
print(pcaplot)

#### Collapse/combine samples
HCddsColl <- collapseReplicates(HCdds, HCdds$SampleName, HCdds$Well)
CombinedObj <- HCddsColl
HCddsColl=estimateSizeFactors(HCddsColl)
HCddsColl=DESeq(HCddsColl,test = "LRT", reduced = ~ 1)
plotDispEsts(HCddsColl)
HCgroupComp <- results(HCddsColl)
HCgroupComp <- HCgroupComp[order(HCgroupComp$pvalue),]
write.csv(HCgroupComp,'HCgroupComp.csv')

### Time point Analysis

HCPhenoData$SampleName <- droplevels(HCPhenoData$SampleName)
#HCPhenoData$Status:SampleName <- factor(paste(HCPhenoData$Status,HCPhenoData$SampleName,sep=":"))

countHC <- assay(CombinedObj)
pData <- colData(CombinedObj)
pData$SampleName <- droplevels(pData$SampleName)
print(data.frame(pData[c(1,2,6)]))
HCddsColl = DESeqDataSetFromMatrix(countData = countHC, colData = pData,design= ~ Status + Status:SampleName)
#HCddsColl <- collapseReplicates(HCdds, HCdds$SampleName, HCdds$Well)
#HCddsColl=estimateSizeFactors(HCddsColl)
#HCddsColl=DESeq(HCddsColl,test="LRT", reduced = ~ Status + SampleName)

#print(HCddsColl)






