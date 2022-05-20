#latest adaptations were performed on 01/08/2022.
#please also refer to the WGCNA base analysis

#==========================================================================================================#
#                                              required packages                                           #
#==========================================================================================================#
library(DESeq2)
library(xlsx)
library(Hmisc)
library(dplyr)
library(viridis)
library(pheatmap)
library(variancePartition)
library(RColorBrewer)
library(pcaExplorer)
library(vsn)



#==========================================================================================================#
#                                                    functions                                             #
#==========================================================================================================#
IQRoutlier.detector <- function(df, barplot=T){
  d <- rcorr(as.matrix(df), type="pearson")
  
  barplot.df <- data.frame(Intersample.correlation=c(rowMeans(d$r), rowMedians(d$r)), 
                           Metric=c(rep("Correlation Means",
                                        length(colnames(d$r))), rep("Correlation Medians", length(colnames(d$r)))))
  if(barplot){
    #700 x 470
    outlier.barplot <- 
      barplot.df %>% 
      ggplot( aes(x=Metric, y=Intersample.correlation, fill=Metric)) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      geom_text(aes(label=c(rownames(d$r), rownames(d$r))))+
      theme(
        legend.position="Type",
        plot.title = element_text(size=11)
      ) +
      #  ggtitle("Boxplots showing outliers based on inter-sample correlation means and medians before data correction") +
      theme_bw()+
      theme(plot.title = element_text(hjust=0.5))+
      xlab("")+
      #geom_hline(yintercept=median(d_core)-3*mad(d_core),  linetype="dashed")+
      geom_segment(aes(x=1.7,xend=2.3,y=ifelse(median(d$r)-3*mad(d$r) < -1, -1, median(d$r)-3*mad(d$r)),yend=ifelse(median(d$r)-3*mad(d$r) < -1, -1, median(d$r)-3*mad(d$r))),color="red", linetype="dashed")+
      geom_segment(aes(x=1.7,xend=2.3,y=ifelse(median(d$r)+3*mad(d$r) > 1, 1, median(d$r)+3*mad(d$r)),yend=ifelse(median(d$r)+3*mad(d$r) > 1, 1, median(d$r)+3*mad(d$r))),color="red", linetype="dashed")+
      geom_segment(aes(x=0.7,xend=1.3,y=ifelse(mean(d$r)-3*sd(d$r) < -1, -1, mean(d$r)-3*sd(d$r)),yend=ifelse(mean(d$r)-3*sd(d$r) < -1, -1, mean(d$r)-3*sd(d$r))),color="red", linetype="dashed")+
      geom_segment(aes(x=0.7,xend=1.3,y=ifelse(mean(d$r)+3*sd(d$r) > 1, 1, mean(d$r)+3*sd(d$r)),yend=ifelse(mean(d$r)+3*sd(d$r) > 1, 1, mean(d$r)+3*sd(d$r))),color="red", linetype="dashed")
    return(outlier.barplot)}
  else {
    outlier.heatmap <- pheatmap(d$r)
    return(outlier.heatmap)
  }
}

plotPCA.custom <- function (object, intgroup = "condition", ntop = 44099, returnData = TRUE, metadata = metadata, pc.1=1, pc.2=2) 
{
  rv <- rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(object[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% colnames(metadata))) {
    stop("the argument 'intgroup' should specify columns of metadata")
  }
  intgroup.df <- as.data.frame(metadata[intgroup])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    metadata[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc.1], PC2 = pca$x[, pc.2], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = paste("PC",pc.1), y = paste("PC",pc.2), 
                              color = "group")) + geom_point(size = 3) + xlab(paste0(paste("PC", pc.1, ": ", sep=""),
                                                                                     round(percentVar[1] * 100), "% variance")) + ylab(paste0(paste("PC", pc2, ": ", sep=""), 
                                                                                                                                              round(percentVar[2] * 100), "% variance")) + coord_fixed()
}



#==========================================================================================================#
#                                   RNA-seq analysis start: counts files                                   #
#==========================================================================================================#
#set working directory to counts here
setwd(...)

b1 <- read.table("UMC-AB-b002 - Batch1.csv", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
b2 <- read.table(file = "UMC-AB-b002.ReadCounts.tsv", sep = '\t', header = TRUE, row.names=1, check.names = F)

#Merge dataframes but drop NAs
df <- merge(b1, b2, by=0, all=FALSE)
#Set gene names as index
rownames(df) <- df[,1]
#Remove gene names from df
df <- df[,-1]



#==========================================================================================================#
#                                      RNA-seq analysis start: covariate file                              #
#set working directory to counts here
setwd(...)
#Reading in metadata
metadata.b1 <- read.xlsx("Metadata_MKMD_RNAseq_Batch1.xlsx", sheetIndex=1, as.data.frame=T)
metadata.b2 <- read.xlsx("Metadata_MKMD_RNAseq_Batch2.xlsx", sheetIndex=1, as.data.frame=T)
colnames(metadata.b1) <- colnames(metadata.b2)
meta.df <- rbind(metadata.b1, metadata.b2)
meta.df <- meta.df[colnames(meta.df) != c("QC", "Sequencing.lane")]
meta.df$Sample.Name. <- as.character(meta.df$Sample.Name.)
meta.df$stimulation <- NA
meta.df[rownames(meta.df[meta.df$Concentration == "CTR",]),]$stimulation <- 1
meta.df[rownames(meta.df[meta.df$Concentration != "CTR",]),]$stimulation <- 2
meta.df$stimulation <- as.factor(meta.df$stimulation)

for (i in rownames(meta.df)){
  if (meta.df[i,]$Concentration == "10 mM")
    meta.df[i,]$Concentration <- "10mM"
  if (meta.df[i,]$Concentration == "50 mM")
    meta.df[i,]$Concentration <- "50mM"
}

meta.df$Concentration <- as.factor(meta.df$Concentration)


#Creating a new metadata file with additional covariates
new.meta <- data.frame("librarysize" = colSums(df))
rownames(new.meta) <- colnames(df)
new.meta$batch <- NA
new.meta[rownames(new.meta) %in% colnames(b1),]$batch <- "1"
new.meta[rownames(new.meta) %in% colnames(b2),]$batch <- "2"
new.meta$percent.ERCC <- colSums(x = df[grep(pattern = "^ERCC-", x = rownames(df), value = TRUE), 
                                     drop = FALSE,])/colSums(df) * 100
new.meta$percent.MT <- colSums(x = df[grep(pattern = "^MT-", x = rownames(df), value = TRUE), 
                                      drop = FALSE,])/colSums(df) * 100

#Parsing sample names
meta.df$Sample.Name.[is.na(meta.df$Sample.Name.)] <- as.character(rownames(new.meta[new.meta$batch == 2,]))

meta.df.all <- merge(meta.df, new.meta, by.x=1, by.y=0)
rownames(meta.df.all) <- meta.df.all$Sample.Name.



#==========================================================================================================#
#                             RNA-seq analysis start: doublet comparison & removal                         #
#==========================================================================================================#
#Doublets are: 24, 25, 31, 32

#Subset:
meta.doub <- meta.df.all[meta.df.all$Sample.NR. %in% c(24, 25, 31, 32),]

write.table(meta.doub, "D20 sample doublets.txt")

summary(meta.df.all$DIV)



#==========================================================================================================#
#            RNA-seq analysis start: failed sample removal & new df & covariate file creation              #
#==========================================================================================================#
remove.samples <- c(rownames(meta.df.all[meta.df.all$Sample.Name. %in% rownames(meta.doub) & meta.df.all$batch == 2,]), 
                    rownames(meta.df.all[meta.df.all$Sample.NR. == "41" & meta.df.all$batch == 1,]))
df.clean <- df[!colnames(df) %in% remove.samples]
meta.clean <- meta.df.all[!rownames(meta.df.all) %in% remove.samples,]



#==========================================================================================================#
#                                           RNA-seq analysis start: QC on all data                           #
#==========================================================================================================#
#Library size all samples
hist(log10(meta.clean$librarysize), xlab="Library sizes (log10 counts)", 
     main=paste("Counts together \n\n", "mean = ", round(mean(meta.clean$librarysize),2), ";",
                "median = ", round(median(meta.clean$librarysize),2)), 
     breaks=20, col="grey80", ylab="Number of samples")
#Library size b1 vs b2
hist(log10(meta.clean[meta.clean$batch == 1,]$librarysize), xlab="Library sizes (log10 counts)", 
     main=paste("Counts batch 1\n\n", "mean = ", round(mean(meta.clean[meta.clean$batch == 1,]$librarysize),2), ";",
                "median = ", round(median(meta.clean[meta.clean$batch == 1,]$librarysize),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

hist(log10(meta.clean[meta.clean$batch == 2,]$librarysize), xlab="Library sizes (log10 counts)", 
     main=paste("Counts batch 2\n\n", "mean = ", round(mean(meta.clean[meta.clean$batch == 2,]$librarysize),2), ";",
                "median = ", round(median(meta.clean[meta.clean$batch == 2,]$librarysize),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

filt.meta <- meta.clean[meta.clean$librarysize > 500000 & meta.clean$librarysize < 15000000,]
filt.df <- df.clean[colnames(df.clean) %in% rownames(filt.meta)]

#Mapping the library size distributions
hist(log10(filt.meta$librarysize), xlab="Library sizes (log10 counts)", 
     main=paste("Counts (filtered; both batches together) \n\n", "mean = ", round(mean(filt.meta$librarysize),2), ";",
                "median = ", round(median(filt.meta$librarysize),2)), 
     breaks=20, col="grey80", ylab="Number of samples")
#Library sizes b1 vs b2 (after sample filtering)
hist(log10(filt.meta[filt.meta$batch == 1,]$librarysize), xlab="Library sizes (log10 counts)", 
     main=paste("Counts batch 1\n\n", "mean = ", round(mean(filt.meta[filt.meta$batch == 1,]$librarysize),2), ";",
                "median = ", round(median(filt.meta[filt.meta$batch == 1,]$librarysize),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

hist(log10(filt.meta[filt.meta$batch == 2,]$librarysize), xlab="Library sizes (log10 counts)", 
     main=paste("Counts batch 2\n\n", "mean = ", round(mean(filt.meta[filt.meta$batch == 2,]$librarysize),2), ";",
                "median = ", round(median(filt.meta[filt.meta$batch == 2,]$librarysize),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

#mtRNA percentage
hist(log10(filt.meta$percent.MT), xlab="Mitochondrial proportion", 
     main=paste("Mitochondrial percentage \n\n", "mean = ", round(mean(filt.meta$percent.MT ),2), ";",
                "median = ", round(median(filt.meta$percent.MT ),2)), 
     breaks=20, col="grey80", ylab="Number of Samples")

filt.meta <- filt.meta[filt.meta$percent.MT < 37,]
filt.df <- df.clean[colnames(df.clean) %in% rownames(filt.meta)]

#Post-filter distribution of mtRNA percentage
hist(log10(filt.meta$percent.MT), xlab="Mitochondrial proportion", 
     main=paste("Mitochondrial percentage \n\n", "mean = ", round(mean(filt.meta$percent.MT ),2), ";",
                "median = ", round(median(filt.meta$percent.MT ),2)), 
     breaks=20, col="grey80", ylab="Number of Samples")



#Outlier detection based on IQR correlation
sOutliers <- IQRoutlier.detector(filt.df,barplot=T)
#gOutliers <- IQRoutlier.detector(t(filt.df))

filt.meta$Cell.line <- as.factor(as.numeric(filt.meta$Cell.line))

pca <- plotPCA.custom(as.matrix(filt.df), intgroup=c("Cell.line", "Concentration"), ntop = 20000, returnData=TRUE, 
                      metadata=filt.meta)
PoV <- round(100 * attr(pca, "percentVar"))
ggplot(pca, aes(PC1, PC2, color=pca[["Cell.line"]], label=name)) +
  geom_point(size=2) +
  geom_text(aes(label=name),hjust=0, vjust=0)+
  #  geom_label(data=pca, aes(label=all_covs_sub$cell_type_advanced))+
  labs(x = paste0("PC1: ",PoV[1],"% variance"), y = paste0("PC2: ",PoV[2],"% variance"), 
       color = "Cell.line") +  
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))



#==========================================================================================================#
#                                     RNA-seq analysis start: filtering                                    #
#==========================================================================================================#
gfilt.df <- filt.df[rowSums(filt.df) >= length(colnames(filt.df)),] #Removes ERCC SIs



#==========================================================================================================#
#                              RNA-seq analysis start: post-filtering QC                                   #
#==========================================================================================================#
sOutliers.filt <- IQRoutlier.detector(gfilt.df,barplot=T)
#gOutliers.filt <- IQRoutlier.detector(t(gfilt.df), barplot=T)

gfilt.df <- gfilt.df[sort(colnames(gfilt.df),decreasing=T)]
filt.meta <- filt.meta[sort(rownames(filt.meta), decreasing=T),]
filt.meta$librarysize.postfilter <- colSums(gfilt.df)

#Distribution
hist(log10(filt.meta$librarysize.postfilter), xlab="Library sizes (log10 counts)", 
     main=paste("Counts \n\n", "mean = ", round(mean(filt.meta$librarysize.postfilter),2), ";",
                "median = ", round(median(filt.meta$librarysize.postfilter),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

#Library sizes b1 vs b2 (after gene filtering)
hist(log10(filt.meta[filt.meta$batch == 1,]$librarysize.postfilter), xlab="Library sizes (log10 counts)", 
     main=paste("Counts batch 1 post-filter \n\n", "mean = ", round(mean(filt.meta[filt.meta$batch == 1,]$librarysize),2), ";",
                "median = ", round(median(filt.meta[filt.meta$batch == 1,]$librarysize.postfilter),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

hist(log10(filt.meta[filt.meta$batch == 2,]$librarysize.postfilter), xlab="Library sizes (log10 counts)", 
     main=paste("Counts batch 2 post-filter \n\n", "mean = ", round(mean(filt.meta[filt.meta$batch == 2,]$librarysize),2), ";",
                "median = ", round(median(filt.meta[filt.meta$batch == 2,]$librarysize.postfilter),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

#PCA
pca.filt <- plotPCA.custom(as.matrix(gfilt.df), intgroup=c("Concentration"), 
                           ntop = 50000, returnData=TRUE, metadata=filt.meta)
PoV.filt <- round(100 * attr(pca.filt, "percentVar"))
ggplot(pca.filt, aes(PC1, PC2, color=pca.filt[["Concentration"]], label=name)) +
  geom_point(size=2) +
  geom_text(aes(label=name),hjust=0, vjust=0)+
  #  geom_label(data=pca, aes(label=all_covs_sub$cell_type_advanced))+
  labs(x = paste0("PC1: ",PoV[1],"% variance"), y = paste0("PC2: ",PoV[2],"% variance"), 
       color = "Concentration") +  
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))



#==========================================================================================================#
#                                RNA-seq analysis start: Covariate analysis                                #
#==========================================================================================================#
#Changing factors to factors
filt.meta$Cell.line <- as.factor(filt.meta$Cell.line)
filt.meta$Duplicate <- as.factor(filt.meta$Duplicate)

filt.meta_clean <- filt.meta[!colnames(filt.meta) %in% c("batch", "QC", "Sequencing.lane","Timepoint", "Sample.NR.")]
filt.meta_clean$Experiment <- as.factor(filt.meta_clean$Experiment)
filt.meta_clean$Cell.line <- as.factor(filt.meta_clean$Cell.line)
filt.meta_clean$DIV <- as.factor(filt.meta_clean$DIV)
filt.meta_clean$Concentration <- as.factor(filt.meta_clean$Concentration) #Condition is the Concentration as factor
filt.meta_clean$Sample.Name. <- as.numeric(as.factor(filt.meta_clean$Sample.Name.))
filt.meta_clean$librarysize <- log(filt.meta_clean$librarysize)

#Check which samples were removed along the way
samples.removed <- rownames(new.meta)[!rownames(new.meta) %in% rownames(filt.meta_clean)]
write.table(samples.removed, file="Samples removed.txt")


#====Colinearity analysis of covariates====#
heatmap.color.code <- rev(brewer.pal(11,"RdYlBu"))
#Check co-linearity of covariates for subsetted samples
cov_for.cor <- filt.meta_clean
cov_for.cor$Concentration <- as.factor(as.numeric(cov_for.cor$Condition))
cov_for.cor$Experiment <- as.factor(as.numeric(cov_for.cor$Experiment))
cov_for.cor$DIV <- as.factor(as.numeric(cov_for.cor$DIV))
cov_for.cor$Cell.line <- as.factor(as.numeric(cov_for.cor$Cell.line))
covar.cor.all <- rcorr(as.matrix(cov_for.cor), type="spearman")
covar.cor <- covar.cor.all$r
pheatmap(covar.cor, color = heatmap.color.code)


#====PC loadings for covariates====#
pca.cov <- prcomp(t(filt.df))
PC.covs <- pca.cov$x[,1:25]
PC_cov_correlation <- rcorr(as.matrix(PC.covs), as.matrix(cov_for.cor), type="spearman")
PC_variance.explained <- PC_cov_correlation$r[1:25,26:length(rownames(PC_cov_correlation$r))]
PC_cov.cor.p <- PC_cov_correlation$P[1:25,26:length(rownames(PC_cov_correlation$r))]

pheatmap(PC_variance.explained, cluster_rows = F, cluster_cols = F, show_colnames = T,
         color = heatmap.color.code)

pheatmap(PC_cov.cor.p, cluster_rows = F, cluster_cols = F, show_colnames = T,
         color = colorRampPalette(c("red", "white"))(50),
         display_numbers = T)


#====Variance partitioning====#
#Regression formula for variance partitioning
#regform <- ~ Sample.Name. + (1|Condition) + (1|Concentration) + (1|Cell.line) + librarysize + (1|Experiment) + (1|Duplicate)
#regform <- ~ (1|Concentration) + (1|Cell.line) + librarysize + (1|Experiment) + (1|Duplicate)
#regform <- ~ (1|Cell.line) + (1|Concentration)
regform <- ~ percent.ERCC + percent.MT + (1|Concentration) + (1|Cell.line) + (1|Experiment) + (1|Duplicate)
varPart <- fitExtractVarPartModel(filt.df, regform, filt.meta_clean)

vp <- sortCols(varPart)
#Plotting
plotVarPart(vp) +
  ggeasy::easy_rotate_x_labels(angle = 70, side =c("right")) #Indeed the variance partition plot shows that the SVs explain a large amount of the variation



#==========================================================================================================#
#                                   RNA-seq analysis start: normalization                                  #
#==========================================================================================================#
#Sort columns
filt.meta_clean <- filt.meta_clean[sort(rownames(filt.meta_clean)),]
gfilt.df <- gfilt.df[sort(colnames(gfilt.df))]

all(rownames(filt.meta_clean) == colnames(gfilt.df))


filt.meta_clean$log.librarysize <- log(filt.meta_clean$librarysize)
#Create a DESeq2 matrix
dds.base <- DESeqDataSetFromMatrix(countData = as.matrix(gfilt.df),
                                  colData = filt.meta_clean,
                                  design = ~ stimulation)

dds.complex <- DESeqDataSetFromMatrix(countData = as.matrix(gfilt.df),
                                   colData = filt.meta_clean,
                                   design = ~ log.librarysize + Cell.line + stimulation)



#version 2 based on dds object
# Estimate library size correction scaling factors
dds.base <- estimateSizeFactors(dds.base)
dds.complex <- estimateSizeFactors(dds.complex)

vsd.base <- vst(dds.base)
vsd.complex <- vst(dds.complex)

#Doing QC on the vst-transformed values
meanSdPlot(assay(vsd.complex))

sampleDists <- dist(t(assay(vsd.complex)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd.complex$stimulation, vsd.complex$Concentration, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)



#==========================================================================================================#
#                              RNA-seq analysis start: post-normalization QC                               #
#==========================================================================================================#
#QC QC QC
norm.df <- assay(vsd.complex)
sOutliers.norm <- IQRoutlier.detector(norm.df,barplot=T)
#gOutliers.norm <- IQRoutlier.detector(t(norm.df), barplot=T)

#PCA
pca.norm <- plotPCA.custom(as.matrix(norm.df), intgroup=c("Cell.line", "stimulation","Concentration","DIV"), 
                           ntop = 50000, returnData=TRUE, metadata=filt.meta_clean)
PoV.norm <- round(100 * attr(pca.norm, "percentVar"))
ggplot(pca.norm, aes(PC1, PC2, color=pca.norm[["DIV"]], label=name)) +
  geom_point(size=2) +
#  geom_text(aes(label=name),hjust=0, vjust=0)+
  #  geom_label(data=pca, aes(label=all_covs_sub$cell_type_advanced))+
  labs(x = paste0("PC1: ",PoV.norm[1],"% variance"), y = paste0("PC2: ",PoV.norm[2],"% variance"), 
       color = "DIV") +  
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))

DESeq2:::plotPCA.DESeqTransform
plotPCA(vsd.complex, intgroup=c("Concentration","stimulation"), ntop = 500, returnData=F)

#pcaExplorer(dds=dds.complex)

#Remove batch effects
#Create a meta df that corresponds to the sorted columns of the vsd assay but has covariates as numeric factors
covs_forremBatch <- filt.meta_clean
covs_forremBatch$Concentration <- as.factor(as.numeric(covs_forremBatch$Condition))
covs_forremBatch$Experiment <- as.factor(as.numeric(covs_forremBatch$Experiment))
covs_forremBatch$DIV <- as.factor(as.numeric(covs_forremBatch$DIV))
covs_forremBatch$Cell.line <- as.factor(as.numeric(covs_forremBatch$Cell.line))

batch.rem <- removeBatchEffect(as.matrix(assay(vsd.complex)), batch=covs_forremBatch$Cell.line, 
                               covariates=as.matrix(cbind(covs_forremBatch$log.librarysize)),
                               design=model.matrix(~ covs_forremBatch$Concentration))

#Create a variable that reflects the interaction effect.

#PCA (700 x 550)
pca.cor <- plotPCA.custom(as.matrix(batch.rem), intgroup=c("Concentration","stimulation","Cell.line","DIV"),
                          ntop = 50000, returnData=TRUE, metadata=filt.meta_clean, pc.1=1,pc.2=2)
PoV.cor <- round(100 * attr(pca.cor, "percentVar"))
ggplot(pca.cor, aes(PC1, PC2, color=pca.cor[[c("Cell.line")]], label=name)) +
  geom_point(size=2) +
#  geom_text(aes(label=name),hjust=0, vjust=0)+
  #  geom_label(data=pca, aes(label=all_covs_sub$cell_type_advanced))+
  labs(x = paste0("PC1: ",PoV.cor[1],"% variance"), y = paste0("PC2: ",PoV.cor[2],"% variance"), 
       color = "Cell line") +  
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))

ggplot(pca.cor, aes(PC1, PC2, color=pca.cor[[c("Concentration")]], label=name)) +
  geom_point(size=2) +
#  geom_text(aes(label=name),hjust=0, vjust=0)+
  #  geom_label(data=pca, aes(label=all_covs_sub$cell_type_advanced))+
  labs(x = paste0("PC1: ",PoV.cor[1],"% variance"), y = paste0("PC2: ",PoV.cor[2],"% variance"), 
       color = "Concentration") +  
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))

ggplot(pca.cor, aes(PC1, PC2, color=pca.cor[[c("DIV")]], label=name)) +
  geom_point(size=2) +
  #  geom_text(aes(label=name),hjust=0, vjust=0)+
  #  geom_label(data=pca, aes(label=all_covs_sub$cell_type_advanced))+
  labs(x = paste0("PC1: ",PoV.cor[1],"% variance"), y = paste0("PC2: ",PoV.cor[2],"% variance"), 
       color = "DIV") +  
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))


#====PC loadings for covariates====#
PC.covs.postcor.all <- prcomp(t(filt.df))
PC.covs.postcor <- PC.covs.postcor.all$x[,1:25]
PC_cov_correlation.postcor <- rcorr(as.matrix(PC.covs.postcor), as.matrix(covs_forremBatch), type="spearman")
PC_variance.explained.postcor <- PC_cov_correlation.postcor$r[1:25,26:length(rownames(PC_cov_correlation.postcor$r))]

PC_cov.cor.p.postcor <- PC_cov_correlation.postcor$P[1:25,26:length(rownames(PC_cov_correlation.postcor$r))]

pheatmap(PC_variance.explained.postcor, cluster_rows = F, cluster_cols = F, show_colnames = T,
         color = heatmap.color.code)

pheatmap(PC_cov.cor.p.postcor, cluster_rows = F, cluster_cols = F, show_colnames = T,
         color = colorRampPalette(c("red", "white"))(50),
         display_numbers = T)



#==========================================================================================================#
#                                            RNA-seq analysis: DESeq2                                      #
#==========================================================================================================#
dds.base <- DESeq(dds.base)
dds.complex <- DESeq(dds.complex)

substr(names(mcols(dds.complex)),1,10) 



res.base <- results(dds.base)
res.complex <- results(dds.complex)

summary(res.base)
summary(res.complex)
sum(res.base$padj < 0.05, na.rm=TRUE)
sum(res.complex$padj < 0.05, na.rm=TRUE)

sum(res.complex$log2FoldChange < -2, res.complex$log2FoldChange > 2, na.rm=TRUE)

top20upreg.base <- data.frame(data.frame(res.base)[order(res.base$log2FoldChange, decreasing=T),])[1:20,]
top20downreg.base <- data.frame(data.frame(res.base)[order(res.base$log2FoldChange, decreasing=F),])[1:20,]

top20upreg.complex <- data.frame(data.frame(res.complex)[order(res.complex$log2FoldChange, decreasing=T),])[1:20,]
top20downreg.complex <- data.frame(data.frame(res.complex)[order(res.complex$log2FoldChange, decreasing=F),])[1:20,]



write.table(top20upreg.base, "upreg.genes_base.txt")
write.table(top20downreg.base, "downreg.genes_base.txt")

write.table(top20upreg.complex, "upreg.genes_complex.txt")
write.table(top20downreg.complex, "downreg.genes_complex.txt")

pheatmap(assay(dea.dds)[rownames(assay(dea.dds)) %in% rownames(top50upreg),])

#Look at LRT
dds.complex.2 <- DESeq(dds.complex, test="LRT", reduced=~log.librarysize + Cell.line)

res.complex.2 <- results(dds.complex.2)
sum(res.complex.2$padj < 0.05, na.rm=TRUE)


#LFC shrinkage
resLFC <- lfcShrink(dea.dds, coef="Concentration_3_vs_1", type="normal")

DESeq2::plotMA(resLFC, ylim=c(-2,2))
colData(dea.dds)

summary(resLFC$padj < 0.05)
sum(res$padj < 0.05, na.rm=TRUE)

plotCounts(dea.dds, gene=which.min(res$padj), intgroup="Concentration")



#==========================================================================================================#
#                                              RNA-seq analysis: old code                                  #
#==========================================================================================================#
#Mapping the library size distributions
hist(log10(meta.df.all$librarysize[meta.df.all$batch == 2]), xlab="Library sizes (log10 counts)", 
     breaks=20, col="grey80", ylab="Number of samples")

hist(log10(meta.df.all$librarysize[meta.df.all$batch == 1]), xlab="Library sizes (log10 counts)", 
     main=paste("Counts \n\n", "mean = ", round(mean(meta.df.all$librarysize[meta.df.all$batch == 1]),2), ";",
                "median = ", round(median(meta.df.all$librarysize[meta.df.all$batch == 1]),2)), 
     breaks=20, col="grey80", ylab="Number of cells")

hist(log10(meta.df.all$librarysize[meta.df.all$DIV == c("Day 20", "Day 24") & meta.df.all$batch == 2]), xlab="Library sizes (log10 counts)", 
     breaks=20, col="grey80", ylab="Number of samples")
hist(log10(meta.df.all$librarysize[meta.df.all$DIV == c("Day 104", "Day 106", "Day 107") & meta.df.all$batch == 2]), xlab="Library sizes (log10 counts)", 
     breaks=20, col="grey80", ylab="Number of samples")