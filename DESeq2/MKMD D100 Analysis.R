#latest changes were performed on 04/01/2021 (extraction of data for manuscript)

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

options(check.names = FALSE)


#==========================================================================================================#
#                                   RNA-seq analysis start: normalization                                  #
#==========================================================================================================#
#Subsetting a df to samples of interest (D100)
rownames(gfilt.df) <- make.unique(gsub("__.*","",rownames(gfilt.df)))

subset.df <- gfilt.df[colnames(gfilt.df) %in% rownames(filt.meta[filt.meta$DIV %in% c("W15"),])]
subset.df <- subset.df[rowSums(subset.df) >= length(colnames(subset.df)),] #Also removes ERCC SIs

subset.meta <- filt.meta_clean[rownames(filt.meta_clean) %in% colnames(subset.df),]

#Sort columns
subset.meta <- subset.meta[sort(rownames(subset.meta)),]
subset.df <- subset.df[sort(colnames(subset.df))]

all(rownames(subset.meta) == colnames(subset.df))

#Quick density QC
#550 x 450
plot.density(subset.meta,to.color="Concentration",to.plot="log10(`Library size`)")

#Quick covariate analysis
cov_for.cor <- cov_for.cor[sort(rownames(cov_for.cor)),]
cov_for.cor.sub <- cov_for.cor[cov_for.cor$DIV %in% c("1"),]
cov_for.cor.sub <- cov_for.cor.sub[!colnames(cov_for.cor.sub) %in% c("DIV","Timepoint","Sequencing batch")]
covar.cor.sub <- rcorr(as.matrix(cov_for.cor.sub), type="pearson")
covar.cor.sub_1 <- covar.cor.sub$r

#500x450
pheatmap(covar.cor.sub_1, color = heatmap.color.code)

#D100 PCA-covariate analysis
#500x500
D100.PCAcovar <- PCA.covar(subset.df,cov_for.cor.sub,data.r=T,type="pearson")

#For correcting the data
subset.meta$'Log(library size)' <- scale(subset.meta$`Library size`)
subset.meta$`Scaled percent. MT` <- scale(subset.meta$"Percent. MT")

#Create a DESeq2 matrix
dds.complex <- DESeqDataSetFromMatrix(countData = as.matrix(subset.df),
                                      colData = data.frame(subset.meta),
                                      design = ~ Log.library.size. + Experiment + Percent..ERCC + Stimulation)


#version 2 based on dds object
# Estimate library size correction scaling factors
dds.complex <- estimateSizeFactors(dds.complex)

vsd.complex <- vst(dds.complex)

#Doing QC on the vst-transformed values
meanSdPlot(assay(vsd.complex))

sampleDists <- dist(t(assay(vsd.complex)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd.complex)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

annot_df <- data.frame(colData(dds.complex)$Concentration,colData(dds.complex)$Cell.line)
colnames(annot_df) <- c("Concentration","Cell line")
rownames(annot_df) <- colnames(dds.complex)

annot_colors <- vector('list', length=1)

annot_colors$Concentration[["CTR"]] <- "#2E7281"
annot_colors$Concentration[["10mM"]] <- "black"
annot_colors$Concentration[["50mM"]] <- "#50B8CF"

annot_colors$`Cell line`[["1.5"]] <- "#E7B800"
annot_colors$`Cell line`[["2.6"]] <- "darkblue"
annot_colors$`Cell line`[["4.6"]] <- "#FC4E07"

#Heatmap intersample-distance
pheatmap(sampleDistMatrix, 
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,annotation_row = annot_df, annotation_colors = annot_colors)



#==========================================================================================================#
#                              RNA-seq analysis start: post-normalization QC                               #
#==========================================================================================================#
#QC QC QC
norm.df <- assay(vsd.complex)
sOutliers.norm <- IQRoutlier.detector(norm.df,barplot=T)
#gOutliers.norm <- IQRoutlier.detector(t(norm.df), barplot=T)

#PCA
#550x450
pca.norm <- plotPCA.custom(as.matrix(norm.df), intgroup=c("Cell line", "Stimulation","Concentration","DIV","Experiment"), 
                           ntop = 50000, returnData=TRUE, metadata=subset.meta)
PoV.norm <- round(100 * attr(pca.norm, "percentVar"))
PCAplot(pca.norm, Condition="Concentration", Shape="Cell line", PoV.df=PoV.norm, colors=c('#2E7281','black','#50B8CF'))


PCAplot(pca.norm, "Experiment", PoV.df=PoV.norm, colors=viridis(5))

#plotPCA(vsd.complex, intgroup=c("Cell.line","stimulation"), ntop = 50000, returnData=F)

#pcaExplorer(dds=dds.complex)

#Remove batch effects
#Create a meta df that corresponds to the sorted columns of the vsd assay but has covariates as numeric factors
covs_forremBatch <- subset.meta
covs_forremBatch$Concentration <- as.factor(as.numeric(covs_forremBatch$Concentration))
covs_forremBatch$Experiment <- as.factor(as.numeric(covs_forremBatch$Experiment))
covs_forremBatch$DIV <- as.factor(as.numeric(covs_forremBatch$DIV))
covs_forremBatch$`Cell line` <- as.factor(as.numeric(covs_forremBatch$`Cell line`))
covs_forremBatch$`Replicate` <- as.factor(as.numeric(covs_forremBatch$`Replicate`))
covs_forremBatch$`Batch` <- as.factor(as.numeric(covs_forremBatch$`Batch`))
covs_forremBatch <- covs_forremBatch[!colnames(covs_forremBatch) %in% c("Sequencing batch","Timepoint","DIV")]

batch.rem <- removeBatchEffect(as.matrix(assay(vsd.complex)), batch=covs_forremBatch$`Experiment`, 
                               covariates=as.matrix(cbind(covs_forremBatch$`Log(library size)`,
                                                          covs_forremBatch$`Percent. ERCC`
                                                          )),
                               design=model.matrix(~ covs_forremBatch$Stimulation))

#====PC loadings for covariates====#
#500x500
covs_forremBatch <- covs_forremBatch[colnames(covs_forremBatch) %in% colnames(cov_for.cor.sub)]
D100.PCAcovar.postcorrection <- PCA.covar(batch.rem,covs_forremBatch,data.r=T,type="pearson",heatmap.r=T)

#PCA (550 x 450)
pca.cor <- plotPCA.custom(as.matrix(batch.rem), intgroup=c("Concentration","Stimulation","Cell line","Experiment"),
                          ntop = 50000, returnData=TRUE, metadata=subset.meta, pc.1=1,pc.2=2)
PoV.cor <- round(100 * attr(pca.cor, "percentVar"))

#600x450
PCAplot(pca.cor, "Concentration", Shape="Cell line",PoV.df=PoV.cor, colors=c('#2E7281','black','#50B8CF'))
PCAplot(pca.cor, "Experiment", PoV.df=PoV.cor, colors=viridis(5))


#Inter sample distances post correction
sampleDists.postcor <- dist(t(batch.rem))
sampleDistMatrix.postcor <- as.matrix(sampleDists.postcor)
rownames(sampleDistMatrix.postcor) <- colnames(vsd.complex)
colnames(sampleDistMatrix.postcor) <- NULL

#Heatmap intersample-distance
#650x500
pheatmap(sampleDistMatrix.postcor, show_rownames = F,
         clustering_distance_rows=sampleDists.postcor,
         clustering_distance_cols=sampleDists.postcor,
         col=colors, annotation_row = annot_df, annotation_colors = annot_colors)



#==========================================================================================================#
#                                            RNA-seq analysis: DESeq2                                      #
#==========================================================================================================#
#dds.base <- DESeq(dds.base)
dds.complex <- DESeq(dds.complex)

resultsNames(dds.complex) # lists the coefficients

#res.base <- results(dds.base)
res.complex <- results(dds.complex, name="Stimulation_2_vs_1")

#summary(res.base)
summary(res.complex)
#sum(res.base$padj < 0.05, na.rm=TRUE)
sum(res.complex$padj < 0.01, na.rm=TRUE)
sum(res.complex$log2FoldChange < -2 & res.complex$padj < 0.05, res.complex$log2FoldChange > 2 & res.complex$padj < 0.05, na.rm=TRUE)

#top20upreg.base <- data.frame(data.frame(res.base)[order(res.base$log2FoldChange, decreasing=T),])[1:20,]
#top20downreg.base <- data.frame(data.frame(res.base)[order(res.base$log2FoldChange, decreasing=F),])[1:20,]


#Grab significant downreg & highest lfc genes
res.complex[is.na(res.complex$padj),]$padj <- 1
downreg.genes <- rownames(res.complex[res.complex$log2FoldChange < -2 & res.complex$padj < 0.05,])
upreg.genes <- rownames(res.complex[res.complex$log2FoldChange > 2 & res.complex$padj < 0.05,])


#Only most significant genes
write.xlsx2(res.complex[res.complex$log2FoldChange < 0 & res.complex$padj < 0.05,], 
            file = "04012021_D100_DEGs_down.xlsx", sheetName = "Down",
            col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx2(res.complex[res.complex$log2FoldChange > 0 & res.complex$padj < 0.05,], 
            file = "04012021_D100_DEGs_up.xlsx", sheetName = "Up",
            col.names = TRUE, row.names = TRUE, append = FALSE)


#Plot the genes
#heatmap based on sorted samples for concentration:
fact.to.sort <- as.numeric(annot_df$Concentration)
sorted <- batch.rem[,order(fact.to.sort, decreasing=F)]
sub <- sorted[rownames(sorted) %in% c(upreg.genes,downreg.genes),]
#550x1300
pdf(file="04022021_D100_heatmap_all DEGs_new res.pdf", width=5, height=12)

pheatmap(sub, 
         cluster_rows = T, 
         cluster_cols = F, annotation_legend = T, show_colnames = F,
         annotation_col=annot_df[-2], fontsize_row = 2,
         color = heatmap.color.code, scale='row', annotation_colors = annot_colors)

dev.off()





res.sorted.up <- res.complex[order(res.complex$log2FoldChange, decreasing = T),]
res.sorted.up <- res.sorted.up[res.sorted.up$log2FoldChange > 2 & res.sorted.up$padj < 0.05,]


res.sorted.down <- res.complex[order(res.complex$log2FoldChange, decreasing = F),]
res.sorted.down <- res.sorted.down[res.sorted.down$log2FoldChange < -2 & res.sorted.down$padj < 0.05,]
top.downreg.genes <- rownames(res.sorted.down)[1:15]
top.upreg.genes <- rownames(res.sorted.up)[1:15]


top.sub <- sorted[rownames(sorted) %in% c(top.upreg.genes,top.downreg.genes),]
#550x600
pheatmap(top.sub, 
         cluster_rows = T, 
         cluster_cols = F, annotation_legend = T, show_colnames = F,
         annotation_col=annot_df[-2], fontsize_row = 7,
         color = heatmap.color.code, scale='row', annotation_colors = annot_colors)



#Volcano plot of the results
library(EnhancedVolcano)
#700x750
#EnhancedVolcano(res.complex,
#                lab = rownames(res.complex),
#                x = 'log2FoldChange',
#                y = 'pvalue',FCcutoff = 2)

EnhancedVolcano(res.complex,
                lab = rownames(res.complex),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 2,
                labSize = 6,
                legendPosition = 'right')


res.complex[rownames(res.complex) %in% "SERPINA3",]



#Boxplots of DEGs

#Or manually:

#Dotplot
box <- data.frame(scale(t(sub)), Concentration=annot_df[order(fact.to.sort, decreasing=F),][1])
#box.chr <- data.frame(t(sub.ac), Stimulation=sub.meta.ac$Stimulation)
box.up <- box[colnames(box) %in% c(upreg.genes,"Concentration")]
box.up$Concentration <- factor(box.up$Concentration, levels=c("CTR","10mM","50mM"))
#Boxplots
library(ggplot2)

all.p <- NULL

for (i in upreg.genes){
  n <- as.factor(i)
  p = 
    ggplot(box.up, aes_string(x="Concentration", y=i, fill="Concentration")) +
    geom_boxplot() +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
    #  theme_ipsum() +
    
    theme(
      legend.position="Type",
      plot.title = element_text(size=11)
    ) +
    #  ggtitle("Boxplots showing log2CPM values of SPP1 across cell-type samples") +
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5))+
    #  xlab("Cell type")+
    ylab(paste(i,"(scaled counts)"))+
    #  theme(axis.text.x = element_text(angle = 90))+
    theme(axis.text.x=element_blank(), axis.title.x=element_blank())
  #  theme(legend.position = "none")
  all.p[[i]] <- p}


library(ggpubr)
#700x550
all.up <- ggarrange(plotlist=all.p[1:length(all.p)], common.legend = T)



#==========================================================================================================#
#                                              RNA-seq analysis: boxplots                                  #
#==========================================================================================================#
#Boxplot for single gene per stimulation

gene.box <- data.frame(t(batch.rem), Stimulation=subset.meta$Stimulation, `Concentration`=annot_df$Concentration, `Cell line`=annot_df$`Cell line`, check.names=F)

ggplot(gene.box[colnames(gene.box) %in% c("SERPINA3","Concentration")], aes_string(x="Concentration", y="SERPINA3", fill="Concentration")) +
  geom_boxplot() +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  #  theme_ipsum() +
  
  theme(
    legend.position="Type",
    plot.title = element_text(size=11)
  ) +
  #  ggtitle("Boxplots showing log2CPM values of SPP1 across cell-type samples") +
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))+
  #  xlab("Cell type")+
  ylab(paste("SERPINA3","(scaled corrected counts)"))+
  #  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())



#From all time points and per cell line:
box <- data.frame(scale(t(gfilt.df)),Stimulation = filt.meta$Stimulation ,Concentration = filt.meta$Concentration,`DIV`=filt.meta$DIV, `Cell line`=filt.meta$`Cell line`, check.names=F)
D100.box <- box[rownames(box) %in% c(colnames(subset.df),"Concentration"),]
gene1 <- D100.box[colnames(D100.box) %in% c("SERPINA3","Concentration","Cell line")]
ggplot(gene1, aes_string(x="Concentration", y="SERPINA3", fill="Concentration")) +
  geom_boxplot() +
    scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  #  theme_ipsum() +
  geom_point(aes(shape=`Cell line`))+
  theme(
    legend.position="Type",
    plot.title = element_text(size=11)
  ) +
  #  ggtitle("Boxplots showing log2CPM values of SPP1 across cell-type samples") +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  #  xlab("Cell type")+
  ylab(paste("SERPINA3","(raw counts)"))+
  #  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())


#And with classic theme
ggplot(MTOR, aes_string(x="DIV", y="SERPINA3")) +
  geom_boxplot() +
  #  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+
  #  theme_ipsum() +
  geom_point(aes(shape=`Cell line`))+
  theme(
    legend.position="Type",
    plot.title = element_text(size=11)
  ) +
  #  ggtitle("Boxplots showing log2CPM values of SPP1 across cell-type samples") +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("mTOR","(scaled counts)"))+
  #  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_blank())



#==========================================================================================================#
#                                         RNA-seq analysis: extracting data                                #
#==========================================================================================================#
#Extracting the results for all genes:
write.xlsx2(res.complex, 
            file = "04012021_D100_all_DEGs.xlsx", sheetName = "D100 all",
            col.names = TRUE, row.names = TRUE, append = FALSE, check.names=F)

options(java.parameters = "-Xmx8000m")
#Extracting the normalized values
write.xlsx2(data.frame(batch.rem,check.names=F), 
            file = "04012021_D100_vst_corrected.xlsx", sheetName = "vst-normalized & corrected",
            col.names = TRUE, row.names = TRUE, append = FALSE)


#For mTOR & P70S6K (RPS6KB1)
#Raw
write.xlsx2(data.frame(batch.rem,check.names=F)[rownames(batch.rem) %in% c("MTOR", "RPS6KB1"),], 
            file = "04012021_D100_mTOR_vst_corrected.xlsx", sheetName = "mTOR & P70S6K vst-norm and corrected",
            col.names = TRUE, row.names = TRUE, append = FALSE)

#Filtered
write.xlsx2(data.frame(gfilt.df,check.names=F)[rownames(gfilt.df) %in% c("MTOR", "RPS6KB1"),], 
            file = "04012021_D100_mTOR_filtered_raw.xlsx", sheetName = "mTOR & P70S6K raw",
            col.names = TRUE, row.names = TRUE, append = FALSE)


setwd()
dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", 
                 "Internship", "Internship_Start", "Side projects", 
                 "Amber bulk RNA-seq of stimulated organoids", "04022021_updates","D35")
setwd(dir)
norm.saved <- read.table("04012021_D35_vst normalized_corrected.xlsx", check.names=F)


dir <- file.path("C:", "Users", "rapha", "Dokumente", "UM FN", 
                 "Internship", "Internship_Start", "Side projects", 
                 "Amber bulk RNA-seq of stimulated organoids", "01042021_manuscript figures","Data","D100")

setwd(dir)
norm.saved <- read.xlsx("04012021_D100_vst_corrected.xlsx",sheetIndex = 1, as.data.frame=T)



