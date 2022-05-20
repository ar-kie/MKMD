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

#PCA
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
  intgroup.df <- data.frame(metadata[intgroup],check.names=F)
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    metadata[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc.1], PC2 = pca$x[, pc.2], group = group, 
                  intgroup.df, name = colnames(object), check.names=F)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(pc.1,pc.2)]
    return(d)
  }
  ggplot(data = d, aes_string(x = paste("PC",pc.1), y = paste("PC",pc.2), 
                              color = "group")) + geom_point(size = 3) + xlab(paste0(paste("PC", pc.1, ": ", sep=""),
                                                                                     round(percentVar[1] * 100), "% variance")) + ylab(paste0(paste("PC", pc2, ": ", sep=""), 
                                                                                                                                              round(percentVar[2] * 100), "% variance")) + coord_fixed()
}

#Plot PCA (ggplot)
PCAplot <- function(pca.df, Condition, Shape=NULL, pc.1=1, pc.2=2, PoV.df, colors=c("#00AFBB", "#E7B800", "#FC4E07")){
  if(!is.null(Shape))
    ggplot(pca.df, aes(PC1, PC2, color=pca.df[[c(Condition)]], label=name)) +
      geom_point(size=2, aes(shape=pca.df[[Shape]]))+
      labs(x = paste0("PC", paste(pc.1), ": ", PoV.df[1], "% variance"), y = paste0("PC", paste(pc.2), ": ",PoV.df[2],"% variance"), 
         color = Condition, shape=Shape) +  
      scale_color_manual(values = colors)+
      theme_bw()+
      theme(plot.title = element_text(hjust=0.5))
  else
    ggplot(pca.df, aes(PC1, PC2, color=pca.df[[c(Condition)]], label=name)) +
      geom_point(size=2)+
      labs(x = paste0("PC", paste(pc.1), ": ", PoV.df[1], "% variance"), y = paste0("PC", paste(pc.2), ": ",PoV.df[2],"% variance"), 
         color = Condition) +  
      scale_color_manual(values = colors)+
      theme_bw()+
      theme(plot.title = element_text(hjust=0.5))
}



#Density plot
plot.density <- function(meta.df,to.color="Condition",to.plot="Library.size",ylab.title="Density",xlab.title="Library size"){
  meta.df %>% 
    ggplot(aes_string(color=to.color, x=to.plot, fill=to.color)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab(ylab.title) +
    xlab(xlab.title)
}

#Histogram
plot.hist <- function(data.to.plot,log10ofdat=T,xlab.title="Data",
                      main.title="Counts of X",ylab.title="Number of samples",color="grey80",
                      thresh1=1,thresh2=10){
  if (log10ofdat)
    log.data <-log10(data.to.plot)
  hist(log.data, xlab=xlab.title, 
       main=paste(main.title, "mean = ", round(mean(data.to.plot),2), ";",
                  "median = ", round(median(data.to.plot),2)), 
       breaks=20, col=color, ylab=ylab.title)
  abline(v = thresh1, col="red", lwd=3, lty=2)
  abline(v = thresh2, col="red", lwd=3, lty=2)
}

#Function for covar-PC correlation
#function for PCA
PCA.covar <- function(df, meta, heatmap.r=F, data.r=F, heatmap.p=T,type="spearman"){
  l <- length(colnames(meta))
  l2 <- l+1
  pca.cov <- prcomp(t(df))
  PC.covs <- pca.cov$x[,1:l]
  PC_cov_correlation <- rcorr(as.matrix(PC.covs), as.matrix(meta), type=type)
  PC_variance.explained <- PC_cov_correlation$r[1:l,l2:length(rownames(PC_cov_correlation$r))]
  PC_cov.cor.p <- PC_cov_correlation$P[1:l,l2:length(rownames(PC_cov_correlation$r))]
  if (heatmap.r)
    pheatmap(PC_variance.explained, cluster_rows = F, cluster_cols = F, show_colnames = T,
             color = heatmap.color.code)
  if (heatmap.p)
    pheatmap(PC_cov.cor.p, cluster_rows = F, cluster_cols = F, show_colnames = T,
             color = colorRampPalette(c("red", "white"))(50),
             display_numbers = T)
  if (data.r)
    return (PC.covs)
}

createDT <- function(DF, caption="", scrollY=500){
  data <- DT::datatable(DF, caption=caption,
                        extensions = 'Buttons',
                        options = list( dom = 'Bfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        scrollY = scrollY, scrollX=T, scrollCollapse = T, paging = F,
                                        columnDefs = list(list(className = 'dt-center', targets = "_all"))
                        )
  )
  return(data)
}

PCA_cov_cor_R <- function(cov.df,df){
  l <- length(colnames(df))
  covariates <- colnames(cov.df)
  pca.cov <- prcomp(df)
  var <- get_pca_var(pca.cov) # description of PCs
  ind <- get_pca_ind(pca.cov) # PCs for individuals
  matrix_rsquared = matrix(NA, nrow = length(cov.df), ncol = l) #Number of factors
  matrix_pvalue = matrix(NA, nrow = length(cov.df), ncol = l)
  
  for (x in 1:length(covariates)){
    for (y in 1:l){
      matrix_rsquared[x,y] <- summary( lm(var$cor[,y] ~ cov.df[,covariates[x]]) )$adj.r.squared
      matrix_pvalue[x,y] <- glance( lm(var$cor[,y] ~ cov.df[,covariates[x]]) )$p.value #To insert pvalues in the heatmap
    }
  }
  
  rownames(matrix_rsquared) = covariates
  rownames(matrix_pvalue) = covariates 
  colnames(matrix_rsquared) = paste("PC",1:l, sep="")
  
  matrix_pvalue = matrix(p.adjust(as.vector(as.matrix(matrix_pvalue)), method='bonferroni'),ncol=ncol(matrix_pvalue))
  matrix_pvalue = formatC(matrix_pvalue, format = "e", digits = 2)
  
  # png(paste0(work_dir, "LinearReg_15pcs_covariates.png"), width = 10, height = 10, res = 600, units = "in")
  # pdf(paste0(work_dir, "LinearReg_15pcs_covariates.pdf"), width = 10, height = 10)
  pheatmap(matrix_rsquared, cluster_rows = T, cluster_cols = F, show_colnames = T, show_rownames = T,
           color = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30),
           display_numbers = F)
}



#==========================================================================================================#
#                                   RNA-seq analysis start: counts files                                   #
#==========================================================================================================#
dir <- file.path("D:", "Mount Sinai", "de Witte lab", "Side projects", "Amber bulk RNA-seq of stimulated organoids", "Counts")
setwd(dir)

b1 <- read.table("UMC-AB-b002 - Batch1.csv", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
b2 <- read.table(file = "UMC-AB-b002.ReadCounts.tsv", sep = '\t', header = TRUE, row.names=1, check.names = F)

#Merge dataframes but drop NAs
df <- merge(b1, b2, by=0, all=FALSE)
#Set gene names as index
rownames(df) <- df[,1]
#Remove gene names from df
df <- df[,-1]

rownames(df) <- make.unique(gsub("__.*","",rownames(df)))


#==========================================================================================================#
#                                      RNA-seq analysis start: covariate file                              #
#==========================================================================================================#
dir <- file.path("D:", "Mount Sinai", "de Witte lab", "Side projects", "Amber bulk RNA-seq of stimulated organoids", "Metadata")
setwd(dir)
#Reading in metadata
metadata.b1 <- read.xlsx("Metadata_MKMD_RNAseq_Batch1.xlsx", sheetIndex=1, as.data.frame=T)
metadata.b2 <- read.xlsx("Metadata_MKMD_RNAseq_Batch2.xlsx", sheetIndex=1, as.data.frame=T)

colnames(metadata.b1) <- colnames(metadata.b2)
meta.df <- rbind(metadata.b1, metadata.b2)
colnames(meta.df) <- sub("\\.", " ",colnames(meta.df))
colnames(meta.df) <- sub("\\.", "",colnames(meta.df))
meta.df <- meta.df[!colnames(meta.df) %in% c("QC", "Sequencing lane", "Low concentration")]
meta.df$`Sample Name` <- as.character(meta.df$`Sample Name`)
meta.df$Stimulation <- NA
meta.df[rownames(meta.df[meta.df$Concentration == "CTR",]),]$Stimulation <- 1
meta.df[rownames(meta.df[meta.df$Concentration != "CTR",]),]$Stimulation <- 2
meta.df$Stimulation <- as.factor(meta.df$Stimulation)
meta.df$Replicate <- as.factor(meta.df$Duplicate)
meta.df$Batch <- as.factor(meta.df$Experiment)

#colnames(meta.df[colnames(meta.df) %in% "Experiment"]) <- "`Sequencing batch`"

for (i in rownames(meta.df)){
  if (meta.df[i,]$Concentration == "10 mM")
    meta.df[i,]$Concentration <- "10mM"
  if (meta.df[i,]$Concentration == "50 mM")
    meta.df[i,]$Concentration <- "50mM"
}

meta.df$Concentration <- as.factor(as.character(meta.df$Concentration))


#Creating a new metadata file with additional covariates
new.meta <- data.frame("Library size" = colSums(df), check.names = F)
rownames(new.meta) <- colnames(df)
new.meta$`Sequencing batch` <- NA
new.meta[rownames(new.meta) %in% colnames(b1),]$`Sequencing batch` <- "1"
new.meta[rownames(new.meta) %in% colnames(b2),]$`Sequencing batch` <- "2"
new.meta$`Percent. ERCC` <- colSums(x = df[grep(pattern = "^ERCC-", x = rownames(df), value = TRUE), 
                                        drop = FALSE,])/colSums(df) * 100
new.meta$`Percent. MT` <- colSums(x = df[grep(pattern = "^MT-", x = rownames(df), value = TRUE), 
                                      drop = FALSE,])/colSums(df) * 100
identical(rownames(new.meta), colnames(df))
#Parsing sample names
meta.df$`Sample Name`[is.na(meta.df$`Sample Name`)] <- as.character(rownames(new.meta[new.meta$`Sequencing batch` == 2,]))

meta.df.all <- merge(meta.df, new.meta, by.x=1, by.y=0)
rownames(meta.df.all) <- meta.df.all$`Sample Name`

meta.df.all$Concentration <- relevel(meta.df.all$Concentration,"CTR")



#==========================================================================================================#
#                             RNA-seq analysis start: doublet comparison & removal                         #
#==========================================================================================================#
#Doublets are: 24, 25, 31, 32

#Subset:
meta.doub <- meta.df.all[meta.df.all$'Sample NR' %in% c(24, 25, 31, 32),]

write.table(meta.doub, "D100 sample doublets.txt")

summary(meta.df.all$DIV)



#==========================================================================================================#
#            RNA-seq analysis start: failed sample removal & new df & covariate file creation              #
#==========================================================================================================#
remove.samples <- c(rownames(meta.df.all[meta.df.all$`Sample Name` %in% rownames(meta.doub) & meta.df.all$`Sequencing batch` == 2,]), 
                    rownames(meta.df.all[meta.df.all$'Sample NR' == "41" & meta.df.all$`Sequencing batch` == 1,]))
df.clean <- df[!colnames(df) %in% remove.samples]
meta.clean <- meta.df.all[!rownames(meta.df.all) %in% remove.samples,]
meta.clean$Experiment <- as.factor(as.numeric(as.factor(paste(meta.clean$Experiment,meta.clean$'Cell line', sep="_"))))

meta.clean <- meta.clean[!colnames(meta.clean) %in% c("Sample NR","Duplicate")]
meta.clean$DIV <- as.character(meta.clean$DIV)
meta.clean[meta.clean$DIV == "Day 24",]$DIV <- "Day 20"
meta.clean[meta.clean$DIV == "Day 37",]$DIV <- "Day 35"
meta.clean[meta.clean$DIV %in% c("Day 104","Day 106","Day 107"),]$DIV <- "Day 100"

#Change to weeks
meta.clean[meta.clean$DIV == "Day 20",]$DIV <- "W4"
meta.clean[meta.clean$DIV == "Day 35",]$DIV <- "W5"
meta.clean[meta.clean$DIV == "Day 100",]$DIV <- "W15"
meta.clean$DIV <- as.factor(meta.clean$DIV)
meta.clean$DIV  <- relevel(meta.clean$DIV, "W5")
meta.clean$DIV  <- relevel(meta.clean$DIV, "W4")



#==========================================================================================================#
#                                           RNA-seq analysis start: QC on all data                           #
#==========================================================================================================#
#Library size all samples
plot.hist(meta.clean$'Library size', xlab.title="Library sizes (log10 counts)",
          main.title="Counts sequencing batch 1 & 2 \n\n",ylab.title="Number of samples",
          thresh1=log10(500000),thresh2=log10(15000000))

plot.density(meta.clean,to.color="Concentration",to.plot="log10(`Library size`)")


#ERCC all samples
hist(log10(meta.clean$`Percent. ERCC`), xlab="ERCC SI percentage (log10)", 
     main=paste("ERCC SIs sequencing batch 1 & 2 \n\n", "mean = ", round(mean(meta.clean$`Percent. ERCC`),2), ";",
                "median = ", round(median(meta.clean$`Percent. ERCC`),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

plot.density(meta.clean,to.color="Concentration",to.plot="`Percent. ERCC`",xlab.title="ERCC proportion")

plot.density(meta.clean,to.color="`Sequencing batch`",to.plot="`Percent. ERCC`",xlab.title="ERCC proportion")


barplot(meta.clean$`Percent. ERCC`,
        main=paste("ERCC SIs sequencing batch 1 & 2 \n\n", "mean = ", round(mean(meta.clean$`Percent. ERCC`),2), ";",
                   "median = ", round(median(meta.clean$`Percent. ERCC`),2)), 
        col="grey80", ylab="ERCC SI percentage")


#==========================================================================================================#
#                                     RNA-seq analysis start: QC on `Sequencing batch` subsets                          #
#==========================================================================================================#


#Library size b1 vs b2
hist(log10(meta.clean[meta.clean$`Sequencing batch` == 1,]$'Library size'), xlab="Library sizes (log10 counts)", 
     main=paste("Counts sequencing batch 1\n\n", "mean = ", round(mean(meta.clean[meta.clean$`Sequencing batch` == 1,]$'Library size'),2), ";",
                "median = ", round(median(meta.clean[meta.clean$`Sequencing batch` == 1,]$'Library size'),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

hist(log10(meta.clean[meta.clean$`Sequencing batch` == 2,]$'Library size'), xlab="Library sizes (log10 counts)", 
     main=paste("Counts sequencing batch 2\n\n", "mean = ", round(mean(meta.clean[meta.clean$`Sequencing batch` == 2,]$'Library size'),2), ";",
                "median = ", round(median(meta.clean[meta.clean$`Sequencing batch` == 2,]$'Library size'),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

#Filtering
#Based on counts
#filt.meta <- meta.clean[meta.clean$'Library size' > 500000 & meta.clean$'Library size' < 15000000,] #remove large library size-samples
filt.meta <- meta.clean[meta.clean$'Library size' > 500000,]
filt.df <- df.clean[colnames(df.clean) %in% rownames(filt.meta)]

#Mapping the library size distributions
hist(log10(filt.meta$'Library size'), xlab="Library sizes (log10 counts)", 
     main=paste("Counts sequencing batch 1 & 2 (filtered) \n\n", "mean = ", round(mean(filt.meta$'Library size'),2), ";",
                "median = ", round(median(filt.meta$'Library size'),2)), 
     breaks=20, col="grey80", ylab="Number of samples")


#Library sizes b1 vs b2 (after sample filtering)
hist(log10(filt.meta[filt.meta$`Sequencing batch` == 1,]$'Library size'), xlab="Library sizes (log10 counts)", 
     main=paste("Counts sequencing batch 1\n\n", "mean = ", round(mean(filt.meta[filt.meta$`Sequencing batch` == 1,]$'Library size'),2), ";",
                "median = ", round(median(filt.meta[filt.meta$`Sequencing batch` == 1,]$'Library size'),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

hist(log10(filt.meta[filt.meta$`Sequencing batch` == 2,]$'Library size'), xlab="Library sizes (log10 counts)", 
     main=paste("Counts sequencing batch 2\n\n", "mean = ", round(mean(filt.meta[filt.meta$`Sequencing batch` == 2,]$'Library size'),2), ";",
                "median = ", round(median(filt.meta[filt.meta$`Sequencing batch` == 2,]$'Library size'),2)), 
     breaks=20, col="grey80", ylab="Number of samples")


#Based on mtRNA
#mtRNA percentage
hist(filt.meta$`Percent. MT`, xlab="Mitochondrial proportion", 
     main=paste("Mitochondrial percentage \n\n", "mean = ", round(mean(filt.meta$`Percent. MT` ),2), ";",
                "median = ", round(median(filt.meta$`Percent. MT` ),2)), 
     breaks=20, col="grey80", ylab="Number of Samples")
abline(v = 37, col="red", lwd=3, lty=2)

plot.density(filt.meta,to.color="Concentration",to.plot="`Percent. MT`",xlab.title="mtRNA percentage")

filt.meta <- filt.meta[filt.meta$`Percent. MT` < 37,]
filt.df <- df.clean[colnames(df.clean) %in% rownames(filt.meta)]

#Post-filter distribution of mtRNA percentage
hist(filt.meta$`Percent. MT`, xlab="Mitochondrial proportion", 
     main=paste("Mitochondrial percentage \n\n", "mean = ", round(mean(filt.meta$`Percent. MT` ),2), ";",
                "median = ", round(median(filt.meta$`Percent. MT` ),2)), 
     breaks=20, col="grey80", ylab="Number of Samples")

#550 x 450
plot.density(filt.meta,to.color="Concentration",to.plot="`Percent. MT`",xlab.title="mtRNA percentage")
plot.density(filt.meta,to.color="Concentration",to.plot="`Percent. ERCC`",xlab.title="ERCC percentage")
plot.density(filt.meta,to.color="Concentration",to.plot="`Library size`",xlab.title="Library size")


#Outlier detection based on IQR correlation
sOutliers <- IQRoutlier.detector(filt.df,barplot=T)
library(outliers)
library(flashClust)
library(CoExpNets)

outlist = NULL
#Traspose data frame and get a distance matrix between pairs of samples
sampdist = as.matrix(dist(CoExpNets::trasposeDataFrame(filt.df,F)))
result = grubbs.test(apply(sampdist,1,sum))
#Is the test significant at 95% confidence level?
plot(flashClust(dist(CoExpNets::trasposeDataFrame(filt.df,F)), method = "average"),cex=0.5,
     main=paste0("Detected outlier ",names(result$p.value)," in MKMD STIM/CTR"))

s.filt <- filt.df[!colnames(filt.df) %in% c("069","077","005")]
filt.meta <- filt.meta[rownames(filt.meta) %in% colnames(s.filt),]


sOutliers <- IQRoutlier.detector(s.filt,barplot=T)

filt.meta$'Cell line' <- as.factor(as.numeric(filt.meta$'Cell line'))

pca <- plotPCA.custom(as.matrix(s.filt), intgroup=c("Cell line", "Concentration", "Experiment","DIV"), ntop = 20000, returnData=TRUE, 
                      metadata=filt.meta)

PoV <- round(100 * attr(pca, "percentVar"))
PCAplot(pca, "DIV", PoV.df=PoV)
PCAplot(pca, "Experiment", PoV.df=PoV, colors=viridis(7))



#==========================================================================================================#
#                                     RNA-seq analysis start: filtering                                    #
#==========================================================================================================#
library(edgeR)
#cpm <- data.frame(cpm(filt.df))
#L <- mean(filt.meta$'Library size') * 1e-6
#M <- median(filt.meta$'Library size') * 1e-6

#One might filter based on cpm
#gfilt.df <- data.frame(filt.df[rowSums(cpm) >= min(summary(meta.clean$Concentration))/L,]) #Filter based on min rowSums(CPM) to ignore library size 

#Or
#gfilt.df <- filt.df[filterByExpr(filt.df,group=filt.meta$Concentration,min.prop=0.5),] #Filters based on CPM that at least X samples have X CPM


#Or the normal procedure:
#gfilt.df <- filt.df[rowSums(filt.df) >= length(colnames(filt.df)),] #Also removes ERCC SIs


#Filter based on read counts minimum per condition
#gfilt.df <- filt.df[rowSums(filt.df) >= length(min(summary(meta.clean$Concentration))),] #Also removes ERCC SIs


#Favorite filtering:
#g.filter <- rowSums(df) < quantile(rowSums(df), 0.95) & 
#  rowSums(df) > quantile(rowSums(df), 0.15)
g.filter <- rowSums(s.filt) > quantile(rowSums(s.filt), 0.15)
gfilt.df <- s.filt[g.filter,]



#==========================================================================================================#
#                              RNA-seq analysis start: post-filtering QC                                   #
#==========================================================================================================#
sOutliers.filt <- IQRoutlier.detector(gfilt.df,barplot=T)
#gOutliers.filt <- IQRoutlier.detector(t(gfilt.df), barplot=T)

gfilt.df <- gfilt.df[sort(colnames(gfilt.df),decreasing=T)]

filt.meta <- filt.meta[sort(rownames(filt.meta), decreasing=T),]
#Sanity check:
identical(colnames(gfilt.df),rownames(filt.meta))

#'Library size'
filt.meta$'Library size post filter' <- colSums(gfilt.df)

#Distribution
hist(log10(filt.meta$'Library size post filter'), xlab="Library sizes (log10 counts)", 
     main=paste("Counts \n\n", "mean = ", round(mean(filt.meta$'Library size post filter'),2), ";",
                "median = ", round(median(filt.meta$'Library size post filter'),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

#Library sizes b1 vs b2 (after gene filtering)
hist(log10(filt.meta[filt.meta$`Sequencing batch` == 1,]$'Library size post filter'), xlab="Library sizes (log10 counts)", 
     main=paste("Counts sequencing batch 1 post-filter \n\n", "mean = ", round(mean(filt.meta[filt.meta$`Sequencing batch` == 1,]$'Library size'),2), ";",
                "median = ", round(median(filt.meta[filt.meta$`Sequencing batch` == 1,]$'Library size post filter'),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

hist(log10(filt.meta[filt.meta$`Sequencing batch` == 2,]$'Library size post filter'), xlab="Library sizes (log10 counts)", 
     main=paste("Counts sequencing batch 2 post-filter \n\n", "mean = ", round(mean(filt.meta[filt.meta$`Sequencing batch` == 2,]$'Library size'),2), ";",
                "median = ", round(median(filt.meta[filt.meta$`Sequencing batch` == 2,]$'Library size post filter'),2)), 
     breaks=20, col="grey80", ylab="Number of samples")

plot.density(filt.meta,to.color="Concentration",to.plot="log10(`Library size`)")

#PCA
pca.filt <- plotPCA.custom(as.matrix(gfilt.df), intgroup=c("Concentration","Experiment","DIV"), 
                           ntop = 50000, returnData=TRUE, metadata=filt.meta)
PoV.filt <- round(100 * attr(pca.filt, "percentVar"))
PCAplot(pca.filt, "Concentration", PoV.df=PoV.filt, colors=c('#2E7281','black','#50B8CF'))



#==========================================================================================================#
#                                RNA-seq analysis start: Covariate analysis                                #
#==========================================================================================================#
#Changing factors to factors
filt.meta$Replicate <- as.factor(filt.meta$Replicate)
filt.meta$Condition <- as.factor(filt.meta$Condition)
filt.meta_clean <- filt.meta[!colnames(filt.meta) %in% c("Library size post filter","Condition","Sample Name")]

filt.meta_clean$'Library size' <- log(filt.meta_clean$'Library size')

#Check which samples were removed along the way
samples.removed <- rownames(new.meta)[!rownames(new.meta) %in% rownames(filt.meta_clean)]
write.table(samples.removed, file="Samples removed.txt")


#====Colinearity analysis of covariates====#
heatmap.color.code <- rev(brewer.pal(11,"RdYlBu"))
#Check co-linearity of covariates for subsetted samples
cov_for.cor <- filt.meta_clean
cov_for.cor$DIV <- as.character(cov_for.cor$DIV)

cov_for.cor$Experiment <- as.factor(as.numeric(cov_for.cor$Experiment))
cov_for.cor$Concentration <- as.factor(as.numeric(cov_for.cor$Concentration))
cov_for.cor$DIV <- as.factor(as.numeric(as.factor(cov_for.cor$DIV)))
cov_for.cor$Replicate <- as.factor(as.numeric(as.factor(cov_for.cor$Replicate)))
cov_for.cor$Batch <- as.factor(as.numeric(as.factor(cov_for.cor$Batch)))
cov_for.cor$'Cell line' <- as.factor(as.numeric(cov_for.cor$'Cell line'))
covar.cor.all <- rcorr(as.matrix(cov_for.cor), type="pearson")
covar.cor <- covar.cor.all$r
pheatmap(covar.cor, color = heatmap.color.code)


#====PC loadings for covariates====#
identical(colnames(gfilt.df), rownames(cov_for.cor))
All.PCAcovar <- PCA.covar(gfilt.df,cov_for.cor,data.r=T)

library(factoextra)
library(broom)
PCA_cov_cor_R(cov_for.cor,gfilt.df)


#==========================================================================================================#
#                                RNA-seq analysis start: data correction                                   #
#==========================================================================================================#
#For correcting the data
filt.meta_clean$'Log(library size)' <- scale(filt.meta_clean$`Library size`)
filt.meta_clean$`Scaled percent. MT` <- scale(filt.meta_clean$"Percent. MT")

#Create a DESeq2 matrix
dds.complex <- DESeqDataSetFromMatrix(countData = as.matrix(gfilt.df),
                                      colData = data.frame(filt.meta_clean),
                                      design = ~ Log.library.size. + Cell.line + Batch + Sequencing.batch + Percent..ERCC + Stimulation)


#version 2 based on dds object
# Estimate library size correction scaling factors
dds.complex <- estimateSizeFactors(dds.complex)

vsd.complex <- vst(dds.complex)

#PCA (550 x 450)
pca.norm1 <- plotPCA.custom(as.matrix(assay(vsd.complex)), intgroup=c("Stimulation","Cell line","Experiment","Replicate"),
                            ntop = 50000, returnData=TRUE, metadata=filt.meta_clean, pc.1=1,pc.2=2)
PoV.norm1 <- round(100 * attr(pca.norm1, "percentVar"))


#600x450
PCAplot(pca.norm1, "Stimulation", Shape="Experiment",PoV.df=PoV.norm1)



#Regression without any SVs (data correction)

covs.rs_norm.norm <- cov_for.cor[colnames(cov_for.cor) %in% c("Log(library size)","Percent. ERCC","Cell line","Sequencing batch","Batch")]

resids.rembatch_norm.norm <- apply(assay(vsd.complex), 1, function(y){
  lm( y ~ . , data=cbind(covs.rs_norm.norm))$residuals
})


pca.resRem_norm.norm <- plotPCA.custom(as.matrix(t(resids.rembatch_norm.norm)), intgroup=c("Concentration","Cell line","DIV","Stimulation"),
                                  ntop = 50000, returnData=TRUE, metadata=filt.meta_clean, pc.1=1,pc.2=2)
PoV.resRem_norm.norm <- round(100 * attr(pca.resRem_norm.norm, "percentVar"))
PCAplot(pca.resRem_norm.norm, "DIV", Shape="Stimulation",PoV.df=PoV.resRem_norm.norm)
PCAplot(pca.resRem_norm.norm, "Stimulation",PoV.df=PoV.resRem_norm.norm)

PCA.covar(t(resids.rembatch_norm.norm),cov_for.cor,data.r=T,type="pearson",heatmap.r=T)
PCA_cov_cor_R(cov_for.cor,t(resids.rembatch_norm.norm))

residual.expr <- t(resids.rembatch_norm.norm)


save(residual.expr,file="MKMD_residual_expr_all.Rds")

