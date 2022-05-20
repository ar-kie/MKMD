#latest change: 02/24/2022 (V10.1 and figure updates)


#==========================================================================================================#
#                                              required packages                                           #
#==========================================================================================================#
library(xlsx)
library(Hmisc)
library(dplyr)
library(viridis)
library(pheatmap)
library(variancePartition)
library(RColorBrewer)
library(pcaExplorer)
library(broom)
library(CoExpNets)
library(clusterProfiler)
library(biomaRt)
library(readr)
library(DESeq2)
library(tidyr)
library(readxl)
library(ggplot2)
library(enrichplot)

#run multiple cores
library(doParallel)
library(foreach)
cluster = makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)
#installr::install.Rtools(check = TRUE, check_r_update = TRUE, GUI = TRUE)

dir <- file.path("D:", "R", "R workspace")
setwd(dir)

load("05212021_MKMD_base and correction all data.RData")
MKMD.networks <- readRDS("netMKMD.WGCNA.5.it.20.rds")

#What's the aim now?
#Check for mTOR / Pwhatever?
#Predict eigengenes based on mTOR expression?
#Check eigengene differential expression?
#What else? Prioritize stimulation modules and construct a figure story for it?
#What is the story we want to tell?
#Answer: this is the network we created. These are potentially the stimulation clusters: a) aassociated with stimulation, b) enriched for DEGs (cyan, tan)
#These modules carry these GO terms
#These modules are associated with these cell functions
#These modules fluctuate with expression mTOR
#Look at hub genes of these modules (i.e., cyan has RACK1 which regulates neurodevelopment)
#They interact with modules associated with XX (you now also have more gene lists)


#So the plan:
#Fig. 1A: ME-covariate correlations of modules (restrict to those that are significant),
#Fig. 1B: GSEA enrichment of DEGs (compound up & down for figure; barplot)
#Fig. 1C: barplot of cyan and tan MEs between CTR and stimulation (or 50mM and CTR)
#Fig. 1D: GO enrichment plots of tan, signature enrichment plot of cyan
#Fig. 1E: mTOR expression predicting MEs of relevant modules



#==========================================================================================================#
#                                                  new functions                                           #
#==========================================================================================================#

MEs_lin <- function(MEs,cov.df, return.p=F, return.r=F){
  l <- length(colnames(MEs))
  covariates <- colnames(cov.df)
  matrix_rsquared = matrix(NA, nrow = length(covariates), ncol = l) #Number of factors
  matrix_pvalue = matrix(NA, nrow = length(covariates), ncol = l)
  
  
  for (x in 1:length(covariates)){
    for (y in 1:l){
      matrix_rsquared[x,y] <- summary( lm(MEs[,y] ~ cov.df[,covariates[x]]) )$adj.r.squared
      matrix_pvalue[x,y] <- glance( lm(MEs[,y] ~ cov.df[,covariates[x]]) )$p.value #To insert pvalues in the heatmap
    }
  }
  
  rownames(matrix_rsquared) = covariates
  rownames(matrix_pvalue) = covariates 
  colnames(matrix_rsquared) = colnames(MEs)
  
  matrix_pvalue = matrix(p.adjust(as.vector(as.matrix(matrix_pvalue)), method='bonferroni'),ncol=ncol(matrix_pvalue))
  matrix_pvalue = formatC(matrix_pvalue, format = "e", digits = 2)
  
  matrix_rsquared[is.na(matrix_rsquared)] <- 0
  # png(paste0(work_dir, "LinearReg_15pcs_covariates.png"), width = 10, height = 10, res = 600, units = "in")
  # pdf(paste0(work_dir, "LinearReg_15pcs_covariates.pdf"), width = 10, height = 10)
  ComplexHeatmap::pheatmap(matrix_rsquared, cluster_rows = T, cluster_cols = T, annotation_legend = F, show_colnames = T,
                           legend_breaks= c(0,0.5,1,1), 
                           legend = T, color = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30),
                           legend_labels=c("0","0.5","1","Adjusted R? \n\n\n"),
           display_numbers = matrix_pvalue)
  if (return.p)
    return(matrix_pvalue)
  if (return.r)
    return(matrix_rsquared)
}

setEnrichment <- function(set1, set2, universe = 18000){
  a = sum(set1 %in% set2)
  c = length(set1) - a
  
  b = length(set2) - a
  d = universe - length(set2) - c
  
  contingency_table = matrix(c(a, c, b, d), nrow = 2)
  # one-tailed test for enrichment only
  fisher_results = fisher.test(contingency_table, alternative = "greater")
  # returns data.frame containing the lengths of the two sets, the overlap, the enrichment ratio (odds ratio) and P value
  df <- tibble::tibble( set1_length = length(set1), set2_length = length(set2), overlap = a, ratio = fisher_results$estimate, p.value = fisher_results$p.value)
  return(df)
}

GSEA.byMod <- function(mod.gl, gl.of.int., universe=18000){ #mod.gl = list of genes per module, gl.of.int = gene list to test enrichment for
  df <- NULL
  for (i in names(mod.gl)){
    df <- data.frame(rbind(df, setEnrichment(mod.gl[[i]], gl.of.int.)))}
  df$q.val <- p.adjust(as.vector(as.matrix(df$p.value)), method='fdr')
  rownames(df) <- names(mod.gl)
  return(df)
}


#==========================================================================================================#
#                                RNA-seq analysis start: data correction                                   #
#==========================================================================================================#
#Remove batch effects
#Create a meta df that corresponds to the sorted columns of the vsd assay but has covariates as numeric factors
covs_forremBatch <- filt.meta_clean
covs_forremBatch$Stimulation <- as.factor(as.numeric(covs_forremBatch$Stimulation))
covs_forremBatch$Concentration <- as.factor(as.numeric(covs_forremBatch$Concentration))
covs_forremBatch$Experiment <- as.factor(as.numeric(covs_forremBatch$Experiment))
covs_forremBatch$`Cell line` <- as.factor(as.numeric(covs_forremBatch$`Cell line`))
covs_forremBatch$`Batch` <- as.factor(as.numeric(covs_forremBatch$`Batch`))
covs_forremBatch$Replicate <- as.factor(as.numeric(as.factor(covs_forremBatch$Replicate)))
covs_forremBatch <- covs_forremBatch[!colnames(covs_forremBatch) %in% c("DIV","Library size","Percent. MT")]

batch.rem <- removeBatchEffect(as.matrix(assay(vsd.complex)), batch = covs_forremBatch$Batch, batch2 = covs_forremBatch$`Cell line`,
                               covariates=as.matrix(cbind(covs_forremBatch$`Log(library size)`,
                                                          covs_forremBatch$`Percent. ERCC`,
                                                          as.numeric(covs_forremBatch$`Sequencing batch`))),
                               design=model.matrix(~ covs_forremBatch$Timepoint + covs_forremBatch$`Scaled percent. MT` + covs_forremBatch$Stimulation))

#Subsetting
rv <- rowVars(batch.rem)
select <- order(rv, decreasing = TRUE)[seq_len(min(18000, 
                                                   length(rv)))]
sub.df <- t(batch.rem[select,])

#Create dummy variables for CTR-concentration effect:
covs_forremBatch$"10mM" <- 0
covs_forremBatch$"50mM" <- 0
covs_forremBatch$"CTR" <- 0
covs_forremBatch[covs_forremBatch$Concentration == 1,]$"CTR" <- 1
covs_forremBatch[covs_forremBatch$Concentration == 2,]$"10mM" <- 1
covs_forremBatch[covs_forremBatch$Concentration == 3,]$"50mM" <- 1


covs_forremBatch$`10mM` <- as.factor(covs_forremBatch$`10mM`)
covs_forremBatch$`50mM` <- as.factor(covs_forremBatch$`50mM`)
covs_forremBatch$`CTR` <- as.factor(covs_forremBatch$`CTR`)
#Create dummy variables for CTR-concentration effect:

covs_forremBatch$"W4" <- 0
covs_forremBatch$"W5" <- 0
covs_forremBatch$"W15" <- 0
covs_forremBatch[covs_forremBatch$Timepoint == 1,]$"W4" <- 1
covs_forremBatch[covs_forremBatch$Timepoint == 2,]$"W5" <- 1
covs_forremBatch[covs_forremBatch$Timepoint == 3,]$"W15" <- 1


covs_forremBatch$"CTR_W4" <- 0
covs_forremBatch$"CTR_W5" <- 0
covs_forremBatch$"CTR_W15" <- 0

covs_forremBatch$"10mM_W4" <- 0
covs_forremBatch$"10mM_W5" <- 0
covs_forremBatch$"10mM_W15" <- 0

covs_forremBatch$"50mM_W4" <- 0
covs_forremBatch$"50mM_W5" <- 0
covs_forremBatch$"50mM_W15" <- 0


covs_forremBatch[covs_forremBatch$Concentration == 1 & covs_forremBatch$Timepoint == 1,]$"CTR_W4" <- 1
covs_forremBatch[covs_forremBatch$Concentration == 1 & covs_forremBatch$Timepoint == 2,]$"CTR_W5" <- 1
covs_forremBatch[covs_forremBatch$Concentration == 1 & covs_forremBatch$Timepoint == 3,]$"CTR_W15" <- 1

covs_forremBatch[covs_forremBatch$Concentration == 2 & covs_forremBatch$Timepoint == 1,]$"10mM_W4" <- 1
covs_forremBatch[covs_forremBatch$Concentration == 2 & covs_forremBatch$Timepoint == 2,]$"10mM_W5" <- 1
covs_forremBatch[covs_forremBatch$Concentration == 2 & covs_forremBatch$Timepoint == 3,]$"10mM_W15" <- 1

covs_forremBatch[covs_forremBatch$Concentration == 3 & covs_forremBatch$Timepoint == 1,]$"50mM_W4" <- 1
covs_forremBatch[covs_forremBatch$Concentration == 3 & covs_forremBatch$Timepoint == 2,]$"50mM_W5" <- 1
covs_forremBatch[covs_forremBatch$Concentration == 3 & covs_forremBatch$Timepoint == 3,]$"50mM_W15" <- 1

for (i in colnames(covs_forremBatch))
  if (!i %in% c("Percent. ERCC", "Log(library size)", "Scaled percent. MT"))
    covs_forremBatch[[i]] <- as.factor(covs_forremBatch[[i]])



#==========================================================================================================#
#                                    Annotation: Loading signature lists                                   #
#==========================================================================================================#
#MKMD DEGs

setwd("D:/Mount Sinai/de Witte lab/Side projects/Amber bulk RNA-seq of stimulated organoids/02232021_DEGs")
file_list <- list.files(path="D:/Mount Sinai/de Witte lab/Side projects/Amber bulk RNA-seq of stimulated organoids/02232021_DEGs")
DEGs <- list()

for (i in file_list){
  temp_data <- read_excel(i)[1] #each file will be read in, specify which columns you need read in to avoid any errors
  temp_data$Class <- gsub(".xlsx", "", i) #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
  DEGs[[i]] <- make.unique(gsub("__.*","",temp_data$...1)) #for each iteration, bind the new data to the building dataset
}

names(DEGs) <- gsub(".xlsx*", "", file_list)
names(DEGs) <- gsub("02232021_*", "", names(DEGs))


#SCZ markers
setwd("D:/Mount Sinai/de Witte lab/Side projects/Amber bulk RNA-seq of stimulated organoids/WGCNA/gene lists")
SCZ_TWAS <- read_excel("TWAS_SCZ.xlsx")
SCZ.genes <- SCZ_TWAS$gene_name
SCZ.genes <- SCZ.genes[SCZ.genes %in% colnames(sub.df)]

#Cell-type markers
cell.markers <- data.frame(read_csv("ALL_Fidelity.csv"))
cell.markers <- cell.markers[cell.markers$Gene %in% colnames(sub.df),]
rownames(cell.markers) <- cell.markers$Gene
cell.markers <- cell.markers[,c(4:7)]

top.markers <- list()
for (i in names(cell.markers))
   top.markers[[i]] <- rownames(cell.markers[order(cell.markers[i], decreasing = T),][1:75,])

names(top.markers) <- c("AST", "OLG", "MIC", "NEU")

select_cell.markers <- list(microglia = c("P2RY12", "CX3CR1","AIF1", "IRF5", "IRF8","TMEM119", "SPI1","CD74",
                                     "CSF1R","C3","C1QC"),synaptic = c("SYT1",
                                                                       "SNAP25",
                                                                       "GRIN",
                                                                       "MAP2SST",
                                                                       "RELN"),exc.neuron = c("SLC17A7",
                                                                                              "CAMK2A",
                                                                                              "NRGN"),inh.neuron = c("GAD1","GAD2"),oligo = c("MOBP",
                                                                                                                                              "PLP1",
                                                                                                                                              "OLIG2",
                                                                                                                                              "MAG"),
                       opc = c("PDGFRA",
                               "VCAM","CSPG4"),
                       astro = c("AQP4",
                                 "GFAP",
                                 "SLC1A2", 
                                 "VIM"),endo = c("FL1",
                                                 "CLDN5"))

#Gandal markers
setwd("D:/Mount Sinai/de Witte lab/Side projects/Amber bulk RNA-seq of stimulated organoids/WGCNA/gene lists")
Gandal.clusters <- read_excel("Gandal (2018) - Gene and isoform co-expression module annotation.xlsx", sheet=2)
#Restrict to expressed genes only
Gandal.clusters.annot <- Gandal.clusters[colnames(Gandal.clusters) %in% c("ensembl_gene_id", "Module", "gene_name")]
Gandal.clusters.annot <- Gandal.clusters.annot[Gandal.clusters.annot$gene_name %in% colnames(sub.df),]

Gandal.clusters.annot[Gandal.clusters.annot$Module %in% "geneM6",]$Module <- "Microglia"
Gandal.clusters.annot[Gandal.clusters.annot$Module %in% "geneM3",]$Module <- "Astrocytes"
Gandal.clusters.annot[Gandal.clusters.annot$Module %in% "geneM5",]$Module <- "NfKB"
Gandal.clusters.annot[Gandal.clusters.annot$Module %in% "geneM32",]$Module <- "IFNy"

#Quadrato organoid markers
setwd("D:/Mount Sinai/de Witte lab/Side projects/Amber bulk RNA-seq of stimulated organoids/WGCNA/gene lists")
Quad.markers <- read_excel("Quadrato cluster genes.xlsx")


#Patir microglia
setwd("D:/Mount Sinai/de Witte lab/Side projects/PCA microglia-iPSC comparison/Genelists")
Patir <- read_excel("PatirEtAl_coreMicrogliaGenes.xlsx", col_names = "Patir_microglia")



#Read in ASD gene list
dir <- file.path("D:", "Mount Sinai", "de Witte lab", "Project 1 - 10x and VASAseq", "10x", "Analysis", "10x ASD")
setwd(dir)
ASD_sign <- read_excel("Genes ASD_syndromes.xlsx")




#==========================================================================================================#

#==========================================================================================================#
#                                         Analysis start: ME-covariate cor                                 #
#==========================================================================================================#
#Question: where is the effect?

#Adding some covariates
ME.covR <- rcorr(MEs, as.matrix(cbind("Timepoint"=cov_for.cor$Timepoint,"Concentration"=cov_for.cor$Concentration)), type="spearman")$r
l <- length(colnames(MEs))+1

#Two heatmaps: #1 ME-covariate correlations, #2 ME-ME correlations
ComplexHeatmap::pheatmap(ME.covR[1:length(colnames(MEs)), l:length(colnames(ME.covR))], cluster_cols = F,
                         heatmap.color.code)

identical(rownames(MEs), rownames(cov_for.cor))

MEs_lin(MEs, cbind("Concentration"=covs_forremBatch$Concentration, "CTR"=covs_forremBatch$CTR, "10mM"=covs_forremBatch$`10mM`, "50mM"=covs_forremBatch$`50mM`))

res.check <- MEs_lin(MEs, cbind("Concentration"=covs_forremBatch$Concentration, "CTR"=covs_forremBatch$CTR, "10mM"=covs_forremBatch$`10mM`, "50mM"=covs_forremBatch$`50mM`), return.r = T)


#Check differential expression of module eigengene (plot boxplot)

dat.df <- NULL
dat.df$"CTR" <- t(sub.df[rownames(sub.df) %in% rownames(covs_forremBatch[covs_forremBatch$Stimulation == 1,]),])
dat.df$"STIM" <- t(sub.df[rownames(sub.df) %in% rownames(covs_forremBatch[covs_forremBatch$Stimulation == 2,]),])

dat.df <- data.frame(MEs, "Condition" = as.numeric(covs_forremBatch$Stimulation))
dat.df$Condition[dat.df$Condition == 1] <- "CTR"
dat.df$Condition[dat.df$Condition == 2] <- "STIM"

ggplot(dat.df, aes_string(x="Condition", y="MEtan", fill="Condition")) +
  geom_boxplot() +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("ME expression of tan module"))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())

ggplot(dat.df, aes_string(x="Condition", y="MEcyan", fill="Condition")) +
  geom_boxplot() +
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("ME expression of cyan module"))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())


#==========================================================================================================#
#                             Different approach: Only DEG-enriched modules                                #
#==========================================================================================================#
MEs <- as.matrix(MKMD.networks$MEs)

mod.cols <- data.frame(MKMD.networks$moduleColors)
mod.cols$genes <- rownames(mod.cols)
mod.cols <- mod.cols[order(mod.cols$MKMD.networks.moduleColors),]

ModGenes <- list()
for (i in unique(mod.cols[,1])){
  ModGenes[[i]] <- rownames(mod.cols[mod.cols[,1] %in% i,])
}

#=========================For DEGs==================================#
DEG.enrich <- matrix(nrow=length(names(ModGenes)), ncol=length(names(DEGs)))
rownames(DEG.enrich) <- names(ModGenes)
colnames(DEG.enrich) <- names(DEGs)
for (x in names(DEGs)){
  DEG.enrich[,x] <- -log(GSEA.byMod(mod.gl=ModGenes, DEGs[[x]], universe=18000)$q.val)
}

#Taking all results into one matrix
DEG.enrich.res <- matrix(nrow=27*6, ncol=6)
DEG.rep <- NULL
for (i in names(DEGs))
  DEG.rep <- c(DEG.rep, rep(i, 27))
rownames(DEG.enrich.res) <- DEG.rep
colnames(DEG.enrich.res ) <- colnames(GSEA.byMod(mod.gl=ModGenes, DEGs[[x]], universe=18000))

for (x in names(DEGs)){
  DEG.enrich.res[rownames(DEG.enrich.res) %in% x,] <- as.matrix(GSEA.byMod(mod.gl=ModGenes, DEGs[[x]], universe=18000))
}
DEG.enrich.res <- data.frame(DEG.enrich.res, check.names=F)
DEG.enrich.res$Module <- rep(names(ModGenes),6)


for (x in c(DEGs[[1]], DEGs[[2]], DEGs[[3]], DEGs[[4]], DEGs[[5]], DEGs[[6]])){
  DEG.enrich.res[rownames(DEG.enrich.res) %in% x,] <- as.matrix(GSEA.byMod(mod.gl=ModGenes, c(DEGs[[1]], DEGs[[2]], DEGs[[3]], DEGs[[4]], DEGs[[5]], DEGs[[6]]), universe=18000))
}

DEG.modules <- DEG.enrich[rowMeans(DEG.enrich) > 2.75/length(colnames(DEG.enrich)),]

ComplexHeatmap::pheatmap(DEG.modules, cluster_rows = F, cluster_cols = F, annotation_legend = F, show_colnames = T,
                         legend_breaks= c(0,-log(0.05),10, 20, 30, 40), 
                         legend = T, color = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30), 
                         legend_labels=c("0","2.995", "10", "20", "30", "-log(q) \n\n\n"))

ComplexHeatmap::pheatmap(t(DEG.enrich), cluster_rows = F, cluster_cols = F, annotation_legend = F, show_colnames = T,
                         legend_breaks= c(0,-log(0.05),10, 20, 30, 40), 
                         legend = T, color = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30), 
                         legend_labels=c("0","2.995", "10", "20", "30", "-log(q) \n\n\n"))


DEG.mods <- rownames(DEG.modules)



#Compounding up and down-reg genes into one
DEGs.compound <- list("W10"=c(DEGs$D100_DEGs_down,DEGs$D100_DEGs_up), "W5"=c(DEGs$D35_DEGs_down,DEGs$D35_DEGs_up), "W4"=c(DEGs$D20_DEGs_down, DEGs$D20_DEGs_up))
DEG.compounds.enrich<- matrix(nrow=27*3, ncol=6)
DEG.compound.rep <- NULL
for (i in names(DEGs.compound))
  DEG.compound.rep <- c(DEG.compound.rep, rep(i, 27))

rownames(DEG.compounds.enrich) <- DEG.compound.rep
colnames(DEG.compounds.enrich ) <- colnames(GSEA.byMod(mod.gl=ModGenes, DEGs.compound[[1]], universe=18000))

for (x in names(DEGs.compound)){
  DEG.compounds.enrich[rownames(DEG.compounds.enrich) %in% x,] <- as.matrix(GSEA.byMod(mod.gl=ModGenes, DEGs.compound[[x]], universe=18000))
}
DEG.compounds.enrich <- data.frame(DEG.compounds.enrich, check.names=F)
DEG.compounds.enrich$Module <- rep(names(ModGenes),3)

DEG.compounds.short <- matrix(nrow=length(names(ModGenes)), ncol=length(names(DEGs.compound)))
rownames(DEG.compounds.short) <- names(ModGenes)
colnames(DEG.compounds.short) <- names(DEGs.compound)

for (x in names(DEGs.compound)){
  DEG.compounds.short[,x] <- -log(GSEA.byMod(mod.gl=ModGenes, DEGs.compound[[x]], universe=18000)$q.val)
}
ComplexHeatmap::pheatmap(t(DEG.compounds.short), cluster_rows = T, cluster_cols = T, annotation_legend = F, show_colnames = T,
                         legend_breaks= c(0,-log(0.05),10, 30, 60, 98, 98), 
                         legend = T, color = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30), 
                         legend_labels=c("0","2.995", "10", "30", "60", "98", "-log(q) \n\n\n"))

template <- GO.all
names(template) <- c("W4", "W5")
template$W4@result <- template$W4@result[1:27,]
template$W4@result$Description <- rownames(DEG.compounds.short)
template$W4@result$p.adjust <- DEG.compounds.short[,1]

barplot(template$W4, showCategory=27)



#==========================================================================================================#
#                                             Annotate DEG modules                                         #
#==========================================================================================================#

#===============================For cell-type markers====================================#
CELL.enrich <- matrix(nrow=length(names(ModGenes)), ncol=4) #for visualization (heatmap)
rownames(CELL.enrich) <- names(ModGenes)
colnames(CELL.enrich) <- names(top.markers)
for (x in names(top.markers)){
  CELL.enrich[,x] <- GSEA.byMod(mod.gl=ModGenes, top.markers[[x]], universe=18000)$q.val
}
#Taking all results into one matrix
CELL.enrich.res <- matrix(nrow=27*4, ncol=6)
top.markers.rep <- NULL
for (i in names(top.markers))
  top.markers.rep <- c(top.markers.rep, rep(i, 27))
rownames(CELL.enrich.res) <- top.markers.rep
colnames(CELL.enrich.res ) <- colnames(GSEA.byMod(mod.gl=ModGenes, top.markers[[x]], universe=18000))

for (x in names(top.markers)){
  CELL.enrich.res[rownames(CELL.enrich.res) %in% x,] <- as.matrix(GSEA.byMod(mod.gl=ModGenes, top.markers[[x]], universe=18000))
}
CELL.enrich.res <- data.frame(CELL.enrich.res, check.names=F)
CELL.enrich.res$Module <- rep(names(ModGenes),4)


#Based on selected cell-type markers
manual.CELL.enrich <- matrix(nrow=length(names(ModGenes)), ncol=length(names(select_cell.markers)))
rownames(manual.CELL.enrich) <- names(ModGenes)
colnames(manual.CELL.enrich) <- names(select_cell.markers)
for (x in colnames(manual.CELL.enrich)){
  manual.CELL.enrich[,x] <- -log(GSEA.byMod(mod.gl=ModGenes, select_cell.markers[[x]], universe=18000)$q.val)
}

#Taking all results into one matrix
manual.CELL.enrich.res <- matrix(nrow=27*length(names(select_cell.markers)), ncol=6)
man.CELL.rep <- NULL
for (i in names(select_cell.markers))
  man.CELL.rep <- c(man.CELL.rep, rep(i, 27))
rownames(manual.CELL.enrich.res) <- man.CELL.rep
colnames(manual.CELL.enrich.res) <- colnames(GSEA.byMod(mod.gl=ModGenes, select_cell.markers[[x]], universe=18000))
for (x in names(select_cell.markers)){
  manual.CELL.enrich.res[rownames(manual.CELL.enrich.res) %in% x,] <- as.matrix(GSEA.byMod(mod.gl=ModGenes, select_cell.markers[[x]], universe=18000))
}
manual.CELL.enrich.res <- data.frame(manual.CELL.enrich.res)
manual.CELL.enrich.res$Module <- rep(names(ModGenes),length(names(select_cell.markers)))




#Exemplary plot
ComplexHeatmap::pheatmap(CELL.enrich[rownames(CELL.enrich) %in% DEG.mods,], cluster_rows = T, cluster_cols = T, annotation_legend = F, show_colnames = T,
                         legend_breaks= c(0,-log(0.05),10, 20, 30, 40),
                         legend = T, color = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30), 
                         legend_labels=c("0","2.995", "10", "20", "30", "-log(q) \n\n\n"))

#=====Gandal enrich=================#

#For visualization: only q values
Cluster.enrich <- matrix(nrow=length(names(ModGenes)), ncol=4)
rownames(Cluster.enrich) <- names(ModGenes)
colnames(Cluster.enrich) <- c("Microglia", "Astrocytes", "NfKB", "IFNy")
for (i in colnames(Cluster.enrich))
  Cluster.enrich[,i] <- -log(GSEA.byMod(ModGenes, Gandal.clusters.annot[Gandal.clusters.annot$Module %in% i,]$gene_name, universe=18000)$q.val)

#All results in a matrix
Cluster.enrich.res <- matrix(nrow=27*4, ncol=6)
Cluster.enrich.rep <- NULL
for (i in colnames(Cluster.enrich))
  Cluster.enrich.rep <- c(Cluster.enrich.rep, rep(i, 27))
rownames(Cluster.enrich.res) <- Cluster.enrich.rep
colnames(Cluster.enrich.res) <- colnames(GSEA.byMod(mod.gl=ModGenes, select_cell.markers[[x]], universe=18000))
for (x in colnames(Cluster.enrich)){
  Cluster.enrich.res[rownames(Cluster.enrich.res) %in% x,] <- as.matrix(GSEA.byMod(mod.gl=ModGenes, Gandal.clusters.annot[Gandal.clusters.annot$Module %in% x,]$gene_name, universe=18000))
}
Cluster.enrich.res <- data.frame(Cluster.enrich.res)
Cluster.enrich.res$Module <- rep(names(ModGenes),length(colnames(Cluster.enrich)))

#=========================For Quadrato organoid markers==================================#
Quad.markers <- Quad.markers[,!colnames(Quad.markers) %in% c("UNKNOWN1", "UNKNOWN2", "UNKNOWN3")]
Quad.enrich <- matrix(nrow=27*7, ncol=6)
Quad.enrich.rep <- NULL
Quad.enrich.rep <- NULL
for (i in colnames(Quad.markers))
  Quad.enrich.rep <- c(Quad.enrich.rep, rep(i, 27))
rownames(Quad.enrich) <- Quad.enrich.rep
colnames(Quad.enrich) <- colnames(GSEA.byMod(mod.gl=ModGenes, select_cell.markers[[x]], universe=18000))
for (x in colnames(Quad.markers)){
  Quad.enrich[rownames(Quad.enrich) %in% x,] <- as.matrix(GSEA.byMod(mod.gl=ModGenes, na.omit(Quad.markers[[x]]), universe=18000))
}
Quad.enrich <- data.frame(Quad.enrich, check.rows = F, check.names = F)
Quad.enrich$Module <- rep(names(ModGenes),length(colnames(Quad.markers)))

Quad.enrich.q <- matrix(nrow=length(names(ModGenes)), ncol=7)
rownames(Quad.enrich.q) <- names(ModGenes)
colnames(Quad.enrich.q) <- colnames(Quad.markers)
for (i in colnames(Quad.markers))
  Quad.enrich.q[,i] <- -log(GSEA.byMod(ModGenes, na.omit(Quad.markers[[i]]), universe=18000)$q.val)


ComplexHeatmap::pheatmap(Quad.enrich.q, cluster_rows = T, cluster_cols = T, annotation_legend = F, show_colnames = T,
                         legend_breaks= c(0,-log(0.05),10, 20, 30, 40),
                         legend = T, color = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30), 
                         legend_labels=c("0","2.995", "10", "20", "30", "-log(q) \n\n\n"))


#=========================For disease markers==================================#
#For ASD
ASD <- list("ASD"=ASD_sign$`ASD Gene`, "Rad"=ASD_sign$`Radboud genes`[!as.character(ASD_sign$`Radboud genes`) %in% NA])

#Restrict to expressed genes
ASD$ASD <- ASD$ASD[ASD$ASD %in% colnames(sub.df)]
ASD$Rad <- ASD$Rad[ASD$Rad %in% colnames(sub.df)]

ASD.enrich <- matrix(nrow=length(names(ModGenes)), ncol=2)
rownames(ASD.enrich) <- names(ModGenes)
colnames(ASD.enrich) <- names(ASD)
for (x in names(ASD)){
  ASD.enrich[,x] <- -log(GSEA.byMod(mod.gl=ModGenes, ASD[[x]], universe=18000)$q.val)
}


#For schizophrenia
SCZ.genes
disease.enrich <- cbind(ASD.enrich, "SCZ"=-log(GSEA.byMod(mod.gl=ModGenes, SCZ.genes, universe=18000)$q.val))


#All results in a matrix
disease.enrich.res <- matrix(nrow=27*3, ncol=6)
disease.enrich.rep <- NULL
for (i in colnames(disease.enrich))
  disease.enrich.rep <- c(disease.enrich.rep, rep(i, 27))
rownames(disease.enrich.res) <- disease.enrich.rep
colnames(disease.enrich.res) <- colnames(GSEA.byMod(mod.gl=ModGenes, select_cell.markers[[x]], universe=18000))
for (x in 1:3){
  disease.enrich.res[rownames(disease.enrich.res) %in% unique(rownames(disease.enrich.res))[x],] <- as.matrix(GSEA.byMod(mod.gl=ModGenes, c(ASD,list(SCZ.genes))[[x]], universe=18000))
}
disease.enrich.res <- data.frame(disease.enrich.res)
disease.enrich.res$Module <- rep(names(ModGenes),length(colnames(disease.enrich)))



#==========================Summarizing===============================#
all.enrich <- cbind(disease.enrich, Cluster.enrich, CELL.enrich, manual.CELL.enrich, Quad.enrich.q)
all.enrich <- all.enrich[rowSums(all.enrich) > 2.995,]
all.enrich <- all.enrich[,colSums(all.enrich) > 2.995]
ComplexHeatmap::pheatmap(all.enrich[rownames(all.enrich) %in% DEG.mods,], cluster_rows = T, cluster_cols = T, annotation_legend = F, show_colnames = T,
                         legend_breaks= c(0,-log(0.05),10, 20, 30, 60, 100, 110), 
                         legend = T, color = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30),
                         legend_labels=c("0","2.995", "10", "20", "30", "60", "100", "-log(q) \n\n\n"))

ComplexHeatmap::pheatmap(t(all.enrich[rownames(all.enrich) %in% c("cyan", "tan"),]), cluster_rows = T, cluster_cols = T, annotation_legend = F, show_colnames = T,
                         legend_breaks= c(0,-log(0.05),10, 20, 30, 60, 100, 110), 
                         legend = T, color = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30),
                         legend_labels=c("0","2.995", "10", "20", "30", "60", "100", "-log(q) \n\n\n"))



#==========================================================================================================#
#                                              GO DEG module annotation                                    #
#==========================================================================================================#
#Module annotation of top 100 genes (correlated with hub gene)
cor <- WGCNA::cor
MM <- getMM(
  net = MKMD.networks,
  expr.data.file = sub.df,
  tissue=tissue,
  genes = NULL,
  silent = F,
  keep.grey = F,
  identicalNames = T,
  alt.gene.index = NULL,
  dupAware = T
)

hub.genes <- NULL
mod.names <- NULL
for (i in names(ModGenes)){
  mod <- MM[MM$module %in% i,]
  ordered <- mod[order(mod$mm, decreasing=T),]
  hub.genes <- c(hub.genes,ordered[1,]$name)
  mod.names <- c(mod.names,i)
  names(hub.genes) <- mod.names
  }


hub.genes <- data.frame(hub.genes)

MM <- MM[!MM$name %in% grep(pattern = "^ERCC-", x = MM$name, value = TRUE),]

top25.hubgenes <- NULL
for (i in names(ModGenes)){
  mod <- MM[MM$module %in% i,]
  ordered <- mod[order(mod$mm, decreasing=T),]
  top25.hubgenes <- rbind(top25.hubgenes,na.omit(ordered[1:25,]))
}
top25.hubgenes <- top25.hubgenes[,-1]
  
top25.cor <- NULL
for (i in mod.names){
  cgene <- hub.genes[rownames(hub.genes) %in% i,]
  cgene.expr <- sub.df[,colnames(sub.df) %in% cgene]
  mod.expr <- sub.df[,colnames(sub.df) %in% ModGenes[[i]]]
  cor.expr <- cor(mod.expr,cgene.expr)^2
  top25.cor <- rbind(top25.cor,cbind(names(cor.expr[order(cor.expr[,1], decreasing=T),][1:26]),rep(i,26),cor.expr[order(cor.expr[,1], decreasing=T),][1:26]))}
colnames(top25.cor) <- c("Gene","Module","R^2")
top25.cor <- data.frame(top25.cor)
rownames(top25.cor) <- NULL

write.xlsx2(top25.hubgenes, "09172021_MKMD WGCNA V9_HubGenes.xlsx", sheetName = "25 hub genes (R^2)") #to annotate cell-type with disease-module interactions
write.xlsx2(top25.cor, "09172021_MKMD WGCNA V9_HubGenes.xlsx", sheetName = "25 top correlated genes (with top hub gene, R^2)", append=T) #to annotate cell-type with disease-module interactions



top100 <- list()
for (i in unique(mod.cols$MKMD.networks.moduleColors)){
  cgene <- hub.genes[rownames(hub.genes) %in% i,]
  cgene.expr <- sub.df[,colnames(sub.df) %in% cgene]
  mod.expr <- sub.df[,colnames(sub.df) %in% ModGenes[[i]]]
  cor.expr <- cor(mod.expr,cgene.expr)^2
  top100[[i]] <- names(cor.expr[order(cor.expr[,1], decreasing=T),][1:100])}


library(biomaRt)
genenames <- colnames(sub.df)
listMarts(host="www.ensembl.org")
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
dfs.ensembl <- listDatasets(ensembl)
genemap.mic <- getBM(values = genenames,
                     filters = "external_gene_name",
                     mart = ensembl,
                     attributes = c("ensembl_gene_id", "entrezgene_id",
                                    "hgnc_symbol", "external_gene_name",
                                    "description", "chromosome_name",
                                    "strand"))

GO <- list()
for (i in names(top100)){
  GO[[i]] <- enrichGO(gene    = genemap.mic[genemap.mic$external_gene_name %in% top100[[i]],]$ensembl_gene_id,
                universe      = genemap.mic[genemap.mic$external_gene_name %in% colnames(sub.df),]$ensembl_gene_id,
                OrgDb         = "org.Hs.eg.db",
                keyType       = 'ENSEMBL',
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE, pool=T)}


for (i in DEG.mods){
  print(paste("Module", i, ":"))
  print(head(GO[[i]])[colnames(head(GO[[i]])) %in% c("Description","p.adjust")])}


DEG.GOplots <- list()
for (i in DEG.mods){
  DEG.GOplots[[i]] <- dotplot(GO.all$cyan, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
}

remove(GO)

DEG.GOplots$cyan
DEG.GOplots$tan


DEG.GO.plot <- ggpubr::ggarrange(plotlist=DEG.GOplots, common.legend = F)




#=========================GO enrichment based on all module genes================================#
GO.all <- list()
for (i in c("cyan","tan")){
  GO.all[[i]] <- enrichGO(gene  = genemap.mic[genemap.mic$external_gene_name %in% ModGenes[[i]],]$ensembl_gene_id,
                      universe      = names(genemap.mic[genemap.mic$external_gene_name %in% colnames(sub.df),]$ensembl_gene_id),
                      OrgDb         = "org.Hs.eg.db",
                      keyType       = 'ENSEMBL',
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE, pool=T)}

for (i in c("cyan","tan")){
  print(paste("Module", i, ":"))
  print(head(GO.all[[i]])[colnames(head(GO.all[[i]])) %in% c("Description","p.adjust")])}


DEG.GOplots.all <- list()
for (i in c("cyan","tan")){
  DEG.GOplots.all[[i]] <- dotplot(GO.all$cyan, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")}

#Immune response in tan module
DEG.GOplots.all$tan


#==========================================================================================================#
#                                     Analysis: DEG-cell.type/disease cor                                  #
#==========================================================================================================#
#Oh yes, correlate the PCs of certain cell-type, disease modules etc. with the expression of DEGs (kME) to get an idea how they influence each other!
disease.PC <- data.frame(ASD = prcomp(sub.df[,colnames(sub.df) %in% ASD_sign$`ASD Gene`])$x[,1])
disease.PC$Radboud <- prcomp(sub.df[,colnames(sub.df) %in% ASD_sign$`Radboud genes`])$x[,1]
disease.PC$SCZ <- prcomp(sub.df[,colnames(sub.df) %in% SCZ.genes])$x[,1]


DEG.disease_cor <- data.frame(cor(DEG.expr.df,disease.PC))


ComplexHeatmap::pheatmap(DEG.disease_cor, cluster_rows = T, cluster_cols = F, annotation_legend = T, 
                         show_colnames = T, show_rownames = F, legend = T,
                         color = heatmap.color.code)

cor(DEG.disease_cor$ASD, DEG.disease_cor$Rad)

library(ggpmisc)
library(ggplot2)
formula <- y ~ x
ggplot(DEG.disease_cor, aes_string(x = as.matrix(DEG.disease_cor$ASD), y = as.matrix(DEG.disease_cor$SCZ))) + 
  stat_smooth(method = lm) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + geom_smooth(method = lm, linetype = "dashed") + 
  labs(y = "SCZ disease signature", x = "ASD disease signature") + 
  stat_poly_eq(aes(label = c(paste(..adj.rr.label..))), label.x = "right", label.y = 0.37, 
               formula = formula, parse = TRUE, size = 3) + # geom_cor(method='pearson')+
  stat_fit_glance(method = "lm", method.args = list(formula = formula), geom = "text", 
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")), 
                  label.x = "right", label.y = 0.3, size = 3)



#===========================for all DEGs========================#
allDEG.expr <- NULL
for (i in names(DEGs))
  allDEG.expr[[i]] <- sub.df[,colnames(sub.df) %in% DEGs[[i]]]

allDEG.expr.df <- NULL
for (i in names(allDEG.expr)){
  allDEG.expr.df <- cbind(allDEG.expr.df, allDEG.expr[[i]])}

allDEG.disease_cor <- data.frame(cor(allDEG.expr.df,disease.PC))

ggplot(DEG.disease_cor, aes_string(x = as.matrix(DEG.disease_cor$ASD), y = as.matrix(DEG.disease_cor$SCZ))) + 
  stat_smooth(method = lm) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + geom_smooth(method = lm, linetype = "dashed") + 
  labs(y = "SCZ disease signature", x = "ASD disease signature") + 
  stat_poly_eq(aes(label = c(paste(..adj.rr.label..))), label.x = "right", label.y = 0.37, 
               formula = formula, parse = TRUE, size = 3) + # geom_cor(method='pearson')+
  stat_fit_glance(method = "lm", method.args = list(formula = formula), geom = "text", 
                  aes(label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")), 
                  label.x = "right", label.y = 0.3, size = 3)



#==========================================================================================================#
#                 This warrants further investigation: ASD vs SCZ expression modules                       #
#==========================================================================================================#
# that means... they should be extremely similar (because negative expression to ^5 is positive)
#From enrichment analysis: midnightblue is an ASD module
MEs_lin(MEs, covs_forremBatch)
#Midnightblue is enriched for time-effect

disease.ME <- MEs_lin(MEs, disease.PC, return.r=T)


disease.ME <- disease.ME[,colMeans(disease.ME) > 0.11/3]

#How are disease modules co-regulated with DEG modules?
cor(MEs[,colnames(disease.ME)], MEs[,paste("ME",DEG.mods,sep="")])
cor(MEs[,c("MEmidnightblue", "MEdarkturquoise", "MEbrown")], MEs[,paste("ME",DEG.mods,sep="")])

ComplexHeatmap::pheatmap(cor(MEs[,c("MEmidnightblue", "MEdarkturquoise", "MEbrown")], MEs[,paste("ME",DEG.mods,sep="")]), 
                         cluster_rows = T, cluster_cols = T, annotation_legend = T, 
                         show_colnames = T, show_rownames = T, legend = T,
                         color = heatmap.color.code
)

module.interactions <- MEs_lin(MEs[,c("MEmidnightblue", "MEdarkturquoise", "MEbrown")], MEs[,paste("ME",DEG.mods,sep="")], return.r=T)
head(GO.all)[colnames(head(GO.all$lightcyan)) %in% c("Description","p.adjust")]

module.interactions.all <- MEs_lin(MEs, MEs, return.r=T)

write.xlsx2(module.interactions.all, "09132021_MKMD WGCNA V9_Module Interactions.xlsx", sheetName = "adjusted R^2") #to annotate cell-type with disease-module interactions


dotplot(GO.all$cyan, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
head(GO$midnightblue)[colnames(head(GO$midnightblue)) %in% c("Description","p.adjust")]
head(GO.all$midnightblue)[colnames(head(GO.all$midnightblue)) %in% c("Description","p.adjust")]
head(GO$darkturquoise)[colnames(head(GO$darkturquoise)) %in% c("Description","p.adjust")]
head(GO.all$brown)[colnames(head(GO.all$brown)) %in% c("Description","p.adjust")]
head(GO.all$darkred)[colnames(head(GO.all$darkred)) %in% c("Description","p.adjust")]
#Darkturquoise (SCZ): AMPA receptor and postsynapcit processes (and catabolic)
#Midnightblue (ASD): Synaptic membrane, glutamatergic synapse
#Brown (ASD & SCZ): Synaptic membrane, postsynaptic density/specialization. 


#Next step: check for interaction effect between time-stimulation?


#More important question: development of these modules over time?
#Can we plot them over time (W4 - W5 - W15)

#Predict ME-covariate correlation based on mTOR expression; check mTOR expression correlation with covariates

sub.df[,colnames(sub.df) %in% "MTOR"]

sub.df[,colnames(sub.df) %in% "MTOR"]

glance( lm(MEs[,colnames(MEs) %in% "MEtan"] ~ covs_forremBatch$`50mM` + sub.df[,colnames(sub.df) %in% "MTOR"]) )
glance( lm(MEs[,colnames(MEs) %in% "MEtan"] ~ covs_forremBatch$`50mM`) )

glance( lm(MEs[,colnames(MEs) %in% "MEcyan"] ~ covs_forremBatch$`50mM` + sub.df[,colnames(sub.df) %in% "MTOR"]) )
glance( lm(MEs[,colnames(MEs) %in% "MEcyan"] ~ covs_forremBatch$`50mM`) )


#For concentration:
glance( lm(MEs[,colnames(MEs) %in% "MEtan"] ~ covs_forremBatch$`Concentration` + sub.df[,colnames(sub.df) %in% "MTOR"]) )

glance( lm(MEs[,colnames(MEs) %in% "MEtan"] ~ covs_forremBatch$`Concentration`) )

glance( lm(MEs[,colnames(MEs) %in% "MEcyan"] ~ covs_forremBatch$`Concentration` + sub.df[,colnames(sub.df) %in% "MTOR"]) )

glance( lm(MEs[,colnames(MEs) %in% "MEcyan"] ~ covs_forremBatch$`Concentration`) )



#==========================================================================================================#
#                                    Conclusion: Write results into excel                                  #
#==========================================================================================================#
setwd("D:/Mount Sinai/de Witte lab/Side projects/Amber bulk RNA-seq of stimulated organoids/WGCNA/")

#For GO done on all module genes
df.results.GO <- matrix(ncol = length(names(GO.all$black[1:100])), nrow = 2700)
mod.names.rep <- NULL

for (i in names(GO.all))
  mod.names.rep <- c(mod.names.rep, rep(i, 100))
rownames(df.results.GO) <- mod.names.rep
colnames(df.results.GO) <- names(GO.all$black[1:100])


for (i in names(GO.all))
  df.results.GO[rownames(df.results.GO) %in% i,]  <- as.matrix(GO.all[[i]][1:100])
df.results.GO <- data.frame(df.results.GO)
df.results.GO <- drop_na(df.results.GO)

#For GO done on all high-confidence module genes
df.GO.highconf <- matrix(ncol = length(names(GO$black[1:20])), nrow = 540)
mod.names.rep2 <- NULL
for (i in names(GO.all))
  mod.names.rep2 <- c(mod.names.rep2, rep(i, 20))
rownames(df.GO.highconf) <- mod.names.rep2
colnames(df.GO.highconf) <- names(GO$black[1:20])

for (i in names(GO.all))
  df.GO.highconf[rownames(df.GO.highconf) %in% i,]  <- as.matrix(GO[[i]][1:20])
df.GO.highconf <- data.frame(df.GO.highconf)
df.GO.highconf <- drop_na(df.GO.highconf)


#Adding whether they are enriched for DEGs
df.GO.highconf$DEGmod <- NA
for (i in DEG.mods)
  df.GO.highconf[rownames(df.GO.highconf) %in% c(i, paste(i, 1:99, sep=".")),]$DEGmod <- "YES"
df.GO.highconf$DEGmod[is.na(df.GO.highconf$DEGmod)] <- "NO"

#RACK1 (hub gene of cyan module) sits in ribosomal pathways and is associated with neurodevelopment (termed a neurodevelopmental regulator)

write.xlsx2(df.results.GO, "09072021_MKMD WGCNA V9_GO.xlsx", sheetName = "all_genes")
write.xlsx2(df.GO.highconf, "09072021_MKMD WGCNA V9_GO.xlsx", append=TRUE, sheetName = "high_confidence_genes")


#================grab beta for status / concentration predicting MEs=====================#
library(lm.beta)
glance( lm(MEs[,'MEtan'] ~ covs_forremBatch[,'Stimulation']) )$r.squared
summary( lm(MEs[,'MEtan'] ~ covs_forremBatch[,'Concentration']) )$adj.r.squared
lm(MEs[,'MEtan'] ~ covs_forremBatch[,'Stimulation'])
lm(MEs[,'MEtan'] ~ covs_forremBatch[,'CTR'])
lm(MEs[,'MEtan'] ~ covs_forremBatch[,'10mM'])
lm(MEs[,'MEtan'] ~ covs_forremBatch[,'50mM'])
lm(MEs[,'MEtan'] ~ covs_forremBatch[,'Concentration'])

summary( lm(MEs[,'MEtan'] ~ covs_forremBatch[,'Concentration']) )

df.coef.MEs <- data.frame(matrix(nrow=2, ncol=4))
colnames(df.coef.MEs) <- c("Stimulation","Concentration", "sd_stim", "sd_con")

df.coef.MEs$Stimulation <- c(coef(lm(MEs[,'MEtan'] ~ covs_forremBatch[,'Stimulation']))[2], coef(lm(MEs[,'MEcyan'] ~ covs_forremBatch[,'Stimulation']))[2])
df.coef.MEs$Concentration <- c(coef(lm(MEs[,'MEtan'] ~ covs_forremBatch[,'Concentration']))[3], coef(lm(MEs[,'MEcyan'] ~ covs_forremBatch[,'Concentration']))[3])
df.coef.MEs$sd_stim <- c(coef(summary(lm(MEs[,'MEtan'] ~ covs_forremBatch[,'Stimulation'])))[,2][2], coef(summary(lm(MEs[,'MEcyan'] ~ covs_forremBatch[,'Stimulation'])))[,2][2])
df.coef.MEs$sd_con <- c(coef(summary(lm(MEs[,'MEtan'] ~ covs_forremBatch[,'Concentration'])))[,2][3], coef(summary(lm(MEs[,'MEcyan'] ~ covs_forremBatch[,'Concentration'])))[,2][3])

df.coef.MEs$Module <- c("Tan", "Cyan")

#Plot the results
ggplot(df.coef.MEs, aes_string(x="Module", y="Stimulation", fill = "Module")) +
  geom_bar(stat="identity") +
  theme_classic()+
  scale_y_continuous(expand = c(0,0),
                     limits = c(-0.2,0.2)) +
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("Beta"))+
  ggtitle("ME-Stimulation")+
  geom_errorbar( aes(x=Module, ymin=Stimulation-sd_stim, ymax=Stimulation+sd_stim), width=0.2, colour="black", alpha=0.4, size=0.8)

ggplot(df.coef.MEs, aes_string(x="Module", y="Concentration", fill = "Module")) +
  geom_bar(stat="identity") +
  theme_classic()+
  scale_y_continuous(expand = c(0,0),
                     limits = c(-0.2,0.2)) +
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("Beta"))+
  ggtitle("ME-Concentration")+
  geom_errorbar( aes(x=Module, ymin=Concentration-sd_con, ymax=Concentration+sd_con), width=0.2, colour="black", alpha=0.4, size=0.8)






#============================manually selected gene lists=========#
write.xlsx2(CELL.enrich.res, "09072021_MKMD WGCNA V9_FisherEnrichment.xlsx", sheetName = "unb_cell_type_enrichment")
write.xlsx2(manual.CELL.enrich.res, "09072021_MKMD WGCNA V9_FisherEnrichment.xlsx", append=T, sheetName = "select_cell_type_enrichment")
write.xlsx2(Cluster.enrich.res, "09072021_MKMD WGCNA V9_FisherEnrichment.xlsx", append=T, sheetName = "Gandal_cluster_enrichment")
write.xlsx2(disease.enrich.res, "09072021_MKMD WGCNA V9_FisherEnrichment.xlsx", append=T, sheetName = "diseaseSignature_enrichment")
write.xlsx2(DEG.enrich.res, "09072021_MKMD WGCNA V9_FisherEnrichment.xlsx", append=T, sheetName = "DEG_enrichment")
write.xlsx2(Quad.enrich, "09072021_MKMD WGCNA V9_FisherEnrichment.xlsx", append=T, sheetName = "Quadrato_cluster_enrichment")


#Patir microglia
Patir.enrich <- data.frame(GSEA.byMod(mod.gl=ModGenes, Patir$Patir_microglia, universe=18000))
Patir.enrich <- cbind(Patir.enrich, rownames(Patir.enrich))
rownames(Patir.enrich) <- c("Microglia",paste("Microglia", 1:26, sep="."))
colnames(Patir.enrich) <- colnames(Quad.enrich)
cell.typesQuadplusMic <- rbind(Quad.enrich, Patir.enrich)


write.xlsx2(cell.typesQuadplusMic, "12082021_MKMD WGCNA V10.1_FisherEnrichment_CellTypes.xlsx", append=T, sheetName = "Quadrato+Microglia_cluster_enrichment")




#=============follow-up investigation===============#
write.xlsx2(DEG.disease_cor, "09072021_MKMD WGCNA V9_DEG coReg Structure.xlsx", sheetName = "DEG-dis_cor (Pearson's r)")
write.xlsx2(allDEG.disease_cor, "09072021_MKMD WGCNA V9_DEG coReg Structure.xlsx", append=T, sheetName = "all DEGs-dis_cor (Pearson's r)")
write.xlsx2(disease.ME, "09072021_MKMD WGCNA V9_DEG coReg Structure.xlsx", append=T, sheetName = "disease ME-disease cor (adjusted R^2)")
write.xlsx2(module.interactions, "09072021_MKMD WGCNA V9_DEG coReg Structure.xlsx", append=T, sheetName = "module-module interactions (adjusted R^2)") #to annotate cell-type with disease-module interactions


#===========Mod gene as excel=========================#
library(plyr)
library(xlsx)
ModGenes.table <- NULL
for (i in names(ModGenes))
  ModGenes.table <- rbind.fill(ModGenes.table,data.frame(t(ModGenes[[i]]), check.names = F))
ModGenes.table <- t(ModGenes.table)
colnames(ModGenes.table) <- names(ModGenes)

write.xlsx2(ModGenes.table, "09132021_MKMD WGCNA V9_ModGenes.xlsx", sheetName = "Genes by module") #to annotate cell-type with disease-module interactions










#============Figure updates 12/17/2021==================#
#1. Top 10 upregulated genes Top 10 downregulated genes based on FC (no fold change cut of - as we did over here)
#2. Top 10 upregulated Top 10 downregulated genes based on p-Value
#3. Top 20 regulated DEG based on FC (no cut off - not up and down seperate) and sorted for up and downregulated genes (could be that a greater proportion of the highest affected genes at one timepoint are downregulated for example).


ComplexHeatmap::pheatmap(cor(MEs[,c("MEmidnightblue", "MEdarkturquoise", "MEbrown")], MEs[,paste("ME",DEG.mods,sep="")]), 
                         cluster_rows = T, cluster_cols = T, annotation_legend = T, 
                         show_colnames = T, show_rownames = T, legend = T,
                         color = heatmap.color.code
)


#For D20
setwd(file.path("D:", "R", "R workspace"))
load(file="04012021_MKMD D20_correct model.RData")
d20.res <- res.complex
d20.annot <- annot_df
d20sub.df <- sorted


load(file="04012021_MKMD D35_correct model.RData")
d35.res <- res.complex
d35.annot <- annot_df
d35sub.df <- sorted


load(file="04012021_MKMD D100_correct model.RData")
d100.res <- res.complex
d100.annot <- annot_df
d100sub.df <- sorted


annot_df.tog <- rbind(d20.annot,d35.annot,d100.annot)
library(plyr)
#df.together <- rbind.fill(data.frame(t(d20sub.df)),data.frame(t(d35sub.df)),data.frame(t(d100sub.df)))
#df.together <- data.frame(t(df.together))
#colnames(df.together) <- c(colnames(d20sub.df),colnames(d35sub.df),colnames(d100sub.df))

object.list <- list(d20.res,d35.res,d100.res)
df.list <- list(d20sub.df,d35sub.df,d100sub.df)

#lFC combined, only up and down
#Based on fold change 10 most up 10 most down regulated for each timepoint seperately
#Based on p-adj 10 most significant up and 10 most significant down regulated for each timepoint seperately
#Based on fold change all genes together 20 most differentially expressed genes - then sort them by up and down regulated genes just for the visual ( I made the example for day 24) 

#1)Based on lFC
for (i in 1:3){
  object <- object.list[[i]]
  res.sorted.up <- object[order(object$log2FoldChange, decreasing = T),]
  res.sorted.up <- res.sorted.up[res.sorted.up$log2FoldChange > 0 & res.sorted.up$padj < 0.05,]
  res.sorted.down <- object[order(object$log2FoldChange, decreasing = F),]
  res.sorted.down <- res.sorted.down[res.sorted.down$log2FoldChange < 0 & res.sorted.down$padj < 0.05,]
  top.downreg.genes <- rownames(res.sorted.down)[1:10]
  top.upreg.genes <- rownames(res.sorted.up)[1:10]
  df <- df.list[[i]]
  top.sub.df <- df[rownames(df) %in% c(top.upreg.genes,top.downreg.genes),]
  top.sub.df <- top.sub.df[match(c(top.upreg.genes,top.downreg.genes), rownames(top.sub.df)),]
  annot <- annot_df.tog[rownames(annot_df.tog) %in% colnames(df),]
  name <- c("12272021_D20_lFC heatmap.pdf", "12272021_D35_lFC heatmap.pdf", "12272021_D100_lFC heatmap.pdf")[i]
  pdf(file=name)
  ComplexHeatmap::pheatmap(top.sub.df, 
                           cluster_rows = F, 
                           cluster_cols = F, annotation_legend = T, show_colnames = F,
                           annotation_col=annot[-2], fontsize_row = 7,
                           color = heatmap.color.code, scale='row', annotation_colors = annot_colors)
  dev.off()
  
}

#2)Based on pval
for (i in 1:3){
  object <- object.list[[i]]
  sort.p <- object[order(object$padj, decreasing=F),]
  top.p <- rbind(sort.p[sort.p$log2FoldChange > 0,][1:10,],sort.p[sort.p$log2FoldChange < 0,][1:10,])
  df <- df.list[[i]]
  top.p_df <- df[rownames(df) %in% rownames(top.p),]
  top.p_df <- top.p_df[match(rownames(top.p), rownames(top.p_df)),]
  annot <- annot_df.tog[rownames(annot_df.tog) %in% colnames(df),]
  name <- c("12272021_D20_padj heatmap.pdf", "12272021_D35_padj heatmap.pdf", "12272021_D100_padj heatmap.pdf")[i]
  pdf(file=name)
  ComplexHeatmap::pheatmap(top.p_df, 
                           cluster_rows = F, 
                           cluster_cols = F, annotation_legend = T, show_colnames = F,
                           annotation_col=annot[-2], fontsize_row = 7,
                           color = heatmap.color.code, scale='row', annotation_colors = annot_colors)
  dev.off()
  
}


#3) Top 20 based on lFC, sorted after
for (i in 1:3){
  object <- object.list[[i]]
  top.lfc <- object[object$padj < 0.05,]
  top.lfc <- top.lfc[order(top.lfc$log2FoldChange^2, decreasing = T),][1:20,]
  top.lfc <- top.lfc[order(top.lfc$log2FoldChange, decreasing=T),]
  df <- df.list[[i]]
  top20lfc_df <- df[rownames(df) %in% rownames(top.lfc),]
  top20lfc_df <- top20lfc_df[match(rownames(top.lfc), rownames(top20lfc_df)),]
  annot <- annot_df.tog[rownames(annot_df.tog) %in% colnames(df),]
  name <- c("12272021_D20_lFC top20 heatmap.pdf", "12272021_D35_lFC top20 heatmap.pdf", "12272021_D100_lFC top20 heatmap.pdf")[i]
  pdf(file=name)
  ComplexHeatmap::pheatmap(top20lfc_df, 
                           cluster_rows = F, 
                           cluster_cols = F, annotation_legend = T, show_colnames = F,
                           annotation_col=annot[-2], fontsize_row = 7,
                           color = heatmap.color.code, scale='row', annotation_colors = annot_colors)
  dev.off()
  
}


###############Module structures#################
#Genecluster tree (don't have the TOM structure)
geneTree = hclust(as.dist(MKMD.networks$adjacency), method = "average");
# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


#PCA / dim reduction plot
module.plot <- data.frame("MEs"=rowMeans(MEs))
PoV.norm <- round(100 * attr(pca.norm, "percentVar"))


ME.covs <- data.frame(Color=1:27)
rownames(ME.covs) <- colnames(MEs)
ME.covs$Color <- gsub(".*ME","", colnames(MEs))

module.dat <- plotPCA.custom(MEs, intgroup=c("Color"), 
               ntop = 50000, returnData=TRUE, metadata=ME.covs)
PoV.mods <- round(100 * attr(module.dat, "percentVar"))

PCAplot(module.dat, PoV.df = PoV.mods, Condition = "Color", colors=gsub(".*ME","", colnames(MEs)))


#Dim reduction
distPC <- 1-abs(WGCNA::bicor(MEs,use="p", maxPOutliers = 0.05))
distPC <- ifelse(is.na(distPC), 0, distPC)
MDS <- cmdscale(as.dist(distPC),2)
modNames = substring(names(MEs), 3)
col <- gsub(".*ME","", colnames(MEs))
plot(MDS, col=col, main="Multi-Dimensional Scaling (MDS) of Module Eigengenes (MEs)", cex=2, pch=19, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")


#Correlation plot of MEs
ComplexHeatmap::pheatmap(WGCNA::bicor(MEs,use="p", maxPOutliers = 0.05))





#=============Module-status cor===================#
#ME.covs <- cov_for.cor[rownames(cov_for.cor) %in% rownames(MEs),]
#ME.dat <- MEs_lin(MEs,ME.covs, return.r = T)
#ME.stim.R2 <- t(ME.dat["Concentration",])
#rownames(ME.stim.R2) <- "Concentration"
#ComplexHeatmap::pheatmap(ME.stim.R2, 
#                         cluster_rows = F,
#                         cluster_cols = F, annotation_legend = T, show_colnames = T, show_rownames = T,
#                         color = heatmap.color.code, scale='none')



ME.status.r <- t(rcorr(MEs, covs_forremBatch$Concentration, type="spearman")$r["y",][colnames(MEs)])
rownames(ME.status.r) <- "Concentration"
ComplexHeatmap::pheatmap(ME.status.r, 
                         cluster_rows = F,
                         cluster_cols = F, annotation_legend = T, show_colnames = T, show_rownames = T,
                         color = heatmap.color.code, scale='none')

identical(rownames(MEs), rownames(covs_forremBatch))

ME.status <- MEs_lin(MEs,cov.df=covs_forremBatch,return.r=T)
ComplexHeatmap::pheatmap(data.frame(t(ME.status["Concentration",])), 
                         cluster_rows = F,
                         cluster_cols = T, annotation_legend = T, show_colnames = T, show_rownames = T,
                         color = heatmap.color.code, scale='none')


#Using WGCNA bicor (similar to Spearman)
ME.concnetration.bi <- t(WGCNA::bicor(MEs, covs_forremBatch$Concentration,use="p", maxPOutliers = 0.05))
rownames(ME.concnetration.bi) <- "Concentration"
ComplexHeatmap::pheatmap(ME.concnetration.bi, 
                         cluster_rows = F,
                         cluster_cols = T, annotation_legend = T, show_colnames = T, show_rownames = T,
                         color = heatmap.color.code, scale="none")



#=============Network plots for cyan===================#
library(GGally)
library(igraph)

#https://kateto.net/networks-r-igraph

scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
cyan.hubs <- top25.hubgenes[top25.hubgenes$module == "cyan",]$name
tan.hubs <- top25.hubgenes[top25.hubgenes$module == "tan",]$name

cyan.expr <- sub.df[,colnames(sub.df) %in% cyan.hubs]
mod.coreg <- WGCNA::bicor(cyan.expr,use="p", maxPOutliers = 0.05)
mat <- as.matrix(as.dist(mod.coreg))
mat[is.na(mat)] <- 0

g <- graph.adjacency( #Generating a graph object of the genes within a module of interest
  mat,
  mode="directed",
  weighted=TRUE,
  diag=FALSE)

#g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE) #Removing multiple and looped edges

E(g)$weight <- abs(E(g)$weight)

E(g)$width <- E(g)$weight*4

#Set node sizes
colrs <- c("gray50", "tomato", "gold")
vSizes <- (scale01(apply(t(tan.expr), 1, mean)) + 1.0) * 10
V(g)$size <- vSizes*0.7
V(g)[c(cyan.hubs[cyan.hubs %in% c(DEGs.compound$W4,DEGs.compound$W5,DEGs.compound$W10)])]$Colors <- "tomato"

plot(g,edge.arrow.size=.3)

colrs[V(g)$DEGtype %in% "W10"]


V(g)[c(cyan.hubs[cyan.hubs %in% c(DEGs.compound$W4,DEGs.compound$W5,DEGs.compound$W10)])]$DEGtype <- "W10"
V(g)[c(cyan.hubs[!cyan.hubs %in% c(DEGs.compound$W4,DEGs.compound$W5,DEGs.compound$W10)])]$DEGtype  <- "Other"



#=============Supplementary plots for variance analysis===================#
#VarPar before and after
#prior to correction
identical(rownames(covs_forremBatch),colnames(vsd.complex))
regform <- ~ `Batch`+`Cell line`+`Timepoint`+`Stimulation`+`Percent. ERCC` + `Scaled percent. MT` + `Log(library size)` + `Sequencing batch`
varPart <- fitExtractVarPartModel(as.matrix(assay(vsd.complex)), regform, covs_forremBatch)
vp <- sortCols(varPart)
plotVarPart(vp)

#post correction
var.expl.cor <- fitExtractVarPartModel(batch.rem, regform, covs_forremBatch)
vp2 <- sortCols(var.expl.cor, decreasing = T)
plotVarPart(vp2)

#Same for PC-covariate linear correlation
library(factoextra)
PCA_cov_cor_R(cov.df=covs_forremBatch, df=as.matrix(assay(vsd.complex)))#pre

PCA_cov_cor_R(cov.df=covs_forremBatch, df=batch.rem)#post
