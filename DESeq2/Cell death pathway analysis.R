setwd("/Users/kubler01/Documents/R projects/Amber bulk RNA-seq of stimulated organoids")

load("12272021_MKMD WGCNA V10_1_updated figures for paper correct data.RData")


cd_markers <- c('BCL2', 'BAX', 'CASP3', 'CASP6')

for (i in c('GO:0008219','GO:0070265'))
cell_death <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
      filters = 'go', values = 'GO:0008219', mart = ensembl)

necrotic_cell_death <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
      filters = 'go', values = 'GO:0070265', mart = ensembl)

autophagic_cell_death <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
      filters = 'go', values = 'GO:0048102', mart = ensembl)

programmed_cell_death <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
      filters = 'go', values = 'GO:0012501', mart = ensembl)

programmed_necrotic_cell_death <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id','go_id'),
                               filters = "go", values = 'GO:0097300', mart = ensembl)

all.expr <- cbind(df.list[[1]],df.list[[2]],df.list[[3]])

cell.death_expr <- df.together[rownames(df.together) %in% unique(c(cell_death$hgnc_symbol,necrotic_cell_death$hgnc_symbol,autophagic_cell_death$hgnc_symbol,
                                                                          programmed_cell_death$hgnc_symbol,programmed_necrotic_cell_death$hgnc_symbol)),]

colors = colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(20)
annot_df.tog <- annot_df.tog[match(colnames(df.together),rownames(annot_df.tog)),]

ha <- columnAnnotation('Concentration'=annot_df.tog[,1],'Cell line'=annot_df.tog[,2],col=annot_colors)

annot_df.tog[rownames(annot_df.tog) %in% colnames(df.together),]

ComplexHeatmap::Heatmap(t(scale(t(cell.death_expr))), 
                        name='row Z score',
                        show_row_names = T,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        top_annotation = ha,
                        col = colors)

mM10 <- df.together[,colnames(df.together) %in% rownames(annot_df.tog[annot_df.tog$Concentration %in% '10mM',])]
mM50 <- df.together[,colnames(df.together) %in% rownames(annot_df.tog[annot_df.tog$Concentration %in% '50mM',])]

identical(rownames(meta.df.all),rownames(covs.required))

meta.new <- meta.df.all[rownames(meta.df.all) %in% rownames(covs.required),]
meta.new <- meta.new[match(rownames(covs.required),rownames(meta.new)),]
meta.new$DIV <- as.character(covs.required$DIV)
meta.new[meta.new$DIV %in% 1,]$DIV <- 'W4'
meta.new[meta.new$DIV %in% 2,]$DIV <- 'W5'
meta.new[meta.new$DIV %in% 3,]$DIV <- 'W10'
meta.new$DIV <- as.factor(meta.new$DIV)

meta.50mM <- meta.new[meta.new$`Sample Name` %in% colnames(mM50),]
meta.10mM <- meta.new[meta.new$`Sample Name` %in% colnames(mM10),]
meta.10mM <- meta.10mM[match(colnames(mM10),rownames(meta.10mM)),]
meta.50mM <- meta.50mM[match(colnames(mM50),rownames(meta.50mM)),]


annot_colors$DIV <- c('W4'='grey','W5'='#5ebe5a','W15'='darkgreen')
ha_10mM <- columnAnnotation('Concentration'=meta.10mM$Concentration,
                            'Cell line'=meta.10mM$`Cell line`,
                            'DIV'=meta.10mM$DIV, col=annot_colors)
ha_50mM <- columnAnnotation('Concentration'=meta.50mM$Concentration,
                            'Cell line'=meta.50mM$`Cell line`,
                            'DIV'=meta.50mM$DIV, col=annot_colors)
ComplexHeatmap::Heatmap(t(scale(t(mM10[rownames(mM10) %in% cd_markers,]))), 
                        name='row Z score',
                        show_row_names = T,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = T,
                        show_row_dend = F,
                        top_annotation = ha_10mM,
                        col = colors)

ComplexHeatmap::Heatmap(t(scale(t(mM50[rownames(mM50) %in% cd_markers,]))), 
                        name='row Z score',
                        show_row_names = T,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = T,
                        show_row_dend = F,
                        top_annotation = ha_50mM,
                        col = colors)

meta_1050 <- rbind(meta.10mM,meta.50mM)
mM1050 <- cbind(mM10,mM50)
library(tidyr)
meta_1050 <- drop_na(meta_1050)
meta_1050$DIV <- relevel(meta_1050$DIV,'W5')
meta_1050$DIV <- relevel(meta_1050$DIV,'W4')

mM1050 <- mM1050[,colnames(mM1050) %in% rownames(meta_1050)]
meta_1050 <- meta_1050[order(meta_1050$Concentration,meta_1050$DIV),]
mM1050 <- mM1050[,match(rownames(meta_1050),colnames(mM1050))]

all.df <- df.together[,colnames(df.together) %in% rownames(filt.meta)]
meta.all <- filt.meta[order(filt.meta$Concentration,filt.meta$DIV,filt.meta$`Cell line`),]

ha_1050<- columnAnnotation('Concentration'=meta.all$Concentration,
                                      'Cell line'=meta.all$`Cell line`,
                                      'DIV'=meta.all$DIV, col=annot_colors)

tiff("05062022_manual cell death genes_heatmap.tiff", width=1300, height=800, res=300)
ComplexHeatmap::Heatmap(t(scale(t(all.df[rownames(all.df) %in% cd_markers,]))), 
                        name='row Z score',
                        show_row_names = T,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        top_annotation = ha_1050,
                        col = colors)
dev.off()

tiff("05062022_all GO cell death genes_heatmap.tiff", width=1300, height=1600, res=300)
ComplexHeatmap::Heatmap(t(scale(t(drop_na(all.df[rownames(all.df) %in% rownames(cell.death_expr),])))), 
                        name='row Z score',
                        show_row_names = T,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        top_annotation = ha_1050,
                        row_names_gp = gpar(fontsize=6),
                        col = colors)
dev.off()


DEG.celldeat <- unique(c(DEGs.compound$W4,DEGs.compound$W5,DEGs.compound$W10)[c(DEGs.compound$W4,DEGs.compound$W5,DEGs.compound$W10) %in% rownames(cell.death_expr)])
tiff("05062022_GO cell death genes among DEGs_heatmap.tiff", width=1300, height=800, res=300)
ComplexHeatmap::Heatmap(t(scale(t(drop_na(all.df[rownames(all.df) %in% DEG.celldeat,])))), 
                        name='row Z score',
                        show_row_names = T,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        top_annotation = ha_1050,
                        col = colors)
dev.off()


