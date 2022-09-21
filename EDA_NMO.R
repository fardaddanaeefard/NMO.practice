################################################################################
#Installing required libraries
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("BiocParallel", "DESeq2", "BiocParallel",
#                       "apeglm", "org.Hs.eg.db", "enrichR"))
#install.packages(c("tidyverse", "gridExtra"))


################################################################################
#Loading required libraries
library(DESeq2)
library(tidyverse)
library(BiocParallel)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(vsn)
library(apeglm)
library(pheatmap)
library(gridExtra)
library(clusterProfiler)


################################################################################
#setting cpu parameter
register(SnowParam(15)) # set your cpu threads


#first we get a view on sample types
practice_data <- read.csv("raw_data.csv", header = F)
practice_info <- as.data.frame(t(practice_data[1:2, 1:57]))
colnames(practice_info) <- c("Treatment", "Antibody")
practice_info <- practice_info[-1,]
rownames(practice_info) <- practice_data[3, 2:57]
practice_info[is.na(practice_info)] <- "Untreated"
write.table(practice_info, file = "practice_info.csv",
            sep = ',', col.names = T, row.names = T, quote = F)



#later we process our gene expression data to become a proper table
practice_data <- practice_data[-1:-2,]
colnames(practice_data) <- practice_data[1,]
practice_data <- practice_data[-1,]
rownames(practice_data) <- practice_data$Ensembl_Gene_ID
practice_data <- practice_data[,-1]
practice_data <- data.matrix(practice_data)
write.table(practice_data, file = "practice_data.csv",
            sep = ',', col.names = T, row.names = T, quote = F)


################################################################################
#checking data
#are all rownames in info data matching column names in count data>
all(colnames(practice_data) %in% rownames(practice_info))

#are they in same order?
all(colnames(practice_data) == rownames(practice_info))


################################################################################
#change gene names
ens <- rownames(practice_data)
symbols <- mapIds(org.Hs.eg.db, keys = ens, column = c("SYMBOL")
                  , keytype = 'ENSEMBL')
symbols <- symbols[!is.na(ens)]
symbols <- symbols[match(rownames(practice_data), names(symbols))]
rownames(practice_data) <- symbols
keep <- !is.na(rownames(practice_data))
practice_data <- practice_data[keep,]


################################################################################
#factoring treatments and abs
practice_info$Treatment <- factor(practice_info$Treatment,
                                  levels = c("Untreated", "Ritux", "Other"))
practice_info$Antibody <- factor(practice_info$Antibody,
                                 levels = c("HC", "MOG", "AQP4"))


################################################################################
#now let's make a DESeqDataSet object for ab
dds_ab <- DESeqDataSetFromMatrix(countData = practice_data,
                                 colData = practice_info,
                                 design = ~ Antibody)


#another one for treatment
dds_trt <- DESeqDataSetFromMatrix(countData = practice_data,
                                  colData = practice_info,
                                  design = ~ Treatment)

################################################################################
#filtering genes with less than 10 count
dds_ab <- dds_ab[rowSums(counts(dds_ab) > 1) < 20,]
dds_trt <- dds_trt[rowSums(counts(dds_trt) > 1) < 20,]


################################################################################
#DEG extraction for Antibody factor
dds_ab <- DESeq(dds_ab, fitType = "local", parallel = T)

#QC of Antibody condition
vst_ab <- vst(dds_ab, blind = F)
plotPCA(vst_ab, intgroup = "Antibody")


#getting results for antibody condition
res_ab_aqp <- results(dds_ab, name = "Antibody_AQP4_vs_HC", alpha = 0.05)

res_ab_aqp <- lfcShrink(dds_ab, coef = "Antibody_AQP4_vs_HC", type = "apeglm")

res_ab_mog <- results(dds_ab, name = "Antibody_MOG_vs_HC", alpha = 0.05)

res_ab_mog <- lfcShrink(dds_ab, coef = "Antibody_MOG_vs_HC", type = "apeglm")


################################################################################
#DEG extraction for Treatment factor
dds_trt <- DESeq(dds_trt, fitType = "local", parallel = T)

#QC of treatment condition
vst_trt <- vst(dds_trt, blind = F)
plotPCA(vst_trt, intgroup = "Treatment")

#getting results for treatment condition
res_trt_rt <- results(dds_trt, name = "Treatment_Ritux_vs_Untreated",
                      alpha = 0.05)

res_trt_rt <- lfcShrink(dds_trt, coef = "Treatment_Ritux_vs_Untreated",
                        type = "apeglm")

res_trt_oth <- results(dds_trt, name = "Treatment_Other_vs_Untreated",
                      alpha = 0.05)

res_trt_oth <- lfcShrink(dds_trt, coef = "Treatment_Other_vs_Untreated",
                        type = "apeglm")

################################################################################
#getting significant DEGs
sig_trt_rt <- res_trt_rt[res_trt_rt$padj < 0.05 & res_trt_rt$baseMean > 50,]
sig_trt_oth <- res_trt_oth[res_trt_oth$padj < 0.05 & res_trt_oth$baseMean > 50,]
sig_ab_aqp <- res_ab_aqp[res_ab_aqp$padj < 0.05 & res_ab_aqp$baseMean > 50,]
sig_ab_mog <- res_ab_mog[res_ab_mog$padj < 0.05 & res_ab_mog$baseMean > 50,]



#volcano plots
volcano_plot_ab_aqp <- EnhancedVolcano(res_ab_aqp, x = "log2FoldChange",
              y = "pvalue", lab = rownames(res_ab_aqp), pCutoff = 0.05,
              FCcutoff = 0.57, legendPosition = "", title = "Anti-AQP4 VS HC",
              col=c('black', 'black', 'black', 'red3'), caption = "",
              subtitle = "") + coord_cartesian(xlim=c(-10, 10)
              ) + scale_x_continuous(breaks=seq(-10,10, 2))


volcano_plot_ab_mog <- EnhancedVolcano(res_ab_mog, x = "log2FoldChange",
              y = "pvalue", lab = rownames(res_ab_aqp), pCutoff = 0.05,
              FCcutoff = 0.57, legendPosition = "", title = "Anti-MOG VS HC",
              col=c('black', 'black', 'black', 'red3'), caption = "",
              subtitle = "") + coord_cartesian(xlim=c(-10, 10)
              ) + scale_x_continuous(breaks=seq(-10,10, 2))


volcano_plot_trt_rt <- EnhancedVolcano(res_trt_oth, x = "log2FoldChange",
              y = "pvalue", lab = rownames(res_trt_oth), pCutoff = 0.05,
              FCcutoff = 0.57, legendPosition = "", title = "Rituximab VS HC",
              col=c('black', 'black', 'black', 'red3'), caption = "",
              subtitle = "") + coord_cartesian(xlim=c(-10, 10)
              ) + scale_x_continuous(breaks=seq(-10,10, 2))


volcano_plot_trt_oth <- EnhancedVolcano(res_trt_oth, x = "log2FoldChange",
              y = "pvalue", lab = rownames(res_ab_aqp), pCutoff = 0.05,
              FCcutoff = 0.57, legendPosition = "", subtitle = "",
              title = "Other treatments VS HC", , caption = "",
              col=c('black', 'black', 'black', 'red3')) + coord_cartesian(
              xlim=c(-10, 10)) + scale_x_continuous(breaks=seq(-10,10, 2))


#combile all plots
volcano_plot_all <- grid.arrange(volcano_plot_ab_aqp, volcano_plot_ab_mog,
                                 volcano_plot_trt_rt, volcano_plot_trt_oth)

#save it as pdf
ggsave(plot = volcano_plot_all, filename = "Volcano Plot.pdf", width = 35,
       height = 20)


#Enrichment analysis
#GO enrichment
GO_results_aqp_all <- enrichGO(gene = rownames(sig_ab_aqp), OrgDb = org.Hs.eg.db,
                               keyType = "SYMBOL", ont = "ALL")
plot(barplot(GO_results_aqp_all, showCategory = 20, title = "AQP_VS_HC"))


GO_results_mog_all <- enrichGO(gene = rownames(sig_ab_mog), OrgDb = org.Hs.eg.db,
                               keyType = "SYMBOL", ont = "ALL")
plot(barplot(GO_results_mog_all, showCategory = 20, title = "MOG_VS_HC"))


GO_results_rit_all <- enrichGO(gene = rownames(sig_trt_rt), OrgDb = org.Hs.eg.db,
                               keyType = "SYMBOL", ont = "ALL")
plot(barplot(GO_results_rit_all, showCategory = 20, title = "RIT_VS_HC"))


GO_results_oth_all <- enrichGO(gene = rownames(sig_trt_oth), OrgDb = org.Hs.eg.db,
                               keyType = "SYMBOL", ont = "ALL")
plot(barplot(GO_results_oth_all, showCategory = 20, title = "OTH_VS_HC"))



################################################################################
#AQP
ranks_aqp <- res_ab_aqp$log2FoldChange
names(ranks_aqp) <- rownames(res_ab_aqp)
title(main = "AQP", barplot(sort(ranks_aqp, decreasing = T)))

#MOG
ranks_mog <- res_ab_mog$log2FoldChange
names(ranks_mog) <- rownames(res_ab_mog)
title(main = "MOG", barplot(sort(ranks_mog, decreasing = T)))

#Ritux
ranks_rit <- res_trt_rt$log2FoldChange
names(ranks_rit) <- rownames(res_trt_rt)
title(main = "Rituximab", barplot(sort(ranks_rit, decreasing = T)))

#Other
ranks_oth <- res_trt_oth$log2FoldChange
names(ranks_oth) <- rownames(res_trt_oth)
title(main = "Other Treatments", barplot(sort(ranks_oth, decreasing = T)))

################################################################################
#heatmap plotting
ntd <- normTransform(dds_ab)
vsd <- vst(dds_ab, blind = F)
rlg <- rlog(dds_ab, blind = F)
select <- order(rowMeans(counts(dds_ab,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_ab)[,c("Treatment", "Antibody")])

ntd_heatmap <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE,
                        cluster_cols=FALSE, annotation_col=df)

vsd_heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE,
                        cluster_cols=FALSE, annotation_col=df)

rlg_heatmap <- pheatmap(assay(rlg)[select,], cluster_rows=FALSE,
                        cluster_cols=FALSE, annotation_col=df)

#combination PCA
plotPCA(vsd, intgroup = c("Antibody", "Treatment"))
#save R image
save.image(".RData")
################################################################################
