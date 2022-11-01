#setting work directory
setwd("")

################################################################################
#Loading required libraries
library(vsn)
library(DESeq2)
library(apeglm)
library(biomaRt)
library(pheatmap)
library(tidyverse)
library(gridExtra)
library(enrichplot)
library(BiocParallel)
library(org.Hs.eg.db)
library(clusterProfiler)
library(EnhancedVolcano)
library(EnsDb.Hsapiens.v86)

################################################################################
#setting CPU parameter
register(SnowParam(10))


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
write.table(practice_data, file = "extdata/practice_data.csv",
            sep = ',', col.names = T, row.names = T, quote = F)


################################################################################
#checking data
#are all rownames in info data matching column names in count data>
all(colnames(practice_data) %in% rownames(practice_info))

#are they in same order?
all(colnames(practice_data) == rownames(practice_info))


################################################################################
#factoring treatments and abs
practice_info$Treatment <- factor(practice_info$Treatment,
                levels = c("Untreated", "Ritux", "Other"))
practice_info$Antibody <- factor(practice_info$Antibody,
                 levels = c("HC", "MOG", "AQP4"))


################################################################################
#rearranging datasets
practice_info <- practice_info[order(rownames(practice_info)),]
practice_data <- practice_data[,rownames(practice_info)]


################################################################################
#now let's make a DESeqDataSet object for ab
raw_dds_ab <- DESeqDataSetFromMatrix(countData = practice_data,
                                 colData = practice_info,
                                 design = ~ Antibody)


#another one for treatment
raw_dds_trt <- DESeqDataSetFromMatrix(countData = practice_data,
                                  colData = practice_info,
                                  design = ~ Treatment)



################################################################################
#DEG extraction for Antibody factor
dds_ab <- DESeq(raw_dds_ab, sfType = "poscounts", fitType = "local", parallel = T)

#QC of Antibody condition
vst_ab <- vst(dds_ab, fitType = "local")
PCA_ab <- plotPCA(vst_ab, intgroup = "Antibody")
plotMA(dds_ab, ylim = c(-10, 10))

#getting results for antibody condition
res_ab_aqp <- results(dds_ab, name = "Antibody_AQP4_vs_HC")

res_ab_aqp <- lfcShrink(dds_ab, coef = "Antibody_AQP4_vs_HC",
                      type = "apeglm")

res_ab_mog <- results(dds_ab, name = "Antibody_MOG_vs_HC")

res_ab_mog <- lfcShrink(dds_ab, coef = "Antibody_MOG_vs_HC",
                      type = "apeglm")

summary(res_ab_aqp)
summary(res_ab_mog)

################################################################################
#DEG extraction for Treatment factor
dds_trt <- DESeq(raw_dds_trt, sfType = "poscounts", fitType = "local", parallel = T)

#QC of treatment condition
vst_trt <- vst(dds_trt, fitType = "local")
PCA_trt <- plotPCA(vst_trt, intgroup = "Treatment")

#PCA of both factors
PCA_all <- grid.arrange(PCA_ab, PCA_trt, nrow = 1)

#getting results for treatment condition
res_trt_rt <- results(dds_trt, name = "Treatment_Ritux_vs_Untreated")

res_trt_rt <- lfcShrink(dds_trt, coef = "Treatment_Ritux_vs_Untreated",
                        type = "apeglm")


res_trt_oth <- results(dds_trt, name = "Treatment_Other_vs_Untreated")

res_trt_oth <- lfcShrink(dds_trt, coef = "Treatment_Other_vs_Untreated",
                        type = "apeglm")

summary(res_trt_rt)
summary(res_trt_oth)

################################################################################
#gene name conversion
geneIDs_aqp <- ensembldb::select(EnsDb.Hsapiens.v86,
                                 keys = rownames(res_ab_aqp), keytype = "GENEID",
                                 columns = c("SYMBOL", "GENEID"))

res_ab_aqp <- res_ab_aqp[rownames(res_ab_aqp) %in% geneIDs_aqp$GENEID,]
res_ab_aqp$symbol <- geneIDs_aqp$SYMBOL

geneIDs_mog <- ensembldb::select(EnsDb.Hsapiens.v86,
                                 keys = rownames(res_ab_mog), keytype = "GENEID",
                                 columns = c("SYMBOL", "GENEID"))

res_ab_mog <- res_ab_mog[rownames(res_ab_mog) %in% geneIDs_mog$GENEID,]
res_ab_mog$symbol <- geneIDs_mog$SYMBOL


geneIDs_rt <- ensembldb::select(EnsDb.Hsapiens.v86,
                                keys = rownames(res_trt_rt), keytype = "GENEID",
                                columns = c("SYMBOL", "GENEID"))

res_trt_rt <- res_trt_rt[rownames(res_trt_rt) %in% geneIDs_rt$GENEID,]
res_trt_rt$symbol <- geneIDs_rt$SYMBOL

geneIDs_oth <- ensembldb::select(EnsDb.Hsapiens.v86,
                                 keys = rownames(res_trt_oth), keytype = "GENEID",
                                 columns = c("SYMBOL", "GENEID"))

res_trt_oth <- res_trt_oth[rownames(res_trt_oth) %in% geneIDs_oth$GENEID,]
res_trt_oth$symbol <- geneIDs_oth$SYMBOL

################################################################################
#getting significant DEGs
sig_trt_rt <- res_trt_rt[!is.na(res_trt_rt$padj) & res_trt_rt$padj < .1 &
            abs(res_trt_rt$log2FoldChange) > 1 & res_trt_rt$baseMean > 5,]
sig_trt_rt <- sig_trt_rt[order(sig_trt_rt$log2FoldChange, decreasing = T),]


sig_trt_oth <- res_trt_oth[!is.na(res_trt_oth$padj) & res_trt_oth$padj < .1 &
            abs(res_trt_oth$log2FoldChange) > 1 & res_trt_oth$baseMean > 5,]
sig_trt_oth <- sig_trt_oth[order(sig_trt_oth$log2FoldChange, decreasing = T),]


sig_ab_aqp <- res_ab_aqp[!is.na(res_ab_aqp$padj) & res_ab_aqp$padj < .1 &
            abs(res_ab_aqp$log2FoldChange) > 1 & res_ab_aqp$baseMean > 5,]
sig_ab_aqp <- sig_ab_aqp[order(sig_ab_aqp$log2FoldChange, decreasing = T),]


sig_ab_mog <- res_ab_mog[!is.na(res_ab_mog$padj) & res_ab_mog$padj < .1 &
            abs(res_ab_mog$log2FoldChange) > 1 & res_ab_mog$baseMean > 5,]
sig_ab_mog <- sig_ab_mog[order(sig_ab_mog$log2FoldChange, decreasing = T),]


################################################################################
#volcano plots
volcano_plot_ab_aqp <- EnhancedVolcano(res_ab_aqp, x = "log2FoldChange",
              y = "pvalue", lab = res_ab_aqp$symbol, pCutoff = 0.05,
              FCcutoff = 0.57, legendPosition = "", title = "Anti-AQP4 VS HC",
              col=c('black', 'black', 'black', 'red3'), caption = "",
              subtitle = "") + coord_cartesian(xlim=c(-10, 10)
              ) + scale_x_continuous(breaks=seq(-10,10, 2))


volcano_plot_ab_mog <- EnhancedVolcano(res_ab_mog, x = "log2FoldChange",
              y = "pvalue", lab = res_ab_mog$symbol, pCutoff = 0.05,
              FCcutoff = 0.57, legendPosition = "", title = "Anti-MOG VS HC",
              col=c('black', 'black', 'black', 'red3'), caption = "",
              subtitle = "") + coord_cartesian(xlim=c(-10, 10)
              ) + scale_x_continuous(breaks=seq(-10,10, 2))


volcano_plot_trt_rt <- EnhancedVolcano(res_trt_rt, x = "log2FoldChange",
              y = "pvalue", lab = res_trt_rt$symbol, pCutoff = 0.05,
              FCcutoff = 0.57, legendPosition = "", title = "Rituximab VS HC",
              col=c('black', 'black', 'black', 'red3'), caption = "",
              subtitle = "") + coord_cartesian(xlim=c(-10, 10)
              ) + scale_x_continuous(breaks=seq(-10,10, 2))


volcano_plot_trt_oth <- EnhancedVolcano(res_trt_oth, x = "log2FoldChange",
              y = "pvalue", lab = res_trt_oth$symbol, pCutoff = 0.05,
              FCcutoff = 0.57, legendPosition = "", subtitle = "",
              title = "Other treatments VS HC", caption = "",
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
d_aqp <- dotplot(enrichGO(sort(sig_ab_aqp[sig_ab_aqp$log2FoldChange < 0,"symbol"],
                  decreasing = T), OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL", ont = "ALL"),
                  title = "Downregulated genes in anti-AQP4 patients")
#u_aqp <- dotplot(enrichGO(sort(sig_ab_aqp[sig_ab_aqp$log2FoldChange < -1,"symbol"],
#                         decreasing = T), OrgDb = org.Hs.eg.db,
#                    keyType = "SYMBOL", ont = "ALL"),
#           title = "Upregulated genes in anti-AQP4 patients")
##
d_mog<- dotplot(enrichGO(sort(sig_ab_mog[sig_ab_mog$log2FoldChange < 0,"symbol"],
                         decreasing = T), OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL", ont = "ALL"),
           title = "Downregulated genes in anti-MOG patients")
#u_mog<- dotplot(enrichGO(sort(sig_ab_mog[sig_ab_mog$log2FoldChange > 1,"symbol"],
#                         decreasing = T), OrgDb = org.Hs.eg.db,
#                    keyType = "SYMBOL", ont = "ALL"),
#           title = "Upregulated genes in anti-MOG patients")
##
d_rt <- dotplot(enrichGO(sort(sig_trt_rt[sig_trt_rt$log2FoldChange < 0,"symbol"],
                         decreasing = T), OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL", ont = "ALL"),
           title = "Downregulated genes in Rituximab patients")

#u_rt <- dotplot(enrichGO(sort(sig_trt_rt[sig_trt_rt$log2FoldChange > 1,"symbol"],
#                         decreasing = T), OrgDb = org.Hs.eg.db,
#                    keyType = "SYMBOL", ont = "ALL"),
#           title = "Upregulated genes in Rituximab patients")
##

d_oth <- dotplot(enrichGO(sort(sig_trt_oth[sig_trt_oth$log2FoldChange < 0,"symbol"],
                         decreasing = T), OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL", ont = "ALL"),
           title = "Downregulated genes in patients with other treatments")

u_oth <- dotplot(enrichGO(sort(sig_trt_oth[sig_trt_oth$log2FoldChange > 0,"symbol"],
                               decreasing = T), OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL", ont = "ALL"),
                 title = "Upregulated genes in patients with other treatments")

downs <- grid.arrange(d_aqp, d_mog, d_rt, d_oth)
################################################################################
#AQP
ranks_aqp <- sig_ab_aqp$log2FoldChange
names(ranks_aqp) <- sig_ab_aqp$symbol
title(main = "Anti-AQP", barplot(sort(ranks_aqp, decreasing = T)))

#MOG
ranks_mog <- sig_ab_mog$log2FoldChange
names(ranks_mog) <- sig_ab_mog$symbol
title(main = "Anti-MOG", barplot(sort(ranks_mog, decreasing = T)))

#Ritux
ranks_rit <- sig_trt_rt$log2FoldChange
names(ranks_rit) <- sig_trt_rt$symbol
title(main = "Rituximab", barplot(sort(ranks_rit, decreasing = T)))

#Other
ranks_oth <- sig_trt_oth$log2FoldChange
names(ranks_oth) <- sig_trt_oth$symbol
title(main = "Other Treatments", barplot(sort(ranks_oth, decreasing = T)))

################################################################################
#heatmap plotting
ntd <- normTransform(dds_ab)
vsd <- vst(dds_ab, blind = F)
rlg <- rlog(dds_ab, blind = F)
select <- order(rowMeans(counts(dds_ab,normalized=TRUE)),
                decreasing=TRUE)[-1:-20]
df <- as.data.frame(colData(dds_ab)[,c("Treatment", "Antibody")])

ntd_heatmap <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE,
                        cluster_cols=FALSE, annotation_col=df)

vsd_heatmap <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE,
                        cluster_cols=FALSE, annotation_col=df)

rlg_heatmap <- pheatmap(assay(rlg)[select,], cluster_rows=FALSE,
                        cluster_cols=FALSE, annotation_col=df)

#combination PCA
plotPCA(vsd, intgroup = c("Antibody", "Treatment"))

#dendrogram of samples
normCounts_ab <- counts(dds_ab, normalized = T)
normCounts_trt <- counts(dds_trt, normalized = T)

plot(hclust(dist(t(normCounts_ab))),
            main = "Clustering of samples based on antibodies",
            labels = colData(dds_ab)$Antibody, sub = "", xlab = "")
plot(hclust(dist(t(normCounts_trt))),
           main = "Clustering of samples based on treatments",
           labels = colData(dds_ab)$Treatment, sub = "", xlab = "")
#save R image
save.image(".RData")
################################################################################
#plotting heatmap for the 25 genes of interest
goi <- c("IL1RN", "IFIT5", "IFIT3", "IFIT2", "IFI6", "HERC5", "CMPK2", "IFIT1",
         "MX1", "ISG15", "OASL", "RSAD2", "IFI44", "XAF1", "EIF2AK2", "LAP3",
         "LY6E", "OAS2", "OAS1", "USP18", "HERC6", "SIGLEC1", "IFI27", "OAS3",
         "IFI44L")


geneIds <- ensembldb::select(EnsDb.Hsapiens.v86,
                  keys = rownames(practice_data), keytype = "GENEID",
                  columns = c("SYMBOL", "GENEID"))

practice_data2 <- practice_data[rownames(practice_data) %in% geneIds$GENEID,]
rownames(practice_data2) <- geneIds$SYMBOL
practice_data3 <- practice_data2[rownames(practice_data2) %in% goi,]
v <- vst(practice_data3, fitType = "local", nsub = 25)
pheatmap(v, cluster_rows = F,
         annotation = as.data.frame(colData(dds_ab)[,c("Treatment", "Antibody")]))

