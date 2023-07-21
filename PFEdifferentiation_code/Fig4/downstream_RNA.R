# install packages

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("DESeq2")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("tximport")


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("tximeta")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("EnsDb.Hsapiens.v86")

# install.packages("extrafont")


if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)

library(EnsDb.Hsapiens.v86)
# library(extrafont)
# library(tximport)
library(tximeta)
library(SummarizedExperiment)
library("magrittr")
library("DESeq2")
library(impute)
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("glmpca")
library("ggbeeswarm")
library("genefilter")
library(MASS)
library(ensembldb)

library("AnnotationDbi")
library("org.Hs.eg.db")
library(UpSetR)


# install.packages('magick')
library(magick)

dir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA"
samples <- read.table(file.path(dir,"final_files_no_tep"), header=TRUE)
samples

files <- file.path(dir, "rsem_output",  samples$filename)
all(file.exists(files))


gsemat <- tximeta(files, type = "rsem", txIn = FALSE, txOut = FALSE)


gsemat$sample = factor(samples$sample)
gsemat$replicate = factor(samples$rep)

round( colSums(assay(gsemat)) / 1e6, 1 )

# Filtering
nrow(gsemat)
gsemat = gsemat[rowSums(gsemat@assays@data$counts) > 0,]
nrow(gsemat)
gsemat = gsemat[rowSums(gsemat@assays@data$length) > 0,]
nrow(gsemat)

assays(gsemat)$length[ assays(gsemat)$length == 0] <- NA # set these as missing
idx <- rowSums(is.na(assays(gsemat)$length)) > 6 #must have counts in at least 2 sample
table(idx)
gsemat <- gsemat[!idx,]
dim(gsemat)
length_imp <- impute.knn(assays(gsemat)$length)
assays(gsemat)$length <- length_imp$data

dds_no_batch <- DESeqDataSet(gsemat, design = ~ sample)

dds_batch <- DESeqDataSet(gsemat, design = ~ replicate + sample)


# filtering DESeq object
keep <- rowSums(counts(dds_no_batch) >= 3) >= 2
dds_no_batch = dds_no_batch[keep,]
dim(dds_no_batch)


dds_batch = dds_batch[keep,]
dim(dds_batch)

# this is for clustering & PCA analysis 
vsd_batch = vst(dds_batch,blind=FALSE)
head(assay(vsd_batch), 3)

vsd_no_batch = vst(dds_no_batch,blind=FALSE)
head(assay(vsd_no_batch), 3)

rld_batch <- rlog(dds_batch, blind = FALSE)
head(assay(rld_batch), 3)

rld_no_batch <- rlog(dds_no_batch, blind = FALSE)
head(assay(rld_no_batch), 3)


# < 30 samples, so use rld

dds_no_batch <- estimateSizeFactors(dds_no_batch)

df <- bind_rows(
  as_data_frame(log2(counts(dds_no_batch, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd_no_batch)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld_no_batch)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


dds_batch <- estimateSizeFactors(dds_batch)

df <- bind_rows(
  as_data_frame(log2(counts(dds_batch, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd_batch)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld_batch)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


dev.off()
dev.off()


# use rlog ( recommended for < 30 samples anyway )

assay_mat = assay(rld_no_batch)
colnames(assay_mat) = paste0(rld_no_batch$sample, rld_no_batch$replicate)
order_fields = c('hESrep1','hESrep2','hDErep1','hDErep2','hAFErep1','hAFErep2','hPFErep1','hPFErep2')
assay_mat = assay_mat[,order_fields]
sampleDists <- dist(t(assay_mat))

sampleDistMatrix <- as.matrix( sampleDists, labels=TRUE )
#rownames(sampleDistMatrix) <- order_fields #paste0(rld_no_batch$sample, rld_no_batch$replicate)
#colnames(sampleDistMatrix) <- order_fields #paste0(rld_no_batch$sample, rld_no_batch$replicate)
# reorder rows
#order_fields = c('hESrep1','hESrep2','hDErep1','hDErep2','hAFErep1','hAFErep2','hPFErep1','hPFErep2')
#sampleDistMatrix = sampleDistMatrix[order_fields,order_fields]
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# font_import()

p = pheatmap(sampleDistMatrix,fontsize=12,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists, # cluster_rows = F,cluster_columns = F,
             col = colors,main='Sample Correlation Distance')


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


savedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA//downstream_figures/"
save_pheatmap_pdf(p, file.path(savedir, "rna_correlation_rld.pdf"))

dev.off()
dev.off()

# PCA

# pcaData <- plotPCA(rld_batch, intgroup = c( "replicate", "sample"), returnData = TRUE)
# pcaData
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# p = ggplot(pcaData, aes(x = PC1, y = PC2, color = sample, shape = replicate)) +
#   geom_point(size =3) +
#   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
#   coord_fixed() +
#   ggtitle("PCA with RLD data")
# ggsave(file.path(savedir, "pca_rld_batch.pdf"),plot=p,width=10, height=10)


# dev.off()
pcaData <- plotPCA(rld_no_batch, intgroup = c( "replicate", "sample"), returnData = TRUE)
pcaData$sample <- factor(pcaData$sample, levels = c("hES","hDE","hAFE","hPFE")) 
  
percentVar <- round(100 * attr(pcaData, "percentVar"))
p = ggplot(pcaData, aes(x = PC1, y = PC2, color = sample, shape = replicate)) + 
  scale_color_manual(values = c("red", "yellow", "skyblue", "orchid")) + 
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA on RLD normalized expression")  + theme_bw()
  
  #theme(# axis.line=element_blank(),axis.text.x=element_blank(),
      # axis.text.y=element_blank(),axis.ticks=element_blank(),
      # axis.title.x=element_blank(),
      # axis.title.y=element_blank(),#legend.position="none",
   #   panel.background=element_blank(),
   #   panel.border=element_blank(),
   #   panel.grid.major=element_blank(),
    #  panel.grid.minor=element_blank())




ggsave(file.path(savedir, "pca_rld.pdf"),plot=p,width=5, height=5)




# dev.off()


# PCA + correlation - done from rld normalization 


# # GLMPCA
# gpca <- glmpca(counts(dds_all), L=2)
# gpca.dat <- gpca$factors
# gpca.dat$batch <- dds_all$batch
# gpca.dat$assay <- dds_all$assay
# gpca.dat$timepoint <- dds_all$timepoint
# gpca.dat$cell_line <- dds_all$cell_line
# 
# p = ggplot(gpca.dat, aes(x = dim1, y = dim2, color = timepoint, shape = cell_line)) +
#   geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")
# 
# ggsave(file.path(savedir, "glmpca_all_samples.pdf"),plot=p,width=10, height=10)
# 
# dev.off()
# dev.off()
# 
# 
# mds <- as.data.frame(colData(vsd))  %>%
#   cbind(cmdscale(sampleDistMatrix))
# p = ggplot(mds, aes(x = `1`, y = `2`,  color = timepoint, shape = cell_line)) +
#   geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
# ggsave(file.path(savedir, "mds_all_samples.pdf"),plot=p,width=10, height=10)
# 
# 
# dev.off()
# dev.off()
# 
# # MDS
# mdsPois <- as.data.frame(colData(dds_all)) %>%
#   cbind(cmdscale(samplePoisDistMatrix))
# ggplot(mdsPois, aes(x = `1`, y = `2`, color = timepoint, shape = cell_line)) +
#   geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")
# ggsave(file.path(savedir, "mds_all_samples_poisson.pdf"),plot=p,width=10, height=10)


# DESEQ2 - with batch correction
dds_all <- DESeq(dds_batch) # size factor estimation is done internally of this 

# dds_no_batch <- DESeq(dds_no_batch)


#res_all <- results(dds_all,alpha=0.01,lfcThreshold = 1) # comparison between TEP and AFE
#res_no_batch <- results(dds_no_batch,alpha=0.01,lfcThreshold = 1)



# pairwise comparisons
res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0.5,contrast = c("sample","hDE","hES")) # comparison between TEP and AFE
resSig_all <- subset(res_all, padj < 0.1)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(resSig_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
resSig_all$geneid = substr(rownames(resSig_all),1,15)
resSig_all@listData = merge(resSig_all@listData,geneIDs1,by.x="geneid",by.y="GENEID")
write.csv(resSig_all@listData,"/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/differential_csvs/DE_over_ESC_differential.csv")


res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0.5,contrast = c("sample","hAFE","hES")) # comparison between TEP and AFE
resSig_all <- subset(res_all, padj < 0.1)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(resSig_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
resSig_all$geneid = substr(rownames(resSig_all),1,15)
resSig_all@listData = merge(resSig_all@listData,geneIDs1,by.x="geneid",by.y="GENEID")
write.csv(resSig_all@listData,"/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/differential_csvs/AFE_over_ESC_differential.csv")


res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0.5,contrast = c("sample","hAFE","hDE")) # comparison between TEP and AFE
resSig_all <- subset(res_all, padj < 0.1)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(resSig_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
resSig_all$geneid = substr(rownames(resSig_all),1,15)
resSig_all@listData = merge(resSig_all@listData,geneIDs1,by.x="geneid",by.y="GENEID")
write.csv(resSig_all@listData,"/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/differential_csvs/AFE_over_DE_differential.csv")



res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0.5,contrast = c("sample","hPFE","hDE")) # comparison between TEP and AFE
resSig_all <- subset(res_all, padj < 0.1)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(resSig_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
resSig_all$geneid = substr(rownames(resSig_all),1,15)
resSig_all@listData = merge(resSig_all@listData,geneIDs1,by.x="geneid",by.y="GENEID")
write.csv(resSig_all@listData,"/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/differential_csvs/PFE_over_DE_differential.csv")


res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0.5,contrast = c("sample","hPFE","hAFE")) # comparison between TEP and AFE
resSig_all <- subset(res_all, padj < 0.1)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(resSig_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
resSig_all$geneid = substr(rownames(resSig_all),1,15)
resSig_all@listData = merge(resSig_all@listData,geneIDs1,by.x="geneid",by.y="GENEID")
write.csv(resSig_all@listData,"/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/differential_csvs/PFE_over_AFE_differential.csv")


res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0.5,contrast = c("sample","hPFE","hES")) # comparison between TEP and AFE
resSig_all <- subset(res_all, padj < 0.1)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(resSig_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
resSig_all$geneid = substr(rownames(resSig_all),1,15)
resSig_all@listData = merge(resSig_all@listData,geneIDs1,by.x="geneid",by.y="GENEID")
write.csv(resSig_all@listData,"/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/differential_csvs/PFE_over_ESC_differential.csv")


library(UpSetR)

# read in differential expression files ( wrote in python )
basedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA//differential_csvs_less_stringent/'
de_over_esc <- read.csv( file.path(basedir, "DE_over_ESC_differential.csv"), header=T, sep="," )
de_over_esc_up = de_over_esc[de_over_esc$log2FoldChange>0,]
de_over_esc_down = de_over_esc[de_over_esc$log2FoldChange<0,]

afe_over_esc <- read.csv( file.path(basedir, "AFE_over_ESC_differential.csv"), header=T, sep="," )
afe_over_esc_up = afe_over_esc[afe_over_esc$log2FoldChange>0,]
afe_over_esc_down = afe_over_esc[afe_over_esc$log2FoldChange<0,]

afe_over_de <- read.csv( file.path(basedir, "AFE_over_DE_differential.csv"), header=T, sep="," )
afe_over_de_up = afe_over_de[afe_over_de$log2FoldChange>0,]
afe_over_de_down = afe_over_de[afe_over_de$log2FoldChange<0,]


pfe_over_esc <- read.csv( file.path(basedir, "PFE_over_ESC_differential.csv"), header=T, sep="," )
pfe_over_esc_up = pfe_over_esc[pfe_over_esc$log2FoldChange>0,]
pfe_over_esc_down = pfe_over_esc[pfe_over_esc$log2FoldChange<0,]


pfe_over_de <- read.csv( file.path(basedir, "PFE_over_DE_differential.csv"), header=T, sep="," )
pfe_over_de_up = pfe_over_de[pfe_over_de$log2FoldChange>0,]
pfe_over_de_down = pfe_over_de[pfe_over_de$log2FoldChange<0,]


pfe_over_afe <- read.csv( file.path(basedir, "PFE_over_AFE_differential.csv"), header=T, sep="," )
pfe_over_afe_up = pfe_over_afe[pfe_over_afe$log2FoldChange>0,]
pfe_over_afe_down = pfe_over_afe[pfe_over_afe$log2FoldChange<0,]


genes_to_label = c('NANOG','POU5F1','SOX2','SOX17','FOXA2','PAX1','PAX9','HOXA1','HOXA2','HOXA3','HOXA4',
                   'HOXB1','HOXB2','HOXB3','HOXB4','HOXC1','HOXC2','HOXC3','HOXC4','HOXD1','HOXD2','HOXD3','HOXD4',
                   'LEFTY1','LEFTY2','CER1','VGLL2','TBX1','HIC1','PBX1','PBX2','MEIS1','MEIS2','MEIS3','SIX1','GBX2','OTX1','OTX2','TFAP4','TFAP2C',
                   'ISL1','ISL2','IRX1','IRX2','IRX3','EYA1','FGF10','HEY1','BMP4','BMP1','BMP2', 'BMP3','BMP5','BMP6','BMP7','HES1','ID1','FN1','EOMES','NODAL','FGF17','SMAD7','FOXE1',
                   'NKX2-1','NKX2-2','NKX2-3','NKX2-5','NKX2-6','RORA','RORB','RORC','GATA3','GATA4','GATA6','TP53','TP63','SHH','DMRT2','FOXP1','FOXP2','WNT5A','ROR2')


# make volcano plots for RNA pairwise
rna_volcano = function(ip,op_file, logfc=1.0,padj=0.05){
  ip$diffexpressed <- "NO"
  # if log2Foldchange > 1.0 and FDR < 0.01, set as "UP" 
  ip$diffexpressed[ip$log2FoldChange > logfc & ip$padj < padj] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  ip$diffexpressed[ip$log2FoldChange < -logfc & ip$padj < padj] <- "DOWN"
  
  genes_to_label_ex = genes_to_label
  
  tmp1 = ip[(ip$log2FoldChange < -logfc & ip$padj < padj),]
  tmp1 = tmp1[order(tmp1$log2FoldChange),]
  genes_to_label_ex = c(genes_to_label_ex,head(tmp1$SYMBOL,2))
  tmp1 = tmp1[order(tmp1$padj),]
  genes_to_label_ex = c(genes_to_label_ex,head(tmp1$SYMBOL,2))
  
  
  tmp2 = ip[(ip$log2FoldChange > logfc & ip$padj < padj),]
  tmp2 = tmp2[order(tmp2$log2FoldChange),]
  genes_to_label_ex = c(genes_to_label_ex,tail(tmp2$SYMBOL,5))
  tmp2 = tmp2[order(tmp2$padj),]
  genes_to_label_ex = c(genes_to_label_ex,head(tmp2$SYMBOL,5))
  
  
  
  
  # ip$label <- "NO"
  ip$SYMBOL[!(ip$SYMBOL %in% genes_to_label_ex)] <- ""
  
  # label top 5 logfc + pval
  
  
  library(ggrepel)
  # plot adding up all layers we have seen so far
  ggplot(data=ip, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=SYMBOL)) + xlim(-30,30) + 
    geom_point() + 
    theme_minimal() + #geom_text_repel(data=subset(anno, label == "YES"),
    #          aes(x=Fold, y=-log10(FDR), col=diffexpressed,label=symbol)) + 
    geom_text_repel(max.overlaps = Inf) +
    scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-logfc, logfc), col="red") +
    geom_hline(yintercept=-log10(padj), col="red")
  
  ggsave(op_file)
}


rna_volcano(de_over_esc,file.path(basedir,"DE_over_ESC.pdf"))
rna_volcano(afe_over_esc,file.path(basedir,"AFE_over_ESC.pdf"))
rna_volcano(afe_over_de,file.path(basedir,"AFE_over_DE.pdf"))
rna_volcano(pfe_over_esc,file.path(basedir,"PFE_over_ESC.pdf"))
rna_volcano(pfe_over_de,file.path(basedir,"PFE_over_DE.pdf"))
rna_volcano(pfe_over_afe,file.path(basedir,"PFE_over_AFE.pdf"))



basedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/downstream_figures"
genesets <- read.csv( file.path(basedir, "upset_ip_less_stringent.csv"), header=T, sep = ",")

to_plot = colnames(genesets)
to_plot = to_plot[to_plot != 'X']


dev.off()

pdf(file.path(basedir,'upset_plot_differential_gene_expr.pdf'))


upset(genesets, sets=to_plot, #sets.bar.color = c("#FF0000","#0000FF","#FF0000","#0000FF","#FF0000","#0000FF","#FF0000","#0000FF","#FF0000","#0000FF","#FF0000","#0000FF"),
      order.by = "freq", empty.intersections = "on")

dev.off()





# use genes differential in at least 1 pair
genes_to_plot = genesets$X
length(genes_to_plot)

library(dbplyr)
library(magrittr)
# install.packages("sets")
library(sets)

dim(rld_no_batch)
# instead of hvgs use all sets of DEGs
rld_no_batch_var = rld_no_batch[substr(rownames(assay(rld_no_batch)),1,15) %in% genes_to_plot]
dim(rld_no_batch_var)


mat <- assay(rld_no_batch)

mat  <- assay(rld_no_batch_var)


tmpcol = colData(rld_no_batch_var)
colnames(mat) = paste(tmpcol$sample,'_',tmpcol$replicate,sep='')

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(mat),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
tmp_rowname = substr(rownames(mat),1,15)
length(tmp_rowname)
length(as.set(tmp_rowname)) # duplicates exist, make unique
rownames(mat) = substr(rownames(mat),1,15)
mat = merge(mat,geneIDs1,by.x='row.names',by.y="GENEID")
length(mat$SYMBOL)
length(as.set(mat$SYMBOL))
unique_genename = make.unique(as.character(mat$SYMBOL), sep = "_")
length(unique_genename)
length(as.set(unique_genename))
rownames(mat) = unique_genename
mat = mat[,!(names(mat) %in% c('Row.names','SYMBOL'))]

write.csv(mat,file.path(savedir,"unscaled_matrix_for_integration.csv"))

# mat = read.csv(file.path(savedir,"RNA_unscaled_matrix_for_integration.csv"))

# mat = mat[,!names(mat) %in% c('hES_rep1','hES_rep2')] # remove ES samples
mat = t(scale(t(mat)))
# mat  <- mat - rowMeans(mat)

#anno <- as.data.frame(colData(rld_no_batch_var)[, c("sample","replicate")])


library(ComplexHeatmap)

# p = pheatmap(mat, scale='row')


mat = read.csv(file.path(savedir,"RNA_scaled_matrix_for_integration.csv"),row.names = 1, header=TRUE)

# col_scale = scale_fill_gradientn(colours = viridis::magma(20)) 



# pick the number of clusters
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 50
data <- mat
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
dev.off()
pdf(file.path(savedir,'elbow_plot_for_optimal_clustering.pdf'))
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()

pushViewport(viewport(gp = gpar(fontfamily = "Arial")))
set.seed(123)
dev.off()
#pdf(file.path(savedir,"hvgs_10_clusters.pdf"))
pdf(file.path(savedir,"hvgs_7_clusters_1.pdf"))

p = ComplexHeatmap::Heatmap(as.matrix(mat), name = "mat", col=viridis::magma(200), 
                            column_order = c('hES_rep1','hES_rep2','hDE_rep1','hDE_rep2','hAFE_rep1','hAFE_rep2','hPFE_rep1','hPFE_rep2'),
                            row_km = 7,cluster_columns = FALSE,show_row_names = FALSE,use_raster = TRUE, row_title = 1:7)
p <- draw(p)
r.dend <- row_dend(p)  #Extract row dendrogram
rcl.list <- row_order(p)  #Extract clusters (output is a list)

for (i in 1:length(row_order(p))){ # number of clusters
  if (i == 1) {
    clu <- t(t(row.names(mat[row_order(p)[[i]],]))) # this is the first cluster in order, not cluster number 1
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
    } else {
      clu <- t(t(row.names(mat[row_order(p)[[i]],])))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)
      }
  }
write.csv(out,file.path(savedir,"7clusters_1.csv"))
dev.off()


write.csv(mat,file.path(savedir,"scaled_matrix_for_integration.csv"))
# rerun from here 
# topVarGenes <- head(order(rownames(assay(rld_no_batch)), decreasing = TRUE), 4000)
# mat  <- assay(vsd_no_batch)[ topVarGenes, ]

# anno <- as.data.frame(colData(vsd_no_batch)[, c("sample","replicate")])
# p = pheatmap(mat, annotation_col = anno,scale='row',show_rownames=F)
# 
# # savedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/cellmatch_rework/processed_rnaseq/downstream_figures/"
# save_pheatmap_pdf(p, file.path(savedir, "top_var_genes_heatmap_scale_vsd_no_batch.pdf"),width=20, height=40)
# 
# dev.off()

# include gene and make results table, improve filtering for top genes

# 
# columns(org.Hs.eg.db)
# ens.str <- substr(rownames(mat), 1, 15)
# resSig$symbol <- mapIds(org.Hs.eg.db,
#                         keys=ens.str,
#                         column="SYMBOL",
#                         keytype="ENSEMBL",
#                         multiVals="first")
# resSig$entrez <- mapIds(org.Hs.eg.db,
#                         keys=ens.str,
#                         column="ENTREZID",
#                         keytype="ENSEMBL",
#                         multiVals="first")
# 
# resOrdered <- resSig[order(resSig$pvalue),]
# 
# resOrdered[!is.na(resOrdered$symbol),]
# 
# resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
# write.csv(resOrderedDF, file = "results.csv")


# Time course experiments 





# use vst counts for inferelator



mat_inferelator  <- assay(vsd_batch)
# mat_inferelator  <- mat_inferelator - rowMeans(mat_inferelator) # this line was added for the scaled version

# add in gene annotations
ens.str <- substr(rownames(mat_inferelator), 1, 15)

genenames = mapIds(org.Hs.eg.db,keys=ens.str,
       column="SYMBOL",
       keytype="ENSEMBL",
       multiVals="first")


write.matrix(mat_inferelator,file="/Users/LoboM/Dropbox/hESC_bulk/RNA/downstream_figures/mat_inferelator_ip_not_scaled.tsv",sep='\t')

write.csv(genenames,file="/Users/LoboM/Dropbox/hESC_bulk/RNA/downstream_figures/gene_names.txt")




# make hox gene plot 
# make "important" gene plot for RNA - this probably doesn't need to go in now, considering we'll label Fig3a with 
# gene names



# make volcano plots for pairiwse DE 


# save normalized matrix for pairwise



# custom plots 
# Hox genes
# misc markers
library(sets)

mat <- assay(rld_no_batch)
tmpcol = colData(rld_no_batch)
colnames(mat) = paste(tmpcol$sample,'_',tmpcol$replicate,sep='')

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(mat),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
tmp_rowname = substr(rownames(mat),1,15)
length(tmp_rowname)
length(as.set(tmp_rowname)) # duplicates exist, make unique
rownames(mat) = substr(rownames(mat),1,15)
mat = merge(mat,geneIDs1,by.x='row.names',by.y="GENEID")
length(mat$SYMBOL)
length(as.set(mat$SYMBOL))
unique_genename = make.unique(as.character(mat$SYMBOL), sep = "_")
length(unique_genename)
length(as.set(unique_genename))
rownames(mat) = unique_genename
mat = mat[,!(names(mat) %in% c('Row.names','SYMBOL'))]

scaled_mat = t(scale(t(mat)))

savedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA//downstream_figures/"

write.csv(scaled_mat,file.path(savedir,"scaled_RNA_matrix_all_genes.csv"))


library(tidyr)
# mat = mat %>% drop_na()
library(ComplexHeatmap)
genes_to_label = c('NANOG','POU5F1','SOX2','SOX17','FOXA2','PAX1','PAX9',#'HOXA1','HOXA2','HOXA3','HOXA4',# 'HOXC1','HOXC2','HOXC3',
                   #'HOXB1','HOXB2','HOXB3','HOXB4','HOXC4','HOXD1','HOXD2','HOXD3','HOXD4',
                   'LEFTY1','LEFTY2','CER1','VGLL2','TBX1','HIC1','PBX1','PBX2','MEIS1','MEIS2','MEIS3','SIX1','GBX2','OTX1','OTX2','TFAP4','TFAP2C',
                   'ISL1','ISL2','IRX1','IRX2','IRX3','EYA1','FGF10','HEY1','BMP4','BMP1','BMP2', 'BMP3','BMP5','BMP6','BMP7','HES1','ID1','FN1','EOMES','NODAL','FGF17','SMAD7','FOXE1',#'NKX2-1','NKX2-2',
                   'NKX2-3','NKX2-5','NKX2-6','RORA','RORB','RORC','GATA3','GATA4','GATA6','TP53','TP63','SHH','DMRT2','FOXP1','FOXP2','WNT5A','ROR2')
# mat = as.numeric(mat)
scaled_mat_tmp = data.frame(scaled_mat[genes_to_label,])  #%>% drop_na()

# average across replicates
scaled_mat_tmp$hESC =   rowMeans(scaled_mat_tmp[,c('hES_rep1','hES_rep2')], na.rm=TRUE)
scaled_mat_tmp$hDE =   rowMeans(scaled_mat_tmp[,c('hDE_rep1', 'hDE_rep2')], na.rm=TRUE)
scaled_mat_tmp$hAFE =   rowMeans(scaled_mat_tmp[,c('hAFE_rep1', 'hAFE_rep2')], na.rm=TRUE)
scaled_mat_tmp$hPFE =   rowMeans(scaled_mat_tmp[,c('hPFE_rep1', 'hPFE_rep2')], na.rm=TRUE)

scaled_mat_tmp <- subset( scaled_mat_tmp, select = c('hESC','hDE','hAFE','hPFE') )


dev.off()
basedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/downstream_figures"



# scale by rows 
pdf(file.path(basedir,"interesting_known_markers.pdf"),width=4, height=10)
p = ComplexHeatmap::Heatmap(as.matrix(scaled_mat_tmp), name = "mat", col=viridis::magma(200), 
                            column_order = c('hESC','hDE','hAFE','hPFE'),
                            cluster_rows=TRUE,cluster_columns = FALSE,show_row_names = TRUE,use_raster = TRUE)
p <- draw(p)
dev.off()

# Hox genes
library(magrittr)
scaled_mat_tmp = data.frame(scaled_mat[grepl("HOX",rownames(scaled_mat)),])
scaled_mat_tmp$hESC =   rowMeans(scaled_mat_tmp[,c('hES_rep1','hES_rep2')], na.rm=TRUE)
scaled_mat_tmp$hDE =   rowMeans(scaled_mat_tmp[,c('hDE_rep1', 'hDE_rep2')], na.rm=TRUE)
scaled_mat_tmp$hAFE =   rowMeans(scaled_mat_tmp[,c('hAFE_rep1', 'hAFE_rep2')], na.rm=TRUE)
scaled_mat_tmp$hPFE =   rowMeans(scaled_mat_tmp[,c('hPFE_rep1', 'hPFE_rep2')], na.rm=TRUE)
scaled_mat_tmp <- subset( scaled_mat_tmp, select = c('hESC','hDE','hAFE','hPFE') )


dev.off()

pdf(file.path(basedir,"hox_genes.pdf"),width=4, height=4)
p = ComplexHeatmap::Heatmap(as.matrix(scaled_mat_tmp), name = "mat", col=viridis::magma(200), 
                            column_order = c('hESC','hDE','hAFE','hPFE'),row_order=sort(rownames(scaled_mat_tmp)),
                            cluster_rows=FALSE,cluster_columns = FALSE,show_row_names = TRUE,use_raster = TRUE)
p <- draw(p)
dev.off()



# play with RNA normalization
tpm = assays(dds_batch)[["abundance"]]
# add sample names
tmpcol = colData(dds_batch)
colnames(tpm) = paste(tmpcol$sample,'_',tmpcol$replicate,sep='')
# add in gene names, instead of ensembl ids
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(tpm),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
tmp_rowname = substr(rownames(tpm),1,15)
length(tmp_rowname)
length(as.set(tmp_rowname)) # duplicates exist, make unique
rownames(tpm) = substr(rownames(tpm),1,15)
tpm = merge(tpm,geneIDs1,by.x='row.names',by.y="GENEID")
length(tpm$SYMBOL)
length(as.set(tpm$SYMBOL))
unique_genename = make.unique(as.character(tpm$SYMBOL), sep = "_")
length(unique_genename)
length(as.set(unique_genename))
rownames(tpm) = unique_genename
tpm = tpm[,!(names(tpm) %in% c('Row.names','SYMBOL'))]
savedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA//downstream_figures/"
write.csv(tpm,file.path(savedir,"TPM_counts_all.csv"))




tpm = assays(dds_batch)[["counts"]]
# add sample names
tmpcol = colData(dds_batch)
colnames(tpm) = paste(tmpcol$sample,'_',tmpcol$replicate,sep='')
# add in gene names, instead of ensembl ids
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(tpm),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
tmp_rowname = substr(rownames(tpm),1,15)
length(tmp_rowname)
length(as.set(tmp_rowname)) # duplicates exist, make unique
rownames(tpm) = substr(rownames(tpm),1,15)
tpm = merge(tpm,geneIDs1,by.x='row.names',by.y="GENEID")
length(tpm$SYMBOL)
length(as.set(tpm$SYMBOL))
unique_genename = make.unique(as.character(tpm$SYMBOL), sep = "_")
length(unique_genename)
library(sets)
length(as.set(unique_genename))
rownames(tpm) = unique_genename
tpm = tpm[,!(names(tpm) %in% c('Row.names','SYMBOL'))]
savedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA//downstream_figures/"
write.csv(tpm,file.path(savedir,"raw_counts_all.csv"))


##### RNA heatmaps  ###### 



savedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA//downstream_figures/"

scaled_mat = read.csv(file.path(savedir,"scaled_RNA_matrix_all_genes.csv"))


library(tidyr)
# mat = mat %>% drop_na()
library(ComplexHeatmap)
genes_to_label = c('NANOG','POU5F1','SOX2','SOX17','FOXA2','PAX1','PAX9','HOXA1','HOXA2','HOXA3','HOXA4',# 'HOXC1','HOXC2','HOXC3',
                   'HOXB1','HOXB2','HOXB3','HOXB4',#'HOXC4','HOXD1','HOXD2','HOXD3','HOXD4',
                   'LEFTY1','LEFTY2','CER1','VGLL2','TBX1','HIC1','PBX1','PBX2','MEIS1','MEIS2','MEIS3','SIX1','GBX2','OTX1','OTX2','TFAP4','TFAP2C',
                   'ISL1','ISL2','IRX1','IRX2','IRX3','EYA1','FGF10','HEY1','BMP4','BMP1','BMP2', 'BMP3','BMP5','BMP6','BMP7','HES1','ID1','FN1','EOMES','NODAL','FGF17','SMAD7','FOXE1',#'NKX2-1','NKX2-2',
                   'NKX2-3','NKX2-5','NKX2-6','RORA','RORB','RORC','GATA3','GATA4','GATA6','TP53','TP63','SHH','DMRT2','FOXP1','FOXP2','WNT5A','ROR2')
# mat = as.numeric(mat)
scaled_mat_tmp = data.frame(scaled_mat[scaled_mat$X %in% genes_to_label,])  #%>% drop_na()

# average across replicates
scaled_mat_tmp$hESC =   rowMeans(scaled_mat_tmp[,c('hES_rep1','hES_rep2')], na.rm=TRUE)
scaled_mat_tmp$hDE =   rowMeans(scaled_mat_tmp[,c('hDE_rep1', 'hDE_rep2')], na.rm=TRUE)
scaled_mat_tmp$hAFE =   rowMeans(scaled_mat_tmp[,c('hAFE_rep1', 'hAFE_rep2')], na.rm=TRUE)
scaled_mat_tmp$hPFE =   rowMeans(scaled_mat_tmp[,c('hPFE_rep1', 'hPFE_rep2')], na.rm=TRUE)

rownames(scaled_mat_tmp) = scaled_mat_tmp$X
scaled_mat_tmp <- subset( scaled_mat_tmp, select = c('hESC','hDE','hAFE','hPFE') )


dev.off()
basedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/downstream_figures"



# scale by rows 
basedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA//'

pdf(file.path(basedir,"interesting_known_markers.pdf"),width=4, height=10)
p = ComplexHeatmap::Heatmap(as.matrix(scaled_mat_tmp), name = "mat", col=viridis::magma(200), 
                            column_order = c('hESC','hDE','hAFE','hPFE'),
                            cluster_rows=TRUE,cluster_columns = FALSE,show_row_names = TRUE,use_raster = TRUE)
p <- draw(p)
dev.off()

# Hox genes
library(magrittr)
scaled_mat_tmp = data.frame(scaled_mat[grepl("HOX",rownames(scaled_mat)),])
scaled_mat_tmp$hESC =   rowMeans(scaled_mat_tmp[,c('hES_rep1','hES_rep2')], na.rm=TRUE)
scaled_mat_tmp$hDE =   rowMeans(scaled_mat_tmp[,c('hDE_rep1', 'hDE_rep2')], na.rm=TRUE)
scaled_mat_tmp$hAFE =   rowMeans(scaled_mat_tmp[,c('hAFE_rep1', 'hAFE_rep2')], na.rm=TRUE)
scaled_mat_tmp$hPFE =   rowMeans(scaled_mat_tmp[,c('hPFE_rep1', 'hPFE_rep2')], na.rm=TRUE)
scaled_mat_tmp <- subset( scaled_mat_tmp, select = c('hESC','hDE','hAFE','hPFE') )


dev.off()

pdf(file.path(basedir,"hox_genes.pdf"),width=4, height=4)
p = ComplexHeatmap::Heatmap(as.matrix(scaled_mat_tmp), name = "mat", col=viridis::magma(200), 
                            column_order = c('hESC','hDE','hAFE','hPFE'),row_order=sort(rownames(scaled_mat_tmp)),
                            cluster_rows=FALSE,cluster_columns = FALSE,show_row_names = TRUE,use_raster = TRUE)
p <- draw(p)
dev.off()




# top genes plus their direct targets heatamp - candidate fig 5
##### volcano plots ##### 

res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0.0,contrast = c("sample","hDE","hES")) # comparison between TEP and AFE
# resSig_all <- subset(res_all, padj < 0.1)
res_all <- lfcShrink(dds_all,
                 contrast = c("sample","hDE","hES"), res=res_all, type = 'normal')
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(res_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
geneIDs1 = geneIDs1[!is.na(geneIDs1$SYMBOL),]
symbols = geneIDs1[match(substr(rownames(res_all),1,15),geneIDs1$GENEID),"SYMBOL"]
symbols[is.na(symbols)] = 'tmp_na'
rownames(res_all) = symbols
keep = rownames(res_all)!='tmp_na'
res_all <- res_all[keep,]
#res_all$geneid = substr(rownames(res_all),1,15)
#res_all@listData = merge(res_all@listData,geneIDs1,by.x="geneid",by.y="GENEID")
# res_all
# res_all@listData
tmp_thresh = max(abs(min(res_all$log2FoldChange )),max(res_all$log2FoldChange ))
p = EnhancedVolcano(res_all,
                lab = rownames(res_all),
                x = 'log2FoldChange',
                y = 'pvalue',
                title='Differential RNA hDE over hESC', pCutoff = 10e-2,
                FCcutoff = 1,pointSize = 3.0,
                labSize = 6.0,colAlpha = 1,boxedLabels = TRUE,drawConnectors = TRUE,
                widthConnectors = 0.75,labFace = 'bold',xlim = c(-tmp_thresh, tmp_thresh),
                selectLab = c('LEFTY1','LEFTY2','EOMES',
                              'CER1','GATA4','SOX17','GATA6','SOX2','POU5F1', 'BMP2','NODAL',
                              'NANOG','GATA3','GATA4','GATA6','OTX2','CER1','FGF17','SOX17'),)
ggsave(file.path(savedir, "DE_over_ESC_volcano.pdf"),plot=p,width=10, height=10)

dev.off()
# write.csv(res_all@listData,"/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/RNA/differential_csvs/DE_over_ESC_differential.csv")


res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0,contrast = c("sample","hAFE","hES")) # comparison between TEP and AFE
res_all <- lfcShrink(dds_all,
                     contrast = c("sample","hAFE","hES"), res=res_all, type = 'normal')
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(res_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
geneIDs1 = geneIDs1[!is.na(geneIDs1$SYMBOL),]
symbols = geneIDs1[match(substr(rownames(res_all),1,15),geneIDs1$GENEID),"SYMBOL"]
symbols[is.na(symbols)] = 'tmp_na'
rownames(res_all) = symbols
keep = rownames(res_all)!='tmp_na'
res_all <- res_all[keep,]
tmp_thresh = max(abs(min(res_all$log2FoldChange )),max(res_all$log2FoldChange ))

p = EnhancedVolcano(res_all,
                    lab = rownames(res_all),
                    x = 'log2FoldChange',
                    y = 'pvalue',pointSize = 3.0,
                    labSize = 6.0,colAlpha = 1,boxedLabels = TRUE,
                    title='Differential RNA hAFE over hES', pCutoff = 10e-2,
                    FCcutoff = 1,drawConnectors = TRUE,
                    widthConnectors = 0.75,labFace = 'bold',xlim = c(-tmp_thresh, tmp_thresh),
                selectLab = c('GATA4','GATA6',
                              'MEIS1','MEIS2','MEIS3','SOX2','HOXB3',
                              'PBX1','HES1','PAX9','PAX1','GATA3','SIX1','VGLL2',
                              'EYA1','ISL1','HOXB2','HOXB4','HOXA2','SOX2','FOXA2',
                              'POU5F1', 'NANOG',
                              'OTX2','RARB'),max.overlaps=Inf)
                    
ggsave(file.path(savedir, "AFE_over_ES_volcano.pdf"),plot=p,width=10, height=10)

# dev.off()







res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0,contrast = c("sample","hAFE","hDE")) # comparison between TEP and AFE
res_all <- lfcShrink(dds_all,
                     contrast = c("sample","hAFE","hDE"), res=res_all, type = 'normal')
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(res_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
geneIDs1 = geneIDs1[!is.na(geneIDs1$SYMBOL),]
symbols = geneIDs1[match(substr(rownames(res_all),1,15),geneIDs1$GENEID),"SYMBOL"]
symbols[is.na(symbols)] = 'tmp_na'
rownames(res_all) = symbols
keep = rownames(res_all)!='tmp_na'
res_all <- res_all[keep,]
tmp_thresh = max(abs(min(res_all$log2FoldChange )),max(res_all$log2FoldChange ))

p = EnhancedVolcano(res_all,
                    lab = rownames(res_all),
                    x = 'log2FoldChange',
                    y = 'pvalue',pointSize = 3.0,
                    labSize = 6.0,colAlpha = 1,boxedLabels = TRUE,
                    title='Differential RNA hAFE over hDE', pCutoff = 10e-2,
                    FCcutoff = 1,drawConnectors = TRUE,
                    widthConnectors = 0.75,labFace = 'bold',xlim = c(-tmp_thresh, tmp_thresh),
                    selectLab = c('GATA4','GATA6',
                                  'MEIS1','MEIS2','MEIS3','SOX2','HOXB3',
                                  'PBX1','HES1','PAX9','PAX1','GATA3','SIX1','VGLL2',
                                  'EYA1','ISL1','HOXB2','HOXB4','HOXA2','SOX2','FOXA2',
                                  'SOX17','CER1','NODAL',
                                  'OTX2','RARB'),max.overlaps=Inf)

ggsave(file.path(savedir, "AFE_over_DE_volcano.pdf"),plot=p,width=10, height=10)








res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0,contrast = c("sample","hPFE","hAFE")) # comparison between TEP and AFE
res_all <- lfcShrink(dds_all,
                     contrast = c("sample","hPFE","hAFE"), res=res_all, type = 'normal')
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(res_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
geneIDs1 = geneIDs1[!is.na(geneIDs1$SYMBOL),]
symbols = geneIDs1[match(substr(rownames(res_all),1,15),geneIDs1$GENEID),"SYMBOL"]
symbols[is.na(symbols)] = 'tmp_na'
rownames(res_all) = symbols
keep = rownames(res_all)!='tmp_na'
res_all <- res_all[keep,]
tmp_thresh = max(abs(min(res_all$log2FoldChange )),max(res_all$log2FoldChange ))

p = EnhancedVolcano(res_all,
                    lab = rownames(res_all),
                    x = 'log2FoldChange',
                    y = 'pvalue', pointSize = 3.0,
                    labSize = 6.0,colAlpha = 1,boxedLabels = TRUE,
                    title='Differential RNA hPFE over hAFE', pCutoff = 10e-2,
                    FCcutoff = 1, drawConnectors = TRUE,
                    widthConnectors = 0.75,labFace = 'bold',xlim = c(-tmp_thresh, tmp_thresh),
                    max.overlaps=Inf,
                selectLab = c('FN1','HOXA4','BMP5','DMRT2','HOXB4','HOXB3','TP63','BMP1','ISL1','HOXA2','HOXB2','IRX1',
                'ISL2','SMAD7',
                'PBX1','HES1','PAX9','PAX1','GATA3','SIX1','VGLL2',
                'EYA1','ISL1','HOXB4','HOXA2','SOX2','FOXA2','HEY1','SHH','FOXE1','PBX1','VGLL2')
                    )

ggsave(file.path(savedir, "PFE_over_AFE_volcano.pdf"),plot=p,width=10, height=10)













res_all <- results(dds_all,alpha=0.1,lfcThreshold = 0,contrast = c("sample","hPFE","hDE")) # comparison between TEP and AFE
res_all <- lfcShrink(dds_all,
                     contrast = c("sample","hPFE","hDE"), res=res_all, type = 'normal')
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= substr(rownames(res_all),1,15), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
geneIDs1 = geneIDs1[!is.na(geneIDs1$SYMBOL),]
symbols = geneIDs1[match(substr(rownames(res_all),1,15),geneIDs1$GENEID),"SYMBOL"]
symbols[is.na(symbols)] = 'tmp_na'
rownames(res_all) = symbols
keep = rownames(res_all)!='tmp_na'
res_all <- res_all[keep,]
tmp_thresh = max(abs(min(res_all$log2FoldChange )),max(res_all$log2FoldChange ))

p = EnhancedVolcano(res_all,
                    lab = rownames(res_all),
                    x = 'log2FoldChange',
                    y = 'pvalue', pointSize = 3.0,
                    labSize = 6.0,colAlpha = 1,boxedLabels = TRUE,
                    title='Differential RNA hPFE over hDE', pCutoff = 10e-2,
                    FCcutoff = 1, drawConnectors = TRUE,
                selectLab = c('FN1','BMP5','HOXB3','TP63','BMP1','ISL1','HOXA2','HOXB2',
                              'MEIS1','MEIS2','CDH1','EOMES','NODAL','LEFTY1','LEFTY2','FGF17',
                              'PBX1','HES1','PAX9','GATA3','SIX1','VGLL2','GATA4','GATA6',
                              'ISL1','HOXA2','HEY1','SHH','FOXE1','PBX1','VGLL2'),
                    widthConnectors = 0.75,labFace = 'bold',xlim = c(-tmp_thresh, tmp_thresh),
                     max.overlaps=Inf
                    
)

ggsave(file.path(savedir, "PFE_over_DE_volcano.pdf"),plot=p,width=10, height=10)


