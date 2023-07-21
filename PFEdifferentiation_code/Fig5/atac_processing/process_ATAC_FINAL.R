library(DESeq2)
library(stringr)
library(lattice) 

#library(devtools)
#install_github("andrelmartins/bigWig", subdir="bigWig", force=TRUE)

library(bigWig)


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("DiffBind")

library(DiffBind)
# library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(dplyr)


samples <- read.csv('/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/samplesheet_FINAL_2.csv')
peaks <- dba(sampleSheet=samples,minOverlap=1)
peaks = dba.blacklist(peaks, blacklist=TRUE, greylist=FALSE) # blacklist already done, step not required


peaks.consensus <- dba.peakset(peaks,consensus = DBA_CONDITION, minOverlap=2)



dba.show(peaks.consensus,peaks.consensus$masks$Consensus)

dev.off()
savedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/'



pdf(file.path(savedir,'venn_peak_overlap_2.pdf'))
dba.plotVenn(peaks.consensus,peaks.consensus$masks$Consensus)
dev.off()

consensus <- dba.peakset(peaks.consensus, bRetrieve=TRUE,
                         peaks=peaks.consensus$masks$Consensus,
                         minOverlap=1) 



binding_matrix <- dba.count(peaks,peaks=consensus,filter=1) # peaks.consensus

# save(binding_matrix,file="/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/peak_count_matrix_post_qc_2.RData")

load("/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/peak_count_matrix_post_qc_2.RData")

savedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/'

# dev.off()
pdf(file.path(savedir,'correlation_heatmap_atac_post_qc_2.pdf'))
plot(binding_matrix,colScheme="Purples")
dev.off()

dev.off()
pdf(file.path(savedir,'pca_plot_atac_post_qc_2.pdf'))

dba.plotPCA(binding_matrix,DBA_CONDITION, label=DBA_ID,vColors = c("red", "yellow", "skyblue", "orchid"))
dev.off()


# normalize
binding_matrix_default <- dba.normalize(binding_matrix)

# scores are the normalized scores
save(binding_matrix_default,file="/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/peak_count_matrix_post_qc_normalized_2.RData")

# load("/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/peak_count_matrix_post_qc_normalized_2.RData")


##### DIFFERENTIAL ANALYSIS - NOT USED IN MANUSCRIPT ####### 
# stringent 0.01, 1
# lax 0.1, 0.5 
fc_thresh = 1
pval_thresh = 0.05

de_over_esc_model <- dba.contrast(binding_matrix_default,design="~Replicate + Condition",contrast=c("Condition","hDE","hESC"))
de_over_esc_model <- dba.analyze(de_over_esc_model)
de_over_esc_gain <- dba.report(de_over_esc_model,bAll=FALSE,bUsePval=FALSE,bGain=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
de_over_esc_loss <- dba.report(de_over_esc_model,bAll=FALSE,bUsePval=FALSE,bLoss=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
write.table( x = data.frame(de_over_esc_gain[[1]]), file = file.path(savedir,'/pairiwse_differential/DE_over_ESC/de_over_esc_gain_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(de_over_esc_loss[[1]]), file = file.path(savedir,'/pairiwse_differential/DE_over_ESC/de_over_esc_loss_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )

#DBA$config$bUsePval

afe_over_esc_model <- dba.contrast(binding_matrix_default,design="~Replicate + Condition",contrast=c("Condition","hAFE","hESC"))
afe_over_esc_model <- dba.analyze(afe_over_esc_model)
afe_over_esc_gain <- dba.report(afe_over_esc_model,bAll=FALSE,bUsePval=FALSE,bGain=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
afe_over_esc_loss <- dba.report(afe_over_esc_model,bAll=FALSE,bUsePval=FALSE,bLoss=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
write.table( x = data.frame(afe_over_esc_gain[[1]]), file = file.path(savedir,'/pairiwse_differential/AFE_over_ESC/afe_over_esc_gain_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(afe_over_esc_loss[[1]]), file = file.path(savedir,'/pairiwse_differential/AFE_over_ESC/afe_over_esc_loss_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )


afe_over_de_model <- dba.contrast(binding_matrix_default,design="~Replicate + Condition",contrast=c("Condition","hAFE","hDE"))
afe_over_de_model <- dba.analyze(afe_over_de_model)
afe_over_de_gain <- dba.report(afe_over_de_model,bAll=FALSE,bUsePval=FALSE,bGain=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
afe_over_de_loss <- dba.report(afe_over_de_model,bAll=FALSE,bUsePval=FALSE,bLoss=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
write.table( x = data.frame(afe_over_de_gain[[1]]), file = file.path(savedir,'/pairiwse_differential/AFE_over_DE/afe_over_de_gain_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(afe_over_de_loss[[1]]), file = file.path(savedir,'/pairiwse_differential/AFE_over_DE/afe_over_de_loss_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )



# for PFE, no replicates, so results not trustworthy, 
# couldnt build a model since no replicates, so specified the contrast explicityly
pfe_over_afe_model <- dba.contrast(binding_matrix_default,
                                   group1=binding_matrix_default$masks$hPFE, 
                                   group2=binding_matrix_default$masks$hAFE,
                                   name1="hPFE", name2="hAFE")
pfe_over_afe_model <- dba.analyze(pfe_over_afe_model)
pfe_over_afe_gain <- dba.report(pfe_over_afe_model,bAll=FALSE,bUsePval=FALSE,bGain=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
pfe_over_afe_loss <- dba.report(pfe_over_afe_model,bAll=FALSE,bUsePval=FALSE,bLoss=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
write.table( x = data.frame(pfe_over_afe_gain[[1]]), file = file.path(savedir,'/pairiwse_differential/PFE_over_AFE/pfe_over_afe_gain_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(pfe_over_afe_loss[[1]]), file = file.path(savedir,'/pairiwse_differential/PFE_over_AFE/pfe_over_afe_loss_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )



pfe_over_de_model <- dba.contrast(binding_matrix_default,
                                  group1=binding_matrix_default$masks$hPFE, 
                                  group2=binding_matrix_default$masks$hDE,
                                  name1="hPFE", name2="hDE")
pfe_over_de_model <- dba.analyze(pfe_over_de_model)
pfe_over_de_gain <- dba.report(pfe_over_de_model,bAll=FALSE,bUsePval=FALSE,bGain=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
pfe_over_de_loss <- dba.report(pfe_over_de_model,bAll=FALSE,bUsePval=FALSE,bLoss=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
write.table( x = data.frame(pfe_over_de_gain[[1]]), file = file.path(savedir,'/pairiwse_differential/PFE_over_DE/pfe_over_de_gain_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(pfe_over_de_loss[[1]]), file = file.path(savedir,'/pairiwse_differential/PFE_over_DE/pfe_over_de_loss_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )


pfe_over_esc_model <- dba.contrast(binding_matrix_default,
                                   group1=binding_matrix_default$masks$hPFE, 
                                   group2=binding_matrix_default$masks$hESC,
                                   name1="hPFE", name2="hESC")
pfe_over_esc_model <- dba.analyze(pfe_over_esc_model)
pfe_over_esc_gain <- dba.report(pfe_over_esc_model,bAll=FALSE,bUsePval=FALSE,bGain=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
pfe_over_esc_loss <- dba.report(pfe_over_esc_model,bAll=FALSE,bUsePval=FALSE,bLoss=TRUE,th=pval_thresh,fold=fc_thresh,DataType =  DBA_DATA_GRANGES)
write.table( x = data.frame(pfe_over_esc_gain[[1]]), file = file.path(savedir,'/pairiwse_differential/PFE_over_ESC/pfe_over_esc_gain_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
write.table( x = data.frame(pfe_over_esc_loss[[1]]), file = file.path(savedir,'/pairiwse_differential/PFE_over_ESC/pfe_over_esc_loss_lax_2.bed'), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )





###### DOWNSTREAM ###### 
all_peaks = dba.peakset(binding_matrix_default,peaks=NULL,bRetrieve = TRUE)

export(all_peaks, file.path(savedir,'consensus_peaks_400bp.bed')) # save as bed


##### annotate peaks - Homer ##### 
# source ~/.bashrc
# convert to 0-based indexing before running Homer
# annotatePeaks.pl consensus_peaks_400bp.bed hg38 -annStats all_DARs_annostats.txt > consensus_peaks_400bp_annotated.bed
# Homer output is 1-based indexing






# make heatmap with the high cv consensus set 
savedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/"
con_heatmap = read.csv(file.path(savedir,"/FINAL/high_cv_consensus_40k_median_plus_1sd.csv"),header=TRUE,sep='\t')
rownames(con_heatmap) = con_heatmap$X
con_heatmap = con_heatmap[,!(names(con_heatmap) %in% c('cv','peak_id','chr','start','end','X'))]
con_heatmap = t(scale(t(con_heatmap)))

library(ComplexHeatmap)


set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 50
data <- con_heatmap
wss <- sapply(1:k.max,
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
dev.off()
pdf(file.path(savedir,'elbow_plot_for_optimal_cluster_num.pdf'))
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()




pushViewport(viewport(gp = gpar(fontfamily = "Arial",fontsize = 8)))
set.seed(123)
dev.off()

pdf(file.path(savedir,"consensus_4_clusters_with_cv_40k.pdf"))

p = ComplexHeatmap::Heatmap(as.matrix(con_heatmap), name = "mat", col=viridis::magma(200), 
                            column_order = c("hESC","hDE","hAFE","hPFE"),
                            row_km = 4,cluster_columns = FALSE,show_row_names = FALSE,use_raster = TRUE, row_title = 1:4)
p <- draw(p)
r.dend <- row_dend(p)  #Extract row dendrogram
rcl.list <- row_order(p)  #Extract clusters (output is a list)

for (i in 1:length(row_order(p))){ # number of clusters
  if (i == 1) {
    clu <- t(t(row.names(con_heatmap[row_order(p)[[i]],]))) # this is the first cluster in order, not cluster number 1
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("PeakID", "Cluster")
  } else {
    clu <- t(t(row.names(con_heatmap[row_order(p)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
write.csv(out,file.path(savedir,"4_clusters_ATAC_consensus_with_cv_cutoff_40k.csv"))
dev.off()






# extract quantile normalized TMMs on 40k set
# subset in python

#######  TMM normalization   ####### 
binding_matrix_tmm <- dba.normalize(binding_matrix,normalize=DBA_NORM_TMM)
matrix_counts_tmm <- dba.peakset(binding_matrix_tmm, bRetrieve=T, DataType=DBA_DATA_FRAME)
write.table(matrix_counts_tmm, file = file.path(savedir,"/FINAL/consensus_count_matrix_tmm.csv"), sep = "\t",
            row.names = TRUE, col.names = TRUE, quote=FALSE)

#######  TMM + quantile normalization   ####### 
library(preprocessCore)
tmp1 = normalize.quantiles(as.matrix(matrix_counts_tmm[4:10]))
rownames(tmp1) = paste(paste(matrix_counts_tmm$CHR,matrix_counts_tmm$START,sep=':'),matrix_counts_tmm$END,sep='-')
colnames(tmp1) = colnames(matrix_counts_tmm)[4:10]
write.table(tmp1,file = file.path(savedir,"/FINAL/consensus_count_matrix_tmm_quan_norm.csv"), sep = "\t",
            row.names = TRUE, col.names = TRUE, quote=FALSE)




##### MAKE HEATMAP - SUPP 5D ####### 
# make heatmap with the high cv consensus set 
savedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/"
# file made in python script get_varying_peaks_40k.py
con_heatmap = read.csv(file.path(savedir,"/FINAL/high_cv_consensus_40k_median_plus_1sd.csv"),header=TRUE,sep='\t')
rownames(con_heatmap) = con_heatmap$X
con_heatmap = con_heatmap[,!(names(con_heatmap) %in% c('cv','peak_id','chr','start','end','X'))]
con_heatmap = t(scale(t(con_heatmap)))

library(ComplexHeatmap)


set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 50
data <- con_heatmap
wss <- sapply(1:k.max,
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
wss
dev.off()
pdf(file.path(savedir,'elbow_plot_for_optimal_cluster_num.pdf'))
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE,
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()

# optimal number of clusters at elbow was 4


pushViewport(viewport(gp = gpar(fontfamily = "Arial",fontsize = 8)))
set.seed(123)
dev.off()

pdf(file.path(savedir,"consensus_4_clusters_with_cv_40k.pdf"))

p = ComplexHeatmap::Heatmap(as.matrix(con_heatmap), name = "mat", col=viridis::magma(200), 
                            column_order = c("hESC","hDE","hAFE","hPFE"),
                            row_km = 4,cluster_columns = FALSE,show_row_names = FALSE,use_raster = TRUE, row_title = 1:4)
p <- draw(p)
r.dend <- row_dend(p)  #Extract row dendrogram
rcl.list <- row_order(p)  #Extract clusters (output is a list)

for (i in 1:length(row_order(p))){ # number of clusters
  if (i == 1) {
    clu <- t(t(row.names(con_heatmap[row_order(p)[[i]],]))) # this is the first cluster in order, not cluster number 1
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("PeakID", "Cluster")
  } else {
    clu <- t(t(row.names(con_heatmap[row_order(p)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
write.csv(out,file.path(savedir,"4_clusters_ATAC_consensus_with_cv_cutoff_40k.csv"))
dev.off()




