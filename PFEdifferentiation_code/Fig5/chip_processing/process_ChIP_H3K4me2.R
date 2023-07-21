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



samples <- read.csv('/Users/LoboM/Dropbox/chip_bams/H3K4me2/samplesheet_2.csv')
peaks <- dba(sampleSheet=samples,minOverlap=1)
peaks = dba.blacklist(peaks, blacklist=TRUE, greylist=FALSE) # blacklist already done, step not required


peaks.consensus <- dba.peakset(peaks,consensus = DBA_CONDITION, minOverlap=2)



dba.show(peaks.consensus,peaks.consensus$masks$Consensus)

# dev.off()
savedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ChIP/downstream_H3K4me2/'



pdf(file.path(savedir,'venn_peak_overlap_2.pdf'))
dba.plotVenn(peaks.consensus,peaks.consensus$masks$Consensus)
dev.off()

consensus <- dba.peakset(peaks.consensus, bRetrieve=TRUE,
                         peaks=peaks.consensus$masks$Consensus,
                         minOverlap=1) 



binding_matrix <- dba.count(peaks,peaks=consensus,filter=1)

save(binding_matrix,file=file.path(savedir,"/peak_count_matrix_post_qc_2.RData"))



# dev.off()
pdf(file.path(savedir,'correlation_heatmap_atac_post_qc_2.pdf'))
plot(binding_matrix)
dev.off()


pdf(file.path(savedir,'pca_plot_atac_post_qc_2.pdf'))
#dev.off()
dba.plotPCA(binding_matrix,DBA_CONDITION, label=DBA_ID)
dev.off()



# save TMM quan norm counts in H3K2ac regions 
con_dir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/per_sample_integration_grn_etc/nov17_2022/input/"
atac_con = import(file.path(con_dir,"final_peakset.bed"))

binding_matrix_tmm_consensus_2 <- dba.count(peaks,peaks=atac_con,summits=0,score=DBA_NORM_TMM) 
matrix_counts_tmm <- dba.peakset(binding_matrix_tmm_consensus_2, bRetrieve=T, DataType=DBA_DATA_FRAME)
library(preprocessCore)
tmp1 = normalize.quantiles(as.matrix(matrix_counts_tmm[4:9]))
rownames(tmp1) = paste(paste(matrix_counts_tmm$CHR,matrix_counts_tmm$START,sep=':'),matrix_counts_tmm$END,sep='-')
colnames(tmp1) = colnames(matrix_counts_tmm)[4:9]
write.table(tmp1, file = file.path(savedir,"consensus_atac_H3K4me2_counts_tmm_quan_norm.csv"), sep = "\t",
            row.names = TRUE, col.names = TRUE, quote=FALSE)
write.table(matrix_counts_tmm, file = file.path(savedir,"consensus_atac_H3K4me2_counts_tmm.csv"), sep = "\t",
            row.names = TRUE, col.names = TRUE, quote=FALSE)



