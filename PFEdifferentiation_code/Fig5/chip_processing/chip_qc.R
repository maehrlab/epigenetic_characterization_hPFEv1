# all minCount = 1 should be minCount = 0


library(DESeq2)
library(stringr)
library(lattice) 

#library(devtools)
#install_github("andrelmartins/bigWig", subdir="bigWig", force=TRUE)

library(bigWig)


# from https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("DiffBind")

library(DiffBind)
# library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(dplyr)


# from https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("DiffBind")

library(DiffBind)
# library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(dplyr)


samples <- read.csv('/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/all_chip.csv')
peaks <- dba(sampleSheet=samples,minOverlap=1)
peaks = dba.blacklist(peaks, blacklist=TRUE, greylist=FALSE) # blacklist already done, step not required


peaks.consensus <- dba.peakset(peaks,consensus = DBA_CONDITION, minOverlap=2)



dba.show(peaks.consensus,peaks.consensus$masks$Consensus)


# dev.off()
savedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ChIP/downstream/'




consensus <- dba.peakset(peaks.consensus, bRetrieve=TRUE,
                         peaks=peaks.consensus$masks$Consensus,
                         minOverlap=1) 



binding_matrix <- dba.count(peaks,peaks=consensus,filter=1)

save(binding_matrix,file=file.path(savedir,"/all_histone_peak_count_matrix_post_qc_2.RData"))

# dev.off()
pdf(file.path(savedir,'correlation_heatmap_all_histone_post_qc_2.pdf'))
plot(binding_matrix,colScheme="YlGnBu")
dev.off()


pdf(file.path(savedir,'pca_plot_allhistone_post_qc_2.pdf'))
#dev.off()
dba.plotPCA(binding_matrix,DBA_CONDITION, label=DBA_ID)# ,vColors = c("red", "yellow", "skyblue", "orchid"))
dev.off()