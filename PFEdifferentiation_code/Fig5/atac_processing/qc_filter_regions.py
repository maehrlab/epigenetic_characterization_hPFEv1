#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# - go with qval = 3, signal = 3 

# - run through a list of peak files, set qval = 3, signal = 3 cutoff

# In[2]:


filenames = """hAFE_ATAC_hAFE_ATAC_rep1_2015NOV14.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak,hAFE_ATAC_hAFE_ATAC_rep2_2015NOV14.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak,hAFE_ATAC_rep.pooled.pval0.01.300K.bfilt.narrowPeak,hDE_ATAC_hDE_ATAC_rep1_2015NOV14.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak,hDE_ATAC_hDE_ATAC_rep2_2015NOV14.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak,hDE_ATAC_hDE_rep3_2017AUG29_ATAC.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak,hDE_ATAC_hDE_rep4_2017AUG29_ATAC.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak,hDE_ATAC_rep.pooled.pval0.01.300K.bfilt.narrowPeak,hES_ATAC_hESC_ATAC_rep1_2015NOV14.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak,hES_ATAC_hESC_ATAC_rep2_2015NOV14.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak,hES_ATAC_hES_rep3_2017AUG29_ATAC.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak,hES_ATAC_rep.pooled.pval0.01.300K.bfilt.narrowPeak,hPFE_ATAC_hPFE_ATAC_rep1_2016JAN07.1.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.bfilt.narrowPeak""".split(',')


# In[4]:


import pandas as pd


# In[5]:


sig_cutoff = 3
q_cutoff = 3
for entry in filenames: 
    tmp_file = pd.read_csv("../peak_pipe/" + entry ,sep='\t',header=None)
    tmp_file = tmp_file[((tmp_file[6]>sig_cutoff) & (tmp_file[8]>q_cutoff))]
    tmp_file.to_csv("../peaks_clean/" + entry,sep='\t',index=False,header=False)

