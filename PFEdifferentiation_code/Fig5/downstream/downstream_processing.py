#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[29]:


import seaborn as sns


# In[70]:


basedir = "/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/FINAL/"


# In[ ]:


# TMM + quantile file made in R code 
consensus_counts = pd.read_csv("/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/FINAL/consensus_count_matrix_tmm_quan_norm.csv",sep='\t')



consensus_counts['hESC'] = 0.5*(consensus_counts['hESC_rep1'] + consensus_counts['hESC_rep3'])
consensus_counts['hDE'] = 0.5*(consensus_counts['hDE_rep1'] + consensus_counts['hDE_rep4'])
consensus_counts['hAFE'] = 0.5*(consensus_counts['hAFE_rep1'] + consensus_counts['hAFE_rep2'])
consensus_counts['hPFE'] = consensus_counts['hPFE_rep1']

new_idx = []
for entry in consensus_counts.index:
    new_idx.append(entry.split(':')[0] + ':' + str(int(entry.split(':')[1].split('-')[0]) -1 ) + '-' + (entry.split(':')[1].split('-')[1] ))
consensus_counts['new_idx'] = new_idx
consensus_counts = consensus_counts.set_index('new_idx')
consensus_counts.index.name = None

for entry in ['hESC_rep1', 'hESC_rep3', 'hDE_rep1', 'hDE_rep4', 'hAFE_rep1',
       'hAFE_rep2', 'hPFE_rep1']:
    del consensus_counts[entry]
    
consensus_counts.to_csv("/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/differential_peaks/FINAL/consensus_count_matrix_tmm_quan_norm_for_maelstrom.bed",sep='\t')


# In[71]:


# read consensus, subset with low cv 
consensus_peaks = pd.read_csv(basedir + "consensus_count_matrix_tmm_quan_norm_for_maelstrom.bed",sep='\t',index_col=0)


# In[72]:


consensus_peaks


# In[73]:


consensus_peaks['cv'] = (consensus_peaks[['hESC','hDE','hAFE','hPFE']].std(axis=1))/(consensus_peaks[['hESC','hDE','hAFE','hPFE']].mean(axis=1))


# In[74]:


sns.distplot(consensus_peaks['cv'])


# In[75]:


print(consensus_peaks['cv'].median(),consensus_peaks['cv'].mean(),consensus_peaks['cv'].std())


# In[76]:


cv_thresh = consensus_peaks['cv'].median() + consensus_peaks['cv'].std()
print(cv_thresh)
consensus_peaks = consensus_peaks[consensus_peaks['cv'] > cv_thresh]


# In[77]:


consensus_peaks.shape


# In[78]:


consensus_peaks.to_csv(basedir + "high_cv_consensus_40k_median_plus_1sd.csv")


# In[79]:


consensus_peaks['peak_id'] = consensus_peaks.index

chr_tmp = []
start_tmp = []
end_tmp = []

for entry in consensus_peaks.index:
    chr_tmp.append(entry.split(':')[0])
    start_tmp.append(int(entry.split(':')[1].split('-')[0]) + 1)
    end_tmp.append(int(entry.split(':')[1].split('-')[1]))
    
consensus_peaks['chr'] = chr_tmp
consensus_peaks['start'] = start_tmp
consensus_peaks['end'] = end_tmp
    


# In[80]:


consensus_peaks.to_csv(basedir + "high_cv_consensus_40k_median_plus_1sd.csv",sep='\t')


# #### Histone mark avg TMM quan-norm files

# #### H3K27ac

# In[ ]:


basedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/per_sample_integration_grn_etc/nov17_2022/input/'

# file made in process_ChIP_H3K27ac.R 

consensus_all = pd.read_csv(basedir + "consensus_atac_H3K27ac_counts_tmm_quan_norm.csv",index_col=0,sep='\t')

new_idx = []
for entry in consensus_all.index:
    new_idx.append(entry.split(':')[0] + ':' + str(int(entry.split(':')[1].split('-')[0]) - 1) + '-' + entry.split(':')[1].split('-')[1])
    
consensus_all['new_idx'] = new_idx
consensus_all = consensus_all.set_index('new_idx')
consensus_all.index.name = None

# average replicates
consensus_all['hESC'] = 0.5*(consensus_all['hESC_rep1'] + consensus_all['hESC_rep2'])
consensus_all['hDE'] = 0.5*(consensus_all['hDE_rep1'] + consensus_all['hDE_rep2'])
consensus_all['hAFE'] = 0.5*(consensus_all['hAFE_rep1'] + consensus_all['hAFE_rep2'])
consensus_all['hPFE'] = 0.5*(consensus_all['hPFE_rep1'] + consensus_all['hPFE_rep2'])


for entry in ['hESC_rep1', 'hESC_rep2', 'hDE_rep1', 'hDE_rep2',
        'hAFE_rep1', 'hAFE_rep2', 'hPFE_rep1', 'hPFE_rep2']:
     del consensus_all[entry]
    
consensus_40k = pd.read_csv(basedir + "consensus_40k.bed",header=None,sep='\t')
consensus_40k['peak_id'] = consensus_40k[0] + ':' + consensus_40k[1].astype(str) + '-' + consensus_40k[2].astype(str)
consensus_40k = consensus_40k.set_index('peak_id')
consensus_40k.index.name = None


consensus_40k_h3k27ac = consensus_all.loc[list(set(consensus_40k.index).intersection(consensus_all.index)),]


consensus_40k_h3k27ac.to_csv(basedir + "consensus_40k_for_maelstrom_ATAC_H3K27ac_with_avg.bed",sep='\t')


# #### H3K4me2

# In[ ]:


basedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/per_sample_integration_grn_etc/nov17_2022/input/'

# file made in process_ChIP_H3K4me2.R 

consensus_all = pd.read_csv(basedir + "consensus_atac_H3K4me2_counts_tmm_quan_norm.csv",index_col=0,sep='\t')

new_idx = []
for entry in consensus_all.index:
    new_idx.append(entry.split(':')[0] + ':' + str(int(entry.split(':')[1].split('-')[0]) - 1) + '-' + entry.split(':')[1].split('-')[1])
    
consensus_all['new_idx'] = new_idx
consensus_all = consensus_all.set_index('new_idx')
consensus_all.index.name = None

# average replicates
consensus_all['hESC'] = 0.5*(consensus_all['hESC_rep1'] + consensus_all['hESC_rep2'])
consensus_all['hDE'] = 0.5*(consensus_all['hDE_rep1'] + consensus_all['hDE_rep2'])
consensus_all['hAFE'] = consensus_all['hAFE_rep1'] # 0.5*(consensus_all['hAFE_rep1'] + consensus_all['hAFE_rep2'])
consensus_all['hPFE'] = consensus_all['hPFE_rep1'] # 0.5*(consensus_all['hPFE_rep1'] + consensus_all['hPFE_rep2'])


for entry in ['hESC_rep1', 'hESC_rep2', 'hDE_rep1', 'hDE_rep2',
       'hAFE_rep1',  'hPFE_rep1']:
    del consensus_all[entry]
    
consensus_40k = pd.read_csv(basedir + "consensus_40k.bed",header=None,sep='\t')
consensus_40k['peak_id'] = consensus_40k[0] + ':' + consensus_40k[1].astype(str) + '-' + consensus_40k[2].astype(str)
consensus_40k = consensus_40k.set_index('peak_id')
consensus_40k.index.name = None


consensus_40k_H3K4me2 = consensus_all.loc[list(set(consensus_40k.index).intersection(consensus_all.index)),]

consensus_40k_H3K4me2.to_csv(basedir + "consensus_40k_for_maelstrom_ATAC_H3K4me2_with_avg.bed",sep='\t')


# #### H3K27me3 

# In[ ]:


basedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/per_sample_integration_grn_etc/nov17_2022/input/'

# file made in process_ChIP_H3K27me3.R 

consensus_all = pd.read_csv(basedir + "consensus_atac_H3K27me3_counts_tmm_quan_norm.csv",index_col=0,sep='\t')

new_idx = []
for entry in consensus_all.index:
    new_idx.append(entry.split(':')[0] + ':' + str(int(entry.split(':')[1].split('-')[0]) - 1) + '-' + entry.split(':')[1].split('-')[1])
    
consensus_all['new_idx'] = new_idx
consensus_all = consensus_all.set_index('new_idx')
consensus_all.index.name = None

# average replicates
consensus_all['hESC'] = 0.5*(consensus_all['hESC_rep1'] + consensus_all['hESC_rep2'])
consensus_all['hDE'] = 0.5*(consensus_all['hDE_rep1'] + consensus_all['hDE_rep2'])
consensus_all['hAFE'] = consensus_all['hAFE_rep1'] # 0.5*(consensus_all['hAFE_rep1'] + consensus_all['hAFE_rep2'])
consensus_all['hPFE'] = 0.5*(consensus_all['hPFE_rep1'] + consensus_all['hPFE_rep2'])


for entry in ['hESC_rep1', 'hESC_rep2', 'hDE_rep1', 'hDE_rep2',
       'hAFE_rep1',  'hPFE_rep1','hPFE_rep2']:
    del consensus_all[entry]
    
consensus_40k = pd.read_csv(basedir + "consensus_40k.bed",header=None,sep='\t')
consensus_40k['peak_id'] = consensus_40k[0] + ':' + consensus_40k[1].astype(str) + '-' + consensus_40k[2].astype(str)
consensus_40k = consensus_40k.set_index('peak_id')
consensus_40k.index.name = None


consensus_40k_H3K27me3 = consensus_all.loc[list(set(consensus_40k.index).intersection(consensus_all.index)),]

consensus_40k_H3K27me3.to_csv(basedir + "consensus_40k_for_maelstrom_ATAC_H3K27me3_with_avg.bed",sep='\t')


# #### motif command

# ##### gimme maelstrom --no-filter bed_file hg38 output_dir -N 1

# In[ ]:




