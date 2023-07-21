#!/usr/bin/env python
# coding: utf-8


import pandas as pd


import gget



all_tfs = pd.read_csv("code_plus_input_final_atac_only/PageRank_scores.csv",index_col=0)


use_tfs = list(all_tfs.max(axis=1).sort_values(ascending=False)[0:150].index)
use_tfs = all_tfs.loc[use_tfs,]

samp_tops = {}
for entry in use_tfs.columns:
    samp_tops[entry] = list(use_tfs.loc[use_tfs.idxmax(axis=1) == entry,].index)



# ## get high confidence GRN for viz



basedir = 'code_plus_input_final_atac_only/'
samp = 'hPFE'
weight_fix = False # fixed weight or use the weights from GRN 

grn = pd.read_csv(basedir + 'Adjacency-' + samp + '.csv',index_col=0) 

grn = grn.loc[samp_tops['PFE'],samp_tops['PFE']]
# subset to tf grn


# subset to require gene to be differentially up at that time point
print(grn.shape)

grn = grn.T # now row = TF, column = gene

# create adjacency matrix
grn_tuples = []
for tf in grn.index: # for every TF
    
    
    for target in grn.columns:
        if grn.loc[tf,target] > 0:
            if weight_fix:
                wt = 1 
            else:
                wt = grn.loc[tf,target]
            # print('{},{},{}'.format(tf,target,wt))
            grn_tuples.append((tf,target,wt))
            
            
# create df from tuples 
grn_df = pd.DataFrame.from_records(grn_tuples,columns=['source','target','weight'])

# save grn 
grn_df.to_csv(basedir + samp + '_wt_edges_tops.csv')

grn_df[grn_df['weight']>10].to_csv(basedir + samp + '_wt_edges_tops_plus_edge_wt_thresh_10.csv')

