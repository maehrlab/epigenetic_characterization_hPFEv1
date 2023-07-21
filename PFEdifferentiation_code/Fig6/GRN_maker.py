#!/usr/bin/env python

# Code By: Ryan Abramowitz

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import pandas as pd
import seaborn as sns
import numpy as np
from datetime import datetime
import random
import gc
import warnings
plt.rcParams['pdf.fonttype'] = 42
plt.xticks(rotation='horizontal')
plt.yticks(rotation='horizontal')
warnings.filterwarnings(action='once')



figfolder = 'topTFs'
os.system(f'rm -r {figfolder}')
os.makedirs(figfolder,exist_ok=True)
os.system(f'cp GRN_maker.py {figfolder}/')

def run_stuff(save_adjacency_files=True):
    np.random.seed(0)
    random.seed(0)
    corrcutoff = 0
    basedir = './GRN_Ryan_ATAC_only/'

    global_rna = pd.read_csv(basedir+'TPM_counts_all.csv',index_col=0)
    global_rna = np.log1p(global_rna)
    
    global_rna = pd.DataFrame({s:global_rna.loc[:,global_rna.columns.str.contains(s)].mean(axis=1) for s in ['hESC','hDE','hAFE','hPFE']})
    global_rna.index = global_rna.index.str.upper() # RNA upper
    global_personalization_rna = global_rna.T
    global_personalization_rna -= global_personalization_rna.mean()
    global_personalization_rna /= global_personalization_rna.std()
    global_personalization_rna = np.exp(global_personalization_rna).fillna(0).T

    def toRank(series):
        series = series.sort_values(ascending=False)
        torank = {}
        for i,s in enumerate(series):
            torank[s] = torank.get(s,i+1)
        return series.apply(torank.get)


    figdir = f'{figfolder}/output'
    print('='*250)
    print(figdir,datetime.now())
    print('='*250)
    os.makedirs(figdir,exist_ok=True)

    corrdir = 'RNA_Motif_correlation/all_motif_rna_correlation_fig4.csv'

    filter_genes = lambda df: df.loc[(df['padj'] < 0.1) * (np.abs(df['log2FoldChange']) > 0.5),'SYMBOL']
    not_housekeeping_global = list(set([j.upper() for i in os.listdir('GRN_Ryan/pairwise_differential(RNA)') if i.endswith('.csv') for j in filter_genes(pd.read_csv(f'GRN_Ryan/pairwise_differential(RNA)/{i}'))])) # not housekeeping upper
    motif2TFs_df = pd.read_csv(corrdir)[['MotifName','TF','rna_motif_correlation']]
    motif2TFs_df['TF'] = motif2TFs_df['TF'].apply(lambda s: s.upper()) # motif2tf upper
    motif2TFs_df['rna_motif_correlation'] = motif2TFs_df['rna_motif_correlation']
    motif2TFs_df = motif2TFs_df.loc[motif2TFs_df['TF'].isin(not_housekeeping_global)]
    dic_motif2TFs_global = {}
    for _,row in motif2TFs_df.iterrows():
        if row['rna_motif_correlation'] >= corrcutoff:
            dic_motif2TFs_global[row['MotifName']] = dic_motif2TFs_global.get(row['MotifName'],[]) + [row['TF']] # weigh by number of correlation
    TFlist_global = list(set([j for i in dic_motif2TFs_global.values() for j in i]))
    TFlist_global.sort()
    print('TFlist_global:',len(TFlist_global),TFlist_global[:10])

    for sample in ['hPFE','hAFE','hDE','hESC',]:
        print(f'on sample {sample}...')

        rna = global_rna.loc[not_housekeeping_global,sample]
        dic_motif2TFs = {k:[i for i in v if i in rna.index] for k,v in dic_motif2TFs_global.items()} # <- ALL SOURCES ARE NOT HOUSEKEEPING
        dic_motif2TFs = {k:[(v0,1/len(v)) for v0 in v] for k,v in dic_motif2TFs_global.items()}  # weigh by number of TFs

        loc2genes = pd.read_csv(basedir + sample + f"/{sample}_loops.bed",sep='\t',header=None)
        loc2genes['peak'] = loc2genes[0].astype('str') + ':' + (loc2genes[1]+1).astype('str') + '-' + loc2genes[2].astype('str')
        loc2genes = loc2genes.loc[loc2genes[6].isin(not_housekeeping_global)].reset_index() # <- ALL TARGETS ARE NOT HOUSEKEEPING
        loc2genes['gene'] = loc2genes[6]
        loc2genes = loc2genes[['peak','gene']]

        loc2motifs = pd.read_csv(basedir + sample + f'/{sample}_activator_regions_with_footprints.bed',sep='\t',header=None)
        loc2motifs['motif'] = loc2motifs[len(loc2motifs.columns)-3]
        loc2motifs['peak'] = loc2motifs[0].astype('string') + ':' + (loc2motifs[1]+1).astype('string') + '-' + loc2motifs[2].astype('string')
        loc2motifs['tfs'] = loc2motifs['motif'].apply(lambda x: dic_motif2TFs.get(x,[]))
        loc2motifs = loc2motifs[['tfs','peak',]].dropna()

        gene2locs = {}
        for _,row in loc2genes.iterrows():
            gene2locs[row['gene']] = gene2locs.get(row['gene'],[]) + [row['peak']]
        loc2tfs = {}
        for _,row in loc2motifs.iterrows():
            for tf_and_normalization in row['tfs']:
                loc2tfs[row['peak']] = loc2tfs.get(row['peak'],[]) + [tf_and_normalization]
        def list2dict(items):
            countdict = {}
            for tf,w in items:
                countdict[tf] = countdict.get(tf, 0) + w
            return countdict
        print('Making Adjacency!!!')
        adjacency = pd.DataFrame({g:list2dict([(tf,rna[tf]*rna[g]*normalization) for p in ps for tf,normalization in loc2tfs.get(p,[])]) for g,ps in gene2locs.items()}).fillna(0)
        for i in adjacency.columns[adjacency.columns.isin(adjacency.index)==False]:
            adjacency.loc[i] = 0
        adjacency = adjacency.T
        for i in adjacency.columns[adjacency.columns.isin(adjacency.index)==False]:
            adjacency.loc[i] = 0
        adjacency = adjacency.loc[adjacency.columns,adjacency.columns]
        if save_adjacency_files:
            adjacency.to_csv(f'{figdir}/Adjacency-{sample}.csv')
        # columns = tf, row = genes
        G = nx.from_numpy_matrix(np.array(adjacency),create_using=nx.DiGraph())
        G = nx.relabel_nodes(G, dict(enumerate(adjacency.columns)))


        print('Running PageRank...')
        pr = nx.pagerank(G,weight='weight',personalization=dict(global_personalization_rna.loc[list(G.nodes),sample]))
        pagerankRNA = pd.DataFrame(pr,index=['pagerank']).T
        pagerankRNA = pagerankRNA.sort_values('pagerank',ascending=False)
        pagerankRNA.to_csv(f'{figdir}/pagerank_with_RNA{sample}.csv')

    print('Done with GRN creations! Time for PageRank result analysis...')

    sns.color_palette('cividis')
    pagerank_df = pd.DataFrame([],columns=[],index=None)
    for i in ['ESC','DE','AFE','PFE']:
        this_df = pd.read_csv(f'{figdir}/pagerank_with_RNAh{i}.csv',index_col=0,)
        this_df.columns=[i]
        pagerank_df = pd.concat([pagerank_df,this_df],axis=1,join='outer')
    pagerank_df = pagerank_df.fillna(0)

    tfs = pagerank_df.loc[TFlist_global]
    ordering = tfs.max(axis=1).sort_values(ascending=False).index
    tfs = tfs.loc[ordering]
    tfs.to_csv(f'{figdir}/PageRank_scores.csv')
    
    def plotheatmap_multiple_ways(toplot,name):
        plt.figure(figsize=[toplot.shape[0]//55 + 4,toplot.shape[0]*12/50])
        type2order = {'ESC':0,'DE':1,'AFE':2,'PFE':3}
        df = toplot.loc[toplot.apply(lambda x: (type2order[x.idxmax()],3*x['PFE']+x['AFE']-x['DE']-3*x['ESC']),axis=1).sort_values().index]
        sns.heatmap(df, cmap='hot', cbar_kws={'label': "pagerank zscore"},xticklabels=True,yticklabels=True,linewidths=1, linecolor='black',square=True)
        plt.savefig(f'{figdir}/{name}_max_ordered.pdf')
        plt.close(); plt.cla(); plt.clf(); gc.collect()

        for celltype in ['PFE','AFE','DE','ESC']:
            df = toplot.loc[toplot.apply(lambda x: x.idxmax() == celltype,axis=1)]
            plt.figure(figsize=[toplot.shape[0]//55 + 4,toplot.shape[0]*12/50])
            plt.tight_layout()
            g = sns.heatmap(df, cmap='hot', xticklabels=True,yticklabels=True,linewidths=1, linecolor='black',square=True, cbar=False)
            sns.set(rc = {'figure.figsize':[toplot.shape[0]//55 + 4,toplot.shape[0]*12/50]})
            g.tick_params(axis='x', rotation=0)
            g.tick_params(axis='y', rotation=0)
            plt.yticks(fontname = 'Arial')
            plt.xticks(fontname = 'Arial')
            plt.tight_layout()
            plt.savefig(f'{figdir}/{name}_max_ordered - {celltype}.pdf')
            plt.close(); plt.cla(); plt.clf(); gc.collect()
    sig_tfs = tfs
    sig_tfs -= sig_tfs.values.mean()
    sig_tfs /= sig_tfs.values.std(ddof=1)
    for n in [100,125,150]:
        plotheatmap_multiple_ways(((tfs.T - tfs.mean(axis=1))/tfs.std(axis=1)).T.loc[ordering[:n]],f'SCOREPLOT-zscore_genes_top{n}')
        # plotheatmap_multiple_ways(tfs.iloc[:n],f'SCOREPLOT-pagerank_scores_top{n}')
    tfranks = pd.DataFrame({c:toRank(pagerank_df[c]) for c in ['ESC','DE','AFE','PFE']})
    tfranks = tfranks.loc[ordering]
    tfranks.to_csv(f'{figdir}/ranks.csv')

    plt.subplots(figsize=(4,60))
    # x = np.log2(tfranks.loc[TFlist_global]).sort_index()
    x = np.log2(tfranks.loc[ordering])
    sns.heatmap(x,xticklabels=-1, cmap='cubehelix_r',square=True,cbar_kws={'label':'log2(Ranking)'},vmin=0,vmax=np.log2(len(ordering)))
    plt.savefig(f'{figdir}/total_tf_ranks_heatmap.pdf')

    print(f'Finished with {figdir} at {datetime.now()}!!!')
    return None

print('Starting at:',datetime.now())
run_stuff(save_adjacency_files=True)
print('Done at:',datetime.now())
