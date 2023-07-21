library(igraph)
library(dplyr)
library(magrittr)
library(leiden)
library(ggplot2)

basedir = '/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/per_sample_integration_grn_etc/notebooks/code_plus_input_final_atac_only/'
samp = 'hPFE'

savedir = basedir 

rna_tpm = read.csv("/Users/LoboM/Dropbox/sharedUMass_Macrina_Rene/hESC_2021/bulk/ATAC/per_sample_integration_GRN_etc/notebooks//rna_avg_tpm.csv",row.names=1)
grn = read.csv(paste(file.path(basedir, samp),'_wt_edges_tops_plus_edge_wt_thresh_10.csv',sep=''),row.names = 1)
# grn contains only those genes pval < 0.01, logfc > 1
pg_rank = read.csv(file.path(basedir,"ranks.csv"))

dim(rna_tpm)
rna_tpm = rna_tpm[rownames(rna_tpm) %in% grn$source,]
dim(rna_tpm)

regulator_df = rna_tpm[rownames(rna_tpm) %in% grn$source,]['hPFE']





set.seed(1)

edge_igraph = graph.data.frame(grn,directed = TRUE)


my_layout = igraph::layout_with_fr(edge_igraph)


regulator_df[igraph::get.vertex.attribute(edge_igraph)[[1]], "fle2"] = my_layout[,2]
regulator_df[igraph::get.vertex.attribute(edge_igraph)[[1]], "fle1"] = my_layout[,1]

regulator_df$gene = rownames(regulator_df)


# add page rank
dim(pg_rank)
pg_rank = pg_rank[pg_rank$X %in% rownames(regulator_df),]
dim(pg_rank)
regulator_df[pg_rank$X,"page_rank"] = pg_rank$PFE


# rank based on expression
regulator_df$rank_rna = rank(-regulator_df$hPFE)



# save grn
regulator_df$label = regulator_df$gene

# add fle's to edge list
regulator_df = data.frame(regulator_df)
rownames(regulator_df) = regulator_df$gene
grn$x_start = regulator_df[grn$source,"fle1"]
grn$y_start = regulator_df[grn$source,"fle2"]
grn$x_end = regulator_df[grn$target,"fle1"]
grn$y_end = regulator_df[grn$target,"fle2"]

#install.packages("viridis")
library(viridis)

ggplot() + 
  geom_segment(data = grn,
               alpha = 0.1,linewidth = 0.5,color = "black",
               mapping = aes(x = x_start, 
                             y = y_start, 
                             xend = x_end, 
                             yend = y_end, ), arrow=arrow(length = unit(0.4, "cm"),angle=10,ends="last", type = "closed"),arrow.fill='black',) + 
  geom_point(
    data = regulator_df, size=5,
    mapping = aes(x = fle1, y = fle2, 
                   colour = page_rank, 
                  label = gene) #, size = page_rank)
  ) + 
  
  
  ggrepel::geom_label_repel(
    data = regulator_df,  # %>% subset(label %in% to_label), #size=10,
    mapping = aes(x = fle1, y = fle2, 
                  # colour = as.character(cluster),
                  label = gene),max.overlaps=Inf
  )  + 
  #scale_colour_manual(values = color_scale) + 
  #scale_fill_manual(values = color_scale)  + 
  scale_size_area(max_size = 2)  + 
  coord_fixed() + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),#legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) + scale_color_viridis(option = "mako",direction = -1)  
ggsave(file.path(savedir,"grn_pfe_tops.pdf"), width = 8, height = 8,limitsize = FALSE)

write.csv(regulator_df,file.path(savedir,"top_pg_rna_pfe.csv"))
