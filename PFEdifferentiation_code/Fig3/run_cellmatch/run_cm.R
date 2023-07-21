#install.packages("Seurat")
library(Seurat)
library(scater)
library(magrittr)

setwd('../run_cellmatch/run_cellmatch/')

source("Robustness.R")
source("Utilities.R")
source("ModelSelection2.R")
source("Transparency.R")
source("DiscreteFunctionTools.R")

# file was converted from h5ad -> rds using anndata2ri library in "reference_with_orthologous genes"


load("../objects/roadmap.rds")

colnames(adata_map) <- gsub('_','.',colnames(adata_map), fixed = TRUE)

LabelsPath = "../objects/labels_curated_annotation.csv"

Labels <- as.matrix(read.csv(LabelsPath))

DataMap <- adata_map
remove(adata_map)

# load all timecourse
Labels <- as.vector(Labels[,2])

# DataMap_Mat <- as.matrix(assays(DataMap,"data")[[1]])
DataMap_Mat <- as.matrix(assays(DataMap,"data",withDimnames = TRUE)[[1]])

# the merged object
load("../objects/differentiation_timecourse_esc_to_pfe.rds")
DataTest_IP <- adata_rds
remove(adata_rds)
# DataTest_Mat <- as.matrix(assays(DataTest_IP,"data")[[1]])
DataTest_Mat <- as.matrix(assays(DataTest_IP,"data",withDimnames = TRUE)[[1]])
LabelTest <- as.vector(DataTest_IP$leiden)

NumGenes = 2000
num_init = 10
niter = 20
verbose = T

results_path <- "results/"
#output_dir <- "results/"

# aggrecate by cluster - train + Select informative genes
DataTrain<-AggregateByCluster(DataMap_Mat, group_labels = Labels, method = "average")
DataTrain = as.matrix(DataTrain)
cgenesA<-SelectInformativeGenes(DataTrain, K = NumGenes)
DataTrain<-as.matrix(DataTrain[cgenesA,])

# aggrecate by cluster - test, subset to informative genes
DataTest<-AggregateByCluster(DataTest_Mat, group_labels = LabelTest, method = "average")
DataTest = as.matrix(DataTest)
DataTest <- DataTest[cgenesA,]

matchmodel = SelectModelThoroughly2(query     = DataTest,
                                    reference = DataTrain,
                                    verbose = verbose, num_init = num_init,
                                    compute_penalty = "correlation_distance",niter = niter)
saveRDS(matchmodel, file.path(results_path, "matchmodel.Robj"))




# query = DataTest
# reference = DataTrain
# equivalents = matchmodel$x
# do_heatmaps = T
# results_path = results_path
# compute_penalty = "correlation_distance"
# 
# # sanitize input
# query %<>% as.matrix
# reference %<>% as.matrix
# if( any(query<0) | any(reference<0) ){
#   stop("EvaluateByGene expects nonnegative input.")
# }
# assertthat::are_equal(length(colnames(query)), length(equivalents))
# pairs = data.frame(query = colnames(query),
#                    reference = equivalents,
#                    stringsAsFactors = F)
# pairs[["query"]] %<>% as.character
# pairs[["reference"]] %<>% as.character
# compute_penalty = CheckInputs(query,
#                               reference,
#                               pairs,
#                               compute_penalty,
#                               do_match_ncol = T)
# 
# if(do_heatmaps){
#   dir.create(results_path, recursive = T, showWarnings = F)
# }
# 
# # functions for handling non one-to-one matches
# compatible_query_clusters =
#   aggregate( pairs["query"], by = pairs["reference"], FUN = c) %>%
#   (function(X) setNames(X[[2]], X[[1]]))
# 
# compatible_ref_clusters =
#   aggregate( pairs["query"], by = pairs["reference"], FUN = c) %>%
#   (function(X) setNames(X[[2]], X[[1]]))
# 
# # print(compatible_query_clusters)
# 
# get_max_match = function(gene, cluster_r){
#   all_matches = compatible_query_clusters[[cluster_r]]
#   if(length(all_matches)==0){
#     return("no_match")
#   } else {
#     return(nwm(query[gene, ][all_matches]))
#   }
# }
# get_query_spec = function(gene, cluster_r){
#   all_matches = compatible_query_clusters[[cluster_r]]
#   if(length(all_matches)==0){
#     return(0)
#   } else {
#     return(get_specificity(x = query[gene, ], i = all_matches))
#   }
# }
# 
# # Get cluster-specific markers shared in reference and query
# df_idx = c()
# ref_name = c()
# query_name = c()
# for (entry in rownames(query)){
#   new_val = paste(paste(pairs$reference,pairs$query,sep='-'),entry,sep='|')
#   # print(new_val)
#   df_idx = c(df_idx,new_val)
#   ref_name = c(ref_name,pairs$reference)
#   query_name = c(query_name,pairs$query)
#   
# }
#   
# X = data.frame(ref_query_gene = df_idx,
#                  gene = rep(rownames(query),each=dim(pairs)[1]),
#                ref_name=ref_name,
#                query_name=query_name,
#                penalty = rep(compute_penalty(query[    ,pairs$query],
#                                          reference[,pairs$reference]),each=dim(pairs[1])),
#                # max_cluster_r = apply(reference[, pairs$reference ], 1, nwm),
#                # change specificity to do ratio in that cluster over bgd
#                
#                # gene is assigned to the reference cluster with maximum specificity
#                #specificity_r = apply(reference[, pairs$reference ], 1, get_specificity),
#                stringsAsFactors = F
# )
# 
# spec_r = c()
# spec_q = c()
# for (entry in rownames(X)){
#   spec_r = c(spec_r,get_specificity(reference[X[entry,]$gene,unique(pairs$reference)],X[entry,]$ref_name))
#   
#   spec_q = c(spec_q,get_specificity(query[X[entry,]$gene,unique(pairs$query)],X[entry,]$query_name))
# }
# 
# X$specificity_r = spec_r
# X$specificity_q = spec_q
# 
# 
# # rownames(X) = X$gene
# #  X$max_compatible_cluster_q = with( X, mapply(get_max_match,  gene, max_cluster_r) )
# # originally specificity on the query cluster is computed 
# # as the avg expression of all query clusters the reference maps to divided by the average of the remaining clusters
# # now I'm changing it to query_map/average over remaining maps
# 
# # originally specificity on the reference cluster was computed 
# # as the all clusters where the max index was attained, 
# # this seems buggy since DE and EPI could both have the same expression of that gene 
# # now I want to do it only for the maximum index
# 
# # max_cluster_r is the maximum reference
# # max_compatible_cluster_q is the maximum query cluster mapping to that reference
# 
# # X$specificity_q            = with( X, mapply(get_query_spec, gene, max_cluster_r) )
# make_stage_ordered_factor = function(x) factor(x, levels = unique(x))
# X$ref_name %<>% make_stage_ordered_factor
# X$query_name %<>% make_stage_ordered_factor
# X$is_shared_cluster_specific = with(
#   X,
#   (specificity_r > fold_change_cutoff) &
#     (specificity_q > fold_change_cutoff)
# )




eval_results = EvaluateByGene2(query = DataTest,
                                         reference = DataTrain,
                                         equivalents = matchmodel$x,
                                         do_heatmaps = T,
                                         results_path = results_path,
                                         compute_penalty = "correlation_distance" )

# only write genes that are "shared cluster specific" - added on March 10, 2023 
eval_results2 = eval_results[eval_results$is_shared_cluster_specific == TRUE,]

write.csv(eval_results2, file.path(results_path, "eval.csv"))

neighbors = MutateModel(query     = DataTest,
                                   reference = DataTrain,
                                   init = matchmodel$x)

p = PlotNeighboringModels(neighbors,
                                     #omit_one_model_scores = omissions,
                                     results_path = results_path)


