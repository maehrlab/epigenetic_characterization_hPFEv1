#' Count genes that support each matching decision.
#'
#' @param gene_eval Results of EvaluateByGene, or RunCellMatch(...)$eval_results .
#'
#' @export
#'
CountSupportingGenes = function(gene_eval){
  subset( gene_eval, is_shared_cluster_specific, select = "max_compatible_cluster_q" ) %>% table %>% sort
}






#' Evaluate matches between conformable count matrices.
#'
#' @param query @param reference Expression matrices to be matched.
#' Each column should be one cluster. Each row should be one gene or genomic feature.
#' The genes should be shared and ordered the same way.
#' @param equivalents a list of columns from reference to be matched with the query.
#' @param do_heatmaps obvious
#' @param results_path where to save heatmaps
#' @param ... Additional parameters passed to ClusterGenes
#'
#' @export
#'
EvaluateByGene2 = function(
    query,
    reference,
    equivalents = colnames(reference),
    compute_penalty = "correlation_distance",
    fold_change_cutoff = 2,
    results_path = NULL,
    do_heatmaps = !is.null(results_path),
    reorder_genes = F,
    pseudocount = 0.1,min_gene_pct=1,
    ...
){


  # sanitize input
  query %<>% as.matrix
  reference %<>% as.matrix
  if( any(query<0) | any(reference<0) ){
    stop("EvaluateByGene expects nonnegative input.")
  }
  assertthat::are_equal(length(colnames(query)), length(equivalents))
  pairs = data.frame(query = colnames(query),
                     reference = equivalents,
                     stringsAsFactors = F)
  pairs[["query"]] %<>% as.character
  pairs[["reference"]] %<>% as.character
  compute_penalty = CheckInputs(query,
                                reference,
                                pairs,
                                compute_penalty,
                                do_match_ncol = T)

  if(do_heatmaps){
    dir.create(results_path, recursive = T, showWarnings = F)
  }

  # functions for handling non one-to-one matches
  compatible_query_clusters =
    aggregate( pairs["query"], by = pairs["reference"], FUN = c) %>%
    (function(X) setNames(X[[2]], X[[1]]))

  compatible_ref_clusters =
    aggregate( pairs["reference"], by = pairs["query"], FUN = c) %>%
    (function(X) setNames(X[[2]], X[[1]]))

  # print(compatible_query_clusters)

  get_max_match = function(gene, cluster_r){
    all_matches = compatible_query_clusters[[cluster_r]]
    if(length(all_matches)==0){
      return("no_match")
    } else {
      return(nwm(query[gene, ][all_matches]))
    }
  }
  get_query_spec = function(gene, cluster_r){
    all_matches = compatible_query_clusters[[cluster_r]]
    if(length(all_matches)==0){
      return(0)
    } else {
      return(get_specificity(x = query[gene, ], i = all_matches))
    }
  }

  # Get cluster-specific markers shared in reference and query
  df_idx = c()
  ref_name = c()
  query_name = c()
  for (entry in rownames(query)){
    new_val = paste(paste(pairs$reference,pairs$query,sep='-'),entry,sep='|')
    # print(new_val)
    df_idx = c(df_idx,new_val)
    ref_name = c(ref_name,pairs$reference)
    query_name = c(query_name,pairs$query)

  }

  X = data.frame(ref_query_gene = df_idx,
                   gene = rep(rownames(query),each=dim(pairs)[1]),
                 ref_name=ref_name,
                 query_name=query_name,
                 penalty = rep(compute_penalty(query[    ,pairs$query],
                                           reference[,pairs$reference]),each=dim(pairs)[1]),
                 # max_cluster_r = apply(reference[, pairs$reference ], 1, nwm),
                 # change specificity to do ratio in that cluster over bgd

                 # gene is assigned to the reference cluster with maximum specificity
                 #specificity_r = apply(reference[, pairs$reference ], 1, get_specificity),
                 stringsAsFactors = F
  )

  spec_r = c()
  spec_q = c()
  for (entry in rownames(X)){
    spec_r = c(spec_r,get_specificity(reference[X[entry,]$gene,unique(pairs$reference)],X[entry,]$ref_name))

    spec_q = c(spec_q,get_specificity(query[X[entry,]$gene,unique(pairs$query)],X[entry,]$query_name))
  }

  X$specificity_r = spec_r
  X$specificity_q = spec_q


  # rownames(X) = X$gene
  #  X$max_compatible_cluster_q = with( X, mapply(get_max_match,  gene, max_cluster_r) )
  # originally specificity on the query cluster is computed
  # as the avg expression of all query clusters the reference maps to divided by the average of the remaining clusters
  # now I'm changing it to query_map/average over remaining maps

  # originally specificity on the reference cluster was computed
  # as the all clusters where the max index was attained,
  # this seems buggy since DE and EPI could both have the same expression of that gene
  # now I want to do it only for the maximum index

  # max_cluster_r is the maximum reference
  # max_compatible_cluster_q is the maximum query cluster mapping to that reference

  # X$specificity_q            = with( X, mapply(get_query_spec, gene, max_cluster_r) )
  make_stage_ordered_factor = function(x) factor(x, levels = unique(x))
  X$ref_name %<>% make_stage_ordered_factor
  X$query_name %<>% make_stage_ordered_factor
  X$is_shared_cluster_specific = with(
    X,
    (specificity_r > fold_change_cutoff) &
      (specificity_q > fold_change_cutoff)
  )
  
  
  # for each query cluster, print warning if match is too general
  for (entry in unique(pairs$query)){
    X_tmp = X[((X$is_shared_cluster_specific == TRUE) & (X$query_name == entry)),]
    # print(dim(X_tmp)[1])
    # print(dim(reference)[1]/(100*min_gene_pct))
    if ((dim(X_tmp)[1]) < (dim(reference)[1]/(100*min_gene_pct))){
      # print("hi")
      print(sprintf("WARNING!!! Classification of %s to %s has few (%d) specific genes", entry,compatible_ref_clusters[[entry]],dim(X_tmp)[1]))
      
    }
  }
  
  # Optional plots
   if( do_heatmaps ){
    
    # Cluster the remaining genes and make heatmaps for them.
    remaining_genes = setdiff(rownames(reference),unique(X[X$is_shared_cluster_specific==TRUE,]$gene))
    Y = ClusterGenes( query[    remaining_genes, ] %>% apply(1, div_by_max) %>% t ,
                      reference[remaining_genes, ] %>% apply(1, div_by_max) %>% t ,
                      ... )
    
    Y$cluster = paste("other_patterns",Y$cluster,sep="_")
    
    matched_genes = unique(X[X$is_shared_cluster_specific==TRUE,]$gene)
    Z = ClusterGenes( query[    matched_genes, ] %>% apply(1, div_by_max) %>% t ,
                      reference[matched_genes, ] %>% apply(1, div_by_max) %>% t ,
                      ... )
    Z$cluster = paste("shared_specific",Z$cluster,sep="_")
    
    gene_cluster = rbind(Y,Z)
    
    #X[["final_cluster"]] = X$max_compatible_cluster_q %>% as.character %>% paste0("cluster_specific_", .)
    #X %<>% set_rownames(X$gene)
    
    # add cluster to eval results 
    X = merge(x = X, y = gene_cluster,all.x=TRUE)
    
    
    # X[Y$gene, "final_cluster"] = Y$cluster %>% as.character %>% paste0("other_patterns_", .)
    X %<>% (dplyr::arrange)(cluster, penalty)
    # rownames(X) = X$gene
    
    # Summary of penalties by cluster
    median_penalty = aggregate(X$penalty, X["cluster"], FUN = median) %>% (dplyr::arrange)(x)
    X$cluster %<>% factor(levels = median_penalty$cluster %>% rev)
    
    
    
    X$cluster_genecount = ave(rep(1, nrow(X)), X$cluster, FUN = sum)/dim(pairs)[1]
    
    X$expr_q = rowMeans(query[X$gene, ])
    X$expr_r = rowMeans(reference[X$gene, ])
    
    # gene_spec tells if the gene is specific somewhere or specific nowhere
    X$gene_spec = TRUE
    X[substr(X$cluster,1,4) == 'othe',]$gene_spec = FALSE
    
    X_tmp = unique(X[,c('gene','penalty','cluster','cluster_genecount','expr_q','expr_r')])
    rownames(X_tmp) = X_tmp$gene
    X_tmp$is_spec_somewhere = TRUE
    X_tmp[unique(X[X$gene_spec == FALSE, ]$gene),]$is_spec_somewhere = FALSE
    # specific in at least one cluster vs not specific anywhere 
    p = ggplot(X_tmp, aes(fill = cluster_genecount, y = penalty, x = cluster)) +
      geom_violin() +
      geom_jitter(size = 0.1) +
      facet_wrap(~X_tmp$is_spec_somewhere, scales = "free_x") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(file.path(results_path, "penalty_vs_expression_pattern.pdf"), p, width = 10, height = 8)
    
    
    p = ggplot(X_tmp, aes(color = penalty, y = log10(expr_q + pseudocount), x = log10(expr_r + pseudocount))) +
      geom_point() +
      geom_smooth(method = lm) +
      facet_wrap(~X_tmp$is_spec_somewhere, scales = "free_x") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(file.path(results_path, "penalty_vs_expr_level.pdf"), p, width = 10, height = 8)
    
    
    # rownames(X) = X$gene
    save_one_plot = function(type, k, gu_k, use_matches = T, X_plot = X_tmp){
      if(length(gu_k)==0){
        return()
      }
      if(use_matches){
        ref = reference[gu_k, pairs$reference, drop = F]
      } else {
        ref = reference[gu_k, , drop = F]
      }
      
      # set row index to gene names 
      if (k=='all'){
        rownames(X_plot) = paste(X_plot$gene,X_plot$query_name,sep='_')
        genes_to_label = 'none'
       }
      else{ # one cluster has no overlap genes 
        rownames(X_plot) = X_plot$gene
        genes_to_label = 'all'
        }
      
      
    
    
    
    
      p = PairedHeatmap2(query = query[     gu_k, pairs$query    , drop = F],
                        reference = ref,
                        genes = gu_k,
                        scores = X_plot$penalty/2,
                        reorder_genes = F, genes_to_label = genes_to_label) +
        scale_fill_viridis_c(option = "magma")
      dir.create(file.path(results_path, type), showWarnings = F, recursive = T)
      ggsave(file.path(results_path, type, paste0("k=", k, ".pdf")),
             p, width = 10, height = pmin(4 + length(gu_k)/4, 45))
      if(k=="all"){
        ggsave(file.path(results_path, type, paste0("k=", k, "_short.pdf")),
               p, width = 10, height = 6)
      }
      p
    }
    
    adderall = function(x) x %>% as.list %>% c(list("all"))
    
    
    # cluster-specific genes
    ref_clusters = pairs$reference
    q_clusters = pairs$query
    
    X_tmp = X[X$is_shared_cluster_specific == TRUE,]
    X_tmp$order_temp = X_tmp$query_name %>% factor(levels = colnames(query))
    X_tmp = X_tmp[order(X_tmp$order_temp),]
    
    
    for(k in q_clusters ){
      # ref_clusters[["all"]] = "placeholder" # prevents error on ref_clusters[["all"]]
      #gu_k = X_tmp %>%
      #  subset(  query_name %in% k, select = "gene", drop = T)
      
      # subset to only genes that are used 
      X_tmp_clus = X_tmp[(X_tmp$query_name == k),]
      # sort above by penalty
      X_tmp_clus %<>% (dplyr::arrange)(penalty)
      gu_k = X_tmp_clus$gene
      save_one_plot(type = "cluster_specific", k, gu_k, use_matches = T, X_plot=X_tmp_clus)
    }
    
    # plot all shared genes
    X_tmp %<>% (dplyr::arrange)(query_name,penalty)
    save_one_plot(type = "cluster_specific", 'all', X_tmp$gene, use_matches = T, X_plot=X_tmp)
    
    
    # Remaining genes - genes that are not shared specific anywhere
    # gene _ spec == FALSE
    
    X_tmp = unique(X[X$gene_spec==FALSE,c('gene','penalty','cluster','cluster_genecount','expr_q','expr_r')])
    
    rownames(X_tmp) = X_tmp$gene
    # Y = Y[order(Y$cluster, X_tmp[Y$gene, "penalty"]),]
    
    for(k in unique(X_tmp$cluster) ){
        
        # subset to only genes that are used 
        X_tmp_clus = X_tmp[(X_tmp$cluster == k),]
        # sort above by penalty
        X_tmp_clus %<>% (dplyr::arrange)(penalty)
        gu_k = X_tmp_clus$gene
        save_one_plot(type = "other_patterns", k, gu_k, use_matches = F, X_plot=X_tmp_clus)
      }
  }
  return(X)
}






#' Evaluate matches between conformable count matrices.
#'
#' @param query @param reference Expression matrices to be matched.
#' Each column should be one cluster. Each row should be one gene or genomic feature.
#' The genes should be shared and ordered the same way.
#' @param equivalents a list of columns from reference to be matched with the query.
#' @param do_heatmaps obvious
#' @param results_path where to save heatmaps
#' @param ... Additional parameters passed to ClusterGenes
#'
#' @export
#'
EvaluateByGene = function(
  query,
  reference,
  equivalents = colnames(reference),
  compute_penalty = "correlation_distance",
  fold_change_cutoff = 2,
  results_path = NULL,
  do_heatmaps = !is.null(results_path),
  reorder_genes = F,
  pseudocount = 0.1,
  ...
){
  
  # sanitize input
  query %<>% as.matrix
  reference %<>% as.matrix
  if( any(query<0) | any(reference<0) ){
    stop("EvaluateByGene expects nonnegative input.")
  }
  assertthat::are_equal(length(colnames(query)), length(equivalents))
  pairs = data.frame(query = colnames(query),
                     reference = equivalents,
                     stringsAsFactors = F)
  pairs[["query"]] %<>% as.character
  pairs[["reference"]] %<>% as.character
  compute_penalty = CheckInputs(query,
                                reference,
                                pairs,
                                compute_penalty,
                                do_match_ncol = T)

  if(do_heatmaps){
    dir.create(results_path, recursive = T, showWarnings = F)
  }

  # functions for handling non one-to-one matches
  compatible_query_clusters =
    aggregate( pairs["query"], by = pairs["reference"], FUN = c) %>%
    (function(X) setNames(X[[2]], X[[1]]))
  
  # print(compatible_query_clusters)
  
  get_max_match = function(gene, cluster_r){
    all_matches = compatible_query_clusters[[cluster_r]]
    if(length(all_matches)==0){
      return("no_match")
    } else {
      return(nwm(query[gene, ][all_matches]))
    }
  }
  get_query_spec = function(gene, cluster_r){
    all_matches = compatible_query_clusters[[cluster_r]]
    if(length(all_matches)==0){
      return(0)
    } else {
      return(get_specificity(x = query[gene, ], i = all_matches))
    }
  }

  # Get cluster-specific markers shared in reference and query
  X = data.frame(
    gene = rownames(query),
    penalty = compute_penalty(query[    ,pairs$query],
                              reference[,pairs$reference]),
    max_cluster_r = apply(reference[, pairs$reference ], 1, nwm),
    specificity_r = apply(reference[, pairs$reference ], 1, get_specificity),
    stringsAsFactors = F
  )
  
  
  rownames(X) = X$gene
  X$max_compatible_cluster_q = with( X, mapply(get_max_match,  gene, max_cluster_r) )
  X$specificity_q            = with( X, mapply(get_query_spec, gene, max_cluster_r) )
  make_stage_ordered_factor = function(x) factor(x, levels = unique(x))
  X$max_cluster_r %<>% make_stage_ordered_factor
  X$max_compatible_cluster_q %<>% make_stage_ordered_factor
  X$is_shared_cluster_specific = with(
    X,
    (specificity_r > fold_change_cutoff) &
      (specificity_q > fold_change_cutoff)
  )

  # Optional plots
  if( do_heatmaps ){

    # Cluster the remaining genes and make heatmaps for them.
    remaining_genes = X %>% subset(!is_shared_cluster_specific ) %>% extract2("gene")
    Y = ClusterGenes( query[    remaining_genes, ] %>% apply(1, div_by_max) %>% t ,
                      reference[remaining_genes, ] %>% apply(1, div_by_max) %>% t ,
                      ... )
    X[["final_cluster"]] = X$max_compatible_cluster_q %>% as.character %>% paste0("cluster_specific_", .)
    X %<>% set_rownames(X$gene)
    X[Y$gene, "final_cluster"] = Y$cluster %>% as.character %>% paste0("other_patterns_", .)
    X %<>% (dplyr::arrange)(final_cluster, penalty)
    rownames(X) = X$gene

    # Summary of penalties by cluster
    median_penalty = aggregate(X$penalty, X["final_cluster"], FUN = median) %>% (dplyr::arrange)(x)
    X$final_cluster %<>% factor(levels = median_penalty$final_cluster %>% rev)
    X$cluster_genecount = ave(rep(1, nrow(X)), X$final_cluster, FUN = sum)
    X$expr_q = rowMeans(query[X$gene, ])
    X$expr_r = rowMeans(reference[X$gene, ])
    p = ggplot(X, aes(fill = cluster_genecount, y = penalty, x = final_cluster)) +
      geom_violin() +
      geom_jitter(size = 0.1) +
      facet_wrap(~X$is_shared_cluster_specific, scales = "free_x") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(file.path(results_path, "penalty_vs_expression_pattern.pdf"), p, width = 10, height = 8)
    p = ggplot(X, aes(color = penalty, y = log10(expr_q + pseudocount), x = log10(expr_r + pseudocount))) +
      geom_point() +
      geom_smooth(method = lm) +
      facet_wrap(~X$is_shared_cluster_specific, scales = "free_x") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(file.path(results_path, "penalty_vs_expr_level.pdf"), p, width = 10, height = 8)

    rownames(X) = X$gene
    save_one_plot = function(type, k, gu_k, use_matches = T){
      if(length(gu_k)==0){
        return()
      }
      if(use_matches){
        ref = reference[gu_k, pairs$reference, drop = F]
      } else {
        ref = reference[gu_k, , drop = F]
      }
      p = PairedHeatmap(query = query[     gu_k, pairs$query    , drop = F],
                        reference = ref,
                        genes = gu_k,
                        scores = X[gu_k, "penalty"]/2,
                        reorder_genes = F) +
        scale_fill_viridis_c(option = "magma")
      dir.create(file.path(results_path, type), showWarnings = F, recursive = T)
      ggsave(file.path(results_path, type, paste0("k=", k, ".pdf")),
             p, width = 10, height = pmin(4 + length(gu_k)/4, 45))
      if(k=="all"){
        ggsave(file.path(results_path, type, paste0("k=", k, "_short.pdf")),
               p, width = 10, height = 6)
      }
      p
    }

    adderall = function(x) x %>% as.list %>% c(list("all"))
    # cluster-specific genes
    ref_clusters = pairs$reference
    X$order_temp = X$max_compatible_cluster_q %>% factor(levels = colnames(query))
    X = X[order(X$order_temp),]
    for(k in seq_along(ref_clusters) %>% adderall ){
      ref_clusters[["all"]] = "placeholder" # prevents error on ref_clusters[["all"]]
      gu_k = X %>%
        subset(is_shared_cluster_specific) %>%
        subset( k=="all" | max_cluster_r %in% ref_clusters[[k]], select = "gene", drop = T)
      save_one_plot(type = "cluster_specific", k, gu_k, use_matches = T)
    }

    # Remaining genes
    Y = Y[order(Y$cluster, X[Y$gene, "penalty"]),]
    for(k in Y$cluster %>% unique %>% adderall ){
      gu_k = Y %>% subset( k=="all" | cluster==k , select = "gene", drop = T)
      save_one_plot(type = "other_patterns", k, gu_k, use_matches = F)
    }

  }
  return(X)
}

#' Helper function for commonly used inputs.
#'
CheckInputs = function(query,
                       reference,
                       pairs = NULL,
                       compute_penalty = "correlation_distance",
                       do_match_ncol ){
  assertthat::are_equal(nrow(query), nrow(reference))
  assertthat::are_equal(rownames(query), rownames(reference))
  assertthat::assert_that( all(pairs$query     %in% colnames(query    )) )
  assertthat::assert_that( all(pairs$reference %in% colnames(reference)) )
  assertthat::assert_that( is.matrix(query))
  assertthat::assert_that( is.matrix(reference))
  assertthat::assert_that( is.numeric(query))
  assertthat::assert_that( is.numeric(reference))

  if( do_match_ncol ){
    assertthat::are_equal( ncol( query ), ncol( reference ) )
  }
  if( is.character(compute_penalty) && compute_penalty %in% names( GetDefaultDistances() ) ){
    compute_penalty = GetDefaultDistances()[[compute_penalty]]
  }
  return(compute_penalty)
}

#' Helper function for EvaluateByGene.
#'
#' @param query @param reference Expression matrices to be matched.
#' Each column should be one cluster. Each row should be one gene or genomic feature.
#' The genes should be shared and ordered the same way.
#' @param K number of clusters
#'
#' @export
#'
ClusterGenes = function(query, reference, K = 20){
  query %<>% as.matrix
  reference %<>% as.matrix
  CheckInputs(query, reference, do_match_ncol = F)
  combined_matrix = cbind(query, reference)
  cluster_mod = kmeans(combined_matrix, centers = K, iter.max = 100, nstart = 20)
  Y = data.frame( gene = rownames(query), cluster = cluster_mod$cluster, stringsAsFactors = F )
  return(Y)
}



PairedHeatmap2 = function(query, reference, genes, reorder_genes = TRUE,
                         genes_to_label = ifelse(length(genes)>50, "none", "all"),
                         query_prefix = "q_", reference_prefix = "r_", include_difference = F, scores = NULL){
query = query %>% as.matrix
reference = reference %>% as.matrix
CheckInputs(query, reference, do_match_ncol = F)
query     %<>% as.matrix
reference %<>% as.matrix
colnames(query)     %<>% make.unique("__dup")
colnames(reference) %<>% make.unique("__dup")
assertthat::assert_that(query_prefix!=reference_prefix)
colnames(query)     %<>% paste0(query_prefix, .)
colnames(reference) %<>% paste0(reference_prefix, .)

if(reorder_genes & length(genes)>1){
  perm_vec = seriation::seriate(dist(cbind(query, reference)[genes,]), "TSP")[[1]] %>% as.numeric
  genes_in_order = genes[perm_vec]
} else {
  genes_in_order = genes
}

if(genes_to_label=="none"){
  genes_to_label = NULL
} else if(genes_to_label=="all"){
  genes_to_label = genes_in_order
}

# Heatmap these genes
make_plot_df = function(values_show, query_ref ){
  cn = c("gene", "cluster", "expression")
  df  = values_show[ genes_in_order, , drop = F] %>%
    as.matrix %>%
    (reshape2::melt) %>%
    set_colnames(cn)
  
  df[["expression"]] %<>% ave(df[["gene"]], FUN = div_by_max )
  df[["gene"]] %<>% as.character  %>% factor(levels = unique(genes_in_order))
  df[['query_ref']] = query_ref %>% factor(levels = c("reference", "query", "difference"))
  df
}

plot_df = rbind(
  make_plot_df(values_show = query,     query_ref = "query"),
  make_plot_df(values_show = reference, query_ref = "reference")
)

if(include_difference){
  difference = apply(query, 1, div_by_max) - apply(reference, 1, div_by_max)
  difference = t(difference)
  difference = (difference + 1)/2 # This avoids wrecking the color scale for the rest of the plot
  plot_df %<>% rbind( make_plot_df(values_show = difference, query_ref = "difference") )
}

plot_df$cluster %<>% factor(levels = c(colnames(query), colnames(reference)))

p = ggplot(plot_df) +
  geom_tile(aes(cluster, gene, fill = expression)) +
  facet_grid(~query_ref, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_discrete(name = "",
                   breaks = genes_in_order,
                   labels = ifelse(genes_in_order %in% genes_to_label, genes_in_order, "")) +
  ggtitle(length(genes_in_order) + " genes")

if(!is.null(scores)){
  assertthat::are_equal(length(scores), length(genes))
  scores_visual_width = ncol(query) / 5
  scores = scores * scores_visual_width
  lateral_offset = ncol(query) + 1
  vertical_offset = nrow(query)
  panel =  if(include_difference) "difference" else "query"
  p = p + geom_tile(aes(y = y, x = x, width = width), fill = "black",
                    data = data.frame( y = genes,
                                       query_ref = panel,
                                       x     = lateral_offset + scores/2,
                                       width = scores ))
  p = p + geom_text(label = "Gene scores",
                    x = lateral_offset + scores_visual_width/2,
                    y = nrow(query) + vertical_offset*0.05,
                    data = data.frame(query_ref = panel))
  p = p + geom_text(label = "0",
                    x = lateral_offset,
                    y = nrow(query) + vertical_offset*0.025,
                    data = data.frame(query_ref = panel))
  p = p + geom_text(label = "2",
                    x = lateral_offset + scores_visual_width,
                    y = nrow(query) + vertical_offset*0.025,
                    data = data.frame(query_ref = panel))
  #this ghost rectangle moves the bounding boxes out
  p = p + geom_tile(aes(x=x,y=y), alpha=0,
                    data = data.frame(#query_ref = panel,
                      x = lateral_offset + scores_visual_width + 1,
                      y = vertical_offset + 2 ))
  p = p + geom_segment(
    aes(y=y,x=x,yend=yend, xend=xend),
    data.frame(
      y = 0 %>% rep(2),
      yend = vertical_offset %>% rep(2),
      x =    lateral_offset + c(0, scores_visual_width),
      xend = lateral_offset + c(0, scores_visual_width),
      query_ref = panel)
  )
}
p
}


#' Depict conformable count matrices.
#'
#' @param query @param reference Expression matrices to be matched.
#' Each column should be one cluster. Each row should be one gene or genomic feature.
#' The dimensions should be equal.
#' @param genes Genes to show.
#' @param reorder_genes If true, use the seriation TSP method to reorder the genes.
#' @param query_prefix @param reference_prefix Prefixes to distinguish column labels. Must be different.
#' @param genes_to_label "none", "all", or a character vector of gene names.
#' @param include_difference Boolean. Include a third panel with query minus reference?
#' @param scores scores. If provided, these are shown as a marginal barplot. The plot is
#' visually optimized for situations where the scores range from 0 to 1.
#'
#' @export
#'

PairedHeatmap = function(query, reference, genes, reorder_genes = TRUE,
                         genes_to_label = ifelse(length(genes)>50, "none", "all"),
                         query_prefix = "q_", reference_prefix = "r_", include_difference = F, scores = NULL){
  query %<>% as.matrix
  reference %<>% as.matrix
  CheckInputs(query, reference, do_match_ncol = F)
  query     %<>% as.matrix
  reference %<>% as.matrix
  colnames(query)     %<>% make.unique("__dup")
  colnames(reference) %<>% make.unique("__dup")
  assertthat::assert_that(query_prefix!=reference_prefix)
  colnames(query)     %<>% paste0(query_prefix, .)
  colnames(reference) %<>% paste0(reference_prefix, .)

  if(reorder_genes & length(genes)>1){
    perm_vec = seriation::seriate(dist(cbind(query, reference)[genes,]), "TSP")[[1]] %>% as.numeric
    genes_in_order = genes[perm_vec]
  } else {
    genes_in_order = genes
  }

  if(genes_to_label=="none"){
    genes_to_label = NULL
  } else if(genes_to_label=="all"){
    genes_to_label = genes_in_order
  }

  # Heatmap these genes
  make_plot_df = function(values_show, query_ref ){
    cn = c("gene", "cluster", "expression")
    df  = values_show[ genes_in_order, , drop = F] %>%
      as.matrix %>%
      (reshape2::melt) %>%
      set_colnames(cn)
    df[["expression"]] %<>% ave(df[["gene"]], FUN = div_by_max )
    df[["gene"]] %<>% as.character  %>% factor(levels = genes_in_order)
    df[['query_ref']] = query_ref %>% factor(levels = c("reference", "query", "difference"))
    df
  }

  plot_df = rbind(
    make_plot_df(values_show = query,     query_ref = "query"),
    make_plot_df(values_show = reference, query_ref = "reference")
  )
  if(include_difference){
    difference = apply(query, 1, div_by_max) - apply(reference, 1, div_by_max)
    difference = t(difference)
    difference = (difference + 1)/2 # This avoids wrecking the color scale for the rest of the plot
    plot_df %<>% rbind( make_plot_df(values_show = difference, query_ref = "difference") )
  }

  plot_df$cluster %<>% factor(levels = c(colnames(query), colnames(reference)))

  p = ggplot(plot_df) +
    geom_tile(aes(cluster, gene, fill = expression)) +
    facet_grid(~query_ref, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_discrete(name = "",
                     breaks = genes_in_order,
                     labels = ifelse(genes_in_order %in% genes_to_label, genes_in_order, "")) +
    ggtitle(length(unique(plot_df$gene)) + " genes")

  if(!is.null(scores)){
    assertthat::are_equal(length(scores), length(genes))
    scores_visual_width = ncol(query) / 5
    scores = scores * scores_visual_width
    lateral_offset = ncol(query) + 1
    vertical_offset = nrow(query)
    panel =  if(include_difference) "difference" else "query"
    p = p + geom_tile(aes(y = y, x = x, width = width), fill = "black",
                      data = data.frame( y = genes,
                                         query_ref = panel,
                                         x     = lateral_offset + scores/2,
                                         width = scores ))
    p = p + geom_text(label = "Gene scores",
                      x = lateral_offset + scores_visual_width/2,
                      y = nrow(query) + vertical_offset*0.05,
                      data = data.frame(query_ref = panel))
    p = p + geom_text(label = "0",
                      x = lateral_offset,
                      y = nrow(query) + vertical_offset*0.025,
                      data = data.frame(query_ref = panel))
    p = p + geom_text(label = "2",
                      x = lateral_offset + scores_visual_width,
                      y = nrow(query) + vertical_offset*0.025,
                      data = data.frame(query_ref = panel))
    #this ghost rectangle moves the bounding boxes out
    p = p + geom_tile(aes(x=x,y=y), alpha=0,
                      data = data.frame(#query_ref = panel,
                                        x = lateral_offset + scores_visual_width + 1,
                                        y = vertical_offset + 2 ))
    p = p + geom_segment(
      aes(y=y,x=x,yend=yend, xend=xend),
      data.frame(
        y = 0 %>% rep(2),
        yend = vertical_offset %>% rep(2),
        x =    lateral_offset + c(0, scores_visual_width),
        xend = lateral_offset + c(0, scores_visual_width),
        query_ref = panel)
    )
  }
  p
}

#' Call transparency functions on high-scoring models.
#'
#' @param neighboring_model_scores Output from MutateModel.
#' @param results_path Where to store plots.
#'
FullComparisonReport = function(neighboring_model_scores, results_path){
  is_successful = standardize(neighboring_model_scores$objective) < -1
  comparisons = list()
  for(ii in colnames(query)){
    alternatives = is_successful[,ii] %>% which %>% names
    if(length(alternatives) <= 1){
      next
    }
    models_to_compare = list()
    for(aa in alternatives){
      aa_full = neighboring_model_scores$init %>% setNames(colnames(query))
      aa_full[[ii]] = aa
      models_to_compare[[aa]] = aa_full
    }
    models_to_compare %<>% as.data.frame
    comparisons[[ii]] = CompareByGene(query,
                               reference,
                               results_path = file.path(results_path, ii),
                               models_to_compare = models_to_compare)
  }
  comparisons
}


#' Display (and save to disk) missing genes from a given cell type.
#'
#' @param query @param reference normalized count matrices sharing rownames (gene names).
#' @param query_baseline @param reference_baseline Character. Baseline for each species or condition: for example, pluripotent cells.
#' Defaults to first column of query or reference respectively.
#' @param query_mature @param reference_mature Character. Clusters to be compared when looking for missing genes.
#' Defaults to last column of query or reference respectively.
#'
#'
#' @export
#'
PlotMissingGenes = function(query,
                            reference,
                            query_baseline     = colnames(query)[[    1]],
                            reference_baseline = colnames(reference)[[1]],
                            query_mature       = colnames(query)[[    ncol(query)]],
                            reference_mature   = colnames(reference)[[ncol(reference)]],
                            results_path,
                            genes = rownames(reference),
                            genes_label = NULL,
                            fc_cutoff = 2,
                            cpm_cutoff = 1,
                            prop_genes_label = 0.01,
                            main = "Change from baseline to mature state"){
  results_path = file.path(results_path, "missing_genes")
  dir.create(results_path, recursive = T)
  assertthat::assert_that(query_baseline     %in% colnames(query))
  assertthat::assert_that(reference_baseline %in% colnames(reference))
  assertthat::assert_that(query_mature     %in% colnames(query))
  assertthat::assert_that(reference_mature %in% colnames(reference))
  genes %<>% intersect(rownames(reference))
  genes %<>% intersect(rownames(query))
  assertthat::assert_that(length(genes)>1)
  X = data.frame(
    expr_r = log2( 1 + reference[genes, reference_mature]),
    expr_q = log2( 1 + query[    genes, query_mature]),
    log2_fc_q      = log2( 1 + reference[genes, reference_mature]) - log2( 1 + reference[genes,reference_baseline]),
    log2_fc_r      = log2( 1 + query[    genes,     query_mature]) - log2( 1 + query[    genes,    query_baseline]),
    gene = genes,
    stringsAsFactors = F
  )
  X$expressed_in_r = X$expr_r > log2(cpm_cutoff)
  X$expressed_in_q = X$expr_q > log2(cpm_cutoff)
  X$abs_log_fc = with(X, abs(log2_fc_r - log2_fc_q) )
  X$designation = NA
  X$designation = with(X, ifelse( expressed_in_r & log2_fc_r - log2_fc_q > log(fc_cutoff), "missing_in_query", designation ) )
  X$designation = with(X, ifelse( expressed_in_q & log2_fc_q - log2_fc_r > log(fc_cutoff), "missing_in_reference", designation ) )

  if(prop_genes_label>0){
    extra_labels = ggrepel::geom_label_repel(data = subset(X, !is.na(designation) & abs_log_fc>quantile(abs_log_fc, 1-prop_genes_label, na.rm = T)))
  } else {
    extra_labels = geom_point()
  }
  p = ggplot(X, aes(log2_fc_r, log2_fc_q, label = gene, color = designation)) +
    geom_point() +
    ggrepel::geom_label_repel(data = subset(X, gene %in% genes_label)) +
    extra_labels +
    coord_fixed() +
    geom_abline(aes(intercept = 0, slope = 1)) +
    xlab("Reference log2 fc (mature state over baseline)") +
    ylab("Query log2 fc (mature state over baseline)") +
    ggtitle(main)
  ggsave(file.path(results_path, "missing_genes.pdf"), p, width = 6, height = 6)
  write.csv(X %>% dplyr::arrange(designation, -abs_log_fc), file.path(results_path, "missing_genes.csv"))
  list(table = X, plot = p)
}
