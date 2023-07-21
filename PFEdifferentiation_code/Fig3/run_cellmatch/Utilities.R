
#' Show equivalents on a tSNE or similar.
#'
#' @param embedding tSNE, UMAP, or embedding of choice. Data frame or matrix with two named columns.
#' @param group_labels Character vector of same size and embedding.
#' @param pairs Small dataframe with columns query, reference, and (optionally) color.
#'
#' This function is for when you have a tSNE for a \emph{reference} dataset, with
#' cluster labels from the \emph{reference} dataset. It will show where the query falls
#' on this reference visualization. If you want to do the opposite, i.e. start with a
#' query dataset and use a reference color-scheme, you'll have to trick the function
#' by swapping the column names in the pairs argument.
#'
#' Reference cells assigned multiple query equivalents will be colored
#' as one or another at random.
#'
#' @export
#' @import ggplot2
#'
PlotEquivalentsTSNE = function(embedding, group_labels, pairs,
                               background_color = "grey",
                               background_name = "background_"){
  assertthat::are_equal(nrow(embedding), length(group_labels))
  assertthat::assert_that(is.character(group_labels))
  assertthat::assert_that(is.vector(group_labels))
  assertthat::assert_that(all(pairs$reference %in% group_labels))
  if(is.null(pairs$color)){
    pairs$color = scales::hue_pal()(nrow(pairs))
  }
  embedding_labels = colnames(embedding)
  X = data.frame(embedding[,1],
                 embedding[,2],
                 group_labels,
                 stringsAsFactors = F) %>%
    set_colnames(c(embedding_labels, "original"))

  # fill in query equivalents
  X$equivalent = background_name
  for(ct in pairs$reference %>% unique){
    relevant_pairs = subset(pairs, reference == ct)
    cells = group_labels %>% make.names %>% is.element(relevant_pairs$reference)
    X$equivalent[ cells ] = sample( x = relevant_pairs$query,
                                    replace = T,
                                    size = sum(cells))
  }
  # This is to get the legend in the right order
  X$equivalent %<>% factor(levels = pairs$query %>% c(background_name))
  # this is to plot non-NA on top, but otherwise evenly dispersed
  X = X[order(X$equivalent!=background_name, runif(length(X$equivalent))), ]

  p = ggplot(X, aes_string(x = embedding_labels[[1]],
                           y = embedding_labels[[2]],
                           colour = "equivalent" )) +
    scale_colour_manual( values =
                           pairs$color %>% c(background_color) %>%
                           setNames(c(pairs$query, background_name))
    ) +
    geom_point() +
    coord_fixed() +
    ggtitle("Equivalents") +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())
  p
}


#' Aggregate a single-cell expression matrix.
#'
#' @param expr Single-cell expression matrix with one cell per column.
#' This function expects untransformed normalized expression, e.g. counts
#' per million.
#' @param group_labels Categories within which to aggregate expression.
#' @param scale_genes Should each gene be scaled to a max of 1?
#' @param method How to aggregate. detection_rate finds the proportion of cells
#' expressing a given gene in a given group_labels. Otherwise, values from expr are
#' simply averaged.
#'
#' @export
#'
#' @import magrittr
#'
AggregateByCluster = function(expr, group_labels, method, scale_genes = F){
  if(!is.character(group_labels)){
    warning("group_labels arg should be a character vector; it will be coerced\n")
  }
  group_labels %<>% as.character
  assertthat::are_equal(length(group_labels), ncol(expr))
  indicator_matrix = model.matrix(~group_labels+0, data.frame(group_labels = group_labels, stringsAsFactors = F))
  cnames = indicator_matrix %>% colnames %>% gsub("^group_labels", "", .)
  indicator_matrix = Matrix::Matrix(indicator_matrix)
  if(method == "sum") {
    # no modification needed to indicators or expression values
  } else if(method == "average"){
    indicator_matrix = indicator_matrix %*% diag(1/Matrix::colSums(indicator_matrix))
  } else if( method=="detection_rate" ){
    expr %<>% is_greater_than(0)
    indicator_matrix = indicator_matrix %*% diag(1/Matrix::colSums(indicator_matrix))
  } else  {
    stop(" 'method' should be 'detection_rate' or 'average'. \n")
  }
  Z = expr %*% indicator_matrix
  colnames(Z) = cnames
  if(scale_genes){
    Z %<>% apply(1, div_by_max) %>% t
  }
  return(Z)
}

#' Select informative genes from aggregated pseudo-bulk profiles.
#'
#' @param reference Tall and skinny matrix of cell type profiles.
#' @param K how many genes to return
#' @param regress_out factors to remove when computing gene variances
#' @param mean_cutoff Genes with mean expression below this value are excluded.
#' @param pseudocount cv is computed as sd / (mean + pseudocount).
#'
#' @export
#'
SelectInformativeGenes = function(reference,
                                  regress_out = NULL,
                                  K = 1000,
                                  mean_cutoff = 1,
                                  pseudocount = 0.01){
  if(!is.null(regress_out)){
    resid = reference - t(predict(lm(t(as.matrix(reference)) ~ regress_out)))
  } else {
    resid = reference
  }
  hvg_stats = data.frame(
    gene_sd = apply(resid, 1, sd),
    gene_mean = apply(reference, 1, mean) + pseudocount,
    gene = rownames(reference),
    stringsAsFactors = F
  )
  hvg_stats$gene_cv = with(hvg_stats, gene_sd / gene_mean )

  hvg_stats$excess_var = with(
    hvg_stats, gene_cv - predict(mgcv::gam(gene_cv~gene_mean))
  )
  p = ggplot( hvg_stats ) +
    geom_point(aes(x=gene_mean, y=gene_cv, color=excess_var), size = 0.1) +
    geom_smooth(aes(x=gene_mean, y=gene_cv), method = "gam", formula = y ~ s(x, bs = "cs")) +
    xlab("Log10 CPM") +
    ylab("cv") +
    scale_x_log10()
  print(p)

  n_genes_passing_filter = hvg_stats %>% subset(gene_mean >= mean_cutoff) %>% nrow
  if(n_genes_passing_filter < K){
    warning("Too few genes pass mean expression cutoff. Setting mean_cutoff to -Inf.\n")
    mean_cutoff = -Inf
  }
  dispersion_cutoff = hvg_stats %>%
    subset(gene_mean >= mean_cutoff) %>%
    extract2("excess_var") %>%
    sort(decreasing = T) %>% extract(K)
  hvg = subset(
    hvg_stats,
    gene_mean >= mean_cutoff & excess_var >= dispersion_cutoff,
    select = "gene", drop = T
  )
  return(hvg)
}




#
# Some miscellaneous utility functions for internal use.
#

"+" = function(x,y) {
  if(is.character(x) || is.character(y)) {
    return(paste(x , y, sep=""))
  } else {
    .Primitive("+")(x,y)
  }
}

ReplaceNA = function(x, filler = 0){
  x[is.na(x)] = filler
  x
}

PadVector = function(x, n, filler = "") {
  if(length(x)>n){
    stop("x is too long.\n")
  } else if(length(x)<n){
    x = c(x, rep(filler, n - length(x)))
  }
  return(x)
}

nwm = function(x) names(which.max(x))

#' Helper function for EvaluateByGenes.
#'
#'
get_specificity = function(x, i = NULL){
  # If index not given, use all indices where max is attained
  if(length(i)==0){
    i = which(max(x)==x) # returns multiple values if they are same, if multiple different clusters have the same expression, 
    # they will also get returned which is a problem
  }
  m1 = mean(x[i]) # in case  of multiple maximum's mean will be same as actual maximum, , otherwise will mean over all i's
  m2 = (sum(x) - sum(x[i])) / (length(x) - length(x[i])) # mean of everything other than m1
  if(m1==0) {
    return(0)
  } else {
    m1/m2 # ratio of max over others
  }
}


div_by_max = function(x) {if(max(x)>0) x/max(x) else 0*x}


#' Extract counts per million, not transformed, from a Seurat version 2 object.
#'
#' @export
#'
get_cpm = function(seurat_object) {
  seurat_object %<>% (Seurat::NormalizeData)
  seurat_object@data@x %<>% exp %>% subtract(1) %>% multiply_by(1e2)
  seurat_object@data
}

#' Pad a matrix with zeroes to facilitate merging of expression data with
#' non-identical gene-sets.
#'
#' @export
#'
PadMatrix = function(X, desired_rownames){
  missing = desired_rownames %>% setdiff(rownames(X))
  if(length(missing)==0){
    return(X)
  }
  Z = matrix(0, ncol = ncol(X), nrow  = length(missing))
  rownames(Z) = missing
  colnames(Z) = colnames(X)
  rbind(X, Z) %>% set_rownames(c(rownames(X), rownames(Z))) %>% extract(desired_rownames, )
}

#' Merge matrices assuming that rownames are overlapping but non-identical sets of genes.
#' For genes present in only some matrices, zeroes are added in the other matrices.
#'
#' @param expr List of matrices
#'
#' @export
#'
MergeMatrixList = function(expr){
  all_genes = lapply(expr, rownames) %>% Reduce(f=union) %>% setdiff("")
  expr %<>% lapply(cellmatch::PadMatrix, all_genes)
  Reduce(expr, f = cbind)
}
