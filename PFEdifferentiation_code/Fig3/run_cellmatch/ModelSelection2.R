
#' Match query and reference RNA profiles and render a full report.
#'
#' @param query @param reference Expression matrices to be matched.
#' Each column should be one cluster. Each row should be one gene or genomic feature.
#' The genes should be shared and ordered the same way.
#' @param results_path Where to save outputted files.
#' @param K How many genes to select.
#' @param num_init Number of random restarts for coordinate descent. Higher means slower and better quality.
#' @param compute_penalty a gene-wise penalty function acting on two matrices of the same size,
#' and producing numeric vector output with one entry per gene. Make sure big is bad.
#' Accepts shorthand: try "correlation_distance" (recommended), "squared", or "absolute".
#'
#' @export
#'
RunCellMatch2 = function(query, reference, results_path, K = 2000, num_init = 5, verbose = T, niter=10){
  query %<>% as.matrix
  reference %<>% as.matrix
  dir.create(results_path, recursive = T)
  
  # gene selection
  shared_genes = intersect(rownames(query), rownames(reference))
  assertthat::assert_that(length(shared_genes)>=K)
  variable_genes = cellmatch::SelectInformativeGenes(reference[shared_genes, ], K = K)
  
  query %<>% as.matrix
  reference %<>% as.matrix
  
  # Run the model. With verbose=T, it will print the matches and the objective value as it goes.
  matchmodel = SelectModelThoroughly2(query     = query[variable_genes,],
                                                reference = reference[variable_genes,],
                                                verbose = verbose, num_init = num_init,
                                                compute_penalty = "correlation_distance",niter = niter)
  saveRDS(matchmodel, file.path(results_path, "matchmodel.Robj"))
  
  # Show plots justifying the model.
  eval_results = cellmatch::EvaluateByGene(query = query[variable_genes,],
                                           reference = reference[variable_genes,],
                                           equivalents = matchmodel$x,
                                           do_heatmaps = T,
                                           results_path = results_path,
                                           compute_penalty = "correlation_distance" )
  
  write.csv(eval_results, file.path(results_path, "eval.csv"))
  
  # Check neighboring models for goodness of fit.
  neighbors = cellmatch::MutateModel(query     = query[variable_genes,],
                                     reference = reference[variable_genes,],
                                     init = matchmodel$x)
  # omissions = cellmatch::MutateModel(query     = query[variable_genes,],
  #                                    reference = reference[variable_genes,],
  #                                    init = matchmodel$x,
  #                                    do_omit = T)
  p = cellmatch::PlotNeighboringModels(neighbors,
                                       #omit_one_model_scores = omissions,
                                       results_path = results_path)
  
  list(
    matchmodel = matchmodel,
    eval_results = eval_results,
    support = CountSupportingGenes(eval_results),
    neighbors = neighbors,
    heatmap = p
  )
}


#' Tools for model selection and comparison, implemented via DiscreteFunctionTools.
#'
#' @param query @param reference Expression matrices to be matched.
#' Each column should be one cluster. Each row should be one gene or genomic feature.
#' The genes should be shared and ordered the same way.
#' @param compute_penalty a gene-wise penalty function acting on two matrices of the same size,
#' and producing numeric vector output with one entry per gene. Make sure big is bad.
#' Accepts shorthand: try "correlation_distance" (recommended), "squared", or "absolute".
#' @param penalty_is_genewise For complete flexibility with the penalty, set this to FALSE and
#' pass in any penalty function that returns a scalar. Make sure big is bad.
#' @param wrap_fun Routine to wrap around the outside.
#' Try CoordinateDescent, ComputeSingleCoordPaths, or SampleUniformly.
#' @param values @param init @param ... Passed to wrap_fun.
#'
#' @export
#'
StudyModel2 = function(
  query,
  reference,
  compute_penalty = "correlation_distance",
  penalty_is_genewise = T,
  wrap_fun,
  values = colnames(reference),
  init = sample(values, ncol(query), replace = T),
  niter = 10,
  ...
){
  query %<>% as.matrix
  reference %<>% as.matrix
  compute_penalty = CheckInputs(query,
                                reference,
                                compute_penalty = compute_penalty,
                                do_match_ncol = T)
  
  EvaluateModel = function(x, index_omit = NULL){
    index_include = setdiff(1:ncol(query), index_omit)
    x = x[index_include]
    compute_penalty(query[,index_include], reference[,x]) %>% mean
  }
  
  wrap_fun(#default wrap_fun = CoordinateDescent
    FUN = EvaluateModel,
    values = values,
    init = init,
    niter = niter,
    ... )
}


#' Form matches between columns of two count matrices.
#'
#' Works like SelectModelFast but with multiple random initializations.
#'
#' @export
#'
SelectModelThoroughly2 = function(  query,
                                   reference,
                                   compute_penalty = "correlation_distance",
                                   penalty_is_genewise = T,
                                   values = colnames(reference),
                                   num_init = 5,
                                   niter = 10,
                                   ...){
  query %<>% as.matrix
  reference %<>% as.matrix
  best_obj = Inf
  best_model = NULL
  for(ii in 1:num_init){
    current = SelectModelFast2(query = query,
                              reference = reference,
                              compute_penalty = compute_penalty,
                              penalty_is_genewise = penalty_is_genewise,
                              values = values, 
                              niter = niter, ...)
    if(current$objective_final < best_obj){
      best_obj = current$objective_final
      best_model = current
    }
  }
  return(best_model)
}

#' Form matches between columns of two count matrices.
#'
#' Use coordinate descent to select a (locally) optimal matching model.
#' See \code{?StudyModel}; this calls that with CoordinateDescent as wrap_fun.
#' See also \code{?SelectModelThoroughly}.
#'
#' @export
#'
SelectModelFast2 = function(  query,
                             reference,
                             compute_penalty = "correlation_distance",
                             penalty_is_genewise = T,
                             values = colnames(reference),
                             init = sample(values, ncol(query), replace = T),
                             niter = 10,
                             ...){
  query %<>% as.matrix
  reference %<>% as.matrix
  StudyModel2( query,
              reference,
              compute_penalty = compute_penalty,
              penalty_is_genewise = T,
              wrap_fun = CoordinateDescent,
              values = colnames(reference),
              init = init,
              niter = niter,
              ...
  )
}


#' Study all neighboring matching models by "mutating" the model at each stage.
#'
#' @param do_omit Boolean. See \code{?StudyModel}. When do_omit is false (default),
#' StudyModel is called with ComputeSingleCoordPaths as wrap_fun. Otherwise,
#' StudyModel is called with OmitSingleCoords as wrap_fun.
#'
#' @export
#'
MutateModel = function(  query,
                         reference,
                         do_omit = F,
                         compute_penalty = "correlation_distance",
                         penalty_is_genewise = T,
                         values = colnames(reference),
                         init,
                         ...){
  query %<>% as.matrix
  reference %<>% as.matrix
  
  if(!do_omit){
    edits = StudyModel2( query,
                        reference,
                        compute_penalty = compute_penalty,
                        penalty_is_genewise = T,
                        wrap_fun = ComputeSingleCoordPaths,
                        values = colnames(reference),
                        init = init,
                        ... )
  } else {
    edits = StudyModel2( query,
                        reference,
                        compute_penalty = compute_penalty,
                        penalty_is_genewise = T,
                        wrap_fun = OmitSingleCoords,
                        init = init,
                        ... )
  }
  colnames(edits$objective) = colnames(query)
  edits
}

#' Randomly draw \code{niter} models and compute their penalty functions.
#'
#' See \code{?StudyModel} and \code{?SampleUniformly}.
#'
#' @export
#'
SampleModels = function( query,
                         reference,
                         compute_penalty = "correlation_distance",
                         penalty_is_genewise = T,
                         values = colnames(reference),
                         init = sample(values, size = ncol(query), replace = T),
                         niter = 100,
                         ...){
  query %<>% as.matrix
  reference %<>% as.matrix
  StudyModel2( query,
              reference,
              compute_penalty = compute_penalty,
              penalty_is_genewise = T,
              wrap_fun = SampleUniformly,
              values = colnames(reference),
              niter = niter,
              init = init,
              ... )
}


#' Display output from MutateModel.
#'
#' @param neighboring_model_scores Output from MutateModel (default mode)
#' @param omit_one_model_scores Output from MutateModel (with do_omit=T)
#' @param omit_one_display_mode How to display scores of omitted models. One of c("sidebar", "main", "both").
#' @param results_path Where to save output.
#' @param plotname filename of plot. Should end in a ggplot-supported file extension
#' such as .png or .pdf.
#'
#' @export
#'
PlotNeighboringModels = function(neighboring_model_scores,
                                 query_order  = NULL,
                                 ref_order  = NULL,
                                 omit_one_model_scores = NULL,
                                 omit_one_display_mode = "both",
                                 results_path,
                                 plotname = "alternatives.pdf"){
  assertthat::assert_that(omit_one_display_mode %in% c("sidebar", "main", "both"))
  dir.create(results_path, showWarnings = F, recursive = T)
  if(!is.null(omit_one_model_scores) & omit_one_display_mode %in% c("main", "both")){
    scores_use = rbind(neighboring_model_scores$objective, omit_one_model_scores$objective)
  } else {
    scores_use = neighboring_model_scores$objective
  }
  if(is.null(query_order)) {query_order = colnames(scores_use)}
  if(is.null(ref_order  )) {ref_order   = rownames(scores_use)}
  
  obj_long = scores_use %>%
    reshape2::melt() %>%
    set_colnames(c("reference", "query", "correlation_distance"))
  obj_long$query     %<>% as.character %>% factor(levels = query_order)
  obj_long$reference %<>% as.character %>% factor(levels = ref_order)
  
  p = ggplot(obj_long) +
    ggtitle( "Alternative models" ) +
    geom_tile(aes(x = query, y = reference, fill = correlation_distance)) +
    # Mark the selected model with a black dot
    annotate(geom = "point",
             x = neighboring_model_scores$objective %>% colnames %>% as.character,
             y = neighboring_model_scores$init %>% as.character)  +
    scale_fill_viridis_c(direction = -1) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
  
  if(!is.null(omit_one_model_scores) & omit_one_display_mode %in% c("sidebar", "both") ){
    changes_in_fit = omit_one_model_scores$objective[1,] - omit_one_model_scores$objective_at_init
    # Manually map colors, red-white-blue with red < 0, white=0, blue > 0
    colormap = colorRampPalette(c("red", "white", "blue"))(101)
    put_in_window = function(x) round(51 + ( 50*x / max(abs(x)) ))
    colors = colormap[put_in_window(changes_in_fit)]
    # Add sidebar to plot
    p = p + annotate( geom = "tile",
                      x = omit_one_model_scores$objective %>% colnames,
                      y = -2,
                      fill = colors )
  }
  
  ggsave(file.path(results_path, plotname), p, height = 8, width = 8)
  return(p)
}


#' Return another matching model that seems just as good.
#'
#' @param neighboring_model_scores Output from MutateModel.
#'
#' @export
#'
GetReasonableAlternative = function(neighboring_model_scores){
  get_something = function(x, what) {
    w = sort(x, decreasing = F)
    if( what == "gap"){
      return((w[1] - w[2])/w[1])
    } else if( what == "top"){
      return(names(w[1]))
    }else if( what == "second"){
      return(names(w[2]))
    }
    
  }
  gaps = neighboring_model_scores$objective %>% apply(2, get_something, "gap")
  ifelse(
    gaps > mean(gaps),
    neighboring_model_scores$objective %>% apply(2, get_something, "second") ,
    neighboring_model_scores$objective %>% apply(2, get_something, "top")
  )
}



PlotNullDistribution = function(null_distribution,
                                optimal_model_score,
                                results_path,
                                plotname = "null_distribution.pdf"){
  dir.create(results_path, showWarnings = F, recursive = T)
  pdf(file.path(results_path, plotname))
  both = c(null_distribution, optimal_model_score)
  hist(null_distribution,
       xlim=c(min(both),max(both)),
       40, probability = T,
       main = "Random pairings versus optimal model",
       xlab = "Penalty",
       ylab = "Null distribution probability density")
  abline(v = optimal_model_score,
         col = "red")
  dev.off()
}
