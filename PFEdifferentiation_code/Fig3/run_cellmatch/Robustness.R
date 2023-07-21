
#' Return a named list of methods for transforming input.
#'
#' @export
#'
GetDefaultTransformations = function(){
  list(
    identity = function(x) x,
    div_by_max  = function(x) {d = max(x);         if( d == 0 ){ x*0 } else {x/d}},
    div_by_sum  = function(x) {d = sum(x);         if( d == 0 ){ x*0 } else {x/d}},
    div_by_norm = function(x) {d = sqrt(sum(x*x)); if( d == 0 ){ x*0 } else {x/d}}
  )
}


#' Return a named list of methods for distance computation.
#'
#' @export
#'
GetDefaultDistances = function(){
   list(
     negative_correlation = function(x, y) {
       sapply(1:nrow(x), function(ii) suppressWarnings(cor(x[ii,], y[ii,]))) %>%
       ReplaceNA %>%
       multiply_by(-1)
     },
     correlation_distance_slow = function(x, y) {
       warning("This is slow and deprecated; use correlation_distance.\n")
       sapply(1:nrow(x), function(ii) suppressWarnings(cor(x[ii,], y[ii,]))) %>%
         ReplaceNA %>%
         subtract(1, .) %>%
         multiply_by(2) %>%
         sqrt
     },
     correlation_distance = function(x, y) {
       center = function(z) sweep(z, 1, STATS =      rowMeans(z  ) , FUN = "-")
       scale  = function(z) sweep(z, 1, STATS = sqrt(rowMeans(z^2)), FUN = "/")
       rowMeans( scale(center(x)) * scale(center(y)) ) %>%
         ReplaceNA %>%
         subtract(1, .) %>%
         multiply_by(2) %>%
         pmax(0) %>%
         sqrt
     },
     squared = function(x, y) {(x-y) %>% raise_to_power(2) %>% rowSums %>% sqrt},
     absolute = function(x, y) {abs(x-y) %>% rowSums},
     is_specific = function(x, y) {
       -EvaluateByGene(x, y, do_heatmaps = F)[["is_shared_cluster_specific"]]
     }
   )
}

#' Check different genesets, distance functions, and gene-wise data transformations.
#'
#' @param genesets Named list of character vectors
#' @param distances_use Named list of functions compatible with the
#' compute_penalty argument of SelectModel. See examples in GetDefaultDistances().
#' @param transforms_use Named list of functions. See examples in GetDefaultTransformations().
#' @param verbose Logical. If TRUE, print progress.
#'
#' @export
#'
ExamineRobustness = function(
  genesets,
  distances_use  = GetDefaultDistances(),
  transforms_use = GetDefaultTransformations(),
  verbose = T
){
  # check input
  assertthat::assert_that(!is.null(names(genesets)))
  assertthat::assert_that(!is.null(names(transforms_use)))
  assertthat::assert_that(!is.null(names(distances_use)))
  # allocate output structure
  dn = list(
    names(transforms_use),
    names(genesets),
    names(distances_use)
  )
  match_model = array(list(NULL), dimnames = dn, dim = sapply(dn, length))
  # do it up
  for(transform_name    in names(transforms_use)){
    for(gu_name         in names(genesets)){
      for(distance_name in names(distances_use)){
        if(verbose){
          cat(transform_name, gu_name, distance_name, sep = "; ")
        }
        gu = genesets[[gu_name]]
        transform = transforms_use[[transform_name]]
        set.seed(0)
        mm = SelectModel(
          query[gu,] %>% apply(2, transform),
          reference[gu,] %>% apply(2, transform),
          compute_penalty = distances_use[[distance_name]]
        )
        match_model[transform_name, gu_name, distance][[1]] = mm
      }
    }
  }
  match_model
}
