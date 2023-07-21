#' Minimize the given function by coordinate descent on a finite domain.
#'
#' @param FUN Function to minimize. Should take a list or vector as the first input.
#' @param init Initialization -- first list to put into FUN.
#' @param verbose If true, print optimization updates.
#' @param values What to iterate over at each coordinate.
#' @param ... Additional args passed to FUN.
#'
CoordinateDescent = function( FUN, values, init = NULL, niter = 10, verbose = T, ... ){
  assertthat::assert_that( all( init %in% values ) )
  if(is.null(init)){
    stop("You must specify 'init'.\n")
  }

  x = init #random init 
  best_in = init
  best_out = FUN(init, ...) #FUN = correlation penalty function
  objective = matrix(nrow = niter, ncol = length(init), data = NA )
  for( ii in 1:niter){
    previous_out = best_out
    for(coord in seq_along(init)){# for each query cluster
      x = best_in
      for(val in values){#for each ref cluster
        x[[coord]] = val
        current_out = FUN(x, ...)
        if(current_out < best_out){
          best_out = current_out
          best_in = x #index
        }
      } 
      objective[ii, coord] = best_out #corresponding value of index
      if(verbose){
        best_in_display = best_in
        best_in_display[coord] %<>% paste0("|", "|")
        print(best_in_display)
        print(best_out)
      }
    }
    if(previous_out <= best_out){
      break
    }
  }
  return(list(
    x = best_in,
    objective = objective,
    objective_final = best_out
  ))
}


#' At each coordinate of init, omit that coordinate and record the output.
#'
#' @param FUN Function to minimize. Should take a list or vector as the first input.
#' Must also accept an optional argument 'index_omit' giving the coordinate to omit.
#' The function itself must know how to omit a coordinate in a manner appropriate to the context.
#'
#' @param values Ignored. This is only present for consistency with similar functions.
#' @param init Initialization -- first list to put into FUN.
#' @param ... Additional args passed to FUN.
#'
OmitSingleCoords = function( FUN, init = NULL, values = NULL, ... ){
  objective = matrix(nrow = 1, ncol = length(init), data = NA )
  for(coord in seq_along(init)){
    objective[1, coord] = FUN( init, index_omit = coord, ...)
  }
  rownames(objective) = "omitted"
  colnames(objective) = names(init)
  return(list( init = init,
               objective = objective,
               objective_at_init = FUN(init, ...)))
}

#' At each coordinate of init, cycle through all possible values and record the output. (while keeping other queyr-ref matches fixed)
#'
#' @param FUN Function to minimize. Should take a list or vector as the first input.
#' @param init perturbations of each coordinate of this input are studied.
#' @param values What to iterate over at each coordinate.
#' @param ... Additional args passed to FUN.
#'
ComputeSingleCoordPaths = function( FUN, values, init, niter = 10, ... ){
  assertthat::assert_that( all( init %in% values ) )
  objective = matrix(nrow = length(values), ncol = length(init), data = NA )
  for(coord in seq_along(init)){#query clusters
    for(ii in seq_along(values)){#ref clusters
      current_in = init
      current_in[[coord]] = values[[ii]]
      objective[ii, coord] = FUN(current_in, ...)
    }
  }
  rownames(objective) = values
  colnames(objective) = names(init)
  return(list( init = init,
               objective = objective))
}


#' Sample inputs uniformly at random.
#'
#' @param FUN Function to minimize. Should take a list or vector as the first input.
#' @param init If coords_fix is given, \code{init[coords_fix]} remains under the null model and only
#' the remaining entries are sampled.
#' Otherwise, only the length matters, as it determines the dimension of the cube to be sampled.
#' @param values What to iterate over at each coordinate.
#' @param niter How many to sample.
#' @param seed input to set.seed (for repeatable results).
#' @param coords_fix Coordinates to leave fixed. This can generate a different null distribution.
#' @param ... Additional args passed to FUN.
#'
SampleUniformly = function( FUN, values, init, niter, coords_fix = NULL, seed = 0, ... ){
  assertthat::assert_that( all( init %in% values ) )
  objective = rep(NA, niter)
  inputs = rep(NA, niter) %>% as.list
  set.seed(seed)
  for(ii in 1:niter){
    x = sample(values, size = length(init), replace = )
    x[coords_fix] = init[coords_fix]
    inputs[[ii]] = x
    objective[[ii]] = FUN(x)
  }
  return(list(inputs = inputs, objective = objective))
}

