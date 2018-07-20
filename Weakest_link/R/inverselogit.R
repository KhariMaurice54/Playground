#' logit and inverselogit
#' 
#' These functions have distinct uses:
#'  (1) in defining phi for WLContinuousdata models, and
#'  (2) for H and Hinv in the quantile stitching models.
#'  
#'  @param z A value in (-Inf, Inf)
#'  @param p A value in [0,1]

H <- pnorm ### could also be inverselogit
inverselogit <-
  function(z){
    exp(z)/ (1 + exp(z))
  }

Hinv <- qnorm   ### could also be logit
logit <-
  function(p) {log(p/(1-p))}
