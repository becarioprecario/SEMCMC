#' @rdname utils
#' @name utils
#' @title Internal functions of package SEMCMC.
#' @description Function to convert ouput from jags and stan to a suitable
#'   format to compute the impacts.
#'
#' @return A 'mcmc.list' with the samples of the variables in the model.
#'
#' @param fit Model fitted with SEMCMC().
#'
#' @importFrom coda mcmc
#' @importFrom coda mcmc.list

#Convert stan output to coda style (used when computing the weights)
stan2coda <- function(fit) {
     mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}


#Convert jags output to coda style (used when computing the weights)
jags2coda <- function(fit) {

  mcmc.list(lapply(fit, function(X) {
    if(length(dim(X) == 4 )) {
      mcmc(X)
    } else {
      mcmc(X)
    }
  }))
}
