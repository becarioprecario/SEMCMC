#Some internal functions to data handling


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
