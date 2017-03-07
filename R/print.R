#' @name print.impacts.SEjags
#' @rdname print.impacts.SEjags
#' @title Display summary of impacts from fitted with SEjags.
#'
#' @description This function will print  summary of the impacts (direct, indirect and total) from
#' a SEjags object. Output is similar to spdep::impacts
#' @param x A SEjags.impacts object.
#' @param ... Extra argument to compute the impacts (not used).
#' @keywords spatial models
#' @export
#' @examples
#' data(columbus)
#'
#' W <- nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' slm.mcmc <- SEjags(m.form, data = columbus, W = W, model = "slm")
#' impacts(slm.mcmc, W)

print.impacts.SEjags <- function(x, ...) {

  #Lots of checks here

  #Summary statistics
  post.values <- lapply(x, function(X) {
    apply(X, 2, function(Y) { c(mean(Y), sd(Y))})
  })

  #Posterior means
  post.mean <- t(do.call(rbind, lapply(post.values, function(X){X[1,]})))
  post.sd <- t(do.call(rbind, lapply(post.values, function(X){X[2,]})))

 cat("Impact measures (posterior mean):\n")
 print(post.mean)
 cat("\nImpact measures (posterior st. dev.):\n")
 print(post.sd)

 return(list(post.mean = post.mean, post.sd = post.sd))
}

