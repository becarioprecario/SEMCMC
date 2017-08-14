#' @name impacts
#' @rdname impacts
#' @title Compute impacts from a Bayesian spatial econometrics model ftited with SEjags.
#'
#' @description This function will compute the impacts (direct, indirect and total) from
#' a SEjags object.
#' @param obj A SEjags object.
#' @param ... Extra argument to compute the impacts.
#' @param W An adjacency matrix, same as used in the call to SEjags().
#' @return A named list with MCMC objects for direct, indirect and 
#' total impacts. 
#' @keywords spatial models
#' @export
#' @examples
#' data(columbus)
#'
#' W <- nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' slm.mcmc <- SEjags(m.form, data = columbus, W = W, model = "slm")
#' impacts(slm.mcmc, W)

impacts <- function(obj, ...) {
    UseMethod("impacts", obj)
}

#' @name impacts.SEjags
#' @rdname impacts
#' @export 
#' @examples
#' data(columbus)
#'
#' W <- nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' slm.mcmc <- SEjags(m.form, data = columbus, W = W, model = "slm")
#' impacts(slm.mcmc, W)

impacts.SEjags <- function(obj, ...) {
  #'...' is simply the W matrix.
  #Lots of checks here

  #Check if we have an intercept
  intercept <- attr(terms(attr(obj, "formula")), "intercept")

  #Define index of variables
  n.var <- dim(obj$b)[1]
  if(intercept) {
    idx.var <- 2:n.var
  } else {
    idx.var <- 1:n.var
  }

  #Variable names
  var.names <- attr(terms(attr(obj, "formula")), "term.labels")

  #Get adj. matrix for SAC model
  if(attr(obj, "model") %in% c("sac", "sacmixed") & class (W) == "list") {
    if(length(W) == 2 & class(W[[1]]) == "matrix") {
      W <- W[[1]]
    } else {
      stop("For SAC model, W must be a matrix or list of length 2.")
    }
  }

  #Obtain variable names

  impacts <- switch(attr(obj, "model"),
    sem = impacts.SEjags.sem(obj, idx.var, var.names),
    slm = impacts.SEjags.slm(obj, W, idx.var, var.names),
    sdm = impacts.SEjags.sdm(obj, W, idx.var, var.names),
    sdem = impacts.SEjags.sdem(obj, W, idx.var, var.names),
    slx = impacts.SEjags.sdem(obj, W, idx.var, var.names),#Same as SDEM
    sac = impacts.SEjags.slm(obj, W, idx.var, var.names),
    sacmixed = impacts.SEjags.sdm(obj, W, idx.var, var.names),
    car = impacts.SEjags.sem(obj, idx.var, var.names),
  )

  return(impacts)
}


#' @rdname impacts.SEjags.xxx
#' @title Compute impacts for different models
#'
#' @description This is an internal function to compute the impacts (direct, indirect 
#' and total)  from a SEjags object.
#' @param obj A SEjags object.
#' @param W Adjacency matrix.
#' @param idx.var An index to subset the covariates coeffiecients
#' @param var.names Vector with variable names
#' @keywords spatial models
#' @examples
#' data(columbus)
#'
#' W <- nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' #SEM model
#' sem.mcmc <- SEjags(m.form, data = columbus, W = W, model = "sem")
#' impacts(sem.mcmc, W)

impacts.SEjags.sem <- function(obj, idx.var, var.names) {

  #No checks here as this is an internal function

  #Av. tot. imp: beta 
  totimp <- t(obj$b[idx.var,1,,1])
  colnames(totimp) <- var.names

  #Av. direct impact: tr((I - rho*W)^{-1} *beta/n
  dirimp <- totimp

  #Indirect impacts
  indimp <- totimp - dirimp

  impacts <- list(direct = dirimp, indirect = indimp,
    total = totimp)
  class(impacts) <- "impacts.SEjags"

  return(impacts)
}

#' @rdname impacts.SEjags.xxx
#' @examples
#' #SLM model
#' slm.mcmc <- SEjags(m.form, data = columbus, W = W, model = "slm")
#' impacts(slm.mcmc, W)
impacts.SEjags.slm <- function(obj, W, idx.var, var.names) {

  #No checks here as this is an internal function

  #Av. tot. imp: beta / (1 - rho)
  rho.sim <- obj$rho[1,,]
  totimp <- apply(obj$b[idx.var,1,,1], 1, function(X) {
   X/(1-rho.sim)
  })
  colnames(totimp) <- var.names

  #Av. direct impact: tr((I - rho*W)^{-1} *beta/n
  nrow.W <- nrow(W)
  tr.aux1 <-  unlist(mclapply(rho.sim, function(X) {
    aux <- mean(diag(solve(diag(nrow.W) - X * W)))
  }))
  dirimp <- apply(obj$b[idx.var,1,,1], 1, function (X) {
   X * tr.aux1
  })
  colnames(dirimp) <- var.names

  #Indirect impacts
  indimp <- totimp - dirimp
  colnames(indimp) <- var.names

  impacts <- list(direct = dirimp, indirect = indimp,
    total = totimp)
  class(impacts) <- "impacts.SEjags"

  return(impacts)
}

#' @rdname impacts.SEjags.xxx
#' @examples
#' sdm.mcmc <- SEjags(m.form, data = columbus, W = W, model = "sdm")
#' impacts(sdm.mcmc, W)

impacts.SEjags.sdm <- function(obj, W, idx.var, var.names) {

  #No checks here as this is an internal function

  #Number of variables: vars. + lagged vars.
  n.var <- length(idx.var)/2

  #Av. tot. imp: beta / (1 - rho)
  rho.sim <- obj$rho[1,,]
  totimp <- sapply(1:n.var, function(X) {
   (obj$b[idx.var[X],1,,1] + obj$b[idx.var[X+n.var],1,,1])/(1-rho.sim)
  })
  colnames(totimp) <- var.names

  #Av. direct impact: tr((I - rho*W)^{-1} *beta/n
  nrow.W <- nrow(W)
  tr.aux1 <-  do.call( rbind, mclapply(rho.sim, function(X) {
    IrhoW.inv <- solve(diag(nrow.W) - X * W)
    aux <- mean(diag(IrhoW.inv))
    auxW <- mean(diag(IrhoW.inv %*% W))
    return(c(aux, auxW))
  }))
  dirimp <- sapply(1:n.var, function(X) {
    tr.aux1[, 1] * obj$b[idx.var[X],1,,1] +
      tr.aux1[, 2] * obj$b[idx.var[X + n.var],1,,1]
  })
  colnames(dirimp) <- var.names

  #Indirect impacts
  indimp <- totimp - dirimp
  colnames(indimp) <- var.names

  impacts <- list(direct = dirimp, indirect = indimp,
    total = totimp)
  class(impacts) <- "impacts.SEjags"

  return(impacts)
}


#' @rdname impacts.SEjags.xxx
#' @examples
#' data(columbus)
#' #SDEM model
#' sdem.mcmc <- SEjags(m.form, data = columbus, W = W, model = "sdem")
#' impacts(sdem.mcmc, W)
impacts.SEjags.sdem <- function(obj, W, idx.var, var.names) {

  #No checks here as this is an internal function

  #Number of variables: vars. + lagged vars.
  n.var <- length(idx.var)/2

  #Av. tot. imp: beta + gamma
  totimp <- t(obj$b[idx.var[1:n.var],1,,1] + obj$b[idx.var[n.var + 1:n.var],1,,1])
  colnames(totimp) <- var.names

  #Av. direct impact: tr((I - rho*W)^{-1} *beta/n
  dirimp <- t(obj$b[idx.var[1:n.var],1,,1])
  colnames(dirimp) <- var.names

  #Indirect impacts
  indimp <- totimp - dirimp

  impacts <- list(direct = dirimp, indirect = indimp,
    total = totimp)
  class(impacts) <- "impacts.SEjags"

  return(impacts)
}

