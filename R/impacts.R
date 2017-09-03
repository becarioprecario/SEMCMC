#' @name impacts
#' @aliases impacts
#' @rdname impacts
#' @title Compute impacts from a Bayesian spatial econometrics model ftited with SEMCMC.
#'
#' @description This function will compute the impacts (direct, indirect and total) from
#' a SEMCMC object.
#' @param obj A SEMCMC object.
#' @param ... Extra argument to compute the impacts.
#' @return A named list with MCMC objects for direct, indirect and 
#' total impacts. 
#' @keywords spatial models
#' @export
#' @examples
#' data(columbus, package = "spdep")
#'
#' W <- spdep::nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' slm.mcmc <- SEMCMC(m.form, data = columbus, W = W, model = "slm")
#' impacts(slm.mcmc, W)

impacts <- function(obj, ...) {
    UseMethod("impacts", obj)
}

#' @name impacts.SEMCMC
#' @aliases impacts.SEMCMC
#' @rdname impacts
#' @export 
#' @examples
#' data(columbus, package = "spdep")
#'
#' W <- spdep::nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' slm.mcmc <- SEMCMC(m.form, data = columbus, W = W, model = "slm")
#' impacts(slm.mcmc, W)

impacts.SEMCMC <- function(obj, ...) {
  #'...' is simply the W matrix.
  #Lots of checks here

  #So far, impacts only implemented for Gaussian responses
  link <- attr(obj, "link")
  if(link %in% c("probit", "logit") | link != "indentity") {
    stop("Impacts not implemented for this type of link.")
  }

  #Check if we have an intercept
  intercept <- attr(terms(attr(obj, "formula")), "intercept")


  nvar <- attr(obj, "nvar")
  #Define index of variables
  n.var <- dim(obj$b)[1]
  if(intercept) {
    idx.var <- 2:nvar
  } else {
    idx.var <- 1:nvar
  }

  #Variable names
  var.names <- attr(terms(attr(obj, "formula")), "term.labels")

  #Sampler
  if(attr(obj, "sampler") == "jags") {
     objres <- as.data.frame(do.call(rbind, obj$results))
  } else {
     objres <- as.data.frame(do.call(rbind, stan2coda(obj$results)))
  }

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
    sem = impacts.SEMCMC.sem(objres, idx.var, var.names),
    slm = impacts.SEMCMC.slm(objres, W, idx.var, var.names),
    sdm = impacts.SEMCMC.sdm(objres, W, idx.var, var.names),
    sdem = impacts.SEMCMC.sdem(objres, W, idx.var, var.names),
    slx = impacts.SEMCMC.sdem(objres, W, idx.var, var.names),#Same as SDEM
    sac = impacts.SEMCMC.slm(objres, W, idx.var, var.names),
    sacmixed = impacts.SEMCMC.sdm(objres, W, idx.var, var.names),
    car = impacts.SEMCMC.sem(objres, idx.var, var.names),
  )

  return(impacts)
}


#' @rdname impacts.SEMCMC.xxx
#' @aliases impacts.SEMCMC.sem
#' @title Compute impacts for different models
#'
#' @description This is an internal function to compute the impacts (direct, indirect 
#' and total)  from a SEMCMC object.
#' @param obj A SEMCMC object.
#' @param idx.var An index to subset the covariates coeffiecients
#' @param var.names Vector with variable names
#' @keywords spatial models
#' @examples
#' data(columbus, package = "spdep")
#'
#' W <- spdep::nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' #SEM model
#' sem.mcmc <- SEMCMC(m.form, data = columbus, W = W, model = "sem")
#' impacts(sem.mcmc, W)

impacts.SEMCMC.sem <- function(obj, idx.var, var.names) {

  #No checks here as this is an internal function

  #Av. tot. imp: beta 
  totimp <- obj[,idx.var]
  colnames(totimp) <- var.names

  #Av. direct impact: tr((I - rho*W)^{-1} *beta/n
  dirimp <- totimp

  #Indirect impacts
  indimp <- totimp - dirimp

  impacts <- list(direct = dirimp, indirect = indimp,
    total = totimp)
  class(impacts) <- "impacts.SEMCMC"

  return(impacts)
}

#' @rdname impacts.SEMCMC.xxx
#' @aliases impacts.SEMCMC.slm
#'
#' @importFrom parallel mclapply
#'
#' @examples
#' data(columbus, package = "spdep")
#'
#' W <- spdep::nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' #SLM model
#' slm.mcmc <- SEMCMC(m.form, data = columbus, W = W, model = "slm")
#' impacts(slm.mcmc, W)
impacts.SEMCMC.slm <- function(obj, W, idx.var, var.names) {

  #No checks here as this is an internal function

  #Av. tot. imp: beta / (1 - rho)
  rho.sim <- obj$rho
  totimp <- apply(obj[, idx.var], 2, function(X) {
   X/(1-rho.sim)
  })
  colnames(totimp) <- var.names

  #Av. direct impact: tr((I - rho*W)^{-1} *beta/n
  nrow.W <- nrow(W)
  tr.aux1 <-  unlist(mclapply(rho.sim, function(X) {
    aux <- mean(diag(solve(diag(nrow.W) - X * W)))
  }))
  dirimp <- apply(obj[, idx.var], 2, function (X) {
   X * tr.aux1
  })
  colnames(dirimp) <- var.names

  #Indirect impacts
  indimp <- totimp - dirimp
  colnames(indimp) <- var.names

  impacts <- list(direct = dirimp, indirect = indimp,
    total = totimp)
  class(impacts) <- "impacts.SEMCMC"

  return(impacts)
}

#' @rdname impacts.SEMCMC.xxx
#' @param W An adjacency matrix, same as used in the call to SEMCMC().
#' @aliases impacts.SEMCMC.sdm
#' @examples
#' data(columbus, package = "spdep")
#'
#' W <- spdep::nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' #SDM model
#' sdm.mcmc <- SEMCMC(m.form, data = columbus, W = W, model = "sdm")
#' impacts(sdm.mcmc, W)

impacts.SEMCMC.sdm <- function(obj, W, idx.var, var.names) {

  #No checks here as this is an internal function

  #Number of variables: vars. + lagged vars.
  n.var <- length(idx.var)/2

  #Av. tot. imp: beta / (1 - rho)
  rho.sim <- obj$rho
  totimp <- sapply(1:n.var, function(X) {
   (obj[, idx.var[X]] + obj[, idx.var[X + n.var]])/(1 - rho.sim)
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
    tr.aux1[, 1] * obj[, idx.var[X]] +
      tr.aux1[, 2] * obj[, idx.var[X + n.var]]
  })
  colnames(dirimp) <- var.names

  #Indirect impacts
  indimp <- totimp - dirimp
  colnames(indimp) <- var.names

  impacts <- list(direct = dirimp, indirect = indimp,
    total = totimp)
  class(impacts) <- "impacts.SEMCMC"

  return(impacts)
}


#' @rdname impacts.SEMCMC.xxx
#' @examples
#' data(columbus, package = "spdep")
#'
#' W <- spdep::nb2mat(col.gal.nb, style = "W")
#' m.form <-  CRIME ~ INC + HOVAL
#' #SDEM model
#' sdem.mcmc <- SEMCMC(m.form, data = columbus, W = W, model = "sdem")
#' impacts(sdem.mcmc, W)
#'

impacts.SEMCMC.sdem <- function(obj, W, idx.var, var.names) {

  #No checks here as this is an internal function

  #Number of variables: vars. + lagged vars.
  n.var <- length(idx.var)/2

  #Av. tot. imp: beta + gamma
  totimp <- obj[, idx.var[1:n.var]] + obj[, idx.var[n.var + 1:n.var]]
  colnames(totimp) <- var.names

  #Av. direct impact: tr((I - rho*W)^{-1} *beta/n
  dirimp <- obj[, idx.var[1:n.var]]
  colnames(dirimp) <- var.names

  #Indirect impacts
  indimp <- totimp - dirimp

  impacts <- list(direct = dirimp, indirect = indimp,
    total = totimp)
  class(impacts) <- "impacts.SEMCMC"

  return(impacts)
}

