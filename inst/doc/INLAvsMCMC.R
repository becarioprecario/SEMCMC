## ----results = "hide"----------------------------------------------------
library(SEjags)
library(INLA)


## ----results = "hide"----------------------------------------------------
data(columbus)
d <- columbus
W <- nb2mat(col.gal.nb, style = "W")
m.form <-  CRIME ~ INC + HOVAL
     
#Fit models with SEjags
if(!file.exists("INLAvsMCMC-MCMC.Rdata") ) {
 sem.mcmc <- SEjags(m.form, data = d, W = W, model = "sem",
   n.burnin = 5000, n.iter = 10000, n.thin = 20, linear.predictor = TRUE)
 slm.mcmc <- SEjags(m.form, data = d, W = W, model = "slm",
   n.burnin = 5000, n.iter = 10000, n.thin = 20, linear.predictor = TRUE)
 sdm.mcmc <- SEjags(m.form, data = d, W = W, model = "sdm",
   n.burnin = 5000, n.iter = 10000, n.thin = 20, linear.predictor = TRUE)

 save(file = "INLAvsMCMC-MCMC.Rdata",
   list = c("sem.mcmc", "slm.mcmc", "sdm.mcmc"))
} else {
  load("INLAvsMCMC-MCMC.Rdata")
}

## ----eval = FALSE, echo = FALSE------------------------------------------
#  #sdem.mcmc <- SEjags(m.form, data = d, W = W, model = "sdem")
#  #slx.mcmc <- SEjags(m.form, data = d, W = W, model = "slx")
#  #sac.mcmc <- SEjags(m.form, data = d, W = W, model = "sac")
#  #sacmixed.mcmc <- SEjags(m.form, data = d, W = W, model = "sacmixed")
#  #car.mcmc <- SEjags(m.form, data = d, W = W, model = "car")
#  
#  #Compute impacts
#  #impacts(sem.mcmc, W)
#  #impacts(slm.mcmc, W)
#  #impacts(sdm.mcmc, W)
#  #impacts(sdem.mcmc, W)
#  #impacts(slx.mcmc, W)
#  #impacts(sac.mcmc, W)
#  #impacts(sacmixed.mcmc, W)
#  #impacts(car.mcmc, W)

## ----results = "hide"----------------------------------------------------
#Area index
columbus$idx <- 1:nrow(columbus)
#Adjacency matrix as sparse matrix
W.inla <- as(W, "CsparseMatrix")

#Model matrix for SLM models
mmatrix <- model.matrix(m.form, columbus)
mmatrix2 <-  cBind(mmatrix, W.inla %*% mmatrix[,-1])
colnames(mmatrix2)[4:5]<- paste("lag", colnames(mmatrix2)[2:3],sep="")

## ----results = "hide"----------------------------------------------------
#Zero-variance for error term
#Zero-variance to remove effect in linear predictor: DOES NOT WORK
zero.variance = list(prec=list(initial = 25, fixed=TRUE))
#Large variance to allow for an uncnstrained estimation of the other parameters
#zero.variance = list(prec=list(initial = 1/100, fixed=TRUE))

## ----results = "hide"----------------------------------------------------
#Compute eigenvalues for SLM model (as in Havard's code)
e = eigen(W.inla)$values
re.idx = which(abs(Im(e)) < 1e-6)
rho.max = 1/max(Re(e[re.idx]))
rho.min = 1/min(Re(e[re.idx]))
rho = mean(c(rho.min, rho.max))

## ----results = "hide"----------------------------------------------------
#
#Variance-covarinace matrix for beta coeffients' prior
#
betaprec <- .001
#Standard regression model
Q.beta = Diagonal(n = ncol(mmatrix), x = 1)
Q.beta = betaprec * Q.beta
#Regression model with lagged covariates
Q.beta2 = Diagonal(n = ncol(mmatrix2), x = 1)
Q.beta2 = betaprec * Q.beta2

## ----results = "hide"----------------------------------------------------
#Arguments for slm latent model
args.slm = list(
   rho.min = rho.min,
   rho.max = rho.max,
   W = W.inla,#as(W.inla,"dgTMatrix"),
   X = matrix(0, nrow(mmatrix),0),
   Q.beta = matrix(1,0,0)
)

#Priors of Hyperparameters
hyper.slm = list(
   prec = list(prior = "loggamma", param = c(0.01, 0.01)),
      rho = list(initial = 0, prior = "logitbeta", param = c(1,1))
)

## ----results = "hide"----------------------------------------------------


#Control fixed
c.fixed <- list(prec = 0.001, prec.intercept = 0.001)

## ----results = "hide"----------------------------------------------------

#SEM model
hyper.sem <- hyper.slm
hyper.sem$rho$initial <- 0.85 #Fixed to posterior mode from MCMC
hyper.sem$rho$fixed <- TRUE

#Change zero variance
zero.variance$prec$initial <- 1/100

#Control inla
c.inla <- list(strategy = "laplace", fast = FALSE,
  tolerance = 0.001,
     int.strategy = "ccd", h = 0.001, dz = 0.05, stencil = 9) 
#  int.strategy = 'grid', diff.logdens = 0.1, h = 0.001, dz = 0.01, stencil = 9)

# Create linear combinations on the covariates to estimate
# linear predictor (and fitted values).

n <- nrow(columbus)

#Test how the structure of the linear combinatios should be
lc1 <- inla.make.lincomb(list("(Intercept)" = 1, 
  INC = columbus$INC[1], HOVAL = columbus$HOVAL[1], 
  idx = c(1, rep(NA, n-1))
))

lc.linpred <- lapply(1:n, function(X) {
  idx.lc <- rep(NA, 49)
  idx.lc[X] <- 1
  aux <- as.list(mmatrix[X, ])
  aux$idx = idx.lc
  inla.make.lincomb(aux)
})


lc.linpred <- do.call(c, lc.linpred)
names(lc.linpred) <- paste("lc.linpred", 1:n, sep = "")

#Linear combination of the fixed effects
lc.fixed <- inla.make.lincombs(INLA:::inla.uncbind(mmatrix))
names(lc.fixed) <- paste("lc.fixed", 1:n, sep = "")


#Linear combinations for SLM and SDM
lc.linpred2 <- lapply(1:n, function(X) {
  idx.lc <- rep(NA, 49)
  idx.lc[X] <- 1
  inla.make.lincomb(list(idx = idx.lc))
})
lc.linpred2 <- do.call(c, lc.linpred2)
names(lc.linpred2) <- paste("lc", 1:n, sep = "")

#SEM model
sem.inla<-inla(CRIME ~ INC + HOVAL +
   f(idx, model = "slm", args.slm = args.slm, hyper = hyper.slm),
   data = as.data.frame(columbus), family = "gaussian",
   lincomb = c(lc.linpred, lc.fixed), control.predictor = list(compute=TRUE),
   control.fixed = c.fixed,
   control.inla = c.inla,
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE)
)


#SLM model
slm.inla<-inla( CRIME ~ -1 +
   f(idx, model="slm",
      args.slm=list(rho.min = rho.min, rho.max = rho.max, W = W.inla, X=mmatrix,
         Q.beta = Q.beta),
      hyper=hyper.slm),
   data=as.data.frame(columbus), family="gaussian",
   lincomb = lc.linpred2, control.predictor = list(compute=TRUE),
   control.fixed = c.fixed,
   control.inla = c.inla,
   control.family = list(hyper=zero.variance),
   control.compute=list(dic=TRUE, cpo=TRUE)
)

#SDM model
hyper.sdm <- hyper.slm

#Fix rho
hyper.sdm$rho$initial <- 0.77 #Fixed to posterior mode from MCMC
#hyper.sdm$rho$fixed <- TRUE

#Use stronger prior
hyper.sdm$rho$param <- c(140, 60)


sdm.inla <- inla( CRIME ~ -1 +
   f(idx, model = "slm",
      args.slm = list(rho.min = rho.min, rho.max = rho.max, W = W.inla, 
        X = mmatrix2, Q.beta = Q.beta2),
      hyper = hyper.sdm),
   data = as.data.frame(columbus), family = "gaussian",
   lincomb = lc.linpred2, control.predictor = list(compute=TRUE),
   control.fixed = c.fixed,
   control.inla = c.inla,
   control.family = list(hyper = zero.variance),
   control.compute = list(dic = TRUE, cpo = TRUE)
)




## ------------------------------------------------------------------------
  ff <- function(z){z * (rho.max - rho.min) + rho.min}
  semmarg <- inla.tmarginal(ff, sem.inla$marginals.hyperpar[[2]])
  slmmarg <- inla.tmarginal(ff, slm.inla$marginals.hyperpar[[2]])
  sdmmarg <- inla.tmarginal(ff, sdm.inla$marginals.hyperpar[[2]])

## ----fig = TRUE, echo = FALSE--------------------------------------------
#Display marginals
par(mfrow = c(2,3))
var.names <- c("Intercept", "INC", "HOVAL", "lagINC", "lagHOVAL")
#Fixed effects
for(v in 1:3) {
  plot(density(sem.mcmc$b[v,,,]), main = var.names[v])
  lines(sem.inla$marginals.fixed[[v]], col = "red")
}

#Precision
plot(density(sem.mcmc$tau[1,,]), main = "Precision")
lines(sem.inla$marginals.hyper[[1]], col = "red")


#Spatial autocorrelation
plot(density(sem.mcmc$lambda[1,,]), main = "Spat. Autocor.")
lines(semmarg, col = "red")

## ----fig = TRUE, echo = FALSE, width = 5, height = 5---------------------
# Display first 9 linear predictors
par(mfrow = c(3, 3))
for(v in 1:9) {
  plot(density(sem.mcmc$mu[v,,,]), main = paste0("Fixed effect, obs. ", v))
  lines(sem.inla$marginals.lincomb.derived[[n + v]], col = "red")
}

## ------------------------------------------------------------------------
sem.raneff.mcmc <- apply(sem.mcmc$mu, 2, function(X) {columbus$CRIME - X})
sem.raneff.mcmc <- matrix(sem.raneff.mcmc, ncol = 49, byrow = TRUE)

sem.raneff.mean <- apply(sem.raneff.mcmc, 2, mean)
sem.raneff.sd <- apply(sem.raneff.mcmc, 2, sd)

## ----fig = TRUE, echo = FALSE--------------------------------------------
par(mfrow = c(1, 2))
plot(sem.raneff.mean, sem.inla$summary.random$idx[, "mean"], xlab = "MCMC",
  ylab = "INLA", main = "Random effects (post. mean)")
  abline(0, 1)
plot(sem.raneff.sd, sem.inla$summary.random$idx[, "sd"], xlab = "MCMC",
  ylab = "INLA", main = "Random effects (post. s.d.)")
  abline(0, 1)

## ----fig = TRUE, echo = FALSE--------------------------------------------
# Display first 9 linear predictors
par(mfrow = c(3, 3))
for(v in 1:9) {
  plot(density(sem.raneff.mcmc[, v]), 
    main = paste0("Random effect, obs. ", v))
  lines(sem.inla$marginals.random$idx[[v]], col = "red")
}

## ----fig = TRUE, echo = FALSE--------------------------------------------
plot(sem.inla$marginals.random$idx[[1]], type = "l")
lines(density(sem.raneff.mcmc[, 1]), lty = 2)

legend ("topleft", lty = 1:2, legend = c("INLA", "MCMC"))


## ----fig = TRUE, echo = TRUE---------------------------------------------
#As reported by INLA
sem.inla$summary.random$idx[1,]
#As computed using the posterior marginal
inla.zmarginal(sem.inla$marginals.random$idx[[1]], FALSE)
#Using the MCMC output
mean(sem.raneff.mcmc[, 1])
sd(sem.raneff.mcmc[, 1])
quantile(sem.raneff.mcmc[, 1], c(0.025, 0.25, 0.5, 0.75, 0.975))

## ----fig = TRUE, echo = FALSE--------------------------------------------
plot(sem.raneff.mcmc[, 1], type = "l") 

## ---- eval = FALSE, echo = FALSE-----------------------------------------
#  #Remove large points
#  idx.x1 <- sem.raneff.mcmc[, 1] < 30
#  #Lets fit a spline
#  dens.x1 <- density(sem.raneff.mcmc[idx.x1, 1])
#  ff.x1 <- splinefun(dens.x1$x, dens.x1$y)
#  
#  integrate(ff.x1, -30, 30)
#  
#  #Summary statistics
#  inla.zmarginal(cbind(dens.x1$x, dens.x1$y))
#  

## ----fig = TRUE, echo = FALSE--------------------------------------------
#Display marginals
par(mfrow = c(2,3))
#Fixed effects
for(v in 1:3) {
  plot(density(slm.mcmc$b[v,,,]), main = var.names[v])
  lines(slm.inla$marginals.random$idx[[nrow(columbus) + v]], col = "red")
}

#Precision
plot(density(slm.mcmc$tau[1,,]), main = "Precision")
lines(slm.inla$marginals.hyper[[1]], col = "red")


#Spatial autocorrelation
plot(density(slm.mcmc$rho[1,,]), main = "Spat. Autocor.")
lines(slmmarg, col = "red")

## ------------------------------------------------------------------------
slm.raneff.mcmc <- apply(slm.mcmc$mu, 2, function(X) {columbus$CRIME - X})
slm.raneff.mcmc <- matrix(slm.raneff.mcmc, ncol = 49, byrow = TRUE)

slm.raneff.mean <- apply(slm.raneff.mcmc, 2, mean)
slm.raneff.sd <- apply(slm.raneff.mcmc, 2, sd)

## ----fig = TRUE, echo = FALSE--------------------------------------------
par(mfrow = c(1, 2))
plot(apply(slm.mcmc$mu, 1, mean), slm.inla$summary.random$idx[1:n, "mean"], 
  xlab = "MCMC",
  ylab = "INLA", main = "Random effects (post. mean)")
  abline(0, 1)
plot(apply(slm.mcmc$mu, 1, sd), slm.inla$summary.random$idx[1:n, "sd"],
  xlab = "MCMC",
  ylab = "INLA", main = "Random effects (post. s.d.)")
  abline(0, 1)

## ----fig = TRUE, echo = FALSE, eval = FALSE------------------------------
#  # Display first 9 linear predictors
#  par(mfrow = c(3, 3))
#  for(v in 1:9) {
#    plot(density(slm.raneff.mcmc[, v]),
#      main = paste0("Random effect, obs. ", v))
#    lines(slm.inla$marginals.random$idx[[v]], col = "red")
#  }

## ----fig = TRUE, echo = FALSE--------------------------------------------
#Display marginals
par(mfrow = c(3,3))
#Fixed effects
for(v in 1:5) {
  plot(density(sdm.mcmc$b[v,,,]), main = var.names[v])
  lines(sdm.inla$marginals.random$idx[[nrow(columbus) + v]], col = "red")
}

#Precision
plot(density(sdm.mcmc$tau[1,,]), main = "Precision")
lines(sdm.inla$marginals.hyper[[1]], col = "red")


#Spatial autocorrelation
plot(density(sdm.mcmc$rho[1,,]), main = "Spat. Autocor.")
lines(sdmmarg, col = "red")

