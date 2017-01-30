#
#Compare how models are fitted with INLA versus a multivariate Normal
#
#Test implemented models using a small dataset and jags
library(spdep)
library(rjags)

#Load data
data(columbus)

#Fit models using spdep
lw <- nb2listw(col.gal.nb, style = "W")

#Model formula
m.form <-  CRIME ~ INC + HOVAL
#SEM
m.sem <- errorsarlm(m.form, data = columbus, listw = lw)

#Adjacency matrix
W <- nb2mat(col.gal.nb, style="W")
W.eigen <- eigen(W)$values

#Data for jags
d.jags <- list (y = columbus$CRIME)
d.jags$X <- model.matrix (CRIME ~ INC +  HOVAL, columbus)
#d.jags$X <- cbind(d.jags$X, W%*%d.jags$X[, -1])
d.jags$W <- W
d.jags$N <- nrow(columbus)
d.jags$n.var <- ncol(d.jags$X)
d.jags$I <- diag (d.jags$N)
#Bound for rho
d.jags$rho.min <- 1/min(W.eigen)
d.jags$rho.max <- 1/max(W.eigen)
#Pprec for response, should be very small
d.jags$prec <- .45

#Inits for jags
d.inits <- list(b = rep(0, d.jags$n.var), tau = 1, rho = 0)
#d.inits$prec <- 1

#Run jags
jm1 <- jags.model('semINLA.bug', data = d.jags, inits = d.inits, n.chains = 1,
  n.adapt = 100)

update(jm1, 1000)

jm1.samp <- jags.samples(jm1, c('b', 'rho', 'prec', 'tau'), 1000)

#MCMC
print(jm1.samp)
#Max. lik.
coef(m.sem)
m.sem$s2 


#NO CONCLUSIVE results!!
