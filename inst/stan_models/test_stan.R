devtools::install("../..")
library(SEjags)

data(columbus)
     
W <- nb2mat(col.gal.nb, style = "W")
m.form <-  CRIME ~ INC + HOVAL
     
#Fit models with SEjags
#sem.mcmc <- SEjags(m.form, data = columbus, W = W, model = "sem", sampler ="stan")

#slm.mcmc <- SEjags(m.form, data = columbus, W = W, model = "slm", sampler ="stan")

#sdm.mcmc <- SEjags(m.form, data = columbus, W = W, model = "sdm", sampler ="stan")

#sac.mcmc <- SEjags(m.form, data = columbus, W = W, model = "sac", sampler = "stan")

#sacmixed.mcmc <- SEjags(m.form, data = columbus, W = W, model = "sacmixed",
#  sampler = "stan")
 

car.mcmc <- SEjags(m.form, data = columbus, W = W, model = "car",
  sampler = "stan")
