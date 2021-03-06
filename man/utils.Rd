% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{utils}
\alias{utils}
\alias{stan2coda}
\title{Internal functions of package SEMCMC.}
\usage{
stan2coda(fit)
}
\arguments{
\item{fit}{Model fitted with SEMCMC().}
}
\value{
A 'mcmc.list' with the samples of the variables in the model.
}
\description{
Function to convert ouput from jags and stan to a suitable
  format to compute the impacts.
}
