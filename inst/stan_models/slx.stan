data {
  int<lower = 1> N; //number of areas
  int<lower = 1> nvar; //Number of covariates

  vector[N] y;  //Response
  matrix[N, nvar] X; //Covariates
}

parameters {
  vector[nvar] b; //Coefficients

  real <lower = 0> tau;  //Precision
}

transformed parameters {
  real sigma;

  sigma = sqrt(1/tau);

}

model {

  y ~ normal(X * b, sigma);

  b ~ normal(0, sqrt(1000));

  tau ~ gamma(0.01, 0.01);
}

