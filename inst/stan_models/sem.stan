data {
  int<lower = 1> N; //number of areas
  int<lower = 1> nvar; //Number of covariates

  vector[N] y;  //Response
  matrix[N, nvar] X; //Covariates
  
  matrix<lower = 0, upper = 1>[N, N] W; //Adjacency matrix

  real lambda_min;
  real lambda_max;
}

parameters {
  vector[nvar] b; //Coefficients

  real <lower = 0> tau;  //Precision

  real<lower = -1, upper = 1> lambda; //Spatial autocorrelation
}

transformed parameters {
  matrix[N, N] PREC;
  matrix[N, N] IlambdaW;

  IlambdaW = diag_matrix(rep_vector(1.0, N)) - lambda * W;
  PREC = tau * ((IlambdaW') * IlambdaW);

}

model {

  y ~ multi_normal_prec(X * b, PREC);

  b ~ normal(0, sqrt(1000));

  lambda ~ uniform(lambda_min, lambda_max);

  tau ~ gamma(0.01, 0.01);
}

