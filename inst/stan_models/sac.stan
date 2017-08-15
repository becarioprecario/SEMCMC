data {
  int<lower = 1> N; //number of areas
  int<lower = 1> nvar; //Number of covariates

  vector[N] y;  //Response
  matrix[N, nvar] X; //Covariates
  
  matrix<lower = 0, upper = 1>[N, N] W; //Adjacency matrix

  real rho_max;
  real rho_min;
  real lambda_max;
  real lambda_min;

}

parameters {
  vector[nvar] b; //Coefficients

  real <lower = 0> tau;  //Precision

  real<lower = -1, upper = 1> rho; //Spatial autocorrelation LINEAR PREDICTOR
  real<lower = -1, upper = 1> lambda; //Spatial autocorrelation ERROR TERM
}

transformed parameters {
  matrix[N, N] PREC;
  matrix[N, N] IrhoW;
  matrix[N, N] IlambdaW;
  matrix[N, N] IrlW;
  vector[N] mu;

  IrhoW = diag_matrix(rep_vector(1.0, N)) - rho * W;
  IlambdaW = diag_matrix(rep_vector(1.0, N)) - lambda * W;
  IrlW = IlambdaW * IrhoW;
  PREC = tau * (IrlW' * IrlW);

  mu = IrhoW\ (X * b);
}

model {

  y ~ multi_normal_prec(mu, PREC);

  b ~ normal(0, 1000);

  rho ~ uniform(rho_min, rho_max);

  lambda ~ uniform(lambda_min, lambda_max);

  tau ~ gamma(0.01, 0.01);
}
