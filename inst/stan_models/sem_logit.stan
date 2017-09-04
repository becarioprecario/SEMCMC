data {
  int<lower = 1> N; //number of areas
  int<lower = 1> nvar; //Number of covariates

  int<lower=0, upper=1> y[N];  //Response (binary)

  matrix[N, nvar] X; //Covariates
  
  matrix<lower = 0, upper = 1>[N, N] W; //Adjacency matrix

  real lambda_min;
  real lambda_max;
}

parameters {
  vector[nvar] b; //Coefficients

  //real <lower = 0> tau;  //Precision

  real<lower = -1, upper = 1> lambda; //Spatial autocorrelation

  vector[N] ylatent;  //Latent variable (real)
}

transformed parameters {
  matrix[N, N] PREC;
  matrix[N, N] IlambdaW;

  real<lower = 0, upper = 1> prob[N];

  IlambdaW = diag_matrix(rep_vector(1.0, N)) - lambda * W;
  PREC =  ((IlambdaW') * IlambdaW);


  for(i in 1:N) {
    prob[i] = exp(ylatent[i])/(1 + exp(ylatent[i]));
  }

}

model {

  y ~ bernoulli(prob);

  ylatent ~ multi_normal_prec(X * b, PREC);

  b ~ normal(0, sqrt(1000));

  lambda ~ uniform(lambda_min, lambda_max);

  //tau ~ gamma(0.01, 0.01);
}

