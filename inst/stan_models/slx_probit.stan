data {
  int<lower = 1> N; //number of areas
  int<lower = 1> nvar; //Number of covariates

  int<lower=0, upper=1> y[N];  //Response (binary)

  matrix[N, nvar] X; //Covariates
}

parameters {
  vector[nvar] b; //Coefficients

  //real <lower = 0> tau;  //Precision

  vector[N] ylatent;  //Latent variable (real)
}

transformed parameters {
  real sigma;
  real<lower = 0, upper = 1> prob[N];

  //sigma = sqrt(1/tau);

  for(i in 1:N) {
    prob[i] = Phi(ylatent[i]);
  }

}

model {
  y ~ bernoulli(prob);

  ylatent ~ normal(X * b, 1);

  b ~ normal(0, sqrt(1000));

  //tau ~ gamma(0.01, 0.01);
}

