#include gev.stan

data {
  int<lower = 0> n;
  vector[n] y;
}

parameters {
  real mu_0;
  real<lower = 0> sigma;
  real<lower = 0, upper = 1> xi_0;
  real delta;
}

transformed parameters {
  real xi = xi_0 - 0.5;
}

model {
  mu_0 ~ normal(0, 10);
  xi_0 ~ beta(5, 5);
  delta ~ normal(0, 10);
  sigma ~ exponential(3);

  y ~ gevt(mu_0, sigma, xi, delta);
}

generated quantities {
  real log_lik;
  log_lik = gevt_lpdf(y | mu_0, sigma, xi, delta);
}

