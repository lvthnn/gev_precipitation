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
  sigma ~ exponential(3);
  xi_0 ~ beta(4, 4);
  delta ~ normal(0, 10);

  y ~ gevt(mu_0, sigma, xi, delta);
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n)
    log_lik[i] = gevt_lpdf([y[i]]' | mu_0, sigma, xi, delta);
}
