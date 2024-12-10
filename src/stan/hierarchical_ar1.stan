#include gev.stan

data {
  int<lower = 0> n;
  int<lower = 0> m;
  vector[n] y;
  vector[m] t_m;
}

parameters {
  vector[m] beta;
  real<lower = 0> sigma;
  real<lower = 0, upper = 1> xi_0;
  real<lower = -1, upper = 1> phi;
  real mu_0;
}

transformed parameters {
  real xi = xi_0 - 0.5;
}

model {
  sigma ~ exponential(3);
  phi ~ uniform(-1, 1);
  mu_0 ~ normal(0, 10);
  xi_0 ~ beta(4, 4);

  beta[1] ~ normal(0, 5);
  for (j in 2:m)
    beta[j] ~ normal(phi * beta[j - 1], 5);

  y ~ gevh(mu_0, beta, t_m, sigma, xi);
}

generated quantities {
  real log_lik;
  log_lik = gevh_lpdf(y | mu_0, beta, t_m, sigma, xi);
}

