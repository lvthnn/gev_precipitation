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
  real<lower = 0> tau;
  real<lower = 0, upper = 1> xi_0;
}

transformed parameters {
  real xi = xi_0 - 0.5;
}

model {
  tau ~ exponential(0.5);
  sigma ~ exponential(3);
  xi_0 ~ beta(4, 4);

  beta[1] ~ normal(0, 5);
  for (j in 2:m)
    beta[j] ~ normal(beta[j - 1], tau);

  y ~ gevh(beta, t_m, sigma, xi);
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n)
    log_lik[i] = gevh_lpdf([y[i]]' | beta, t_m, sigma, xi);
}
