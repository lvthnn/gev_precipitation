#include gev.stan

data {
  int<lower = 0> n;
  vector[n] y;
}

parameters {
  real mu;
  real<lower = 0> sigma;
  real<lower = 0, upper = 1> xi_0;
}

transformed parameters {
   real xi = xi_0 - 0.5;
}

model {
  mu ~ normal(0, 5);
  sigma ~ exponential(0.5);
  xi_0 ~ beta(4, 4);
  
  y ~ gev(mu, sigma, xi);
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = gev_lpdf([y[i]]' | mu, sigma, xi);
  }
}
