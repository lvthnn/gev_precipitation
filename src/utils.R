fit_model <- function(file, data, chains) {
  model <- rstan::stan_model(file = file)

  model_samples <- rstan::sampling(
    model,
    data = data,
    iter = n_iter,
    warmup = warmup,
    init_r = 0.025
  )

  return(model_samples)
}

model_diag <- function(model_samples) {
  loglik <- loo::extract_log_lik(model_samples)

  suppressWarnings({
    waic <- loo::waic(loglik)
    loo <- loo::loo(loglik)
    psis <- loo::psis(-loglik)
  })

  return(list(waic, loo, psis))
}

get_tex_table <- function(model_samples, pars) {
  post_sum <- rstan::summary(
    model_samples,
    regex_pars = pars,
    probs = c(0.025, 0.5, 0.975)
  )

  knitr::kable(post_sum$summary, digits = 4, format = "latex")
}
