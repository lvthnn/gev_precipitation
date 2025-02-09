fit_model <- function(file, data, chains, n_iter, n_warmup) {
  model <- rstan::stan_model(file = file)

  model_samples <- rstan::sampling(
    model,
    data = data,
    iter = n_iter,
    warmup = n_warmup,
    init_r = 0.025
  )

  return(model_samples)
}

gev_mean_null <- function(fit) {
  fit_sum <- rstan::summary(fit)$summary

  gev_mean <- do.call(cbind, lapply(c(1, 4, 8), function(col) {
    mu <- fit_sum["mu", col]
    sigma <- fit_sum["sigma", col]
    xi <- fit_sum["xi", col]

    res <- mu + sigma * (gamma(1 - xi) - 1) / xi

    return(res)
  })) |> data.frame()

  colnames(gev_mean) <- c("null_mean", "null_ci_l", "null_ci_u")

  return(gev_mean)
}

gev_mean_lin <- function(fit) {
  samples <- rstan::extract(fit, pars = c("mu_0", "sigma", "xi", "delta"))
  mu_0 <- samples$mu_0
  sigma <- samples$sigma
  xi <- samples$xi
  delta <- samples$delta

  gev_mean_samples <- sapply(1:132, function(t) {
    mu_t <- mu_0 * (1 + delta * (t - 66))
    return(mu_t + sigma * (gamma(1 - xi) - 1) / xi)
  })

  gev_mean <- apply(gev_mean_samples, 2, mean)
  ci_l <- apply(gev_mean_samples, 2, quantile, probs = 0.025)
  ci_u <- apply(gev_mean_samples, 2, quantile, probs = 0.975)

  return(data.frame(lin_mean = gev_mean, lin_ci_l = ci_l, lin_ci_u = ci_u))
}

gev_mean_hier <- function(fit) {
  samples <- rstan::extract(fit, pars = c("mu_0", "sigma", "xi", "beta"))
  mu_0 <- samples$mu_0
  sigma <- samples$sigma
  xi <- samples$xi
  beta <- samples$beta

  m <- ncol(beta)
  t_m <- floor(seq(1, 132, length = m + 1))[-c(m + 1)]

  gev_mean_samples <- sapply(1:132, function(t) {
    mu_t <- mu_0

    for (j in seq_along(t_m)) {
      mu_t <- mu_t + ifelse(t > t_m[j], beta[, j] * (t - t_m[j]), 0)
    }

    return(mu_t + sigma * (gamma(1 - xi) - 1) / xi)
  })

  gev_mean <- apply(gev_mean_samples, 2, mean)
  ci_l <- apply(gev_mean_samples, 2, quantile, probs = 0.025)
  ci_u <- apply(gev_mean_samples, 2, quantile, probs = 0.975)

  return(data.frame(hier_mean = gev_mean, hier_ci_l = ci_l, hier_ci_u = ci_u))
}

ic_hierarchical <- function(ms, fits) {
  # m
  models <- c("I.I.D.", "R.W.", "AR(1)")
  do.call(rbind, lapply(1:6, function(mind) {
    # prior
    ics_mind <- lapply(1:3, function(mod) {
      log_lik <- loo::extract_log_lik(fits[[mind]][[mod]])
      data.frame(
        LOOIC = suppressWarnings(loo::loo(log_lik)$estimates[3, 1]),
        WAIC = suppressWarnings(loo::waic(log_lik)$estimates[3, 1]),
        type = models[mod],
        m = ms[mind]
      )
    })

    do.call(rbind, ics_mind) |>
      tibble::as_tibble()
  }))
}
