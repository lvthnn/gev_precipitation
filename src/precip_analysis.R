# Global variables
n_iter <- 13000
warmup <- 3000

seed <- 31415926

# Load in the data
dt_maxprecip <- readr::read_csv(
  file = "data/dt_maxprecip.csv",
  show_col_types = FALSE
)

model_data <- list(
  n = nrow(dt_maxprecip),
  y = dt_maxprecip$max_precip
)

# Set up parallel cores
options(mc.cores = parallel::detectCores())

fit_model <- function(file, data, chains) {
  model <- rstan::stan_model(file = file)

  model_samples <- rstan::sampling(
    model,
    data = data,
    iter = n_iter,
    warmup = warmup,
    init_r = 0.025,
    seed = seed
  )

  return(model_samples)
}

model_diag <- function(model_samples) {
  loglik <- loo::extract_log_lik(model_samples)

  suppressWarnings({
    waic <- loo::waic(loglik)
    loo <- loo::loo(loglik)
    psis <- loo::psis(loglik)
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

# Null model
null_fit <- fit_model(
  file = "src/stan/nullmodel.stan",
  data = model_data
)

# Time trend model
timetrend_fit <- fit_model(
  file = "src/stan/timetrend.stan",
  data = model_data
)

# Hierarchical model -- i.i.d. prior for random effects

m <- 6
t_m <- floor(seq(1, nrow(dt_maxprecip), length = m + 2))[-c(1, m + 2)]

hmodel_data <- list(
  n = nrow(dt_maxprecip),
  m = m,
  y = dt_maxprecip$max_precip,
  t_m = floor(seq(1, nrow(dt_maxprecip), length = m + 2))[-c(1, m + 2)]
)

iid_fit <- fit_model(
  file = "src/stan/hierarchical_iid.stan",
  data = hmodel_data
)

rw_fit <- fit_model(
  file = "src/stan/hierarchical_rw.stan",
  data = hmodel_data
)

ar1_fit <- fit_model(
  file = "src/stan/hierarchical_ar1.stan",
  data = hmodel_data
)

save(
  null_fit, timetrend_fit, iid_fit, rw_fit, ar1_fit,
  file = "fits_2310.RData"
)

# Evaluation
# - PSIS og LOOIC
# - WAIC
# - QQplot (ath. smá handavinna fyrir non-stationary dreifingar)
# - Trend línur
# - Return level? [hentar ekki vel fyrir non-stationary]
#   - Virkar vel fyrir stationary tímaröð
#   - Hvaða stærð fær maður á flóði fyrir gefinn endurkomutíma?
#     — x. t — 1/(1 - q) — q quantile af breytunni [length of return period]
#     — y. quantile of the random variable [return level]

#
