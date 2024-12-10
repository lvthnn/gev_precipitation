source("src/utils.R")

# Global variables
n_iter <- 13000
warmup <- 3000

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

m <- 33
t_m <- floor(seq(1, nrow(dt_maxprecip), length = m + 1))[-c(m + 1)]

hmodel_data <- list(
  n = nrow(dt_maxprecip),
  m = m,
  y = dt_maxprecip$max_precip,
  t_m = t_m
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

# =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
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
# =#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

gev_mean <- function(t_m, fit) {
  m <- length(t_m)
  betas <- rstan::summary(fit)$summary[1:m, 1]
  mu_0 <- rstan::summary(fit)$summary["mu_0", 1]
  sigma <- rstan::summary(fit)$summary["sigma", 1]
  xi <- rstan::summary(fit)$summary["xi", 1]

  mu_t <- c()
  m <- length(betas)

  for (t in 1:132) {
    mu <- mu_0
    for (j in 1:m) {
      if (t > t_m[j]) mu <- mu + betas[j] * (t - t_m[j])
    }

    mu_t[t] <- mu
  }

  gev_mean <- mu_t + sigma * (gamma(1 - xi) - 1) / xi

  return(gev_mean)
}

gm <- gev_mean(t_m, iid_fit)
plot(1:132, gm,
  type = "l", lwd = 2, col = "red", ylim = range(hmodel_data$y),
  xlab = "Year", ylab = expression(mu)
)
lines(1:132, hmodel_data$y)
