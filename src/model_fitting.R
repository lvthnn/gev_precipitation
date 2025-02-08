source("src/utils.R")

set.seed(4922)

# Global variables
n_iter <- 8000
n_warmup <- 1000

# Load in the data
dt_maxprecip <- readr::read_csv(
  file = "data/dt_maxprecip.csv",
  show_col_types = FALSE
)

n <- nrow(dt_maxprecip)

model_data <- list(
  n = n,
  y = dt_maxprecip$max_precip,
  n_iter = n_iter,
  n_warmup = n_warmup
)

y <- model_data$y

# Set up parallel cores
options(mc.cores = parallel::detectCores())

# Null model
null_fit <- fit_model(
  file = "src/stan/nullmodel.stan",
  data = model_data,
  n_iter = n_iter,
  n_warmup = n_warmup
)

yrep_null <- rstan::extract(null_fit, pars = "y_rep")$y_rep


# Time trend model
timetrend_fit <- fit_model(
  file = "src/stan/timetrend.stan",
  data = model_data,
  n_iter = n_iter,
  n_warmup = n_warmup
)

yrep_timetrend <- rstan::extract(timetrend_fit, pars = "y_rep")$y_rep

# Hierarchical model -- i.i.d. prior for random effects

ms <- c(4, 6, 12, 18, 24, 33)

fits_hierarchical <- pbapply::pblapply(ms, function(m) {
  t_m <- floor(seq(1, nrow(dt_maxprecip), length = m + 1))[-c(m + 1)]

  hmodel_data <- list(
    n = nrow(dt_maxprecip),
    m = m,
    y = dt_maxprecip$max_precip,
    t_m = t_m
  )

  iid_fit <- fit_model(
    file = "src/stan/hierarchical_iid.stan",
    data = hmodel_data,
    n_iter = n_iter,
    n_warmup = n_warmup
  )

  rw_fit <- fit_model(
    file = "src/stan/hierarchical_rw.stan",
    data = hmodel_data,
    n_iter = n_iter,
    n_warmup = n_warmup
  )

  ar1_fit <- fit_model(
    file = "src/stan/hierarchical_ar1.stan",
    data = hmodel_data,
    n_iter = n_iter,
    n_warmup = n_warmup
  )

  yrep_iid <- rstan::extract(iid_fit, pars = "y_rep")$y_rep
  yrep_rw <- rstan::extract(rw_fit, pars = "y_rep")$y_rep
  yrep_ar1 <- rstan::extract(ar1_fit, pars = "y_rep")$y_rep


  return(list(iid_fit, rw_fit, ar1_fit, yrep_iid, yrep_rw, yrep_ar1))
})

save(
  null_fit, yrep_null, timetrend_fit, yrep_timetrend,
  fits_hierarchical,
  file = paste0("fits_", Sys.Date() |> format("%d%m"), ".RData")
)
