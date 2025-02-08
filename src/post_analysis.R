# PACKAGES AND IMPORTS
library(bayesplot)
library(ggplot2)
library(tidyverse)

source("src/utils.R")

load("data/environ/fits_1701.RData")

theme_set(bayesplot::theme_default())
bayesplot::color_scheme_set(scheme = "blue")



# LOAD DATA AND CREATE PLOT
dt_maxprecip <- readr::read_csv(
  file = "data/dt_maxprecip.csv",
  show_col_types = FALSE
)

y <- dt_maxprecip$max_precip

go_trendplot <- ggplot(dt_maxprecip, aes(x = 1:132, y = y)) +
  geom_line(linetype = "dashed", colour = "darkblue") +
  geom_point(alpha = 0.5, colour = "darkblue", size = 2) +
  labs(
    x = "Time",
    y = "Max. annual precipitation"
  )

ggsave(
  go_trendplot,
  width = 7,
  height = 3.5,
  device = "pdf",
  file = "pdf/img/trendplot.pdf",
)



# OPTIMISATION OVER m
ms <- c(4, 6, 12, 18, 24, 33)

dt_optim_m <- ic_hierarchical(ms, fits_hierarchical) |>
  tidyr::pivot_longer(cols = c("LOOIC", "WAIC"))

loo_lin <- loo::loo(timetrend_fit)$estimates[3, 1]
waic_lin <- loo::waic(timetrend_fit |> loo::extract_log_lik())$estimates[3, 1]

col_scheme <- color_scheme_get("brewer-Spectral")
ic_cols <- c("LOOIC" = col_scheme[[6]], "WAIC" = col_scheme[[2]])

go_optim_m <- ggplot(dt_optim_m, aes(x = m, y = value, colour = name)) +
  geom_line(linetype = "dashed") +
  geom_point(size = 2) +
  facet_wrap(~type) +
  geom_hline(
    yintercept = loo_lin,
    linetype = "dotted",
    colour = col_scheme[6]
  ) +
  geom_hline(
    yintercept = waic_lin,
    linetype = "dotted",
    colour = col_scheme[2]
  ) +
  labs(
    x = "*m*",
    y = "Estimate",
    colour = "Information criterion"
  ) +
  theme(
    legend.position = "bottom",
    axis.title.x = ggtext::element_markdown()
  ) +
  scale_colour_manual(values = ic_cols)

ggsave(
  go_optim_m,
  width = 7,
  height = 3.5,
  device = "pdf",
  file = "pdf/img/optim_m.pdf"
)

optim_m <- dt_optim_m |>
  arrange(value) |>
  head(1) |>
  select(m, type) |>
  mutate(
    m = which(ms == m),
    type = which(c("I.I.D.", "R.W.", "AR(1)") == type)
  ) |>
  rename(
    mind = m,
    mod = type
  )

rw_fit <- fits_hierarchical[[optim_m$mind]][[optim_m$mod]]
yrep_rw <- fits_hierarchical[[optim_m$mind]][[optim_m$mod + 3]]
m <- ms[optim_m$mind]


# POSTERIOR MEAN OF mu IN MODELS
post_mean_hier <- gev_mean_hier(rw_fit)
post_mean_lin <- gev_mean_lin(timetrend_fit)
post_mean_null <- gev_mean_null(null_fit)

dt_maxprecip_post <- dt_maxprecip |>
  cbind(
    post_mean_hier,
    post_mean_lin,
    post_mean_null
  ) |>
  pivot_longer(
    cols = matches("^(hier|lin|null)_"),
    names_to = c("model", ".value"),
    names_pattern = "^(hier|lin|null)_(.+)$"
  ) |>
  mutate(
    year = sort(rep(1:132, 3)),
    model = fct_recode(
      model,
      "Non-stationary (non-linear)" = "hier",
      "Non-stationary (linear)" = "lin",
      "Stationary" = "null",
    )
  ) |>
  rename(time = year)

go_post_mean <- ggplot(
  dt_maxprecip_post,
  aes(x = time, y = max_precip)
) +
  geom_point(colour = "darkblue", alpha = 0.5, size = 1.25) +
  geom_line(colour = "darkblue", alpha = 0.5, linetype = "dashed") +
  geom_line(aes(x = time, y = mean), linewidth = 0.9) +
  geom_line(aes(x = time, y = ci_l), linetype = "dashed") +
  geom_line(aes(x = time, y = ci_u), linetype = "dashed") +
  labs(x = "Time", y = "Max. annual precipitation") +
  facet_wrap(~model)

ggsave(
  go_post_mean,
  width = 8,
  height = 1 / 3 * 8,
  device = "pdf",
  file = "pdf/img/post_mean.pdf"
)


# NULL MODEL
go_trace_null <- mcmc_trace(null_fit, pars = c("mu", "sigma", "xi")) +
  theme(legend.position = "none")

go_acf_null <- mcmc_acf(null_fit, pars = c("mu", "sigma", "xi"))

ggsave(
  go_trace_null,
  width = 8,
  height = 2,
  device = "pdf",
  file = "pdf/img/trace_null.pdf",
)

ggsave(
  go_acf_null,
  width = 8,
  height = 8,
  device = "pdf",
  file = "pdf/img/acf_null.pdf"
)

# TIMETREND MODEL
go_trace_timetrend <- mcmc_trace(
  timetrend_fit,
  pars = c("mu_0", "sigma", "xi", "delta")
) +
  theme(legend.position = "none")

go_acf_timetrend <- mcmc_acf(
  timetrend_fit,
  pars = c("mu_0", "sigma", "xi", "delta")
)

ggsave(
  go_trace_timetrend,
  width = 8,
  height = 4,
  device = "pdf",
  file = "pdf/img/trace_timetrend.pdf",
)

ggsave(
  go_acf_timetrend,
  width = 8,
  height = 8,
  device = "pdf",
  file = "pdf/img/acf_timetrend.pdf"
)

go_trace_rw <- mcmc_trace(
  rw_fit,
  pars = c(paste0("beta[", 1:m, "]"), "mu_0", "sigma", "xi"),
  facet_args = list(ncol = 4)
) +
  theme(legend.position = "none")

ggsave(
  go_trace_rw,
  width = 12,
  height = 10,
  device = "pdf",
  file = "pdf/img/trace_rw.pdf"
)

go_acf_rw_1 <- mcmc_acf(
  rw_fit,
  pars = c(paste0("beta[", 1:11, "]"))
) +
  labs(x = "", y = "")

go_acf_rw_2 <- mcmc_acf(
  rw_fit,
  pars = c(paste0("beta[", 12:m, "]"), "mu_0", "sigma", "xi")
) +
  labs(x = "", y = "")

ggsave(
  go_acf_rw_1,
  width = 12,
  height = 6,
  device = "pdf",
  file = "pdf/img/acf_rw_1.pdf",
)
ggsave(
  go_acf_rw_2,
  width = 12,
  height = 6,
  device = "pdf",
  file = "pdf/img/acf_rw_2.pdf"
)

# More analysis for our best model, RW

# QQPLOT

samples <- rstan::extract(rw_fit, pars = c("mu_0", "sigma", "xi", "beta"))
mu_0 <- samples$mu_0
sigma <- samples$sigma
xi <- samples$xi
beta <- samples$beta

t_m <- floor(seq(1, 132, length = m + 1))[-c(m + 1)]

mean_rw <- sapply(1:132, function(t) {
  mu_t <- mu_0

  for (j in seq_along(t_m)) {
    mu_t <- mu_t + ifelse(t > t_m[j], beta[, j] * (t - t_m[j]), 0)
  }

  return(mu_t + sigma * (gamma(1 - xi) - 1) / xi)
})

var_rw <- ifelse(
  xi == 0,
  sigma^2 * pi^2 / 6,
  sigma^2 * (gamma(1 - 2 * xi) - gamma(1 - xi)^2) / xi^2
)

sd_rw <- do.call(rbind, lapply(var_rw, function(s) rep(sqrt(s), 132)))
zrep_rw <- (yrep_rw - mean_rw) / sd_rw
z <- (y - mean(y)) / sd(y)

qs_zrep <- zrep_rw |> quantile(seq(0, 0.995, length.out = 125))
qs_z <- z |> quantile(seq(0, 0.995, length.out = 125))
qqs <- tibble(qs_zrep, qs_z)

go_qqplot_yrep <- ggplot(qqs, aes(x = qs_zrep, y = qs_z)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Standardised PPD quantiles",
    y = "Observed standardised quantiles"
  )

ggsave(
  go_qqplot_yrep,
  width = 8,
  height = 3.4,
  device = "pdf",
  file = "pdf/img/qqplot_yrep.pdf"
)

# PPC PLOTS
ppc_plots <- lapply(c("max", "min", "mean", "median"), function(fn_str) {
  fn <- get(fn_str)
  yrep_fn <- apply(zrep_rw, 1, fn)
  ecdf_fn <- ecdf(yrep_fn)
  p95_int <- quantile(yrep_fn, probs = c(0.025, 0.975)) |> unname()

  go_ppc_fn <- ppc_stat(z, zrep_rw, stat = fn_str) +
    theme(legend.position = "none") +
    geom_vline(
      xintercept = p95_int[1],
      colour = "darkblue",
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = p95_int[2],
      colour = "darkblue",
      linetype = "dashed"
    ) +
    labs(title = paste0(fn_str, " (*p* = ", ecdf_fn(fn(z)) |> round(2), ")")) +
    theme(title = ggtext::element_markdown())

  if (fn_str == "max") {
    go_ppc_fn <- go_ppc_fn + xlim(quantile(zrep_rw, probs = c(0, 0.99999)))
  }

  return(go_ppc_fn)
})

go_ppc_rw <- cowplot::plot_grid(plotlist = ppc_plots)

ggsave(
  go_ppc_rw,
  width = 8,
  height = 5,
  device = "pdf",
  file = "pdf/img/ppc_rw.pdf"
)



# DENSITY OVERLAY
go_epdf_overlay <- ppc_dens_overlay(z, zrep_rw[1:25, ]) +
  theme(legend.position = "none")
go_ecdf_overlay <- ppc_ecdf_overlay(z, zrep_rw[1:25, ]) +
  theme(legend.position = "none")

go_edf_overlay <- cowplot::plot_grid(go_epdf_overlay, go_ecdf_overlay)

ggsave(
  go_edf_overlay,
  width = 8,
  height = 4,
  device = "pdf",
  file = "pdf/img/ppc_edf_overlay.pdf"
)



# INDEX PLOT OF POSTERIOR MEANS OF \beta_j WITH 95% CREDIBLE SETS
beta_mean <- posterior::summarise_draws(rw_fit, mean)
beta_ci <- posterior::summarise_draws(
  rw_fit,
  ~ quantile(.x, probs = c(0.025, 0.975))
)

beta_tibble <- cbind(beta_mean, beta_ci)[1:13, -3] |>
  tibble() |>
  rename(cs_l = `2.5%`, cs_u = `97.5%`)

go_betas <- ggplot(
  beta_tibble,
  aes(x = 1:13, y = mean, ymin = cs_l, ymax = cs_u)
) +
  geom_point() +
  geom_errorbar() +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  labs(
    x = "Variable",
    y = "Posterior mean"
  )

ggsave(
  go_betas,
  width = 6,
  height = 2,
  device = "pdf",
  file = "pdf/img/betas.pdf"
)
