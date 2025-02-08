library(tidyverse)
library(bayesplot)

theme_set(bayesplot::theme_default())
bayesplot::color_scheme_set(scheme = "blue")

set.seed(1552)

load("data/environ/fits_1701.RData")

rw_fit_optim <- fits_hierarchical[[4]][[2]]

m <- 18
t_m <- floor(seq(1, 132, length = m + 1))[-c(m + 1)]

beta <- rstan::summary(rw_fit_optim)$summary[1:18, 1]

beta_lines <- lapply(1:m, function(j) {
  mu_j <- sapply(1:132, function(t) {
    beta[j] * ifelse(t >= t_m[j], 1, 0) * (t - t_m[j])
  })
})

dt_mu_j <- do.call(cbind, beta_lines) |>
  as_tibble() |>
  mutate(t = 1:132) |>
  select(t, everything())

colnames(dt_mu_j)[2:(m + 1)] <- paste0("beta_", 1:m)

dt_mu_j_long <- dt_mu_j |>
  pivot_longer(cols = starts_with("beta"))

mu_j_plot <- ggplot(dt_mu_j_long, aes(x = t, y = value, colour = name)) +
  geom_vline(xintercept = c(t_m + 1, 132), linetype = "dotted", alpha = 0.25) +
  geom_line() +
  xlim(1, 132) +
  labs(
    x = "Time",
    y = expression(mu[j](t))
  ) +
  theme(legend.position = "none")

mu_0 <- rstan::summary(rw_fit_optim)$summary["mu_0", 1]

mu_t <- sapply(1:132, function(t) {
  mu_t <- mu_0
  for (j in seq_along(t_m)) {
    mu_t <- mu_t + ifelse(t >= t_m[j], beta[j] * (t - t_m[j]), 0)
  }
  return(mu_t)
}) |>
  as_tibble() |>
  mutate(
    t = 1:132,
    t_m = cut(t, c(0, t_m, 132))
  ) |>
  select(t, t_m, value)

mu_t_plot <- ggplot(mu_t, aes(x = t, y = value)) +
  geom_vline(xintercept = c(t_m + 1, 132), linetype = "dotted", alpha = 0.25) +
  geom_line() +
  geom_line(aes(colour = t_m)) +
  theme(legend.position = "none") +
  labs(
    x = "Time",
    y = expression(mu[t])
  ) +
  xlim(1, 132)

ggsave(
  cowplot::plot_grid(mu_j_plot, mu_t_plot),
  width = 7,
  height = 3.5,
  device = "pdf",
  file = "pdf/img/beta_basis.pdf"
)
