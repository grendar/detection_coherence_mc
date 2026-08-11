#
#
# Aug 11 2026 Tu -- blind spot of the logrank test
#
# marian.grendar@gmail.com
#


################################################################################
#
#  Logrank test
#  Monte Carlo study of power against blind spot
#
#  Demonstrates: power bounded at alpha for all n at the
#  blind spot distribution (theta(P) = 0, S1 != S2)
#  vs power growing with n at a non-blind-spot crossing
#  configuration.
#
################################################################################


################################################################################
#
# libs
#
library(survival)
library(ggplot2)


################################################################################
#
# fncs
#
###############################################################
#
# Piecewise exponential random number generator
#    Group 1: hazard = a for t <= c, b for t > c
#    Group 2: hazard = b for t <= c, a for t > c
#
rpiecexp <- function(n, a, b, c) {
  #
  u        <- runif(n)
  s_c      <- exp(-a * c)           # S1(c) = exp(-a*c)
  before   <- u > s_c               # event before crossing
  t        <- numeric(n)
  t[before]  <- -log(u[before]) / a
  t[!before] <- c + log(s_c / u[!before]) / b
  t
  #
}

###############################################################
#
# Logrank detection functional (numerical integration)
#    theta(c) = integral of h(t)*{lambda1-lambda2} dt
#    h(t) = S1(t)*S2(t) / (2*(S1(t)+S2(t)))  [no-censoring limit]
# 
logrank_theta <- function(c, a = 0.5, b = 2.0,
                          t_max = 30, n_grid = 50000) {
  t  <- seq(0, t_max, length.out = n_grid)
  dt <- t[2] - t[1]
  
  S1 <- ifelse(t <= c,
               exp(-a * t),
               exp(-a * c - b * (t - c)))
  S2 <- ifelse(t <= c,
               exp(-b * t),
               exp(-b * c - a * (t - c)))
  
  h           <- S1 * S2 / (2 * (S1 + S2))
  lambda_diff <- ifelse(t <= c, a - b, b - a)
  
  sum(h * lambda_diff) * dt
}

###############################################################
#
# Power simulation function
#
sim_power_logrank <- function(n_per_group, a, b, c,
                              n_sim = 2000, alpha = 0.05,
                              seed = 42) {
  #
  set.seed(seed)
  reject <- logical(n_sim)
  for (i in seq_len(n_sim)) {
    t1    <- rpiecexp(n_per_group, a, b, c)
    t2    <- rpiecexp(n_per_group, b, a, c)  # groups swapped
    time  <- c(t1, t2)
    event <- rep(1L, 2L * n_per_group)       # no censoring
    grp   <- factor(rep(1:2, each = n_per_group))
    sfit  <- survdiff(Surv(time, event) ~ grp)
    reject[i] <- pchisq(sfit$chisq, df = 1,
                        lower.tail = FALSE) < alpha
  }
  mean(reject)
  #
}



################################################################################
#
# MC
#
################################################################################
#
#
###############################################################
#
# Find c* numerically (blind spot crossing point)
#  theta(c*) = 0
#
a <- 0.5 
b <- 2.0 
#
c_grid    <- seq(0.01, 3, by = 0.001)
theta_grid <- sapply(c_grid, logrank_theta, a = a, b = b)
#
# Locate sign change
sign_changes <- which(diff(sign(theta_grid)) != 0)
c_star <- c_grid[sign_changes[1]] +
  0.001 * abs(theta_grid[sign_changes[1]]) /
  (abs(theta_grid[sign_changes[1]]) +
     abs(theta_grid[sign_changes[1] + 1]))
#
cat(sprintf("Blind spot crossing point: c* = %.4f\n", c_star))
cat(sprintf("theta(c*) = %.2e\n", logrank_theta(c_star)))
#
# Blind spot crossing point: c* = 0.5182
# theta(c*) = 3.40e-05
#

# Choose a non-blind-spot crossing config:
# same a, b but c shifted -- theta(P) != 0
c_alt <- 0.48 # c_star * 0.5   # earlier crossing -> theta(P) != 0
cat(sprintf("Non-blind-spot crossing:   c_alt = %.4f\n", c_alt))
cat(sprintf("theta(c_alt) = %.4f\n", logrank_theta(c_alt)))
#
# Non-blind-spot crossing:   c_alt = 0.4800
# theta(c_alt) = 0.0138


###############################################################
#
# Survival curve plot: blind spot distributions
#
t_plot <- seq(0, 6, length.out = 500)
S1_blind <- ifelse(t_plot <= c_star,
                   exp(-a * t_plot),
                   exp(-a * c_star - b * (t_plot - c_star)))
S2_blind <- ifelse(t_plot <= c_star,
                   exp(-b * t_plot),
                   exp(-b * c_star - a * (t_plot - c_star)))

surv_df <- data.frame(
  t    = rep(t_plot, 2),
  S    = c(S1_blind, S2_blind),
  grp  = rep(c("S1 (a=0.5 then b=2.0)",
               "S2 (b=2.0 then a=0.5)"), each = 500)
)

p_surv <- ggplot(surv_df, aes(x = t, y = S, colour = grp)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c_star, linetype = "dashed",
             colour = "grey40") +
  annotate("text", x = c_star + 0.1, y = 0.9,
           label = sprintf("c* = %.3f", c_star),
           hjust = 0, size = 3.5, colour = "grey30") +
  scale_colour_manual(values = c("firebrick3", "steelblue4")) +
  labs(
    title    = "Blind Spot Distributions: Crossing Survival Functions",
    subtitle = sprintf(
      "S1 \u2260 S2 but logrank detection functional \u03b8(P) = 0  (a=%.1f, b=%.1f)",
      a, b),
    x      = "Time",
    y      = "Survival probability",
    colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())
#
p_surv



###############################################################
#
# MC of power for growing sample size
#
#
n_vals <- c(25, 50, 100, 200, 500, 1000, 2000, 5000) 
#
#
# power against blind spot
power_blind <- sapply(n_vals, function(n) {
  cat(sprintf("  Blind spot:     n = %5d ... ", n))
  p <- sim_power_logrank(n, a, b, c_star)
  cat(sprintf("power = %.3f\n", p))
  p
})
# Blind spot:     n =    25 ... power = 0.064
# Blind spot:     n =    50 ... power = 0.052
# Blind spot:     n =   100 ... power = 0.050
# Blind spot:     n =   200 ... power = 0.040
# Blind spot:     n =   500 ... power = 0.052
# Blind spot:     n =  1000 ... power = 0.042
# Blind spot:     n =  2000 ... power = 0.046
# Blind spot:     n =  5000 ... power = 0.046

# power against non-blind spot
power_alt <- sapply(n_vals, function(n) {
  cat(sprintf("  Non-blind-spot: n = %5d ... ", n))
  p <- sim_power_logrank(n, a, b, c_alt)
  cat(sprintf("power = %.3f\n", p))
  p
})
# Non-blind-spot: n =    25 ... power = 0.055
# Non-blind-spot: n =    50 ... power = 0.051
# Non-blind-spot: n =   100 ... power = 0.076
# Non-blind-spot: n =   200 ... power = 0.081
# Non-blind-spot: n =   500 ... power = 0.140
# Non-blind-spot: n =  1000 ... power = 0.247
# Non-blind-spot: n =  2000 ... power = 0.459
# Non-blind-spot: n =  5000 ... power = 0.859


###############################################################
#
# Results table
#
#
results <- data.frame(
  n             = rep(n_vals, 2),
  power         = c(power_blind, power_alt),
  configuration = rep(c(
    sprintf("Blind spot"),
    sprintf("Non-blind-spot")
  ), each = length(n_vals))
)
#


###############################################################
#
# Plot of power as fnc of sample size
#
p <- ggplot(results, aes(x = n, y = power,
                         shape = configuration,
                         colour = configuration,
                         linetype = configuration)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_shape_manual(name = ' ', values = c(16, 17)) + 
  scale_x_log10(
    breaks = n_vals,
    labels = n_vals
  ) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.1)) +
  scale_colour_manual(name = ' ', values = c("firebrick3", "steelblue4")) +
  scale_linetype_manual(name = ' ', values = c("solid", "solid")) +
  labs(
    title    = "Logrank test", 
    x        = "n per group (log scale)",
    y        = "Power",
    colour   = NULL,
    linetype = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "bottom",  # 'right', #
    legend.direction = 'horizontal', #"vertical",
    plot.subtitle    = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
#
p
# gg_save(p, 'logrank_blind_spot.tiff', 20, 10, 1.2)

################################################################################
#
#
sessioninfo::session_info(pkgs = c('attached'))
# ─ Session info ───────────────────────────────────────────────────────────
# setting  value
# version  R version 4.0.5 (2021-03-31)
# os       CentOS Stream 8
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# date     2026-08-11
# rstudio  2024.12.1+563 Kousa Dogwood (desktop)
# pandoc   2.0.6 @ /usr/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────
# package   * version date (UTC) lib source
# ggplot2   * 3.5.1   2024-04-23 [1] CRAN (R 4.0.5)
# password  * 1.0-0   2016-03-22 [1] CRAN (R 4.0.5)
# patchwork * 1.2.0   2024-01-08 [1] CRAN (R 4.0.5)
# survival  * 3.7-1   2024-07-18 [1] Github (therneau/survival@b79a343)
# 
# [1] /home/mg/R/x86_64-pc-linux-gnu-library/4.0
# [2] /opt/R/4.0.5/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────






