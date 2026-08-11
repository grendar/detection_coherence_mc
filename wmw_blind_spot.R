#
#
# Aug 11 2026 Tu -- blind spot of the WMW test
#
# marian.grendar@gmail.com
#


################################################################################
#
#  WMW Test 
#  Monte Carlo study of power against blind spot
#
#  Blind spot: N(0,1) vs N(0,sigma^2) for any sigma != 1
#  AUC = P(X1 < X2) = 0.5 by symmetry (both centered at 0)
#
################################################################################



################################################################################
#
# libs
#
library(ggplot2)



################################################################################
#
# fncs
#
#
###############################################################
#
# p-value of WMW test 
wmw_pval <- function(x, y) {
  #
  wilcox.test(x, y, exact = FALSE)$p.value
  #
}


###############################################################
#
# Pratt's d3 factor for normal populations 
# (Pratt 1964, Sec 2.3 & pp. 4,10)
#
pratt_d3_normal <- function(theta, lambda = 0.5) {
  #
  # theta = sigma_y / sigma_x
  # alpha = P(min(X1,X2) > Y1), beta = P(X1 > max(Y1,Y2))
  # Formulas use radian arcsin scaled by 1/(2*pi), equivalent to Pratt's
  # "degrees / 360" statement.
  #  
  alpha <- 0.25 + (1 / (2 * pi)) * asin(theta^2 / (theta^2 + 1))
  beta  <- 0.25 + (1 / (2 * pi)) * asin(1 / (theta^2 + 1))   # sub 1/theta
  
  d3 <- (12 * lambda * (alpha - 0.25) + 12 * (1 - lambda) * (beta - 0.25))^(-0.5)
  return(d3)
  #
}

# Asymptotic two-tailed rejection probability from Pratt's d factor
pratt_power <- function(d3, alpha_level = 0.05) {
  #
  K <- qnorm(1 - alpha_level / 2)          # ~1.96 for alpha=0.05
  power <- 2 * (1 - pnorm(K * d3))
  return(power)
  #
}


###############################################################
#
# MC to check van der Vaart-Pratt formula for limit of power 
# for N(0,1) vs N(0,sigma^2)
#
pratt_limit <- function(sigma, n, alpha = 0.05, nsim = 20000) {
  #
  # Simulate large-n rejection rate to get Pratt's c(P*)
  pv <- replicate(nsim, wmw_pval(rnorm(n), rnorm(n, 0, sigma)))
  mean(pv < alpha, na.rm = TRUE)
  #
}


################################################################################
#
# MC
#
################################################################################
#
# Parameters 
#
sigma_bs  <- 10.0   # blind spot: N(0,1) vs N(0,sigma^2)
mu_near   <- 0.5    # near-blind-spot: N(mu,1) vs N(0,1)
alpha     <- 0.05
N_sim     <- 4000

###############################################################
#
# Checking van der Vaart-Pratt limit 
# 
# sigma_y/sigma_x = 10/1 = 10
theta = 10
d3 = pratt_d3_normal(theta)
power_pratt = pratt_power(d3)
power_pratt
# 0.09462627

# MC check
#
set.seed(123L)
cat("Computing theoretical Pratt limit at n=2000...\n")
c_pratt <- pratt_limit(sigma_bs, n = 2000, nsim = 10000)
cat(sprintf("  N(0,1) vs N(0,%.0f): c(P*) ~= %.3f\n\n",
            sigma_bs^2, c_pratt))
#
# N(0,1) vs N(0,sigma^2 = 100): c(P*) ~= 0.095
# => ok


###############################################################
#
# MC of power for growing n
#
set.seed(123L)
n_vals <- c(5, 25, 50, 100, 200, 500, 1000, 2000, 5000)
#
pow <- data.frame(
  n      = n_vals,
  power  = NA_real_,
  lower = NA_real_,
  upper = NA_real_)
#
#
for (j in seq_along(n_vals)) {
  n <- n_vals[j]
  pv_b <- replicate(N_sim,
                    wmw_pval(rnorm(n), rnorm(n, 0, sigma_bs)))
  power <- mean(pv_b < alpha, na.rm = TRUE)
  #
  se    <- sqrt(power * (1 - power) / N_sim)  # MC uncertainty
  pow$lower[j] <- power - 1.96 * se
  pow$upper[j] <- power + 1.96 * se
  #
  pow$power[j] = power
  #
}
#
round(pow, 3)
# n power lower upper
# 1    5 0.064 0.057 0.072
# 2   25 0.086 0.077 0.094
# 3   50 0.093 0.084 0.102
# 4  100 0.098 0.089 0.108
# 5  200 0.100 0.091 0.109
# 6  500 0.099 0.090 0.108
# 7 1000 0.097 0.088 0.106
# 8 2000 0.099 0.090 0.109
# 9 5000 0.094 0.085 0.103
#



###############################################################
#
# Plot of power as fnc of sample size
#
col_blind = 'red'
#
p <- ggplot(pow, aes(x = n, y = power)) +
  geom_hline(yintercept = power_pratt, linetype = "dashed",
             linewidth = 0.5, alpha = 0.6) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = n_vals,
                labels = as.character(n_vals)) +
  scale_y_continuous(limits = c(0, 0.2),
                     breaks = seq(0, 0.2, 0.1)) +
  annotate("text", x = 4500, y = power_pratt - 0.01,
           label = sprintf("van der Vaart - Pratt limit %.3f", power_pratt),
           colour = col_blind, size = 3.0, hjust = 1, alpha = 0.8) +
  labs(x      = expression(italic(n)~"per group (log scale)"),
       y      = "Power",
       colour = NULL, shape = NULL,
       title  = "WMW, equal groups") +
  theme_classic(base_size = 11) +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 8.5),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10, face = "plain"))

#
p
#


################################################################################
#
#
sessioninfo::session_info(pkgs = c('attached'))
#
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
# 
# [1] /home/mg/R/x86_64-pc-linux-gnu-library/4.0
# [2] /opt/R/4.0.5/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────
#





