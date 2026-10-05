#
#
# Oct 5 2026 Mo -- blind spot of the logrank test, given allocation gamma
#
# marian.grendar@gmail.com
#

################################################################################
#
#  Logrank test
# 
# Monte Carlo power of the logrank test at three allocations, γ ∈ {0.3, 0.5, 0.7}. 
# Within each panel the two piecewise exponential hazard pairs
# share the rates a = 0.5, b = 2.0 and differ only in the breakpoint: the blind spot
# pair c = c ∗ (γ), for which θ(P ; γ) = 0, and a perturbed pair c = 1.1 c ∗ (γ), for
# which θ(P ; γ) ≠ 0. Both have S 1 ≠ S 2 with crossing hazards. Power against the
# blind spot pair stays at a fixed level as n grows, a level that itself depends on the
# allocation; power against the perturbed pair grows towards 1 at every allocation.
#
################################################################################


################################################################################
#
# libs
#
library(survival)
library(ggplot2)
library(patchwork)

################################################################################
#
# fncs
#
###############################################################
#

make_S <- function(a, b, c) {
  #
  S1 <- function(t) ifelse(t <= c, exp(-a * t), exp(-a * c - b * (t - c)))
  S2 <- function(t) ifelse(t <= c, exp(-b * t), exp(-b * c - a * (t - c)))
  list(S1 = S1, S2 = S2)
  #
}

theta <- function(c, a, b, gamma) {
  #
  S <- make_S(a, b, c); S1 <- S$S1; S2 <- S$S2
  h <- function(t, gamma) {
    y1 <- gamma * S1(t); y2 <- (1 - gamma) * S2(t); y <- y1 + y2
    ifelse(y <= 0, 0, y1 * y2 / y)
  }
  I1 <- integrate(h, 0, c, gamma = gamma, subdivisions = 200L)$value
  I2 <- integrate(h, c, 80, gamma = gamma, subdivisions = 200L)$value
  (a - b) * (I1 - I2)
  #
}

find_cstar <- function(a, b, gamma)
  uniroot(theta, interval = c(1e-4, 8.0), a = a, b = b, gamma = gamma)$root

sample_piecewise_exp <- function(n, lam1, lam2, c) {
  #
  u <- runif(n); Sc <- exp(-lam1 * c); t <- numeric(n)
  mask <- u <= 1 - Sc
  t[mask]  <- -log1p(-u[mask]) / lam1
  t[!mask] <- c - log((1 - u[!mask]) / Sc) / lam2
  t
  #
}

one_rep <- function(n1, n2, a, b, cc) {
  #
  T1 <- sample_piecewise_exp(n1, a, b, cc)
  T2 <- sample_piecewise_exp(n2, b, a, cc)
  times <- c(T1, T2)
  g1 <- c(rep.int(1L, n1), rep.int(0L, n2))
  g1s <- g1[order(times)]
  is1 <- g1s; is2 <- 1L - g1s
  Y1 <- n1 - c(0, cumsum(is1)[-length(is1)])
  Y2 <- n2 - c(0, cumsum(is2)[-length(is2)])
  Y  <- Y1 + Y2
  O1 <- sum(is1)
  E1 <- sum((Y1 / Y) * (is1 + is2))
  Vp <- sum(Y1 * Y2 / Y^2)
  (O1 - E1) / sqrt(Vp)
  #
}

################################################################################
#
# design
#
################################################################################
#
a <- 0.5
b <- 2.0
gammas <- c(0.3, 0.5, 0.7)

# arms are multiples of c*(gamma), so the perturbation scales with c*
# (c* ranges from ~0.23 to ~0.90 across hazard pairs)
rels <- c(1.00, 1.10, 1.25, 1.60)
arm_lab <- c("c = c* (blind spot)", "c = 1.10 c*", "c = 1.25 c*", "c = 1.60 c*")

ns <- c(250L, 500L, 1000L, 2000L, 4000L, 8000L)

reps <- 2000L          # raise for smoother curves; runtime is ~linear in reps
zcrit <- 1.959963985

############################################################
#
set.seed(123L)
#
oo <- list()
k <- 0L
for (gamma in gammas) {
  #
  cstar <- find_cstar(a, b, gamma)
  #
  for (j in seq_along(rels)) {
    #
    cc <- rels[j] * cstar
    th <- theta(cc, a, b, gamma)
    for (n_total in ns) {
      n1 <- as.integer(round(n_total * gamma)); n2 <- n_total - n1
      z <- numeric(reps)
      for (r in seq_len(reps)) z[r] <- one_rep(n1, n2, a, b, cc)
      rate <- mean(abs(z) > zcrit)
      k <- k + 1L
      oo[[k]] <- data.frame(
        gamma = gamma, arm = arm_lab[j], rel = rels[j],
        cstar = cstar, c = cc, theta = th, n = n_total,
        rate = rate, se = sqrt(rate * (1 - rate) / reps),
        var_z = var(z), mean_z = mean(z)
      )
    #  
    }
    # 
  }
  #
}
#
dd <- do.call(rbind, oo)
dd$arm <- factor(dd$arm, levels = arm_lab)

################################################################################
#
# figure
#
################################################################################
#
# power vs n: flat at c*, growing off it
#
p1 <- ggplot(dd[dd$arm %in% c('c = c* (blind spot)', 'c = 1.10 c*'),],
             aes(n, rate, colour = arm, shape = arm)) +
  geom_hline(yintercept = 0.05, linetype = 2) +
  geom_line() + 
  geom_point(size = 2.5) +
  facet_wrap(~ gamma, labeller = label_bquote(gamma == .(gamma))) +
  scale_x_log10() + ylim(0, 1) +
  labs(x = "Total sample size n (log scale)", y = "Power",
       colour = NULL) +
  theme_classic(base_size = 12) +
  scale_colour_manual(name = NULL, values = c("firebrick3", "steelblue4")) +
  scale_shape_manual(name = NULL, values = c(16, 17)) +
  theme(legend.position = "bottom") + 
  ggtitle('Logrank test')
p1
# save it
# gg_save(p1, 'logrank_blind_spot_gamma.tiff', 20, 10, 1)


################################################################################
#
#
sessioninfo::session_info(pkgs = c('attached'))
# ─ Session info ──────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.0.5 (2021-03-31)
# os       CentOS Stream 8
# system   x86_64, linux-gnu
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/Bratislava
# date     2026-10-05
# rstudio  2024.12.1+563 Kousa Dogwood (desktop)
# pandoc   2.0.6 @ /usr/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────
# package   * version date (UTC) lib source
# ggplot2   * 3.5.1   2024-04-23 [1] CRAN (R 4.0.5)
# patchwork * 1.2.0   2024-01-08 [1] CRAN (R 4.0.5)
# survival  * 3.7-1   2024-07-18 [1] Github (therneau/survival@b79a343)
# 
# [1] /home/mg/R/x86_64-pc-linux-gnu-library/4.0
# [2] /opt/R/4.0.5/lib64/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────
