#
#
# Aug 11 2026 Tu -- Zaremba's example of non-growing power of WMW test 
#                   against a blind spot alternative
#
# marian.grendar@gmail.com
#

################################################################################
#
# Zaremba's example (simplification of Gnedenko's example) -- blind spot of WMW
#
# F1: dirac at 0
# F2: 1/2 mass at -1, 1/2 mass at +1
#
################################################################################
#


################################################################################
#
# --- Sampling functions --------------------------------------
rF1 <- function(n) rep(0, n)
rF2 <- function(n) sample(c(-1, 1), n, replace = TRUE)

# --- Theoretical Var(U) at the blind spot (n1 = n2 = n) ------
varU_true <- function(n) n^2 * n / 4

# --- Tie-corrected null variance, expected tie pattern --------
# t1 = n (all X = 0), t2 = t3 = n/2 (expected split of Y into +1/-1)
Var0_tiecorr <- function(n) {
  N <- 2 * n
  tie_sum <- (n^3 - n) + 2 * ((n / 2)^3 - n / 2)
  (n * n / 12) * ((N + 1) - tie_sum / (N * (N - 1)))
}

Var0_naive <- function(n) n * n * (2 * n + 1) / 12

d2_tiecorr_theory <- function(n) varU_true(n) / Var0_tiecorr(n)
d2_naive_theory   <- function(n) 3 * n / (2 * n + 1)

# --- Direct WMW U-statistic and z (tie-corrected), for MC check ----
wmw_z_tiecorrected <- function(x, y) {
  n1 <- length(x); n2 <- length(y); N <- n1 + n2
  r  <- rank(c(x, y))
  U  <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  EU <- n1 * n2 / 2
  tab <- table(c(x, y))
  tie_term <- sum(tab^3 - tab)
  V0 <- (n1 * n2 / 12) * ((N + 1) - tie_term / (N * (N - 1)))
  (U - EU) / sqrt(V0)
}


################################################################################
#
set.seed(123L)
N_sim <- 3000
alpha <- 0.05

# --- Verify AUC = 1/2 ------------------------------------------------
auc_check <- mean(replicate(5000, mean(outer(rF1(50), rF2(50), ">"))))
cat(sprintf("AUC check: %.4f (should be 0.500)\n", auc_check))
# AUC check: 0.4998 (should be 0.500)

################################################################################
#
# --- Verify d^2 -> 16/9, MC on the tie-corrected statistic -----------
cat("\nVariance ratio check (equal samples, tie-corrected):\n")
cat(sprintf("%-6s %-14s %-14s %-14s\n", "n", "Var(z_tc) MC", "d^2 tie-corr", "d^2 naive"))
cat(strrep("-", 50), "\n")
for (n in c(25, 100, 500, 2000)) {
  zs <- replicate(N_sim, wmw_z_tiecorrected(rF1(n), rF2(n)))
  cat(sprintf("%-6d %-14.4f %-14.4f %-14.4f\n",
              n, var(zs), d2_tiecorr_theory(n), d2_naive_theory(n)))
}
cat(sprintf("\nAsymptotic limits: tie-corrected -> %.4f (16/9), naive -> %.4f (3/2)\n",
            16 / 9, 3 / 2))

################################################################################
#
# --- Power via wilcox.test() (already tie-corrected internally) ------
cat("\nPower MC (equal samples m = n):\n")
cat(sprintf("%-6s %-10s %-14s\n", "n", "power", "predicted (16/9)"))
cat(strrep("-", 32), "\n")
n_vals <- c(25, 50, 100, 200, 500, 1000, 2000, 5000)
pred_power <- 2 * (1 - pnorm(1.96 / sqrt(16 / 9)))
for (n in n_vals) {
  pv <- replicate(N_sim, suppressWarnings(wilcox.test(rF1(n), rF2(n))$p.value))
  cat(sprintf("%-6d %-10.4f %-14.4f\n", n, mean(pv < alpha), pred_power))
}

##########################################################
#
# Variance ratio check (equal samples, tie-corrected):
# n      Var(z_tc) MC   d^2 tie-corr   d^2 naive 
# 25     1.8161         1.7422         1.4706        
# 100    1.8687         1.7689         1.4925        
# 500    1.7953         1.7760         1.4985        
# 2000   1.7577         1.7773         1.4996
# Asymptotic limits: tie-corrected -> 1.7778 (16/9), naive -> 1.5000 (3/2)
#
#

##########################################################
#
# Power MC (equal samples m = n):
# n      power      predicted (16/9)  
# 25     0.1167     0.1416        
# 50     0.1233     0.1416        
# 100    0.1280     0.1416        
# 200    0.1403     0.1416        
# 500    0.1440     0.1416        
# 1000   0.1363     0.1416        
# 2000   0.1450     0.1416
# 5000   0.1390     0.1416
#
#

#─ Session info ───────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 (2024-04-24)
# os       CentOS Stream 8
# system   x86_64, linux-gnu
# ui       RStudio
# language
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/Bratislava
# date     2026-08-11
# rstudio  2024.12.1+563 Kousa Dogwood (desktop)
# pandoc   2.0.6 @ /usr/bin/pandoc
# quarto   1.5.57 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/quarto
#



  
