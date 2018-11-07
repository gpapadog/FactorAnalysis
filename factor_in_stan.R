# Fitting the factor model using Stan.

rm(list = ls())
library(rstan)
setwd('~/Github/FactorAnalysis/')

set.seed(1234)

sample_size <- 100
p <- 10
true_lambda <- matrix(rnorm(2 * p), ncol = 2)
true_sigma0 <- diag(MCMCpack::rinvgamma(p, shape = 300, scale = 100))
true_sigma <- true_lambda %*% t(true_lambda) + true_sigma0

Y <- t(mvnfast::rmvn(sample_size, mu = rep(0, p), sigma = true_sigma))


nu <- 3  # nu parameter in the prior for the phis.
prior_a1 <- c(1, 2)  # alpha, beta parameters in the Gamma prior for a1.
prior_a2 <- c(1, 2)
prior_sigma <- c(0.001, 0.001)
tune_a1 <- 1
tune_a2 <- 4
K <- 5


stan_data <- list(n = sample_size, p = p, y = Y, nu = nu, prior_a1 = prior_a1,
                  prior_a2 = prior_a2, prior_sigma = prior_sigma, K = K)
stan_fit <- stan("stan_code.stan", data = stan_data, iter = 3000, warmup = 1000,
                 thin = 3, chains = 2)

all_pars <- extract(stan_fit)

par(mfrow = c(4, 4), mar = rep(2, 4))
for (ll in 1 : p) {
  plot(density(all_pars$sigma_sq[, ll]), main = ll)
  abline(v = true_sigma0[ll, ll], col = 'red')
}


keep_Nsims <- nrow(all_pars$sigma_sq)

variances <- array(NA, dim = c(p, p, keep_Nsims))
for (ss in 1 : keep_Nsims) {
  variances[, , ss] <- all_pars$Lambda[ss, , ] %*% t(all_pars$Lambda[ss, , ])
  variances[, , ss] <- variances[, , ss] + diag(all_pars$sigma_sq[ss, ])
}


hist(apply(variances, c(1, 2), mean) - true_sigma)
