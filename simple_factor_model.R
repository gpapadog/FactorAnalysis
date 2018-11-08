# xx <- 1 / rgamma(10000, shape = 20, scale = 60)
# xxx <- MCMCpack::rinvgamma(100000, shape = 20, scale = 1 / 60)
# xxxx <- 1 / rgamma(10000, shape = 20, rate = 1 / 60)
# plot(density(xx))
# lines(density(xxx), col = 'red')
# lines(density(xxxx), col = 'green')



# ------- GENERATING DATA -------- #

set.seed(1234)

sample_size <- 500
p <- 6
true_lambda <- matrix(rnorm(2 * p), ncol = 2)
true_sigma0 <- diag(MCMCpack::rinvgamma(p, shape = 300, scale = 100))
true_sigma <- true_lambda %*% t(true_lambda) + true_sigma0

Y <- t(mvnfast::rmvn(sample_size, mu = rep(0, p), sigma = true_sigma))


# ---------- PRIOR SPECIFICATION --------- # 

nu <- 3  # nu parameter in the prior for the phis.
prior_a1 <- c(1, 2)  # alpha, beta parameters in the Gamma prior for a1.
prior_a2 <- c(1, 2)
prior_sigma <- c(0.001, 0.001)
tune_a1 <- 1
tune_a2 <- 4
K <- 5

Nsims <- 5000

# Where to save results.
etas <- array(NA, dim = c(Nsims, sample_size, K))
phis <- array(NA, dim = c(Nsims, p, K))
a1 <- rep(NA, Nsims)
a2 <- rep(NA, Nsims)
deltas <- array(NA, dim = c(Nsims, K))
taus <- array(NA, dim = c(Nsims, K))
lambdas <- array(NA, dim = c(Nsims, p, K))
sigmas <- array(NA, dim = c(Nsims, p))
acc_a1 <- 0
acc_a2 <- 0

# Initiate values.
etas[1, , ] <- rnorm(sample_size * K, mean = 0, sd = 1)
a1[1] <- rgamma(1, shape = prior_a1[1] * 100, rate = prior_a1[2] * 100)
a2[1] <- rgamma(1, shape = prior_a2[1] * 100, rate = prior_a2[2] * 100)
deltas[1, 1] <- rgamma(1, shape = a1[1], rate = 1)
deltas[1, - 1] <- rgamma(1, shape = a2[1], rate = 1)
taus[1, ] <- sapply(1 : K, function(x) prod(deltas[1, 1 : x]))
phis[1, , ] <- rgamma(p * K, shape = nu / 2, rate = nu / 2)
lambdas[1, , ] <- rnorm(p * K, mean = 0,
                        sd = sqrt(phis[1, , ] %*% diag(taus[1, ])))
sigmas[1, ] <- MCMCpack::rinvgamma(p, shape = prior_sigma[1] * 100,
                                   scale = prior_sigma[2] * 100)


# Trying fixing a1, a2
a1[1] <- 1
a2[1] <- 1

# ------- MCMC -------- #

for (ss in 2 : Nsims) {
  
  if (ss %% 1000 == 0) {
    print(paste('Simulation iteration', ss))
  }
  
  curr_etas <- etas[ss - 1, , ]
  curr_a1 <- a1[ss - 1]
  curr_a2 <- a2[ss - 1]
  curr_deltas <- deltas[ss - 1, ]
  curr_taus <- taus[ss - 1, ]
  curr_phis <- phis[ss - 1, , ]
  curr_lambdas <- lambdas[ss - 1, , ]
  curr_sigmas <- sigmas[ss - 1, ]

  
  # --------- Update the sigmas.
  
  curr_means <- curr_lambdas %*% t(curr_etas)
  a_new <- rep(prior_sigma[1] + sample_size / 2, p)
  b_new <- rep(prior_sigma[2], p) + 1 / 2 * apply((Y - curr_means) ^ 2, 1, sum)
  sigmas[ss, ] <- MCMCpack::rinvgamma(p, shape = a_new, scale = b_new)
  curr_sigmas <- sigmas[ss, ]
  
  
  # ---------- Update the lambdas.
  
  prior_variances <- curr_phis %*% diag(curr_taus)
  for (ll in 1 : p) {
    prior_prec <- diag(prior_variances[ll, ])
    data_prec <- t(curr_etas) %*% curr_etas / curr_sigmas[ll]
    post_var <- solve(prior_prec + data_prec)
    post_mean <- post_var %*% t(curr_etas) %*% matrix(Y[ll, ], ncol = 1)
    post_mean <- post_mean / curr_sigmas[ll]
    lambdas[ss, ll, ] <- mvnfast::rmvn(1, mu = post_mean, sigma = post_var)
  }
  curr_lambdas <- lambdas[ss, , ]
  
  
  # ---------- Update the etas.
  prior_prec <- diag(rep(1, K))
  data_prec <- t(curr_lambdas) %*% diag(1 / curr_sigmas) %*% curr_lambdas
  post_var <- solve(prior_prec + data_prec)
  post_mean <- post_var %*% t(curr_lambdas) %*% diag(1 / curr_sigmas)
  for (ii in 1 : sample_size) {
    post_mean_ii <- post_mean %*% matrix(Y[, ii], ncol = 1)
    etas[ss, ii, ] <- mvnfast::rmvn(1, mu = post_mean_ii, sigma = post_var)
  }
  curr_etas <- etas[ss, , ]

  
  # ---------- Update the phis.
  a_new <- matrix(nu + 1 / 2, nrow = p, ncol = K)
  b_new <- (nu + curr_lambdas ^ 2 %*% diag(curr_taus)) / 2
  phis[ss, , ] <- rgamma(p * K, shape = a_new, rate = b_new)
  curr_phis <- phis[ss, , ]
  
  
  # ---------- Update the deltas and taus.
  phi_lambdasq <- apply(curr_phis * curr_lambdas ^ 2, 2, sum)
  for (kk in 1 : K) {
    a_new <- curr_a2 + (K - kk + 1) * p / 2
    weights <- rep(NA, K)
    weights[1 : kk] <- 0
    if (kk < K) {
      weights[- c(1 : kk)] <- sapply((kk + 1) : K, function(x) prod(curr_deltas[(kk + 1) : x]))
    }
    b_new <- 1 + sum(weights * phi_lambdasq) / 2
    deltas[ss, kk] <- rgamma(1, shape = a_new, rate = b_new)
    curr_deltas[kk] <- deltas[ss, kk]
  }
  taus[ss, ] <- sapply(1 : K, function(x) prod(deltas[ss, 1 : x]))
  curr_taus <- taus[ss, ]
  
  
  # ---------- Update a1.
  # Proposing from gamma centered at the current value.
  prop_rate <- curr_a1 / tune_a1
  prop_shape <- prop_rate * curr_a1
  prop_a1 <- rgamma(1, shape = prop_shape, rate = prop_rate)
  
  # Calculating the acceptance probability.
  post_ratio <- dgamma(curr_deltas[1], shape = prop_a1, rate = 1)
  post_ratio <- post_ratio / dgamma(curr_deltas[1], shape = curr_a1, rate = 1)
  post_ratio <- post_ratio * dgamma(prop_a1, shape = prior_a1[1], rate = prior_a1[2])
  post_ratio <- post_ratio / dgamma(curr_a1, shape = prior_a1[1], rate = prior_a1[2])
  
  rev_rate <- prop_a1 / tune_a1
  rev_shape <- rev_rate * prop_a1
  prop_ratio <- dgamma(curr_a1, shape = rev_shape, rate = rev_rate)
  prop_ratio <- prop_ratio / dgamma(prop_a1, shape = prop_shape, rate = prop_rate)
  
  a1[ss] <- curr_a1
  if (runif(1) < post_ratio * prop_ratio) {
    acc_a1 <- acc_a1 + 1
    a1[ss] <- prop_a1
  }
  curr_a1 <- a1[ss]
  
  
  # ---------- Update a2.
  prop_rate <- curr_a2 / tune_a2
  prop_shape <- curr_a2 * prop_rate
  prop_a2 <- rgamma(1, shape = prop_shape, rate = prop_rate)
  
  # The ratio of the posterior.
  post_ratio <- sum(dgamma(curr_deltas[- 1], shape = prop_a2, rate = 1, log = TRUE)) -
    sum(dgamma(curr_deltas[- 1], shape = curr_a2, rate = 1, log = TRUE)) +
    dgamma(prop_a2, shape = prior_a2[1], rate = prior_a2[2], log = TRUE) -
    dgamma(curr_a2, shape = prior_a2[1], rate = prior_a2[2], log = TRUE)
  
  # Proposal ratio.
  rev_rate <- prop_a2 / tune_a2
  rev_shape <- prop_a2 * rev_rate
  prop_ratio <- dgamma(curr_a2, shape = rev_shape, rate = rev_rate, log = TRUE) -
    dgamma(prop_a2, shape = rev_shape, rate = rev_rate, log = TRUE)
  
  a2[ss] <- curr_a2
  if (log(runif(1)) < post_ratio + prop_ratio) {
    acc_a2 <- acc_a2 + 1
    a2[ss] <- prop_a2
  }
  curr_a2 <- a2[ss]
}

c(acc_a1, acc_a2) / (Nsims - 1)


burn <- Nsims / 3
thin <- 10
keep <- seq(burn + 1, Nsims, by = thin)
keep_Nsims <- length(keep)

etas_short <- etas[keep, , ]
phis_short <- phis[keep, , ]
a1_short <- a1[keep]
a2_short <- a2[keep]
deltas_short <- deltas[keep, ]
taus_short <- taus[keep, ]
lambdas_short <- lambdas[keep, , ]
sigmas_short <- sigmas[keep, ]


ll <- 1
kk <- 1
ii <- 1
plot(sigmas_short[, ll], type = 'l')
plot(lambdas_short[, ll, kk], type = 'l')
plot(etas_short[, ii, kk], type = 'l')
plot(a1_short, type = 'l')
plot(a2_short, type = 'l')


llt <- array(apply(lambdas_short, 1, function(x) x %*% t(x)),
             dim = c(p, p, keep_Nsims))
variances <- array(sapply(1 : keep_Nsims, function(x) llt[, , x] +
                            diag(sigmas_short[x, ])), dim = dim(llt))


plot(diag(true_sigma), diag(apply(variances, c(1, 2), mean)))
abline(a = 0, b = 1)


par(mfrow = dim(variances)[c(1, 2)], mar = rep(0, 4))
for (ii in 1 : p) {
  for (jj in 1 : p) {
    plot(density(variances[ii, jj, ]), axes = FALSE, main = '')
    abline(v = true_sigma[ii, jj], col = 'red')
  }
}


par(mfrow = dim(lambdas_short)[c(2, 3)], mar = rep(1, 4))
for (l in 1 : p) {
  for (k in 1 : K) {
    plot(density(lambdas_short[, l, k]), main = '', xlim = range(lambdas_short))
  }
}
