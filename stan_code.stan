// one parameter binomial model
data {
  int<lower=1> n;               // sample size
  int<lower=2> p;               // Dimension of the outcome.
  real y[p,n];                  // data
  real<lower=0> nu;             // Parameter for the prior on phi
  real<lower=0> prior_a1[2];
  real<lower=0> prior_a2[2];
  real<lower=0> prior_sigma[2];
  int<lower=1> K;
}
parameters {
  matrix[K,n] eta;
  row_vector[K] Lambda[p];
  real<lower=0> phi[p,K];
  real<lower=0> delta[K];
  real<lower=0> sigma_sq[p];
  real<lower=0> a1;
  real<lower=0> a2;
}
transformed parameters{
  real<lower=0> sigma[p];
  real<lower=0> tau[K];
  matrix[p,n] mean_of_y;
  for (l in 1 : p) {
    sigma[l] = sqrt(sigma_sq[l]);
  }
  for (k in 1 : K) {
    tau[k] = 1;
    for (k2 in 1 : k) {
      tau[k] = tau[k] * delta[k2];
    }
  }
  for (i in 1 : n) {
    for (l in 1 : p) {
      mean_of_y[l,i] = Lambda[l] * eta[,i];
    }
  }
}
model {
  // likelihood.
  for (i in 1 : n) {
    for (l in 1 : p) {
      y[l,i] ~ normal(mean_of_y[l,i], sigma[l]);
    }
  }
  // prior on etas.
  for (i in 1 : n) {
    for (k in 1 : K) {
      eta[k,i] ~ normal(0,1);
    }
  }
  // prior on lambdas and phis.
  for (l in 1 : p) {
    for (k in 1 : K) {
      Lambda[l,k] ~ normal(0, sqrt(1 / (phi[l,k] * tau[k])));
      phi[l,k] ~ gamma(nu / 2, nu / 2);
    }
  }
  // prior on the deltas.
  delta[1] ~ gamma(a1, 1);
  for (k in 2 : K) {
    delta[k] ~ gamma(a2, 1);
  }
  // prior on the sigmas.
  sigma_sq ~ inv_gamma(prior_sigma[1], prior_sigma[2]);
  // prior on a1, a2.
  a1 ~ gamma(prior_a1[1], prior_a1[2]);
  a2 ~ gamma(prior_a2[1], prior_a2[2]);
}
