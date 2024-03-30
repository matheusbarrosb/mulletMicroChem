data {
  int N;         // NUMBER OF OBSERVATIONS
  int K;         // NUMBER OF HABITATS/MIXTURE COMPONENTS
  vector[N] y;   // Sr/Ca OBSERVATIONS
}
parameters {
  ordered[K] mu;               // MEAN Sr/Ca FOR EACH K-TH HABITAT
  positive_ordered[K] sigma;   // SD of Sr/Ca FOR EACH K-TH HABITAT
  simplex[K] theta;            // MIXTURE PROPORTIONS
  vector<lower=0>[K] alpha;    // DIRICHLET RATE PARAMETER
}
model {
  for (k in 1:K) {
    mu[1] ~ normal(0.9, 0.3);  // DIFFUSE PRIORS BUT CENTERED AT MEAN
    mu[2] ~ normal(1.5, 0.3);
    sigma[1] ~ cauchy(0,1);    // HALF-CAUCHY PRIORS FOR SDs
    sigma[2] ~ cauchy(0,1); 
  }
  theta ~ dirichlet(alpha);
  alpha ~ uniform(1,500);
  for (i in 1:N) {
    for (k in 1:K) { // LIKELIHOOD
      target += log_sum_exp(log(theta[1]) + normal_lpdf(y[i]|mu[1], sigma[1]),
                            log(theta[2]) + normal_lpdf(y[i]|mu[2], sigma[2]));
    }
  }
}
generated quantities {
  vector[N] y_hat;  // POSTERIOR PREDICTIVE
  for (i in 1:N) {
    for (k in 1:K) {
      y_hat[i] = exp(log_mix(theta[k],
                              normal_lpdf(y[i]|mu[1], sigma[1]),
                              normal_lpdf(y[i]|mu[2], sigma[2])));
    }
  }
}











