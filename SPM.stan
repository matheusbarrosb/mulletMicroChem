data {
  int<lower=1> TT;                // TIME SERIES LENGTH
  vector[TT]   y;                 // OBSERVED DATA
  int<lower=1> K;                 // No OF STATES
}
parameters {
  vector<lower=0>[TT] x;          // STATES
  real                x0;         // INITIAL STATE
  simplex[K]          theta;      // MIXTURE COMPONENTS
  vector<lower=0>[K]  mu;         // MEAN OF MIXTURE COMPONENTS
  vector<lower=0>[K]  sigma;      // SD OF MIXTURE COMPONENTS
  real<lower=0>       sigmaState; // STATE SD
  real<lower=0>       sigmaObs;   // OBSERVATION SD
}
model {
  x0 ~ normal(0.04, 0.01);
  x[1] ~ normal(x0, 0.01);
  mu[1] ~ normal(0.03, 0.01);
  mu[2] ~ normal(0.075, 0.01);
  sigma ~ cauchy(0,1);
  sigmaState ~ cauchy(0,1);
  sigmaObs ~ cauchy(0,1);
  theta ~ dirichlet(rep_vector(1,K));
  for (t in 2:TT) {
    x[t] ~ normal(theta[1]*normal_lpdf(x[t-1] | mu[1], sigma[1]) +
                  theta[2]*normal_lpdf(x[t-1] | mu[2], sigma[2]), sigmaState);
  }
  for (t in 1:TT) {
    y[t] ~ normal(x[t], sigmaObs);
  }
}
         










