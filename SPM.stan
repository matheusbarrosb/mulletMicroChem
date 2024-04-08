data {
  int<lower=1> TT;                // TIME SERIES LENGTH
  vector[TT]   y;                 // OBSERVED DATA
  int<lower=1> K;                 // No OF STATES
  real         x0obs;             // OBSERVED VALUE FOR INITIAL STATE              
}
parameters {
  vector<lower=0,upper=0.1>[TT]  x;          // STATES
  real                           x0;         // INITIAL STATE
  simplex[K]                     theta;      // MIXTURE COMPONENTS
  vector[K]                      mu;         // MEAN OF MIXTURE COMPONENTS
  vector<lower=0>[K]             sigmaInit;  // SD OF INITIAL STATE
  vector<lower=0>[K]             sigma;      // SD OF MIXTURE COMPONENTS
  real<lower=0>                  sigmaState; // STATE SD
  real<lower=0>                  sigmaObs;   // OBSERVATION SD
}
model {
  x0          ~ normal(x0obs, 0.01);         // PRIOR FOR INITIAL STATE
  x[1]        ~ normal(x0, sigmaInit);       // INITIAL STATE
  mu[1]       ~ normal(0.075, 0.01);         // MEAN OF FIRST MIXTURE COMPONENT
  mu[2]       ~ normal(0.03, 0.01);          // MEAN OF SECOND MIXTURE COMPONENT
  sigmaInit   ~ cauchy(0,1);                 // THIS BECOMES A HALF-CAUCHY WITH THE LOWER CONSTRAINT DEFINED ABOVE
  sigma       ~ cauchy(0,1);
  sigmaState  ~ cauchy(0,1);
  sigmaObs    ~ cauchy(0,1);
  theta       ~ dirichlet(rep_vector(1,K));
  for (i in 1:TT) {
    for (k in 1:K) {
      target += log_sum_exp(log(theta[1]) + normal_lpdf(y[i] | mu[1], sigma[1]),
                            log(theta[2]) + normal_lpdf(y[i] | mu[2], sigma[2]));
    }
  }
  for (t in 2:TT) { # STATE PROCESS
    x[t] ~ normal(x[t-1], sigmaState);
  }
  for (t in 1:TT) { # OBSERVATION PROCESS
    target += normal_lpdf(y[t] | x[t], sigmaObs);
  }
}
      








