// STATE SPACE MIXTURE MODEL WITH KALMAN FILTER //
// ADAPTED FROM https://github.com/juhokokkala/kalman-stan-randomwalk //
// Differences from the implementation above include the mixture specification
// and that sigmaY is now estimated instead of supplied with the data
data {
  int                  T;        // LENGTH OF TIME SERIES
  int                  NST;      // NUMBER OF STATES
  vector[T]            y;        // OBSERVATIONS
  real                 m0;       // INITIAL OBSERVED STATE
  real<lower=0>        P0;       // VARIANCE OF INITIAL OBSERVED STATE (maybe dont specify as data - give it a normal prior?)
}
parameters {
  vector<lower=0>[T]   sigmaY;  // OBSERVATION SDs
  real<lower=0>        Z;       // STATE VARIANCE, TIME-INVARIANT
  real<lower=0>        Q;       // PROCESS NOISE VARIANCE
  simplex[NST]         theta;   // MIXTURE PROPORTIONS
  vector<lower=0>[NST] mu;      // MIXTURE COMPONENT MEANS
  vector<lower=0>[NST] sigma;   // MIXTURE COMPONENT SDs
}
transformed parameters {
  vector[T]            X;                  // STATE MEAN, FILTERED
  vector[T]            Xhat;               // PREDICTED MEAN
  vector[T]            P;                  // STATE VARIANCE
  vector[T]            Phat;               // PREDICTED VARIANCE   
  vector[T]            R;                  // MEASUREMENT VARIANCE 
  vector[T]            S;                  //
  real                 K;                  // KALMAN GAIN
  R = Z*Z + sigmaY.*sigmaY;
  Xhat[1] = m0;
  Phat[1] = P0;
  for (t in 1:T) {  // INITIALIZING
    if (t > 1) {
      Xhat[t] = X[t-1];
      Phat[t] = P[t-1] + Q*Q;
    }
   S[t] = Phat[t] + R[t];
   K = Phat[t]/(Phat[t] + R[t]);
   X[t] = Xhat[t] + K*(y[t] - Xhat[t]);
   P[t] = Phat[t] - K*S[t]*K;
  }
}
model {
  Z ~ cauchy(0,1);
  Q ~ cauchy(0,1);
  sigmaY ~ cauchy(0,1);
  theta ~ dirichlet(rep_vector(1,NST));
  mu[1] ~ normal(0.08, 0.01);
  mu[2] ~ normal(0.04, 0.01);
  sigma ~ cauchy(0,1);
  for (i in 1:T) {
    for (q in 1:NST) {
      target += log_sum_exp(log(theta[1]) + normal_lpdf(y[i] | mu[1], sigma[1]),
                            log(theta[2]) + normal_lpdf(y[i] | mu[2], sigma[2]));
    }
  }
  for(t in 1:T) {
    y[t] ~ normal(Xhat[t], sqrt(S[t]));
  }
}
generated quantities {
    vector[T] x; //The random-walking signal
    vector[T] z; //x + jitter
    x[T] = normal_rng(X[T], sqrt(P[T]));
    for (i in 1:T-1) {
        int t;
        real varx;
        real meanx;
        t = T-i;
        varx = 1 / (1/P[t] + 1/(Q*Q));
        meanx = (X[t]/P[t] + x[t+1]/(Q*Q))*varx;
        x[t] = normal_rng(meanx,sqrt(varx));
    }
    for (t in 1:T) {
        real meanz;
        real varz;
        varz = 1/ (1/(Z*Z));
        meanz = varz * (x[t]/(Z*Z));
        z[t] = normal_rng(meanz,sqrt(varz));
    }
}




