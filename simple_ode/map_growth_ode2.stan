functions {
  // Resource-competition model
  real[] growth(real t, // time
                real[] B,  // state
                real[] mu, // parameters
                real[] x_r,
                int[] x_i) {
    
    int N = x_i[1]; // Number of equations
    int J = x_i[2]; // Number of resources
    real dB_dt[N];
    
    real theta[J] = {1.3, 1.0}; // Mode of resource division
    
    real c = 0.5;
    real M[J];
    real r[J];
    real B_prop;
    real B_sym[N, J];
    real B_total[J];
    
    real p[J];
    real growth;
    
    // Symmetry of competition
    for(j in 1:J){
      for(i in 1:N){
        B_sym[i, j] = B[i]^theta[j];
      }
      B_total[j] = sum(B_sym[, j]);
    }
    
    // Growth-maintenance correlation
    M[1] = c * mean(mu);
    
    for(i in 1:N){
      M[2] = c * mu[i];
      
      // Resource division
      for(j in 1:J){
        B_prop = B_sym[i, j] / B_total[j];
        r[j] = x_r[j] * B_prop * B[i]^-1 - M[j];
        p[j] = (r[j] / (1 + r[j]));
      }
      
      // Limiting resource
      growth = B[i] * mu[i] * min(p[]);
      
      // Positive growth only
      dB_dt[i] = (growth > 0 ? growth : 0);
    }
    
    return dB_dt;
  }
  
  // Observation function
  vector obs(
    vector phi,
    vector theta,
    real[] x_r,
    int[] x_i){
    
    // Unpack parameters
    int S = x_i[1];
    int J = x_i[2];
    int T = x_i[3];
    
    real mu[S] = to_array_1d(phi[1:S]);
    real sigma = phi[(S+1)];
    
    real y[S] = x_r[1:S];
    real y_hat[T, S] = integrate_ode_bdf(
                        growth, 
                        rep_array(0.5, S), // Initial values
                        0.0, 
                        x_r[(S+J+1):(S+J+T)], 
                        mu, 
                        x_r[(S+1):(S+J)], // Data variables
                        x_i,
                        1e-5, // rel. tol
                        1e-5, // abs. tol
                        1e3); // num. steps); 

    real lp = normal_lpdf(y | y_hat[T, ], sigma);

    return [lp]';
  }
}
data {
  int<lower=1> N;         // Number of obs.
  int<lower=1> S;         // Number of species
  int<lower=1> J;         // Number of resources
  int<lower=1> T;         // Number of timesteps
  real<lower=0> y[N, S];// Observed abundance
  real<lower=0> R[N, J];  // Resource supply rates
  real<lower=0> ts[T];    // Timesteps
}
transformed data{
  real x_r[N, (S + J + T)];   // Package data for threads
  int x_i[N, 3];
  vector[0] theta[N];     // No thread specific parameters
  
  for(i in 1:N){
    x_r[i] = append_array(append_array(y[i], R[i]), ts);
    x_i[i] = {S, J, T};
  }
}
parameters{
  vector<lower=0.001, upper=3>[S] mu; // Relative growth rates
  real<lower=0> sigma;                // Observation error (gaussian)
}
model{
  mu ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  target += sum(map_rect(obs, append_row(mu, sigma), theta, x_r, x_i));
}
