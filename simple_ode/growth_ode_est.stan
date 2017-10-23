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
}
data {
  int<lower=1> S;       // Number of species
  real<lower=0> y[S];   // Data
  int<lower=1> T;       // End point
  real<lower=0> ts[T];  // Time points
  int<lower=1> J;       // Number of resources
  real<lower=0> R[J];   // Resource supply rate
}
transformed data{
  real x_r[J] = R;
  int x_i[2] = {S, J};
}
parameters{
  real<lower=0, upper=2> mu[S];  // Relative growth rates
  real<lower=0> sigma;  // Observation error (gaussian)
}
transformed parameters{
  real<lower=0> y0[S] = rep_array(0.5, S); // Initial values
  real<lower=0> y_hat[T, S]; // Store ODE output
  
  // Simulate system
  y_hat = integrate_ode_rk45(growth, y0, 
                             0.0, ts, mu, x_r, x_i,
                             1e-8, //rel. tol
                             1e-5, //abs. tol
                             1e3); //num. steps
}
model{
  mu ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  // Inference from final state
  y ~ normal(y_hat[T, ], 0.1 * sigma);
}
