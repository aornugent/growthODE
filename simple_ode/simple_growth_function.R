library(tidyverse)
library(deSolve)
theme_set(theme_bw())


# Discrete growth and resource partitioning -------------------------------
# Wrapper with some default values
grow <- function(n = 20, t = 0:100, mu = NA, theta = c(1.3, 1), 
                 R_light = 60, R_soil = 30){

  # deSolve function
  g <- function(time,init,params) {
    with (as.list(c(time,init,params)), {
      
      # Initialise model
      n = length(init)
      B = init
      dB_dt <- vector(length = n)
      
      r <- vector(length = 2)
      M <- vector(length = 2)
      p <- vector(length = 2)
      
      # Coefficient of maintenance
      c = 0.5
      
      # Above ground maintenance (constant)
      M[1] = c * mean(mu)
      
      for(i in 1:n){
        # Below ground maintenance correlated with growth rate
        M[2] = c * mu[i]
        
        for(j in 1:2) {
          # Division of resources by prop. biomass.
          B_prop = B[i]^theta[j] / sum(B^theta[j])
          r[j] = R[j] * B_prop * B[i]^-1 - M[j]
        }
        
        # Rate limiting resource
        p = min(r / (k + r))
        
        # Growth (t + 1)
        B_p1 = B[i] * mu[i] * p
        
        # Positive grwoth only
        dB_dt[i] <- ifelse(B_p1 > 0, B_p1, 0)
      }
      
      return(list(dB_dt))
    })
  }
  
  # Initial biomass of 0.5g, same for all species
  init <- c(
    B = rep(0.5, times = n)
  )
  
  # Relative growth rates
  if(is.na(mu[1])){
    mu <- runif(n, 0.1, 2)
  }
  
  # Package parameters
  pars <- list(
    mu = c(mu),  
    k = 1,
    theta = theta,
    R = c(R_light, R_soil)
  )
  
  # Run for max(t) steps
  df <- ode(func=g, y=init, parms=pars, 
            times=t, method = "bdf") %>%
    as.data.frame() %>%
    gather(species, Biomass, -1) %>%
    mutate(time = time,
           mu = rep(mu, each = length(t)),
           R_light = R_light,
           R_soil = R_soil)
}