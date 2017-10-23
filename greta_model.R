library(greta)
library(tensorflow)
source("greta_ode.R")

S = 9

R_soil = 100
mu_true = seq(0.5, 0.9, length.out = S)

y = c(3.00, 3.83, 4.98, 6.64, 9.13, 13.14, 20.27, 35.75, 116.09)

y = greta_array(y, dim = c(1, S))

mu = normal(0, 1, dim = S, truncation = c(0, 1.5))

y_hat = growth_ode(mu = mu, tmax = 100)

sigma = greta_array(0.01, dim = c(1, S))

# expected vectorised Stan syntax, not explicit dimensions 
distribution(y) = normal(y_hat, sigma, dim = c(1, S))


m <- model(mu)
draws <- mcmc(m, 
              warmup = 100,
              n_samples = 100)
              

summary(draws)

library(bayesplot)
mcmc_intervals(draws)
