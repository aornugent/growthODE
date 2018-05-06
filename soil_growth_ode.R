library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(foreach)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

# Simulate data with deSolve ----------------------------------------------
source("simple_ode/simple_growth_function.R")

N = 10
S = 10
mu = round(seq(0.5, 1, length.out = S), 2)

R_soil <- round((100/N) * 1:N)

sim <- foreach(i = iter(R_soil), 
               .combine='rbind',
               .packages=c('tidyverse', 'deSolve')) %dopar% {
  grow(n = S, mu = c(mu), R_soil = i)
}

g <- ggplot(sim, aes(x = time, y = Biomass)) +
  geom_line(aes(colour = species, linetype = "deSolve")) +
  facet_wrap(~ R_soil)


# # Simulate data in Stan ---------------------------------------------------
# 
# stan_data <- list(
#   N = n,
#   S = S,
#   mu = mu,
#   T = 100, 
#   ts = seq(1, 100, 1),
#   R_light = 60,
#   R_soil = R_soil
# )
# 
# params = c("y_end")
# 
# # Compile model
# test_sim <- stan("simple_ode/growth_ode_soil2.stan",
#                  data = stan_data,
#                  pars = params, 
#                  chains = 1, 
#                  iter = 1,
#                  algorithm = "Fixed_param")
# 
# post <- rstan::extract(test_sim, pars = params)
# 
# sim_stan <- data.frame(
#   time = 100,
#   bind_rows(post),
#   species = rep(c("B1", "B2"), each = n),
#   R_soil = rep(R_soil, times = S)
# )
# 
# # Compare simulations
# g + geom_point(data = sim_stan, 
#               aes(x = time, y = y_end,
#                   shape = "Stan")) +
#   scale_shape_manual(values = c(1)) +
#   labs(shape = "", linetype = "") +
#   guides(colour = F)


# Estimate growth rates with Stan -----------------------------------------

# Extract end points
# end <- spread(sim_stan, species, y_end) %>%
#   select(B1, B2)

# Issue with ordering >10sp
end <- filter(sim, time == 100) %>%
  select(species, Biomass, R_soil) %>%
  spread(species, Biomass)

stan_data <- list(
  N = N,
  S = S,
  y = end[, -1],
  T = 100, 
  ts = seq(1, 100, 1),
  J = 2,
  R = cbind(rep(60), R_soil)
)

params = c("mu")

# Compile model
test_est <- stan("simple_ode/growth_ode_est2.stan",
                 data = stan_data,
                 pars = params, 
                 chains = 1, 
                 iter = 1)

# Run model
est_stan <- stan(fit = test_est,
                 data = stan_data,
                 pars = params, 
                 chains = 3, 
                 iter = 100,
                 control = list(max_treedepth = 7))

# approx

est_stan <- vb(stan_model("simple_ode/growth_ode_est2.stan"),
                 data = stan_data,
                 pars = params, eta = 5)

get_elapsed_time(est_stan)

post <- rstan::extract(est_stan, par = params)

# Summarise posterior
mu_est <- apply(post$mu, 2, median)

# Compare estimates with truth
mu - mu_est
sum(abs(mu - mu_est))
