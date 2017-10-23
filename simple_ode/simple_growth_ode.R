library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Simulate data with deSolve ----------------------------------------------
source("simple_ode/simple_growth_function.R")

# Growth rates
S = 20
mu = runif(S, 0.5, 1)

# Resource supply rates
R_light = 60
R_soil = 100

sim <- grow(n = S, mu = c(mu), 
            R_light = R_light, 
            R_soil = R_soil)

g <- ggplot(sim, aes(x = time, y = Biomass, group = species)) +
  geom_line(aes(color = "deSolve"))


# Simulate data in Stan ---------------------------------------------------

stan_data <- list(
  S = S,
  mu = mu,
  T = 100,
  ts = seq(1, 100, 1),
  J = 2,
  R = c(R_light, R_soil)
)

params = c("y_hat")

# Compile model
test_sim <- stan("simple_ode/growth_ode_sim.stan",
                 data = stan_data,
                 pars = params, 
                 chains = 1, 
                 iter = 1,
                 algorithm = "Fixed_param")

post <- rstan::extract(test_sim)

# Summarise posterior
sim_stan <- apply(post$y_hat, c(2, 3), median) %>%
  as.data.frame() %>%
  mutate(time = 1:100) %>%
  gather(species, Biomass, -time)

# Compare simulations
g + geom_line(data = sim_stan, 
              aes(x = time, y = Biomass, 
                  group = species, color = "Stan")) +
  scale_colour_manual(values = c("black", "red")) +
  labs(color = "")


# Estimate growth rates with Stan -----------------------------------------

# Extract end points
end <- filter(sim_stan, time == 100) 

stan_data <- list(
  S = S,
  y = end$Biomass,
  T = 100,
  ts = seq(1, 100, 1),
  J = 2,
  R = c(R_light, R_soil)
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
                 iter = 40)

get_elapsed_time(est_stan)

post <- rstan::extract(est_stan, par = "mu")

# Summarise posterior
mu_est <- apply(post$mu, 2, median)

# Compare estimates with truth
error = mu - mu_est
mean(abs(error))
