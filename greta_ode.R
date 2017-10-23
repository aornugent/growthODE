op <- greta::.internals$nodes$constructors$op

tf_growth_function <- function(mu, S, tmax){
  
  # uptake rate
  k = 1
  
  # maintenance coef
  c_ = 0.5
  
  # mode of competition
  theta = c(1.3, 1)
  
  # supply rate
  R = c(60, 100)
  
  # maintenance costs, messy
  M1 <- tf$reshape(c_ * tf$reduce_mean(mu), shape(1, 1))
  M2 <- c_ * mu
  
  M <- tf$concat(c(tf$tile(M1, shape(S, 1)), M2), 1L)
  
  # M <- cbind(c_ * mean(mu), c_ * mu)
  
  # growth function
  grow <- function(B, t){
    # map elementwise
    B_sym = tf$map_fn(function(B) tf$pow(B, theta), B)
    B_sum = tf$reduce_sum(B_sym, axis = 0L)
    
    # limiting resources
    r = R * (B_sym/B_sum) * tf$pow(B, -1) - M
    p = r / (1 + r)
    p_min = tf$reduce_min(p, axis = 1L, keep_dims = T)
    
    # positive growth only
    zeros = tf$zeros(shape(S, 1))
    B_delta = B * mu * p_min
    
    dB = tf$where(tf$less(B_delta, 0), zeros, B_delta)
  }
  
  # initial values
  init = tf$constant(rep(0.5, S), shape = c(S, 1))
  
  # timesteps
  t = seq(0, tmax, 1)
  
  # solve
  cleanly(ode <- tf$contrib$integrate$odeint(
    func = grow,
    y0 = init,
    t = t,
    rtol = 1e-8,
    atol = 1e-5,
    method = 'dopri5',
    options = list(max_num_steps = 1000L),
    full_output = F,
    name = 'growth_ode')[tmax,,])
}

growth_ode <- function(mu, tmax = 100){
  
  # number species
  S = length(mu)
  
  # number of timesteps
  tmax = as.integer(tmax)
  
  dimfun <- function(elem_list) {
    # input dimensions
    state_dim <- dim(elem_list[[1]])
    
    if (length(state_dim) != 2 | state_dim[2] != 1)
      stop ('mu must be a column vector greta array',
            call. = FALSE)
    
    # output dimensions
    c(1, S)
  }
  
  op('growth_ode',
     mu = mu,
     operation_args = list(S = S,
                           tmax = tmax),
     tf_operation = tf_growth_function,
     dimfun = dimfun)
}


# patch to check for max_num_steps error
cleanly <- function (expr) {
  
  res <- tryCatch(expr, error = function (e) e)
  
  # if it errored
  if (inherits(res, 'error')) {
    
    numerical_messages <- c("is not invertible",
                            "Cholesky decomposition was not successful",
                            "max_num_steps exceeded")
    
    numerical_errors <- vapply(numerical_messages,
                               grepl,
                               res$message,
                               FUN.VALUE = 0) == 1
    
    # if it was just a numerical error, quietly return a bad value
    if (any(numerical_errors))
      res <- NA
    else
      stop ("greta hit a tensorflow error:\n\n", res, call. = FALSE)
    
  }
  
  res
}
