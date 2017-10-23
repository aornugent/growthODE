# growthODE
Parameter estimation of DeMalach growth model: http://onlinelibrary.wiley.com/doi/10.1111/1365-2745.12557/

- `simple_ode/` contains a deSolve and Stan implementation.
- `greta_ode.R` has a bespoke Tensorflow ODE solver for this model.
- `greta_model.R` attempts to perform inference.
