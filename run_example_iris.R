#' Example of Using the Multilevel Bayesian Non-Parametric Model
#'
#' This script demonstrates how to use the Multilevel Bayesian Non-Parametric Model,
#' as developed and described in the thesis "Modelo Bayesiano No Paramétrico Multinivel."
#' The data used in this example is derived from the iris dataset, showcasing the basic
#' setup and execution of the model.
#'
#' The script first conducts an Ordinary Least Squares (OLS) estimation to set initial
#' hyperparameters. Subsequently, these are employed to fit the model using the R6 class
#' structure developed.


# Load function scripts and the class definition
source('R/MultilevelBNPModel.R')

# Define data
y <- iris$Petal.Length
X <- iris$Sepal.Length

# Perform Ordinary Least Squares (OLS)
X <- as.matrix(X)
n <- nrow(X)
X <- cbind(rep(1, n), X)
p <- ncol(X)

# Estimate OLS coefficients and residual variance
beta_ols <- solve(t(X) %*% X) %*% t(X) %*% y
sig2_ols <- sum((y - X %*% beta_ols)^2) / (n - p)

# Define hyperparameters and process settings
hyperparams <- list(
  beta_0 = beta_ols,
  Sigma_0 = n * sig2_ols * solve(t(X) %*% X),
  nu0 = 1, v0 = sig2_ols,
  alpha_a = 1, alpha_b = 1,
  update_param_prob = 2
)

process_settings <- list(
  n_sams = 5000,
  n_burn = 1000,
  n_skip = 5
)

# Instantiate and fit the model
set.seed(123)  # Ensures reproducibility
model_instance <- MultilevelBNPModel$new(y, X, hyperparams, process_settings)
results <- model_instance$fit(y, X)

# View model summaries
model_instance$summary(analysis_type = "log_likelihood")
model_instance$summary(analysis_type = "alpha")
model_instance$summary(analysis_type = "gamma")
model_instance$summary(analysis_type = "sigma")
model_instance$summary(analysis_type = "beta")

# Optional: save results
# saveRDS(results, file = "results/multilevel_bnp_results.rds")



