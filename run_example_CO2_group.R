#' Example of Fitting the Multilevel Bayesian Non-Parametric Model with Grouping
#'
#' This script provides an example of setting up and fitting a Multilevel Bayesian Non-Parametric Model 
#' using data that includes repeated measurements. It demonstrates the use of a grouping variable, derived 
#' from the dataset, to appropriately cluster data points that belong to similar categories.
#'
#' The script begins by performing an Ordinary Least Squares (OLS) fitting to initialize primary hyperparameters, 
#' which are subsequently used to fit the model.
#' 

# Libraries
source('R/MultilevelBNPModel.R')  # Includes the class for MultilevelBNPModel
library(tidyverse)                # For data manipulation

# Load Data

df <- readRDS("Data/Canada-fuel-consumption-ratings.rds")

# Select response and predictor variables
y <- df$`CO2 emissions (g/km)`
X <- df %>% select(c("Combined (L/100 km)", "Z", "X"))

# Create a grouping variable for repeated measures or categories
group <- as.numeric(as.factor(df$marca_tipo))

# Perform Ordinary Least Squares (OLS) 
X <- as.matrix(X) 
n <- nrow(X) 
X <- cbind(rep(1, n), X) 
p <- ncol(X)
beta_ols <- solve(t(X) %*% X) %*% t(X) %*% y  # OLS estimates
sig2_ols <- sum((y - X %*% beta_ols)^2) / (n - p)  # Residual variance estimate

# Define hyperparameters and process settings
hyperparams <- list(
  beta_0 = beta_ols, 
  Sigma_0 = nrow(X) * sig2_ols * solve(t(X) %*% X),
  nu0 = 1,
  v0 = sig2_ols,
  alpha_a = 1,
  alpha_b = 10,
  update_param_prob = 3
)

process_settings <- list(
  n_sams = 10000,
  n_burn = 50000,
  n_skip = 10
)

# Instantiate and fit the MultilevelBNPModel 
model_instance <- MultilevelBNPModel$new(y, X, hyperparams, process_settings, group = group)

# Fit the model with grouping enabled
results <- model_instance$fit(y, X, grouped = TRUE) 

# View summaries of the model fitting process
model_instance$summary(analysis_type = "log_likelihood")
model_instance$summary(analysis_type = "alpha")
model_instance$summary(analysis_type = "gamma")
model_instance$summary(analysis_type = "sigma")
model_instance$summary(analysis_type = "beta")