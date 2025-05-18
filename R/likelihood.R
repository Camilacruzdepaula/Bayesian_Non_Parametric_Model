#' Log Likelihood and Prior Sampling Functions
#'
#' This file contains functions that are used to calculate log likelihoods
#' and sample from the prior distribution in the context of a Bayesian
#' non-parametric model. It includes functionality for updating model parameters
#' based on grouped or non-grouped data.
#' 

#' Calculate the Log Probability
#'
#' Computes the log probability of the given data and parameters, optionally considering group structures
#' within the data if grouped is set to TRUE.
#' params: y Numeric vector representing the response variable.
#'         X Numeric matrix representing the predictor variables.
#'         dp_object List containing model parameters and group information.
#'         grouped Logical flag indicating whether to use grouped logic.
#' Return: A numeric value representing the log probability.

calculate_LP <- function(y, X, dp_object, grouped) {
  gamma <- dp_object$gamma
  beta <- dp_object$parameters_K$beta
  sigma <- dp_object$parameters_K$sigma
  
  if (grouped) {
    gamma <- gamma[dp_object$group]
  }

  mean_y_estimate <- rowSums(X * beta[gamma, ])
  lp <- sum(dnorm(x = y, mean = mean_y_estimate, sd = sqrt(sigma[gamma]), log = TRUE))
  
  return(lp)
}

#' Log Likelihood for Each Parameter Group
#'
#' Calculates the log likelihood for each group of parameters given the observations.
#' params: y_ Numeric vector of observations.
#'         X_ Matrix of predictor variables, matching y_.
#'         params List containing parameters beta and sigma for each group.
#' Return: A numeric vector of log likelihood values for each parameter group.

log_likelihood_i <- function(y_, X_, params) {
  log_like_y <- c()
  for (i in 1:length(params$sig)) {
    beta <- as.matrix(params$beta[i, ], nrow = ncol(X_))
    sigma <- params$sig[i]
    mean_y_estimate <- X_ %*% beta
    log_like_y[i] <- sum(dnorm(x = y_, mean = mean_y_estimate, sd = sqrt(sigma), log = TRUE))
  }
  return(log_like_y)
}

#' Sample from Prior Distribution
#'
#' Generates samples from the prior distributions for the parameters
#' params:  dp_object List containing model hyperparameters.
#'          p_times Number of samples to draw for each parameter.
#' Return:  A list containing sampled \code{beta} and \code{sigma}.

SamplePrior <- function(dp_object, p_times) {
  hyperparams <- dp_object$hyperparams
  sig <- 1/rgamma(n = p_times, shape = 0.5 * hyperparams$nu0, rate = 0.5 * hyperparams$nu0 * hyperparams$v0)
  
  beta <- mvtnorm::rmvnorm(p_times, hyperparams$beta_0, hyperparams$Sigma_0)
  
  theta <- list(beta = beta, sig = sig)
  return(theta)
}