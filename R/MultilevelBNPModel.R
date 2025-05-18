#' Multilevel Bayesian Non-Parametric Model Class
#'
#' This class implements a Multilevel Bayesian Non-Parametric Model using R6. It provides methods to initialize,
#' fit, and summarize the model, accounting for both grouped and non-grouped data scenarios.
#'
#' Details:
#' - `initialize`: Sets up the model with data, hyperparameters, and process settings.
#' - `fit`: Performs sampling to fit the model parameters based on specified settings.
#' - `summary`: Provides a comprehensive summary of the model fit, analyzing different aspects such as log-likelihood, 
#'    alpha, gamma, sigma, and beta.

source('R/initialize_dp_process.R')
source('R/utils.R')
source('R/sampling_functions.R')
library(R6)

# Define the class
MultilevelBNPModel <- R6Class("MultilevelBNPModel",
                              public = list(
                                hyperparams = NULL,
                                process_settings = NULL,
                                group = NULL,
                                dp_object = NULL,
                                results_fit = NULL,
                                
                                #' Initialize the Multilevel Bayesian Non-Parametric Model
                                #'
                                #' Sets up the model with the provided data, hyperparameters, and process settings. Supports optional grouping.
                                #'
                                #' Params: y Numeric vector or matrix representing the response variable.
                                #'         X Matrix representing predictor variables.
                                #'         hyperparams List of hyperparameters for initialization.
                                #'         process_settings List containing sampling settings.
                                #'         group Optional numeric vector indicating group membership.
                                #'
                                #' Return: An instantiated object of MultilevelBNPModel.

                                initialize = function(y, X, hyperparams, process_settings, group = NULL) {
                                  self$hyperparams <- hyperparams
                                  self$process_settings <- process_settings
                                  self$group <- group
                                  
                                  y <- as.matrix(y, nrow = length(y))
                                  self$dp_object <- list(
                                    hyperparams = hyperparams,
                                    n = dim(y)[1],
                                    p = dim(X)[2],
                                    group = group
                                  )
                                  
                                  self$dp_object <- initialize_dp_process(self$dp_object)
                                  
                                  self$results_fit <- list()
                                },
                                
                                #' Fit the Model
                                #'
                                #' Conducts the model fitting through iterative sampling based on the provided settings. Accounts for grouped data.
                                #'
                                #' Params: y Numeric vector of response values.
                                #'         X Matrix of predictor variables.
                                #'         grouped Logical indicating if data is grouped.
                                #'         progressBar Logical, displays progress if TRUE.
                                #'
                                #' Return: List with the results of model fitting.

                                fit = function(y, X, grouped = FALSE, progressBar = TRUE) {
                                  set.seed(42)
                                  B <- len_objects_mc(self$process_settings)
                                  len_B <- self$process_settings$n_sams
                                  LP_mc <- numeric(len_B)
                                  alpha_mc <- numeric(len_B)
                                  omega_mc <- vector("list", length = len_B)
                                  gamma_mc <- vector("list", length = len_B)
                                  beta_mc <- vector("list", length = len_B)
                                  sigma_mc <- vector("list", length = len_B)
                                  
                                  if (progressBar) {
                                    pb <- txtProgressBar(min = 0, max = B, width = 50, char = "-", style = 3)
                                  }
                                  
                                  for (b in seq_len(B)) {
                                    self$dp_object <- sample_gamma(y, X, self$dp_object, grouped)
                                    self$dp_object <- sample_beta(y, X, self$dp_object, grouped)
                                    self$dp_object <- sample_sigma(y, X, self$dp_object, grouped)
                                    self$dp_object <- sample_alpha(self$dp_object, grouped)
                                    
                                    b_save <- index_save(b, self$process_settings)
                                    
                                    if (!is.null(b_save)) {
                                      if(grouped){
                                        omega_mc[[b_save]] <- self$dp_object$size_K/self$dp_object$m
                                      }
                                      else{
                                        omega_mc[[b_save]] <- self$dp_object$size_K/self$dp_object$n
                                      }
                                      beta_mc[[b_save]] <- self$dp_object$parameters_K$beta
                                      sigma_mc[[b_save]] <- self$dp_object$parameters_K$sigma
                                      gamma_mc[[b_save]] <- self$dp_object$gamma
                                      alpha_mc[b_save] <- self$dp_object$parameters_K$alpha
                                      LP_mc[b_save] <- calculate_LP(y, X, self$dp_object, grouped)
                                    }
                                    
                                    if (progressBar) {
                                      setTxtProgressBar(pb, b)
                                    }
                                  }
                                  
                                  if (progressBar) {
                                    close(pb)
                                  }
                                  
                                  self$results_fit <- list(
                                    beta = beta_mc, 
                                    sigma2 = sigma_mc, 
                                    omega = omega_mc, 
                                    gamma = gamma_mc, 
                                    alpha = alpha_mc,
                                    LP = LP_mc
                                  )
                                  
                                  return(self$results_fit)
                                },
                                
                                #' Summarize the Model Fit
                                #'
                                #' Produces a summary of the model fit, analyzing specified aspects such as log-likelihood, alpha, gamma,
                                #' sigma, and beta.
                                #'
                                #' Params: analysis_type Character string specifying the type of analysis.
                                #'
                                #' Return: Printed output of the summarized results.
  
                                summary = function(analysis_type = NULL) {
                                  if (is.null(self$results_fit)) {
                                    stop("No results available. Please run the model fit first.")
                                  }
                                  
                                  if (analysis_type == "log_likelihood") {
                                    summary_chain(self$results_fit$LP, 'LP')
                                  }
                                  if (analysis_type == "alpha") {
                                    summary_chain(self$results_fit$alpha, 'Alpha')
                                  }
                                  if (analysis_type == "gamma") {
                                    summary_pred_group(self$results_fit$gamma)
                                  }
                                  if (analysis_type == "sigma") {
                                    max_groups = max(apply(do.call(rbind, self$results_fit$gamma), 2, get_mode))
                                    chain <- do.call(rbind, lapply(self$results_fit$sigma2, function(x) x[1:max_groups]))
                                    for (c in 1:ncol(chain)){
                                      cat("\nGROUP", c, '\n')
                                      summary_chain(chain[, c], paste('Sigma', c))
                                    }
                                  }
                                  if (analysis_type == "beta") {
                                    max_groups = max(apply(do.call(rbind, self$results_fit$gamma), 2, get_mode))
                                    chain <- lapply(self$results_fit$beta, function(x) normalize_chain_beta(x, max_groups))
                                    p <- ncol(chain[[1]])
                                    for(g in 1:max_groups){
                                      cat("\nGROUP", g, '\n')
                                      for (beta in 1:p){
                                        cat("\n> Beta", beta, '\n')
                                        beta_chain <- sapply(chain, function(x) x[g, beta])
                                        summary_chain(beta_chain, paste('Group', g, '- Beta', beta))
                                      }
                                    }
                                  }
                                }
                              )
)