#' Initialize the Dirichlet Process Object
#'
#' This function initializes the Dirichlet Process object by setting initial clusters
#' and parameters for `beta`, `sigma`, and `alpha`. It considers whether the data includes 
#' grouping and adjusts initial values accordingly.

#' Initialize the Dirichlet Process Object
#'
#' Prepares the Dirichlet Process object, initializing parameters including clusters and distribution
#' parameters for the model. Adjusts initialization based on the presence of group data.
#'
#' Params: dp_object List containing model hyperparameters and current state information.
#'
#' Retrun: dp_object initialized ready for model fitting.

initialize_dp_process <- function(dp_object) {

  dp_object$K <- 1
  if (!is.null(dp_object$group)) {
    # Adjustments for grouped data
    dp_object$m <- length(unique(dp_object$group))
    dp_object$size_K <- dp_object$m
    dp_object$gamma <- rep(1, dp_object$m)
  }
  else{
    # Set up for non-grouped data
    dp_object$gamma <- rep(1, dp_object$n)
    dp_object$size_K <- dp_object$n
  }
  dp_object$hyperparams$Sigma_0_inv = solve(hyperparams$Sigma_0)

  # Initialize cluster parameters
  dp_object$parameters_K <- list()
  dp_object$parameters_K$beta <- mvtnorm::rmvnorm(1, 
                                                  dp_object$hyperparams$beta_0, 
                                                  dp_object$hyperparams$Sigma_0)
  dp_object$parameters_K$sigma <- 1/rgamma(n = 1, 
                                           shape = 0.5 * dp_object$hyperparams$nu0, 
                                           rate = 0.5 * dp_object$hyperparams$nu0 * dp_object$hyperparams$v0)
  
  dp_object$parameters_K$alpha <- rgamma(n = 1, 
                                         shape = dp_object$hyperparams$alpha_a, 
                                         rate = dp_object$hyperparams$alpha_b)

  return(dp_object)
}