#' Sampling Functions for Bayesian Non-Parametric Model
#'
#' This file includes functions to sample parameters for a Bayesian non-parametric model.
#' These functions perform sampling for gamma, beta, sigma, and alpha with considerations for 
#' grouped or non-grouped data.
#'

# Imports
source('R/likelihood.R')

#' Sample Gamma
#'
#' Updates the gamma assignments based on current parameters in the Dirichlet Process,
#' considering both grouped and non-grouped approaches.
#'
#' Params: y Numeric vector representing the response variable.
#'         X Numeric matrix representing the predictor variables.
#'         dp_object List containing model parameters and hyperparameters.
#'         grouped Logical flag indicating if the data is grouped.
#' Return: dp_object with new gamma assignments.

sample_gamma <- function(y, X, dp_object, grouped = FALSE) {
  m_probs <- dp_object$hyperparams$update_param_prob

  for (i in if (grouped) seq_len(dp_object$m) else seq_len(dp_object$n)) {
    current_group <- dp_object$gamma[i]
    dp_object$size_K[current_group] <- dp_object$size_K[current_group] - 1
    
    if (grouped){
      kj <- which(dp_object$group == i)
      yk <- y[kj]
      Xk <- X[kj, ]
    }
    else{
      yk <- y[i]
      Xk <- X[i, ]
    }
    
    if (dp_object$size_K[current_group] == 0) {
      aux_params <- SamplePrior(dp_object, m_probs - 1)
      aux_params$beta <- rbind(dp_object$parameters_K$beta[current_group, ], aux_params$beta)
      aux_params$sig <- c(dp_object$parameters_K$sigma[current_group], aux_params$sig)
    } else {
      aux_params <- SamplePrior(dp_object, m_probs)
    }
    
    lw <- c(
      log(dp_object$size_K) + log_likelihood_i(yk, Xk, dp_object$parameters_K),
      log((dp_object$parameters_K$alpha / m_probs)) + log_likelihood_i(yk, Xk, aux_params)
    )
    
    new_gamma <- sample(x = 1:(dp_object$K + m_probs), size = 1, replace = FALSE, prob = exp(lw - max(lw)))
    
    if (new_gamma <= dp_object$K) {
      dp_object$size_K[new_gamma] <- dp_object$size_K[new_gamma] + 1
      dp_object$gamma[i] <- new_gamma
      
      if (dp_object$size_K[current_group] == 0) {
        dp_object$K <- dp_object$K - 1
        dp_object$size_K <- dp_object$size_K[-current_group]
        
        dp_object$parameters_K$beta <- as.matrix(dp_object$parameters_K$beta[-current_group,], ncol = 3)
        
        if (dim(dp_object$parameters_K$beta)[2] == 1) {
          dp_object$parameters_K$beta <- t(dp_object$parameters_K$beta)
        }
        dp_object$parameters_K$sigma <- dp_object$parameters_K$sigma[-current_group]
        
        inds_change <- dp_object$gamma > current_group
        dp_object$gamma[inds_change] <- dp_object$gamma[inds_change] - 1
      }
    } else {
      index_aux <- new_gamma - dp_object$K
      
      if (dp_object$size_K[current_group] == 0) {
        dp_object$size_K[current_group] <- 1
        dp_object$parameters_K$beta[current_group, ] <- aux_params$beta[index_aux, ]
        dp_object$parameters_K$sigma[current_group] <- aux_params$sig[index_aux]
      } else {
        new_gamma <- dp_object$K + 1
        dp_object$parameters_K$beta <- rbind(dp_object$parameters_K$beta, aux_params$beta[index_aux, ])
        dp_object$parameters_K$sigma[new_gamma] <- aux_params$sig[index_aux]
        
        dp_object$K <- dp_object$K + 1
        dp_object$size_K[new_gamma] <- 1
        dp_object$gamma[i] <- new_gamma
      }
    }
  }
  return(dp_object)
}

#' Sample Beta Coefficients
#'
#' Samples the beta coefficients for each group using a conjugate multivariate normal distribution.
#'
#' Params: y Numeric vector of response values.
#'         X Matrix of predictor variables.
#'         dp_object List containing model parameters.
#'         grouped Logical indicating whether to account for group structure.
#'
#' Return: An updated \code{dp_object} with sampled beta coefficients.

sample_beta <- function(y, X, dp_object, grouped = FALSE) {
  K <- dp_object$K # Current number of groups
  p <- dp_object$p  # Size of Beta
  beta_0 <- dp_object$hyperparams$beta_0 # Hyperparameter Beta_0
  Sigma_0_inv <- dp_object$hyperparams$Sigma_0_inv # Inverse Sigma_0
  Sigma_0 <- dp_object$hyperparams$Sigma_0 # Sigma_0
  sigma <- dp_object$parameters_K$sigma # Sigma
  gamma <- dp_object$gamma # Group assignment vector
  
  beta_k <- matrix(NA, K, p)
  for (k in 1:K) {
    kj <- which(gamma == k)

    if (grouped) {
      yk = y[dp_object$group %in% kj]
      Xk = X[dp_object$group %in% kj, ]
    }
    else{
      yk <- y[kj]
      Xk <- X[kj, ]
    }
    njk <- length(yk)
    
    if (njk == 0) {
      beta_k[k, ] <- mvtnorm::rmvnorm(1, beta_0, Sigma_0)
    } else {
      if (njk == 1) {
        Xk <- as.matrix(Xk, 1, p)
        XtX_k <- Xk %*% t(Xk)
        Xty_k <- Xk * yk
      } else {
        XtX_k <- t(Xk) %*% Xk
        Xty_k <- t(Xk) %*% yk
      }
      
      beta.sigma_k <- solve(Sigma_0_inv + XtX_k / sigma[k])
      if (!isSymmetric(beta.sigma_k)) {
        beta.sigma_k <- (beta.sigma_k + t(beta.sigma_k)) / 2
      }
      beta.mean_k <- beta.sigma_k %*% (Sigma_0_inv %*% beta_0 + Xty_k / sigma[k])
      
      beta_k[k, ] <- mvtnorm::rmvnorm(1, beta.mean_k, beta.sigma_k)
    }
  }
  beta_k <- as.matrix(beta_k, ncol = ncol(X))
  
  dp_object$parameters_K$beta <- beta_k
  return(dp_object)
}

#' Sample Sigma (Variance)
#'
#' Samples sigma values for each group using an inverse gamma distribution, incorporating observed data variability.
#'
#' Params: y Numeric vector of response values.
#'         X Matrix of predictor variables.
#'         dp_object List containing model parameters and hyperparameters.
#'         grouped Logical indicating whether to consider group structure in sampling.
#'
#' Return: An updated dp_object with sampled sigma values.

sample_sigma <- function(y, X, dp_object, grouped) {
  K <- dp_object$K
  beta <- dp_object$parameters_K$beta
  nu0 <- dp_object$hyperparams$nu0
  v0 <- dp_object$hyperparams$v0
  gamma <- dp_object$gamma
  
  sigma2 <- vector()
  for (k in 1:K) {
    kj <- which(gamma == k)

    if (grouped) {
      yk = y[dp_object$group %in% kj]
      Xk = X[dp_object$group %in% kj, ]
    }
    else{
      yk <- y[kj]
      Xk <- X[kj, ]
    }
    
    njk <- length(yk)
    if (njk == 0) {
      sigma2[k] <- 1/rgamma(n = 1, shape = 0.5 * nu0, rate = 0.5 * nu0 * v0)
    } else {
      # Calculate SSR (Sum of Squared Residuals)
      y_est <- Xk %*% beta[k, ]
      SSR <- sum((yk - y_est)^2)
      a <- 0.5 * (nu0 + njk)
      b <- 0.5 * (nu0 * v0 + SSR)
      sigma2[k] <- 1/rgamma(n = 1, shape = a, rate = b)
    }
  }
  
  dp_object$parameters_K$sigma <- sigma2
  return(dp_object)
}

#' Sample Alpha Parameter
#'
#' Updates the alpha parameter using auxiliary variable sampling,
#' with different logic for grouped and non-grouped cases.
#'
#' Params: dp_object List containing model hyperparameters and current state.
#'         grouped Logical indicating whether to account for group structure.
#'
#' Returns: The updated dp_object with the new alpha value.

sample_alpha <- function(dp_object, grouped) {
  
  if(grouped){
    aux = dp_object$m
  }
  else{
    aux = dp_object$n
  }
  
  a <- dp_object$hyperparams$alpha_a
  b <- dp_object$hyperparams$alpha_b
  
  x <- rbeta(1, dp_object$parameters_K$alpha + 1, aux)
  
  pi1 <- a + dp_object$K - 1
  pi2 <- aux * (b - log(x))
  ratio <- pi1 / (pi1 + pi2)
  
  postParams1 <- a + dp_object$K
  postParams2 <- b - log(x)
  
  if (runif(1) > ratio) {
    postParams1 <- postParams1 - 1
  }
  
  dp_object$parameters_K$alpha <- rgamma(1, postParams1, postParams2)
  return(dp_object)
}


sample_alpha_2 = function(dp_objeto, grouped){
  
  if(grouped){
    aux = dp_object$m
  }
  else{
    aux = dp_object$n
  }

  a = dp_objeto$hiperparams$alpha_a
  b = dp_objeto$hiperparams$alpha_b
  
  x = rbeta(1, dp_objeto$parameters_K$alpha, aux)
  
  dp_objeto$parameters_K$alpha = rgamma(1, a + dp_objeto$K, b - log(x))
  return(dp_objeto)
}
