
calculate_LP = function (y, X, dp_objeto)
{
  gamma = dp_objeto$gamma
  beta = dp_objeto$parameters_K$beta
  sigma = dp_objeto$parameters_K$sigma
  
  mean_y_estimate = rowSums(X*beta[gamma,])
  lp = sum(dnorm(x = y, mean = mean_y_estimate, sd = sqrt(sigma[gamma]), log = TRUE))
  return(lp)
}

log_likelihood_i =  function(y_, X_, params) {
  log_like_y = c()
  for (i in 1:length(params$sig)){
    beta = as.matrix(params$beta[i, ], nrow = ncol(X_))
    sigma=  params$sig[i]
    mean_y_estimate = X_%*%beta
    log_like_y[i] = sum(dnorm(x = y_, mean = mean_y_estimate, sd = sqrt(sigma), log = TRUE)) 
  }
  return(log_like_y)
}


SamplePrior =  function(dp_objeto, p_times) {
  
  hiperparametros = dp_objeto$hiperparams
  sig =  1/rgamma(n = p_times, shape = 0.5*hiperparametros$nu0, 
                  rate = 0.5*hiperparametros$nu0*hiperparametros$v0)
  
  beta =  mvtnorm::rmvnorm(p_times, hiperparametros$beta_0,
                           hiperparametros$Sigma_0)
  
  theta = list(beta = beta, sig = sig)
  return(theta)
}