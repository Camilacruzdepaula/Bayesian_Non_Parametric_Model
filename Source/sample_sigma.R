

sample_sig= function(y, X, dp_objeto){
  
  K = dp_objeto$K
  beta= dp_objeto$parameters_K$beta
  nu0 = dp_objeto$hiperparams$nu0
  v0 = dp_objeto$hiperparams$v0
  gamma =  dp_objeto$gamma
  
  
  sigma2 = vector()
  for (k in 1:K) {
    kj = which(gamma==k)
    yk = y[kj]
    Xk = X[kj, ]
    
    njk = length(yk)
    if (njk == 0) {
      sigma2[k] = 1/rgamma(n = 1, shape = 0.5*nu0, rate = 0.5*nu0*v0)
    } else {
      # Calculate SSR
      y_est = Xk%*%beta[k, ]
      SSR = sum((yk - y_est)^2)
      a = 0.5*(nu0 + njk)
      b = 0.5*(nu0*v0 + SSR)
      sigma2[k] = 1/rgamma(n = 1, shape = a, rate = b)
    }
  }
  
  dp_objeto$parameters_K$sigma = sigma2
  return(dp_objeto)
}
