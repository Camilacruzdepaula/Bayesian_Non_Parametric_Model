

sample_beta = function (y, X, dp_objeto) {
  
  K = dp_objeto$K # Número actual de grupos
  p = dp_objeto$p  # Tamaño Beta
  beta_0 = dp_objeto$hiperparams$beta_0 # Hiperparam Beta_0
  Sigma_0_inv = dp_objeto$hiperparams$Sigma_0_inv # Inversa Sigma_0
  Sigma_0 = dp_objeto$hiperparams$Sigma_0 # Sigma_0
  sigma = dp_objeto$parameters_K$sigma # sigma
  gamma = dp_objeto$gamma # Vector de asignación de grupos

  beta_k = matrix(NA, K, p)
  for (k in 1:K) {
    kj = which(gamma==k)
    yk = y[kj]
    Xk = X[kj, ]
    
    njk = length(yk)
    
    if (njk == 0) {
      beta_k[k, ] = mvtnorm::rmvnorm(1, beta_0, Sigma_0)
    } else {
      
      if (njk==1){
        Xk = as.matrix(Xk, 1, p)
        XtX_k = Xk%*%t(Xk)
        Xty_k = Xk*yk
      }else{
        XtX_k = t(Xk)%*%Xk
        Xty_k = t(Xk)%*%yk
      }
      
      beta.sigma_k = solve(Sigma_0_inv + XtX_k/sigma[k])
      if(!isSymmetric(beta.sigma_k)){
        beta.sigma_k = (beta.sigma_k + t(beta.sigma_k)) / 2
      }
      beta.mean_k = beta.sigma_k%*%(Sigma_0_inv%*%beta_0 + Xty_k/sigma[k])
      
      beta_k[k, ] = mvtnorm::rmvnorm(1, beta.mean_k, beta.sigma_k)
    }
  }
  beta_k = as.matrix(beta_k, ncol = ncol(X))
  
  dp_objeto$parameters_K$beta = beta_k
  return(dp_objeto)
}