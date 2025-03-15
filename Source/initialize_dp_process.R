

initialize_dp_process = function(dp_objeto){
  
  dp_objeto$gamma = rep(1, dp_objeto$n)
  dp_objeto$K = 1
  dp_objeto$tamaño_K = dp_objeto$n
  
  dp_objeto$parameters_K =  list()
  dp_objeto$parameters_K$beta = mvtnorm::rmvnorm(1, 
                                                 dp_objeto$hiperparams$beta_0, 
                                                 dp_objeto$hiperparams$Sigma_0)
  dp_objeto$parameters_K$sigma =  1/rgamma(n = 1, 
                                           shape = 0.5*dp_objeto$hiperparams$nu0, 
                                           rate = 0.5*dp_objeto$hiperparams$nu0*dp_objeto$hiperparams$v0)

  #dp_objeto$parameters_K$alpha = dp_objeto$hiperparams$alpha
  dp_objeto$parameters_K$alpha =  rgamma(n = 1, 
                                         shape = dp_objeto$hiperparams$alpha_a, 
                                         rate = dp_objeto$hiperparams$alpha_b)

  dp_objeto$update_param_prob = dp_objeto$hiperparams$update_param_prob
  return(dp_objeto)
}
