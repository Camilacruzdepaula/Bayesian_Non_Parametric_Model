setwd("/Users/lauracamcruz/Documents/No_parametrica_HM_individuos") 
source('Likelihood.R')

sample_gamma = function(dp_objeto){

  m_probs = dp_objeto$update_param_prob
  
  for (i in seq_len(dp_objeto$n)){
    grupo_actual = dp_objeto$gamma[i]
    dp_objeto$tamaño_K[grupo_actual] = dp_objeto$tamaño_K[grupo_actual] - 1
    
    if(dp_objeto$tamaño_K[grupo_actual] == 0){
      aux_params = SamplePrior(dp_objeto, m_probs - 1)
      
      aux_params$beta = rbind(dp_objeto$parameters_K$beta[grupo_actual, ], 
                              aux_params$beta)
      aux_params$sig = c(dp_objeto$parameters_K$sigma[grupo_actual], 
                         aux_params$sig)
    } else {
      aux_params = SamplePrior(dp_objeto, m_probs)
    }
    
    lw = c(
      log(dp_objeto$tamaño_K) + log_likelihood_i(y[i], X[i, ], dp_objeto$parameters_K),
      log(dp_objeto$parameters_K$alpha/m_probs) + log_likelihood_i(y[i], X[i, ], aux_params))
    
    new_gamma = sample(x = 1:(dp_objeto$K + m_probs), size = 1, replace = FALSE, prob = exp(lw - max(lw)))
    
    if (new_gamma <= dp_objeto$K) {
      dp_objeto$tamaño_K[new_gamma] = dp_objeto$tamaño_K[new_gamma] + 1
      dp_objeto$gamma[i] = new_gamma
      
      if (dp_objeto$tamaño_K[grupo_actual] == 0) {
        dp_objeto$K = dp_objeto$K - 1
        dp_objeto$tamaño_K = dp_objeto$tamaño_K[-grupo_actual]
        
        dp_objeto$parameters_K$beta = as.matrix(dp_objeto$parameters_K$beta[-grupo_actual,], ncol = 3)
        if(dim(dp_objeto$parameters_K$beta)[2] ==1 ){
          dp_objeto$parameters_K$beta = t(dp_objeto$parameters_K$beta)
        }
        dp_objeto$parameters_K$sigma = dp_objeto$parameters_K$sigma[-grupo_actual]
        
        inds_change = dp_objeto$gamma > grupo_actual
        dp_objeto$gamma[inds_change] = dp_objeto$gamma[inds_change] - 1
      }
    } else {
      
      index_aux = new_gamma - dp_objeto$K
      
      if (dp_objeto$tamaño_K[grupo_actual] == 0) {
        
        dp_objeto$tamaño_K[grupo_actual] = 1
        dp_objeto$parameters_K$beta[grupo_actual, ] = aux_params$beta[index_aux, ]
        dp_objeto$parameters_K$sigma[grupo_actual] = aux_params$sig[index_aux]

      } else {

        new_gamma = dp_objeto$K + 1
        dp_objeto$parameters_K$beta = rbind(dp_objeto$parameters_K$beta,
                                            aux_params$beta[index_aux, ])
        dp_objeto$parameters_K$sigma[new_gamma] = aux_params$sig[index_aux]
        
        dp_objeto$K = dp_objeto$K + 1
        dp_objeto$tamaño_K[new_gamma] = 1
        dp_objeto$gamma[i] = new_gamma
      }
    }
  }
  return(dp_objeto)
}
