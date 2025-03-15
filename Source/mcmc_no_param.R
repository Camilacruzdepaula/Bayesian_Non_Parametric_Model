setwd("/Users/lauracamcruz/Documents/No_parametrica_HM_individuos") 
source('initialize_dp_process.R')
source('utils.R')
source('sample_gamma.R')
source('sample_beta.R')
source('sample_sigma.R')
source('sample_alpha.R')

fit_dp_process = function (y, X, hiperparams, process_settings, progressBar = TRUE){
  
  set.seed(27)
  # Set y 
  y = as.matrix(y, nrow = length(y))
  # Create Process Dirichlet Object
  dp_objeto = list()
  dp_objeto$hiperparams = hiperparams
  dp_objeto$hiperparams$Sigma_0_inv = solve(dp_objeto$hiperparams$Sigma_0)
  dp_objeto$n = dim(y)[1]
  dp_objeto$p = dim(X)[2]
  
  # Inicializar Process Dirichlet Object
  
  dp_objeto = initialize_dp_process(dp_objeto)
 
  # Crear objetos de mc
  
  B = len_objects_mc(process_settings)
  len_B = process_settings$n_sams
  LP_mc = numeric(len_B)
  alpha_mc = numeric(len_B)
  omega_mc = vector("list", length = len_B)
  gamma_mc = vector("list", length = len_B)
  beta_mc = vector("list", length = len_B)
  sigma_mc = vector("list", length = len_B)
  
  ### Hacemos el fit
  if (progressBar) {
    pb = txtProgressBar(min = 0, max = B, width = 50, char = "-", style = 3)
  }
  
  for (b in seq_len(B)) {
    # Update gamma
    dp_objeto = sample_gamma(dp_objeto)
    
    # Update Beta
    dp_objeto = sample_beta(y, X, dp_objeto)

    # Update sigma
    dp_objeto = sample_sig(y, X, dp_objeto)
    
    # Update alpha
    # dp_objeto = sample_alpha(dp_objeto)
    dp_objeto = sample_alpha_2(dp_objeto)
    
    b_save = index_save(b, process_settings)
        
    if(!is.null(b_save)){
      omega_mc[[b_save]] = dp_objeto$tamaño_K/dp_objeto$n
      beta_mc[[b_save]] = dp_objeto$parameters_K$beta
      sigma_mc[[b_save]] = dp_objeto$parameters_K$sigma
      gamma_mc[[b_save]] = dp_objeto$gamma
      alpha_mc[b_save] = dp_objeto$parameters_K$alpha
      
      LP_mc[b_save] = calculate_LP(y, X, dp_objeto)

    }
    
    if (progressBar) {
      setTxtProgressBar(pb, b)
    }
  }
  if (progressBar) {
    close(pb)
  }
  return(list(beta = beta_mc, 
              sigma2 = sigma_mc, 
              omega= omega_mc, 
              gamma = gamma_mc, 
              alpha = alpha_mc,
              LP = LP_mc))
}
