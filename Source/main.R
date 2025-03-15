
####################################
# Simulación modelo con grupos
###################################

# Librerías --------------------------------------
setwd("/Users/lauracamcruz/Documents/No_parametrica_HM_individuos") 
source('mcmc_no_param.R')
source('summary_dp_process.R')

#### Subir resultados

### FUNCIONES ------------------------
simulate_task = function(p, grupo, dummies, param_beta, param_sigma){
  set.seed(27)
  n = 500 * grupo
  grupo_vec = rep(c(1:grupo), each = 500)
  
  # Inicializar matriz de diseño
  X = matrix(1, n, 1)  # Incluye el intercepto
  
  if (dummies != 0){
    for (d in 1:dummies) {
      probs = c(0.2 * d, 0.8)/(sum(0.2 * d) + 0.8)
      x_d= sample(c(0, 1), n, prob = probs, replace = TRUE)
      X = cbind(X, x_d)
    }
  }
  
  for (num in 1:(p-dummies-1)){
    x_num = round(rnorm(n, 5*grupo, 10), 2)
    X = cbind(X, x_num)
  }
  
  beta = param_beta[grupo_vec, ]
  epsilon = rnorm(n, mean = 0, sd = sqrt(param_sigma[grupo_vec]))
  y = rowSums(X * beta) + epsilon
  
  list(y = y, X = X, grupo = grupo_vec)
}

p = 3
dummies = 1
param_beta = matrix(c(5, 8, 3,
                      1, -4, -1, 
                      2, 4, -3), 3, 3, byrow = TRUE); param_beta
param_sigma =  c(4, 10, 7)

prueba= simulate_task(p, 3, dummies, param_beta, param_sigma)

y = prueba$y
X = prueba$X

plot(density(y))

# OLS ----------------------------------
p = ncol(X)
n = nrow(X)
beta_ols = solve(t(X)%*%X)%*%t(X)%*%y
sig2_ols = sum((y - X%*%beta_ols)^2)/(n-p)
# -------------------------------------

hiperparams = list(beta_0 = beta_ols,
              Sigma_0 = nrow(X)*sig2_ols*solve(t(X)%*%X),
              nu0 = 1,
              v0 = 100,
              alpha_a = 1,
              alpha_b = 5,
              update_param_prob = 5
              )

process_settings = list(n_sams = 1000,
                        n_burn = 100,
                        n_skip = 5)

start_time = Sys.time()
results = fit_dp_process(y, X, hiperparams, process_settings)
end_time = Sys.time()

end_time-start_time

# Summary
summary_dp(results, 'log_likelihood')
summary_dp(results, 'gamma')
omega = summary_dp(results, 'omega')
sigma = summary_dp(results, 'sigma')
beta = summary_dp(results, 'beta')
alpha = summary_dp(results, 'alpha')

# save(results, file = 'Primera_iteracion_p_3_k_3.Rdata')
load('Somulación_fail.Rdata')


#