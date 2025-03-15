
sample_alpha = function(dp_objeto){
  
  a = dp_objeto$hiperparams$alpha_a
  b = dp_objeto$hiperparams$alpha_b
  
  x = rbeta(1, dp_objeto$parameters_K$alpha + 1, dp_objeto$n)
  
  pi1 = a + dp_objeto$K - 1
  pi2 = dp_objeto$n*(b- log(x))
  ratio = pi1/(pi1 + pi2)
  
  postParams1 = a + dp_objeto$K
  postParams2 = b - log(x)
  
  if (runif(1) > ratio) {
    postParams1 = postParams1 - 1
  }
  
  dp_objeto$parameters_K$alpha = rgamma(1, postParams1, postParams2)
  return(dp_objeto)
}

sample_alpha_2 = function(dp_objeto){
  
  a = dp_objeto$hiperparams$alpha_a
  b = dp_objeto$hiperparams$alpha_b
  
  x = rbeta(1, dp_objeto$parameters_K$alpha, dp_objeto$n)

  dp_objeto$parameters_K$alpha = rgamma(1, a + dp_objeto$K, b - log(x))
  return(dp_objeto)
}