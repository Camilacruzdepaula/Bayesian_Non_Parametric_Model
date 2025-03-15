
len_objects_mc = function(process_settings){

  if (!is.null(process_settings$n_skip)){
    B = process_settings$n_sams*process_settings$n_skip +  process_settings$n_burn
  } else {
    B = process_settings$n_sams + process_settings$n_burn
  }
  return(B)
}


index_save = function(b, process_settings){
      
      
  if (b > process_settings$n_burn){
  
    if (!is.null(process_settings$n_skip)){
    
      if (b%%process_settings$n_skip == 0) {

        b_save = (b - process_settings$n_burn)/process_settings$n_skip
      } else {
  
        b_save = NULL
      }
      
    } else {
      
      b_save = b
    }
  } else {

    b_save = NULL
  }

  return(b_save)
}

calculate_num_groups = function(gamma_list) {
  num_groups = sapply(gamma_list, function(g) length(unique(g)))
  return(num_groups)
}


trace_plot_mc = function(cadena, titulo){
  
  data_cadena = data.frame(Iteracion= seq_len(length(cadena)), Cadena = cadena)
  
  p = ggplot(data_cadena, aes(x = Iteracion, y = Cadena)) +
    geom_point(color = "black", size = 0.5) +  # Línea del trace plot
    labs(title = paste("Plot", titulo),
         x = "Iteración",
         y = titulo) +
    theme_minimal() 
  print(p)
}

normalize_chain = function(cadena, gamma_df, tipo = 'sigma'){
  
  if(tipo == 'sigma'){
    
    groups = table(sapply(cadena, function(x) length(x)))/length(cadena)
    max_groups = as.integer(names(groups[groups == max(groups)]))
  
    max_pred = max(apply(gamma_df, 2, get_mode))
    max_groups = min(max_groups, max_pred)
    
    chain_norm = lapply(cadena, function(x){
                                  if (length(x) < max_groups){
                                    x[(length(x) + 1):max_groups] = NA}
                                  else{
                                    if (length(x) > max_groups){
                                      x = x[1:max_groups] }
                                  }
                                  return(x)
                                })
  } else {
    
    p = ncol(cadena[[1]])
    
    groups = table(sapply(cadena, function(x) nrow(x)))/length(cadena)
    max_groups = as.integer(names(groups[groups == max(groups)]))

    max_pred = max(apply(gamma_df, 2, get_mode))
    max_groups = min(max_groups, max_pred)
    
    chain_norm = lapply(cadena, function(x){
                            if (nrow(x) < max_groups){
                              for(r in seq_len(max_groups - nrow(x))){
                                x = rbind(x, rep(NA, p))
                              }
                            } 
                            else{
                              if (nrow(x) > max_groups){
                                 x = x[1:max_groups, ]
                            }}
                            return(x)
                          })
    }
  return(chain_norm)
}


plot_trace_chain_grouped = function(chain_norm, titulo){
  
  max_groups = length(chain_norm[[1]])
  chain_norm = do.call(rbind, chain_norm)
  chain_norm_stacked = stack(as.data.frame(chain_norm))
  chain_norm_stacked$ind = gsub("V", "G_", chain_norm_stacked$ind)
  chain_norm_stacked$Iteracion = rep(1:nrow(chain_norm), max_groups)
  
  p = ggplot(chain_norm_stacked, aes(x = Iteracion, y = values, color = ind)) +
    geom_line(alpha = 0.6) +
    labs(title = paste0("Plot of", titulo),
         x = "Iteración",
         y = titulo,
         color = "ind") +
    theme_minimal() +
    facet_wrap(~ ind, scales = "free_y")
  print(p)
  
  return(chain_norm_stacked)
}

get_mode = function(v) {
  uniq_vals = unique(v)
  uniq_vals[which.max(tabulate(match(v, uniq_vals)))]
}


summary_mc = function(df){
  summary_df = df %>%
    group_by(ind) %>%
    summarise(
      n = n(),
      Min = min(values, na.rm = TRUE),
      Q1 = quantile(values, 0.25, na.rm = TRUE),
      Mean = mean(values, na.rm = TRUE),
      Median = median(values, na.rm = TRUE),
      Q3 = quantile(values, 0.75, na.rm = TRUE),
      Max = max(values, na.rm = TRUE),
      StdDev = sd(values, na.rm = TRUE)
    )
  return(summary_df)
}


calculate_coef_gelman = function(chain_norm){
  
  chain_norm = do.call(rbind, chain_norm)
  list_chain = list()
  for (i in seq_len(ncol(chain_norm))){
    list_chain[[i]] = chain_norm[, i]
  }
  
  mcmc_chains = lapply(list_chain, coda::mcmc)
  combined_chains = coda::mcmc.list(mcmc_chains)
  coef_gelman = coda::gelman.diag(combined_chains)
  
  return(coef_gelman)
  
}

plot_acf_chain_grouped = function(chain_df, titulo){
  
  acf_results = chain_df%>%
    filter(!is.na(values))%>%
    group_by(ind) %>%
    do(data.frame(acf = acf(.$values, plot = FALSE)$acf, 
                  lag = acf(.$values, plot = FALSE)$lag))
  
  p = ggplot(acf_results, aes(x = lag, y = acf, color = ind)) +
    geom_hline(yintercept = 0, color = "gray") +  # Línea en y=0
    geom_segment(aes(xend = lag, yend = 0), color = "blue", size = 0.7) +
    labs(title = paste0("ACF ", titulo),
         x = "Lag",
         y = "ACF") +
    facet_wrap(~ ind, scales = "free_y") +
    theme_minimal() 
  print(p)
}


muestra_efectiva = function(chain_norm, titulo){
  
  n_chains = length(chain_norm[[1]])
  n_effective = c()
  
  for (i in seq_len(n_chains)){

    print(paste0('Tamaño de muestra Efectiva ', titulo, ' ',i, ':'))
    n_effective[i] = coda::effectiveSize(results$LP)
    print(n_effective[i] )
  }
  return(n_effective)
}


blue_palette = c(
  "#A4C8E1", "#8EB4DE", "#76A9D9", "#3498DB", "#2980B9", 
  "#1F5985", "#153A69", "#0D3E57", "#5A8B9A", "#477C91",
  "#3B6B82", "#2F586F", "#2A5B5E", "#1E4D4E", "#1A3B3F",
  "#195654", "#2C3E50", "#34495E", "#4A6572", "#4F6D75", 
  "#596B7A", "#617482", "#678B8F", "#72A6A8", "#7BB5B2",
  "#8DC3C4", "#99C8CE", "#A3D1D8", "#B1E3E7", "#C4F1F1",
  "#D0E8E8", "#DBE7E7", "#E2F0F0", "#E8F9F9", "#A9CEDB",
  "#B2D8E7", "#C0E7F0", "#D0F3F5", "#B1D2E1",
  "#2C3E50", "#4B778D", "#6C9AB3","#8BA9C4" 
)
