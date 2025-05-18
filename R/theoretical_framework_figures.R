# Theorical framework figures
#' Simulations in Theoretical Framework of Multilevel Bayesian Non-parametric Model Thesis
#'
#' This script contains simulations presented in the theoretical framework of the thesis "Multilevel Bayesian 
#' Non-parametric Model". Each simulation corresponds to specific figures. To ensure reproducibility,
#' set.seed(2024) should be used before executing any figure-specific code.
#' Some simulations are based on examples presented by https://enesdilber.github.io/Sprojects/565p.pdf.
#' This includes simulations for GEM distributions, Dirichlet processes, and various visualizations.
#' Usage of specific set.seed values is recommended for consistency with the documented figures.

##############################
# Libraries
library(MASS)
library(ggplot2) 
library(LaplacesDemon)
library(extraDistr)
library(ggtern)


##############################
# Auxiliary Functions

#' Generate GEM Probability Vector
#' Generates a probability vector using the Stick-breaking construction (GEM distribution).
rGem <- function(alpha) {
  p <- rbeta(1, 1, alpha)
  while (sum(p) < 1) {
    p <- c(p, (1 - sum(p)) * rbeta(1, 1, alpha))
  }
  return(p)
}

#' Simulate GEM Distributions
#' Simulates multiple GEM distributions and averages the resulting samples.
rGEM <- function(n_sim, alpha) {
  print(alpha)
  p_list <- vector("list", n_sim)
  for (i in 1:n_sim) {
    p_list[[i]] <- rGem(alpha)
  }
  max_len <- max(sapply(p_list, length))
  simul <- sapply(p_list, function(x) set_same_len(x, max_len))
  return(rowMeans(simul))
}

#' Normalize Vector Length
#' Ensures all vectors are of `max_len` by padding zeros.
set_same_len <- function(vector, max_len) {
  length(vector) <- max_len
  vector[is.na(vector)] <- 0  # Replace NAs with zeros
  return(vector)
}

#' Log density of the multivariate t-distribution
log_multivariate_t <- function(x, mu, sigma, df) {
  p <- length(mu)
  x_mu <- x - mu
  log_det_Sigma <- log(det(sigma))
  
  term1 <- lgamma((df + p) / 2) - lgamma(df / 2)
  term2 <- -0.5 * (log_det_Sigma + p * log(df) + p * log(pi))
  term3 <- -((df + p) / 2) * log(1 + (1 / df) * t(x_mu) %*% solve(sigma) %*% x_mu)
  
  return(term1 + term2 + term3)
}

#' Function to get the mode of a vector
get_mode <- function(v) {
  uniq_v <- unique(v)
  freq_v <- tabulate(match(v, uniq_v))
  return(uniq_v[which.max(freq_v)])
}

##############################
# Simulations

# Simulation for Figure 3: Dirichlet Distribution with Parameter Vector alpha = (\alpha, \alpha, \alpha)}
set.seed(2024)
x_sim <- as.data.frame(rdirichlet(n = 500, alpha = rep(0.1, 3)))
colnames(x_sim) <- c('X1', 'X2', 'X3')  

# Create ternary plot for Dirichlet distribution
ggtern(data = x_sim, aes(X1, X2, X3)) + 
  stat_density_tern(
    geom = 'polygon',
    aes(fill = ..level..),
    bdl = 0.010,  
    bdl.val = 0,
    bins = 10,  
    color = 'grey'
  ) + 
  geom_point(size = 1.2) + 
  labs(fill = "PDF") +
  theme(legend.position = "right") +   
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "PuBu")) +  
  theme_minimal()

# Simulation for Figure 4: Effect of \eqn{\alpha} on the Number of Unique Components
n <- 1000
k <- 100

x_sim <- numeric(1000)
for (alpha in c(0.1, 1, 10)) {
  for (i in 1:1000) {
    pi_100 <- rdirichlet(n = 1, alpha = rep(alpha, k))
    x_sim[i] <- length(unique(rcat(n, pi_100)))
  }
}

df <- data.frame(alpha = as.factor(rep(c("0.1", "1", "10"), each = n)),
                 x = rep(1:n, times = 3),
                 y = x_sim)

df$alpha <- factor(df$alpha, levels = c("0.1", "1", "10"))

# Create the plot showing effect of alpha on number of unique components
ggplot(df, aes(x, y)) + 
  geom_line(aes(colour = alpha)) +
  xlab("Sample Size") + 
  ylab("Number of Unique Groups") +
  xlim(c(0, n)) + 
  ylim(c(0, k)) +
  scale_colour_manual(name = expression(alpha),  
                      values = c("#A6BDDB", "#3690C0", "#045A8D"),  
                      labels = c(expression(0.1), expression(1), expression(10))) +
  theme_classic()

# Simulation for Figure 7: GEM Distribution Simulation for \eqn{\alpha = (0.1, 1, 10, 100)}

alpha <- c(0.1, 1, 10)
X <- lapply(alpha, function(x) cumsum(rGEM(10, x)))

label <- prob <- c()
for (i in 1:length(X)) {
  label <- c(label, rep(alpha[i], length(X[[i]])))
  prob <- c(prob, X[[i]])
}

df <- data.frame(alpha = as.character(label), probability = prob)

# Cumulative probabilities plot
ggplot(df, aes(alpha, probability)) + 
  theme(legend.position = "none") +
  geom_jitter(width = 0.04) +
  xlab(expression(alpha)) +
  ylab("Cumulative Probability") +
  theme_classic()

# Simulation for Figure 8: Dirichlet Process (DP) with \eqn{G_0 = N(0, 1)} and Various \eqn{\alpha} Values
alpha_vals <- c(0.1, 1, 10)
X <- lapply(alpha_vals, function(a) rGEM(10, a))

labels <- probability <- x_vals <- c()

for (i in 1:length(X)) {
  labels <- c(labels, rep(alpha_vals[i], length(X[[i]])))
  x_vals <- c(x_vals, rnorm(length(X[[i]])))
  probability <- c(probability, X[[i]])
}

counts <- data.frame(theta = x_vals, probability = probability, alpha = labels)

# Visualize Dirichlet Process with G0 = N(0, 1)
ggplot(counts, aes(x = theta, y = probability)) +
  geom_segment(aes(xend = theta, yend = 0), size = 2, lineend = "butt") +
  facet_wrap(~alpha, nrow = 1, scales = 'free', 
             labeller = label_bquote(alpha == .(alpha))) +
  ylab(expression(f(x))) +
  xlab("x") +
  theme_bw()

# Simulation for Figure 10: Non-Parametric GMM Simulation (Stick-Breaking)

n <- 1000  
mu <- c(0, 0)  
sigma <- diag(2)  
alpha_values <- c(0.1, 1, 5, 10)  

X <- matrix(numeric(), nrow = 0, ncol = 2)
Z <- c()

for (alpha in alpha_values) {
  p <- rGem(alpha)  
  z <- sample(seq_along(p), n, p, replace = TRUE)
  m <- mvtnorm::rmvnorm(length(p), mean = mu, sigma = sigma)  
  
  x <- matrix(0, n, 2)
  for (i in 1:n) {
    x[i, ] <- mvtnorm::rmvnorm(1, mean = m[z[i], ], sigma = sigma * 0.01)
  }
  X <- rbind(X, x)
  Z <- c(Z, z)
}

X <- as.data.frame(X)
X$Group <- as.factor(Z)
X$Alpha <- rep(alpha_values, each = n)

blue_palette <- colorRampPalette(c("#A1C9E1", "#102C54"))(46)

# Display the GMM plot
ggplot(X, aes(x = V1, y = V2, color = Group)) +
  geom_point() +
  facet_wrap(~ Alpha, ncol = 2, labeller = label_bquote(alpha == .(alpha))) +
  scale_color_manual(values = blue_palette) + 
  theme_bw() +
  theme(legend.position = "none")

# Table 1: Number of Unique Tables in CRP Across Different Alpha Values

n_customers <- 500  
alpha_values <- rep(c(0.1, 1, 10, 50), each = 10)  

unique_tables_count <- c()
list_array <- list()

for (a in alpha_values) {
  tables <- rep(0, n_customers)  
  tables[1] <- 1  
  
  for (n in 2:n_customers) {
    table_count <- table(tables[1:(n - 1)])
    table_prob <- c(table_count, a) / (n - 1 + a)
    
    unique_table_labels <- seq(1, length(unique(tables)))
    tables[n] <- sample(unique_table_labels, size = 1, prob = table_prob)
  }
  
  list_array[[as.character(a)]] <- sort(table(tables))
  unique_tables_count <- c(unique_tables_count, length(unique(tables)))
}

rbind(alpha_values, unique_tables_count)

# Simulation for Figure 12

# Simulate data X
set.seed(2024)
n <- 300
means <- list(c(0, 0), c(3, 3), c(-3, -3))
covariances <- list(
  matrix(c(0.3, 0.05, 0.05, 0.3), ncol = 2), 
  matrix(c(0.5, -0.08, -0.08, 0.2), ncol = 2),
  matrix(c(0.1, 0.03, 0.03, 0.1), ncol = 2)
)

X <- do.call(rbind, lapply(1:3, function(i) mvrnorm(n = n / 3, mu = means[[i]], Sigma = covariances[[i]])))

# Hyperparameters: Normal prior
mu_0 <- colMeans(X)
k_0 <- 10
S_0 <- diag(2) * 2
v_0 <- 4

alpha_vector <- c(0.001, 1, 20, 50)

maxIterations <- 100
dim <- ncol(X)
cluster_result <- c()

pb <- txtProgressBar(min = 0, max = maxIterations * length(alpha_vector), style = 3)

for (alpha in alpha_vector) {
  z <- rep(1, n) 
  n_cluster <- n
  max_cluster <- 1
  
  result <- matrix(NA, nrow = n, ncol = maxIterations)
  
  for (iter in seq_len(maxIterations)) {
    for (i in seq_len(n)) {
      
      cluster_i <- z[i]
      n_cluster[cluster_i] <- n_cluster[cluster_i] - 1 
      
      if (n_cluster[cluster_i] == 0) {
        max_cluster <- max_cluster - 1
        n_cluster <- n_cluster[-cluster_i]
        
        inds <- z > cluster_i
        z[inds] <- z[inds] - 1
      }
      
      z[i] <- -1
      
      probs <- numeric(max_cluster + 1)
      
      for (c in seq_len(max_cluster)) {
        k_n <- k_0 + n_cluster[c]
        v_n <- v_0 + n_cluster[c]
        
        X_k <- X[z == c,]
        
        if (length(which(z ==  c)) == 1) {
          mu_k <- (k_0 * mu_0 + X_k) / k_n
          S_k <- S_0 + k_0 * (X_k - mu_0) %*% t(X_k - mu_0) / k_n
        } else {
          X_k_mean <- colMeans(X_k)
          mu_k <- (k_0 * mu_0 + n_cluster[c] * X_k_mean) / k_n
          S_k <- S_0 + cov(X_k) + k_0 * n_cluster[c] * (X_k_mean - mu_0) %*% t(X_k_mean - mu_0) / k_n
        }
        
        df_k <- v_n - dim + 1
        sigma_k <- (k_n + 1) * S_k / (k_n * df_k)
        prob_xi <- log_multivariate_t(X[i,], mu_k, sigma_k, df_k)
        
        probs[c] <- log(n_cluster[c]) + prob_xi
      }
      
      df_k <- v_0 - dim + 1
      sigma_k <- (k_0 + 1) * S_0 / (k_0 * df_k)
      prob_xi <- log_multivariate_t(X[i,], mu_0, sigma_k, df_k)
      
      probs[max_cluster + 1] <- log(alpha) + prob_xi
      
      max_probs <- max(probs)
      probs <- exp(probs - max_probs)
      probs <- probs / sum(probs)
      
      new_comp <- sample(seq_len(max_cluster + 1), 1, prob = probs)
      
      if (new_comp == (max_cluster + 1)) {
        n_cluster <- c(n_cluster, 1)
        max_cluster <- max_cluster + 1
      } else {
        n_cluster[new_comp] <- n_cluster[new_comp] + 1
      }
      
      z[i] <- new_comp
    }
    result[, iter] <- z
    
    setTxtProgressBar(pb, iter + (which(alpha_vector == alpha) - 1) * maxIterations)
  }
  
  cluster_result <- c(cluster_result, apply(result[, 10:100], 1, get_mode))
}

close(pb)

col_n <- length(unique(cluster_result))
blue_palette <- colorRampPalette(c("#A1C9E1", "#102C54"))(col_n)

graph_data <- cbind(rbind(X, X, X, X), cluster_result, rep(alpha_vector, each = n))
colnames(graph_data) <- c('X1', 'X2', 'group', 'alpha')
graph_data <- as.data.frame(graph_data)
graph_data$group <- as.factor(graph_data$group)

ggplot(graph_data, aes(x = X1, y = X2, color = group)) +
  geom_point() +
  facet_wrap(~ alpha, ncol = 2, labeller = label_bquote(alpha == .(alpha))) + 
  labs(x = expression(X[1]), y = expression(X[2])) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = blue_palette)