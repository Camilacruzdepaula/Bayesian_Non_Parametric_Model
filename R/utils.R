
#' Utility Functions for MCMC and Chain Analysis
#'
#' This file contains utility functions essential for processing, analyzing, and summarizing
#' Markov Chain Monte Carlo (MCMC) simulation results. Functions include methods for calculating
#' lengths of MCMC objects, indexing with specific interval skips, summarizing chain characteristics,
#' and normalizing data structures for analysis.
#'
#' Details:
#' - `len_objects_mc`: Calculate total iterations needed given process settings.
#' - `index_save`: Determine which iterations to save based on burn-in and skip intervals.
#' - `summary_chain`: Generate a summary plot and descriptive statistics for a given chain.
#' - `get_mode`: Identify the mode of a given vector.
#' - `summary_pred_group`: Summarize prediction group results from chain data.
#' - `normalize_chain_beta`: Normalize beta chains, ensuring consistent dimensions across observations.
#'

library(ggplot2)

#' Calculate Total Iterations for MCMC Process
#' Given process settings, this function calculates the total number of iterations
#' required for the MCMC procedure, accounting for sampling, burn-in, and skip intervals.
#'
#' Params: process_settings A list containing process settings, n_sams, n_burn, n_skip.
#'
#' Return An integer representing the total number of iterations.

len_objects_mc <- function(process_settings) {
  if (!is.null(process_settings$n_skip)) {
    B <- process_settings$n_sams * process_settings$n_skip + process_settings$n_burn
  } else {
    B <- process_settings$n_sams + process_settings$n_burn
  }
  return(B)
}

#' Determine Iterations to Save
#' Calculates which iterations of the MCMC process should be saved, taking into account
#' burn-in period and thinning intervals.
#'
#' Params: b Current iteration number.
#'         process_settings A list of settings that include \code{n_burn} and \code{n_skip}.
#'
#' Return: An integer index indicating which iterations to save, or \code{NULL} if not to be saved.

index_save <- function(b, process_settings) {
  if (b > process_settings$n_burn) {
    if (!is.null(process_settings$n_skip)) {
      if (b %% process_settings$n_skip == 0) {
        b_save <- (b - process_settings$n_burn) / process_settings$n_skip
      } else {
        b_save <- NULL
      }
    } else {
      b_save <- b
    }
  } else {
    b_save <- NULL
  }
  return(b_save)
}

#' Summarize MCMC Chain
#' Provides a summary of a given MCMC chain, including a trace plot, effective sample size,
#' Geweke diagnostic, and descriptive statistics.
#'
#' Params: chain A numeric vector containing the MCMC sample chain data.
#'         title Optional title for the plot.
#'
#' Return A printed summary of the chain's characteristics and a ggplot object.

summary_chain <- function(chain, title = "Chain Summary") {
  
  chain <- na.omit(chain)
  
  data_chain <- data.frame(Iteration = seq_len(length(chain)), Chain = chain)
  # Plotting the trace plot
  print(ggplot(data_chain, aes(x = Iteration, y = Chain)) +
          geom_line(color = "black", linewidth = 0.5) + 
          theme_minimal() +
          labs(title = title, x = "Iteration", y = "Chain Value"))
  
  ess <- coda::effectiveSize(chain)
  geweke_result <- coda::geweke.diag(chain)$z
  summary_stats <- summary(chain)
  
  # Printing statistics
  cat("\nSummary of the Chain\n")
  cat("------------------------------\n")
  
  cat("Effective Sample Size:\n")
  print(format(ess, digits = 4, nsmall = 2))
  
  cat("\nGeweke Diagnostic (z-scores):\n")
  print(format(geweke_result, digits = 4, nsmall = 2))
  
  cat("\nDescriptive Statistics:\n")
  print(format(summary_stats, digits = 4, nsmall = 2))
  cat("------------------------------\n")
}

#' Get Mode of a Vector
#' This function calculates the mode (most frequently occurring value) of a given numeric vector.
#'
#' Params: v A numeric vector.
#'
#' Return: A single number representing the mode of the input vector.

get_mode <- function(v) {
  unique_vals <- unique(v)
  unique_vals[which.max(tabulate(match(v, unique_vals)))]
}

#' Summarize Prediction Groups from Chain Data
#' Summarizes the group predictions from MCMC chain data by calculating the mode for each time point.
#'
#' Param: chain A list where each element is a numeric vector representing an MCMC sample.
#'
#' Return: Printed output showing the distribution of predicted groups.

summary_pred_group <- function(chain) {
  predicted_group <- apply(do.call(rbind, chain), 2, get_mode)
  predicted_group <- table(predicted_group)
  cat("Group Prediction Summary\n")
  for (group in names(predicted_group)) {
    cat("Group", group, ":", predicted_group[group], "\n")
  }
}

#' Normalize Beta Chains
#' Ensures that all beta chains have consistent dimensions by adjusting the number of groups.
#' Fills in missing rows with NA or truncates to align with max_groups.
#'
#' Params: x A matrix representing beta chains.
#'         max_groups Integer indicating the maximum number of groups.
#'
#' Return: A matrix with aligned dimensions across beta chains.

normalize_chain_beta <- function(x, max_groups) {
  p <- ncol(x)
  if (nrow(x) < max_groups) {
    for (r in seq_len(max_groups - nrow(x))) {
      x <- rbind(x, rep(NA, p))
    }
    return(x)
  } else {
    return(x[1:max_groups, ])
  }
}