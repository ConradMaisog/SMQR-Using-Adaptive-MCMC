####################################################################################################################
# >--------------------------------------------------------------------------------------------------------------< #
# >                                Bayesian Quantile Regression with Adaptive MCMC                               < #
# >                                           John Conrad Seg B. Maisog                                          < #
# >                                          Adaptive MCMC Implementation                                        < #
# >--------------------------------------------------------------------------------------------------------------< #
####################################################################################################################







####################################################################################################################
# >--------------------------------------------------------------------------------------------------------------< #
# >                                                   PRELIMINARIES                                              < #
# >--------------------------------------------------------------------------------------------------------------< #

### Load necessary libraries
library(pacman)
p_load(coda, dplyr, invgamma, truncnorm, ald, mvtnorm, MASS, ggplot2, quantreg, moments, bayesQR)

# >--------------------------------------------------------------------------------------------------------------< #
####################################################################################################################




####################################################################################################################
# >--------------------------------------------------------------------------------------------------------------< #
# >                                                  DATA GENERATION                                             < #
# >--------------------------------------------------------------------------------------------------------------< #

### Parameters
sample_sizes <- c(20, 40, 100)
quantiles <- c(0.1, 0.25, 0.5, 0.75, 0.9)
B0 <- c(-2, -1, 0, 1, 2)
B1 <- c(-0.5, -0.25, 0, 0.25, 0.5)
B2 <- c(0.5, 0.25, 0, -0.25, -0.5)

### Function to generate data
generate_data <- function(sample_sizes, quantiles, B0, B1, B2) {
  ### Initialize empty list
  data_list <- list()
  
  ### Loop through sample sizes
  for (n in sample_sizes) {
    ### Loop through error distributions
    for (error_type in c("Standard Normal", "Right-Skewed Gamma", "Heteroscedastic")) { # , 
      ### Loop through quantiles
      for (i in seq_along(quantiles)) {
        q <- quantiles[i]
        beta_0 <- B0[i]
        beta_1 <- B1[i]
        beta_2 <- B2[i]
        
        ### Generate X from standard normal distribution
        X1 <- rnorm(n, mean = 0, sd = 1)
        X2 <- rnorm(n, mean = 0, sd = 1)
        
        ### Error Distributions
        if (error_type == "Standard Normal") {
          e <- rnorm(n, mean = 0, sd = 1)
        } else if (error_type == "Right-Skewed Gamma") {
          e <- rgamma(n, shape = 0.25, scale = 2)
        } else if (error_type == "Heteroscedastic") {
          e <- abs(X1) * rnorm(n, mean = 0, sd = 1)
        }
        
        Y <- beta_0 + beta_1 * X1 + beta_2 * X2 + e
        
        ### Same data into a dataframe
        df <- data.frame(
          Error = error_type,
          Sample_Size = n * 5,
          Quantile = q,
          B0 = beta_0,
          B1 = beta_1,
          B2 = beta_2,
          X1 = X1,
          X2 = X2,
          Y = Y,
          e = e
        )
        ### Append list
        data_list <- append(data_list, list(df))
      }
    }
  }
  ### Return dataframe
  return(do.call(rbind, data_list))
}

set.seed(123)
### Set number of replications
replications <- 100
replicated_data <- list()

### Generate Replicated Data
set.seed(123)
for (rep in 1:replications) {
  print(paste("Generating Data - Replication:", rep))  # Print replication number
  replicated_data[[rep]] <- generate_data(sample_sizes, quantiles, B0, B1)
}

# >--------------------------------------------------------------------------------------------------------------< #
####################################################################################################################






####################################################################################################################
# >--------------------------------------------------------------------------------------------------------------< #
# >                                          FREQUENTIST QUANTILE REGRESSION                                     < #
# >--------------------------------------------------------------------------------------------------------------< #

### Function to perform quantile regression and return results
perform_quantile_regression <- function(data, quantiles) {
  ### Initialize empty list
  results_list <- list()
  ### Loop through error distributions
  for (error_type in unique(data$Error)) {
    subset_data <- data %>%
      filter(Error == error_type, Sample_Size == n)
    ### Loop through sample sizes
    for (n in unique(data$Sample_Size)) {
      for (tau in quantiles) {
        model <- rq(Y ~ X1 + X2, tau = tau, data = subset_data)
        estimates <- coef(model)
        est_B0 <- estimates[1]
        est_B1 <- estimates[2]
        est_B2 <- estimates[3]
        
        results_list <- append(results_list, list(data.frame(
          Error = error_type,
          Sample_Size = n,
          Quantile = tau,
          Est_B0 = round(est_B0, 4),
          Est_B1 = round(est_B1, 4),
          Est_B2 = round(est_B2, 4)
        )))
      }
    }
  }
  
  return(do.call(rbind, results_list))
}

all_results <- list()
quantiles <- c(0.1, 0.25, 0.5, 0.75, 0.9)
replications <- 100
### Replicate the process 100 times
for (rep in 1:1) {
  print(paste("Running Frequentist Quantile Regression - Replication:", rep))  # Print replication number
  
  data <- plot_data #data <- replicated_data[[rep]]  # Retrieve the dataset for this replication
  results_df <- perform_quantile_regression(data, quantiles)
  results_df$Replication <- rep  # Add replication number
  
  all_results <- append(all_results, list(results_df))
}

### Combine all results into a single dataframe
frequentist_results <- do.call(rbind, all_results) 


# >--------------------------------------------------------------------------------------------------------------< #
####################################################################################################################





####################################################################################################################
# >--------------------------------------------------------------------------------------------------------------< #
# >                                                    DIAGNOSTICS                                               < #
# >--------------------------------------------------------------------------------------------------------------< #

### Compute Monte Carlo Standard Errors
compute_mcse <- function(samples) {
  return(sd(samples) / sqrt(length(samples)))
}

### Compute Inefficiency Factor
compute_inefficiency_factor_autocorr <- function(samples, max_lag = 100) {
  # 1) drop missing draws
  x <- na.omit(samples)
  
  # 2) if fewer than 2 points, no autocorr can be estimated
  if (length(x) < 2) {
    warning("Insufficient non-NA samples to compute inefficiency factor")
    return(NA_real_)
  }
  
  # 3) compute ACF up to max_lag, ignoring any remaining NAs
  acf_obj <- acf(x, plot = FALSE, lag.max = max_lag, na.action = na.pass)
  acf_vals <- acf_obj$acf[-1]
  
  # 4) sum up to get IF, skipping any NA acf lags
  IF <- 1 + 2 * sum(acf_vals, na.rm = TRUE)
  return(IF)
}

### Compute Effective Sample Size
compute_ess <- function(n, inefficiency_factor) {
  return(n / inefficiency_factor)
}

### 
compute_mcse_autocorr <- function(samples, inefficiency_factor) {
  n_eff <- compute_ess(length(samples), inefficiency_factor)
  return(sd(samples) / sqrt(n_eff))
}

# >--------------------------------------------------------------------------------------------------------------< #
####################################################################################################################




####################################################################################################################
# >--------------------------------------------------------------------------------------------------------------< #
# >                                      BAYESIAN QUANTILE REGRESSION WITH GIBBS                                 < #
# >--------------------------------------------------------------------------------------------------------------< #

### Function for Gibbs Bayesian Quantile Regression
run_bayesQR_analysis <- function(data, tau, burn_in = 5000, ndraw = 15000, max_lag = 100, display_plots = FALSE) {
  ### Define Prior (adapted to match variable names)
  prior <- prior(Y ~ X1 + X2, data = data, beta0 = rep(0, 3), V0 = 1 * diag(3), shape0 = 6, scale0 = 5)
  
  ### Start time tracking
  start_time <- proc.time()
  
  ### Fit Bayesian Quantile Regression
  capture.output({
    bqr_normal <- bayesQR(Y ~ X1 + X2, quantile = tau, prior = prior, ndraw = ndraw, normal.approx = FALSE, data = data)
  }, file = NULL)
  
  ### End time tracking
  end_time <- proc.time() - start_time
  
  ### Extract MCMC draws after burn-in
  mcmc_draws_normal_betas <- bqr_normal[[1]]$betadraw[(burn_in + 1):nrow(bqr_normal[[1]]$betadraw), ]
  mcmc_draws_normal_sigma <- bqr_normal[[1]]$sigmadraw[(burn_in + 1):length(bqr_normal[[1]]$sigmadraw)]
  
  ### Optional: Display Plots
  if (display_plots) {
    par(mfrow = c(2, 4))  # Adjust layout
    for (i in 1:3) {
      plot(mcmc_draws_normal_betas[, i], type = "l", main = paste("Traceplot of Beta", i-1), ylab = "Value", xlab = "Iteration")
      abline(h = true_beta[i], col = "red", lty = 1)
      acf(mcmc_draws_normal_betas[, i], main = paste("ACF of Beta", i-1))
    }
    plot(mcmc_draws_normal_sigma, type = "l", main = "Traceplot of Sigma", ylab = "Value", xlab = "Iteration")
    acf(mcmc_draws_normal_sigma, main = "ACF of Sigma")
  }
  
  ### Compute Metrics
  mcse_betas <- apply(mcmc_draws_normal_betas, 2, function(x) compute_mcse(x))
  inefficiency_factors <- apply(mcmc_draws_normal_betas, 2, function(x) compute_inefficiency_factor_autocorr(x, max_lag = max_lag))
  ess_betas <- sapply(1:3, function(i) compute_ess(length(mcmc_draws_normal_betas[, i]), inefficiency_factors[i]))
  
  ### Results Dataframe
  results_df <- data.frame(
    Error_Type   = unique(data$Error),
    Sample_Size  = unique(data$Sample_Size),
    Quantile     = tau,
    Est_B0       = round(mean(mcmc_draws_normal_betas[, 1]), 4),
    Est_B1       = round(mean(mcmc_draws_normal_betas[, 2]), 4),
    Est_B2       = round(mean(mcmc_draws_normal_betas[, 3]), 4),
    MCSE_B0      = mcse_betas[1],
    MCSE_B1      = mcse_betas[2],
    MCSE_B2      = mcse_betas[3],
    ESS_B0       = ess_betas[1],
    ESS_B1       = ess_betas[2],
    ESS_B2       = ess_betas[3],
    IF_B0        = inefficiency_factors[1],
    IF_B1        = inefficiency_factors[2],
    IF_B2        = inefficiency_factors[3]
  )
  #print("here")
  
  ### Return results and MCMC draws
  return(list(
    results = results_df,
    mcmc_draws = list(
      betas = mcmc_draws_normal_betas,
      sigma = mcmc_draws_normal_sigma
    )
  ))
}


### Run Bayesian Quantile Regression
all_bayesian_results <- list()
all_mcmc_draws <- list()

set.seed(123)
for (rep in 1:100) {
  print(paste("Running Bayesian Quantile Regression - Replication:", rep))
  
  # Start timing for this replication
  start_time <- Sys.time()
  
  data <- plot_data  # Retrieve dataset for this replication
  
  # Loop through each unique combination of Error Type and Sample Size
  for (error_type in unique(data$Error)) {
    print(paste("Error: ", error_type))
    for (sample_size in unique(data$Sample_Size)) {
      print(paste("Sample Size: ", sample_size))
      subset_data <- data %>%
        filter(Error == error_type, Sample_Size == sample_size)
      
      # Run Bayesian Quantile Regression for each quantile (tau)
      for (tau in quantiles) {
        print(paste("Quantile: ", tau))
        bayes_output <- run_bayesQR_analysis(subset_data, tau, display_plots = FALSE)
        
        # Store results dataframe
        results_df <- bayes_output$results
        results_df$Replication <- rep  # Add replication number
        results_df$Error_Type <- error_type
        results_df$Sample_Size <- sample_size
        results_df$Quantile <- tau
        
        print(results_df)
        
        all_bayesian_results <- append(all_bayesian_results, list(results_df))
        
        # Store MCMC draws separately (optional for further analysis)
        all_mcmc_draws[[paste0("Rep", rep, "_Err", error_type, "_Size", sample_size, "_Tau", tau)]] <- bayes_output$mcmc_draws
      }
    }
  }
  # End timing for this replication
  end_time <- Sys.time()
  
  # Calculate and print the running time for this replication
  running_time <- end_time - start_time
  print(paste("Running time for Replication", rep, ":", running_time))
}

gibbs_results <- do.call(rbind, all_bayesian_results)

# >--------------------------------------------------------------------------------------------------------------< #
####################################################################################################################







####################################################################################################################
# >--------------------------------------------------------------------------------------------------------------< #
# >                                  BAYESIAN QUANTILE REGRESSION WITH ADAPTIVE MCMC                             < #
# >--------------------------------------------------------------------------------------------------------------< #

### Function to compute log of the posterior distribution (target distribution)
log_posterior <- function(y, X1, X2, beta_0, beta_1, beta_2, sigma, tau,
                          prior_mean_beta0, prior_mean_beta1, prior_mean_beta2,
                          prior_type) {
  # Compute the linear predictor (location parameter for ALD)
  mu <- beta_0 + (beta_1 * X1) + (beta_2 * X2)
  
  # Calculate the log-likelihood using the Asymmetric Laplace Distribution (ALD)
  log_likelihood <- sum(sapply(1:length(y),
                               function(i) likALD(y[i], mu = mu[i], sigma = sigma,
                                                  p = tau, loglik = TRUE)))
  
  ### Dynamic Priors for betas
  log_prior_beta_0 <- if (prior_type == "normal") {
    dnorm(beta_0, mean = prior_mean_beta0, sd = 1, log = TRUE)
  } else {
    if (tau < 0.5) {
      if (beta_0 > prior_mean_beta0) -Inf else dexp(prior_mean_beta0 - beta_0,
                                                    rate = 1, log = TRUE)
    } else {
      if (beta_0 < prior_mean_beta0) -Inf else dexp(beta_0 - prior_mean_beta0,
                                                    rate = 1, log = TRUE)
    }
  }
  log_prior_beta_1 <- if (prior_type == "normal") {
    dnorm(beta_1, mean = prior_mean_beta1, sd = 1, log = TRUE)
  } else {
    if (tau < 0.5) {
      if (beta_1 > prior_mean_beta1) -Inf else dexp(prior_mean_beta1 - beta_1,
                                                    rate = 1, log = TRUE)
    } else {
      if (beta_1 < prior_mean_beta1) -Inf else dexp(beta_1 - prior_mean_beta1,
                                                    rate = 1, log = TRUE)
    }
  }
  log_prior_beta_2 <- if (prior_type == "normal") {
    dnorm(beta_2, mean = prior_mean_beta2, sd = 1, log = TRUE)
  } else {
    if (tau > 0.5) {
      if (beta_2 > prior_mean_beta2) -Inf else dexp(prior_mean_beta2 - beta_2,
                                                    rate = 1, log = TRUE)
    } else {
      if (beta_2 < prior_mean_beta2) -Inf else dexp(beta_2 - prior_mean_beta2,
                                                    rate = 1, log = TRUE)
    }
  }
  
  ### Inverse Gamma Prior for sigma
  log_prior_sigma <- dinvgamma(sigma, shape = 6, scale = 5, log = TRUE)
  
  # Combine log-likelihood and log-priors
  log_likelihood + log_prior_beta_0 + log_prior_beta_1 + log_prior_beta_2 + log_prior_sigma
}

### Adaptive Metropolis-Hastings Analysis with Beta2 included
adaptive_mh_analysis <- function(data, tau, burn_in = 5000, max_lag = 100,
                                 true_beta, n_iter = 10000,
                                 step_size_beta, step_size_sigma,
                                 prior_mean_beta0, prior_mean_beta1,
                                 prior_mean_beta2, prior_type) {
  y <- data$Y
  n <- length(y)
  p <- 4  # Number of parameters: 2 betas and 1 sigma
  
  ### Initialize parameters
  beta_0_current <- prior_mean_beta0
  beta_1_current <- prior_mean_beta1
  beta_2_current <- prior_mean_beta2
  sigma_current <- 1
  
  ### Storage for samples and acceptance rate tracking
  samples <- matrix(NA, nrow = n_iter, ncol = p)
  acceptance_count <- 0
  
  ### Start time
  start_time <- proc.time()
  
  ### RW-MH Stage
  for (i in 1:burn_in) {
    beta_0_prop <- beta_0_current + rnorm(1, sd = step_size_beta[1])
    beta_1_prop <- beta_1_current + rnorm(1, sd = step_size_beta[2])
    beta_2_prop <- beta_2_current + rnorm(1, sd = step_size_beta[3])
    sigma_prop  <- sigma_current + rnorm(1, sd = step_size_sigma)
    if (sigma_prop <= 0.01) sigma_prop <- 0.01
    
    log_acc_ratio <- log_posterior(y, data$X1, data$X2,
                                   beta_0_prop, beta_1_prop, beta_2_prop,
                                   sigma_prop, tau,
                                   prior_mean_beta0, prior_mean_beta1,
                                   prior_mean_beta2, prior_type) -
      log_posterior(y, data$X1, data$X2,
                    beta_0_current, beta_1_current,
                    beta_2_current, sigma_current,
                    tau, prior_mean_beta0, prior_mean_beta1,
                    prior_mean_beta2, prior_type)
    accept_prob <- min(1, exp(log_acc_ratio))
    
    if (runif(1) < accept_prob) {
      beta_0_current <- beta_0_prop
      beta_1_current <- beta_1_prop
      beta_2_current <- beta_2_prop
      sigma_current  <- sigma_prop
    }
    samples[i, ] <- c(beta_0_current, beta_1_current,
                      beta_2_current, sigma_current)
  }
  
  ### Calculate adaptive proposal parameters
  burn_samples <- samples[1:burn_in, ]
  mu_alpha   <- colMeans(burn_samples)
  omega_alpha <- cov(burn_samples)
  
  ### IK-MH Stage
  for (i in (burn_in + 1):n_iter) {
    proposal     <- mvrnorm(1, mu = mu_alpha, Sigma = omega_alpha)
    beta_0_prop <- proposal[1]
    beta_1_prop <- proposal[2]
    beta_2_prop <- proposal[3]
    sigma_prop  <- proposal[4]
    if (sigma_prop <= 0.5) sigma_prop <- 0.5
    
    log_acc_ratio <- log_posterior(y, data$X1, data$X2,
                                   beta_0_prop, beta_1_prop,
                                   beta_2_prop, sigma_prop,
                                   tau, prior_mean_beta0,
                                   prior_mean_beta1, prior_mean_beta2,
                                   prior_type) -
      log_posterior(y, data$X1, data$X2,
                    beta_0_current, beta_1_current,
                    beta_2_current, sigma_current,
                    tau, prior_mean_beta0,
                    prior_mean_beta1, prior_mean_beta2,
                    prior_type) +
      dmvnorm(c(beta_0_current, beta_1_current,
                beta_2_current, sigma_current),
              mean = mu_alpha, sigma = omega_alpha, log = TRUE) -
      dmvnorm(c(beta_0_prop, beta_1_prop,
                beta_2_prop, sigma_prop),
              mean = mu_alpha, sigma = omega_alpha, log = TRUE)
    accept_prob <- min(1, exp(log_acc_ratio))
    
    if (runif(1) < accept_prob) {
      beta_0_current <- beta_0_prop
      beta_1_current <- beta_1_prop
      beta_2_current <- beta_2_prop
      sigma_current  <- sigma_prop
      acceptance_count <- acceptance_count + 1
    }
    samples[i, ] <- c(beta_0_current, beta_1_current,
                      beta_2_current, sigma_current)
  }
  
  end_time <- proc.time() - start_time
  post_burn  <- samples[(burn_in + 1):n_iter, ]
  mcmc_betas <- post_burn[, 1:3]
  mcmc_sigma <- post_burn[, 4]
  
  ### Plotting Diagnostics
  overall_title <- paste("Error:", unique(data$Error),
                         "| Sample Size:", unique(data$Sample_Size),
                         "| Tau:", tau)
  par(mfrow = c(2, 3), oma = c(0, 0, 3, 0))
  plot_list <- c("Trace: Beta0", "Trace: Beta1", "Trace: Beta2",
                 "ACF: Beta0", "ACF: Beta1", "ACF: Beta2")
  
  for(j in 1:3) {
    plot(mcmc_betas[, j], type = "l", main = plot_list[j],
         ylab = "Value", xlab = "Iter")
    abline(h = true_beta[j], col = "blue")
    abline(h = mean(mcmc_betas[, j]), col = "red", lty = 2)
  }
  for(j in 1:3) {
    acf(mcmc_betas[, j], main = plot_list[j + 3])
  }
  title(overall_title, outer = TRUE)
  
  ### Compute Metrics
  biases  <- sapply(1:3, function(i) compute_bias(mcmc_betas[, i], true_beta[i]))
  rmses   <- sapply(1:3, function(i) compute_rmse(mcmc_betas[, i], true_beta[i]))
  mcse    <- apply(mcmc_betas, 2, compute_mcse)
  ineff   <- apply(mcmc_betas, 2, compute_inefficiency_factor_autocorr, max_lag)
  ess     <- sapply(1:3, function(i) compute_ess(nrow(mcmc_betas), ineff[i]))
  acc_rate<- acceptance_count / (n_iter - burn_in)
  
  results_df <- data.frame(
    Parameter       = paste0("Beta_", 0:2),
    beta_true       = true_beta,
    beta_estimate   = colMeans(mcmc_betas),
    Bias            = biases,
    RMSE            = rmses,
    MCSE            = mcse,
    ESS             = ess,
    IF              = ineff,
    Acceptance_Rate = rep(acc_rate, 3),
    Time            = rep(end_time["elapsed"], 3),
    Sample_Size     = unique(data$Sample_Size),
    Error_Type      = unique(data$Error),
    Quantile        = tau
  )
  
  return(list(results = results_df,
              mcmc_draws = list(betas = mcmc_betas, sigma = mcmc_sigma)))
}


### Define empty lists
all_aMCMC_results <- list()
all_aMCMC_draws <- list()

quantiles <- c(0.5, 0.25, 0.1, 0.75, 0.9)
set.seed(1234)
for (rep in 1:100) {
  print(paste("Running Bayesian Quantile Regression with Adaptive MCMC - Replication:", rep))
  
  # Start timing for this replication
  start_time <- Sys.time()
  
  data <- plot_data #replicated_data[[rep]]  # Retrieve dataset for this replication
  
  # Loop through each unique combination of Error Type and Sample Size
  for (error_type in unique(data$Error)) {
    for (sample_size in unique(data$Sample_Size)) {
      ### Define empty list to store estimates (Used for dynamic changing prior)
      previous_betas <- list()
      ### Subset data for combination of Error and Sample Size
      subset_data <- data %>%
        filter(Error == error_type, Sample_Size == sample_size)
      
      # Run Bayesian Quantile Regression for each quantile (tau)
      for (tau in quantiles) {
        print(paste("Error:", error_type,
                    "Sample Size:", sample_size,
                    "Quantile:", tau))
        
        ### Change prior
        if (tau == 0.5) {
          prior_type <- "normal"
          prior_mean_beta0 <- 0 
          prior_mean_beta1 <- 0
          # prior_mean_beta2 <- 0
        } else {
          previous_tau <- quantiles[match(tau, quantiles) - 1]
          if (tau == 0.75) {
            previous_tau <- 0.5
          }
          prior_mean_beta0 <- previous_betas[[as.character(previous_tau)]][1]
          prior_mean_beta1 <- previous_betas[[as.character(previous_tau)]][2]
          prior_mean_beta2 <- previous_betas[[as.character(previous_tau)]][3]
          prior_type <- "exponential"
        }
        
        ### True parameters
        true_parameter <- c(
          unique(subset(subset_data, Quantile == tau)$B0),
          unique(subset(subset_data, Quantile == tau)$B1),
          unique(subset(subset_data, Quantile == tau)$B2)
        )
        
        ### Call function with Beta2 included
        aMCMC_output <- adaptive_mh_analysis(
          subset_data, tau,
          true_beta = true_parameter,
          step_size_beta = c(0.4, 0.25, 0.25),
          step_size_sigma = 0.01,
          prior_mean_beta0 = prior_mean_beta0, 
          prior_mean_beta1 = prior_mean_beta1, 
          prior_mean_beta2 = prior_mean_beta2,
          prior_type = prior_type
        )
        
        ### Store results
        results_df <- aMCMC_output$results
        results_df$Error_Type <- error_type
        results_df$Sample_Size <- sample_size
        results_df$Quantile <- tau
        results_df$Rep <- rep
        
        results_df <- results_df %>%
          dplyr::select(Rep, Error_Type, Sample_Size, Quantile, everything())
        
        print(results_df)
        
        ### Store Estimation Results
        all_aMCMC_results <- append(all_aMCMC_results, list(results_df))
        
        ### Store MCMC draws
        all_aMCMC_draws[[paste0(
          "Rep_", rep,
          "_Err_", error_type,
          "_Size_", sample_size,
          "_Tau_", tau
        )]] <- aMCMC_output$mcmc_draws
        
        # Store current beta means (now Beta0, Beta1, Beta2)
        current_betas <- c(
          mean(aMCMC_output$mcmc_draws$betas[, 1]),
          mean(aMCMC_output$mcmc_draws$betas[, 2]),
          mean(aMCMC_output$mcmc_draws$betas[, 3])
        )
        previous_betas[[as.character(tau)]] <- current_betas
      }
    }
  }
  # End timing for this replication
  end_time <- Sys.time()
  
  # Calculate and print the running time for this replication
  running_time <- end_time - start_time
  print(paste("Running time for Replication", rep, ":", running_time))
}

# Initialize an empty list to store data frames
combined_results <- list()

# Loop through the first 100 elements of all_aMCMC_results
for (i in 1:100) {
  
  # Extract the data frame from the list
  df <- all_aMCMC_results[[i]]
  
  # Rename columns to different variable names
  colnames(df) <- c(
    "Replication", "Error", "Sample.Size", "Quantile", "Parameter",
    "True", "Estimate", "Bias", "RMSE",
    "MCSE", "ESS", "IF",
    "Accept.Rate", "Run.Time"
  )
  
  # Append to the list
  combined_results[[i]] <- df
}

# Combine all results into a single dataframe
amcmc_estimate_results <- do.call(rbind, combined_results)


# >--------------------------------------------------------------------------------------------------------------< #
####################################################################################################################