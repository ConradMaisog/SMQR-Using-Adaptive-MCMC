# SMQR Using Adaptive MCMC

Author: John Conrad Seg B. Maisog  
Project: Adaptive MCMC Implementation for Bayesian Quantile Regression  

## Overview

This R project implements Bayesian Quantile Regression using Adaptive MCMC and compares it against Frequentist Quantile Regression and Gibbs Quantile Regression. Simulated datasets are used to evaluate performance across varying quantiles, error distributions, and sample sizes. Diagnostic tools are included to assess the quality of MCMC estimation.

## Structure

The script is structured into five main sections:

1. Preliminaries  
   Installs and loads all required R packages using pacman::p_load().

2. Data Generation  
   Simulates data for different:
   - Sample sizes: 20, 40, 100  
   - Error types: Standard Normal, Right-Skewed Gamma, Heteroscedastic  
   - Quantiles: 0.1, 0.25, 0.5, 0.75, 0.9  
   - True coefficients for each quantile

3. Frequentist Quantile Regression  
   Runs standard quantile regression (rq()) on the generated data to serve as a baseline.

4. Bayesian Quantile Regression with Gibbs Sampling  
   Uses the bayesQR package to estimate quantile regression coefficients using Bayesian methods. For each combination of error type, sample size, and quantile:
   - MCMC draws are collected  
   - MCSE, inefficiency factors, and ESS are computed  
   - Traceplots and autocorrelation plots can optionally be displayed  
   - Output includes posterior estimates, diagnostics, and MCMC draws

5. Bayesian Quantile Regression with Adaptive MCMC  
   Implements a *Metropolis-Hastings algorithm with adaptive tuning* for Bayesian quantile regression under the *Asymmetric Laplace Distribution*.  
   - Custom priors are used depending on the quantile (symmetric at median, asymmetric elsewhere)  
   - Posterior summaries and convergence diagnostics are computed  
   - Includes visual tools and metrics similar to Gibbs output

## How to Run

⚠️ *IMPORTANT*: Run the script *line by line* from top to bottom.  
The code is structured to automatically install and load the required packages, define all functions, and prepare all necessary datasets and objects in sequence.

Please do not skip or rearrange code blocks when running the script.
