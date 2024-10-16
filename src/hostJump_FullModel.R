##############################################################################
##############################################################################
## R script to reproduce all main text and supplemental results and figures ##
## The main runs produce results using the exact analytic solution (when    ##
## possible), numerical integration, and simulation using C++ code.         ##
##############################################################################
##############################################################################

####################
## Load necessary ##
##    packages    ##
####################
# install.packages("ggplot2")
# install.packages("viridis")
# install.packages("tidyverse")
# install.packages("patchwork")
# install.packages("ggpubr")

library(ggplot2)
library(viridis)
library(tidyverse)
library(patchwork)
library(ggpubr)
####################

## set working directory to project root folder. In this case './HostJump_model'
setwd("~/Desktop/Repos/HostJump_Model/") ## ensure this path matches the local path to the project root folder

#################
#################
## Load  model ##
##  functions  ##
#################
#################

## function for probability density function of specified beta mixture distribution
## weights specify the (un)normalized weights associated with each beta distribution
## shape1 and shape2 vectors provide ordered parameters for the beta distributions
betaMix <- function(x, weights = c(), shape1 = c(), shape2 = c() ){
  
  ## ensure vector lengths match
  if(length(weights) != length(shape1) | length(shape2) != length(shape1)){
    print("Error: lengths of input parameter vectors do not match")
    print(paste0("weight: ", length(weights), ", shape1: ", length(shape1), ", shape2: ",  length(shape2)))
    return(0)
  }
  
  weights = weights/sum(weights) ## normalize weights
  dta = data.frame(phi = x, density = weights[1]*dbeta(x = x, shape1 = shape1[1], shape2 = shape2[1]))
  if(length(weights) > 1){
    for(i in 2:length(weights)){
      dta$density = dta$density + weights[i]*dbeta(x = x, shape1 = shape1[i], shape2 = shape2[i])
    }
  }
  
  return(dta)
  
}


## function to generate parameter file for beta mixture models in simulation
betaMix_data <- function(mixtureDistn){
  
  ## error checking 
  if(any(mixtureDistn$a < 0 | mixtureDistn$b < 0 )){ 
    print("error: beta parameters must be positive ")
    return(0)
  }
  if(any(mixtureDistn$weights < 0)){ 
    print("error: weights must all be positive")
    return(0)
  }
  
  ## normalizing weights
  if(sum(mixtureDistn$weights) != 1){
    print('normalizing distribution weights')
    mixtureDistn$weights <- mixtureDistn$weights / sum(mixtureDistn$weights)
  }
  
  ## convert to cumulative weights
  mixtureDistn$weights <- cumsum(mixtureDistn$weights)
  
  ## ensure path matches the location in the cpp code
  write.table(x = mixtureDistn, sep = ',', row.names = F, col.names = F, file = "./src/simulation_files/betaMix_pars.csv")
  
  return(0)
  
}


## generate full output file for model runs over all 3 priors used in the main text
## includes functionality to plot different past-future spillover relationships
## pars should be a data frame containing all model parameters
## mixtureDistn should be a data frame containing the parameters for the prior distribution
## vals should be a data frame of values for the past and future parameters to iterate over (i.e., N, M, lambda, *k*, *theta*)
fullModel_eval <- function(pars, mixtureDistn, vals){
  
  x_lab = c('N','lambda','k','theta')
  y_lab = c('M','lambda_f','k_f','theta_f')
  
  run_out <- vals
  run_out$prob = 0 # create simulation probability column
  run_out$sol = 0 # create analytic probability column
  run_out$exact = 0 # create exact probability column
  
  for(index in 1:dim(vals)[1]){
    
    ## set values being tested
    pars[c(x_lab[pars$past_model+1], y_lab[pars$future_model+1])] <- vals[index, c('past_val', 'future_val')]
    
    ## run simulation model
    if(pars$runSims == TRUE){
      out <- HJ_simulation(n_reps = pars$reps, parameters = pars, batch_name = 'test', output_mode = 0)
      run_out$prob[index] = out$total_prob[1] # store results -- the [1] prevents errors when long output is requested 
    }else{
      run_out$prob[index] = 0
    }
    
    ## evaluate analytic solution via pseudo- beta-Binomial
    val <- model_solution(past = pars$past_model+1, future = pars$future_model+1, pars = pars, a = mixtureDistn$a, b = mixtureDistn$b, weights = mixtureDistn$weights)
    run_out$sol[index] = val
    
    ## exact analytic solution via exp(lgamma()) trick -- only valid when N,M are constants
    ex_val <- model_solution_exact(past = pars$past_model, future = pars$future_model, pars = pars, a = mixtureDistn$a, b = mixtureDistn$b, weights = mixtureDistn$weights)
    run_out$exact[index] = ex_val
    
  }
  
  return(run_out)
  
}


## (mostly) Analytic model solution 
## results are exact in the beta-binomial case, but this function uses numerical integration in all cases
## all solutions are numerically correct for any N >= H >= 0, with the exception of past = 2
## which indicates the case where past spillover follows a gamma distributed rate
## this case is only exact when H = 0.
## More elaborate cases can be evaluated by simulation.
## UD_prior == data frame with columns "nodes" and "weights" to be used as the prior distribution
## nodes vector should contain only values on [0,1], representing host jump probabilities
## weights vector should have the normalized or un-normalized probability of the associated node
model_solution <- function(past, future, pars, a=NA, b=NA, weights = c(), UD_prior = NA, verbose = 0){
  
  value = -1
  
  if( !all(is.na(UD_prior)) & pars$prior_type == 1){ ## check specifications for user defined prior
    
    ##############################
    ##  functions to evaluate   ##
    ##  the re-scaling constant ##
    ## for the posterior of phi ##
    ##############################
    
    ## error checking 
    if(any(UD_prior$nodes > 1 | UD_prior$nodes < 0)){ 
      print("error: nodes must be on the interval [0,1]")
      return(0)
    }
    if(any(UD_prior$weights < 0)){ 
      print("error: weights must all be positive")
      return(0)
    }
    
    ## normalizing weights
    if(sum(UD_prior) != 1){
      if(verbose>0){ print('normalizing distribution weights') }
      UD_prior$weights <- UD_prior$weights / sum(UD_prior$weights)
    }
    
    
    ## normalizing constant in N spillover past with H host jumps
    norm_const1 <- function(par, UD_prior){ 
      N=as.numeric(par['N']); H=as.numeric(par['H_crit']);
      if(N==0){return(1)}
      return(1/sum( as.numeric( UD_prior$weights * choose(n=N, k=H)*(UD_prior$nodes^H)*(1-UD_prior$nodes)^(N-H) ) )) 
    }
    
    ## normalizing constant in poisson spillover past with H host jumps
    norm_const2 <- function(par, UD_prior){ 
      lambda=as.numeric(par['lambda']); H=as.numeric(par['H_crit']); C = 1;
      if(lambda==0){return(1)}
      if(H>0){ 
        C = (lambda^H)/(factorial(H))
        H_norm = seq(0, H, by=1)[-(H+1)] # all values of N that must have been observed
        C1 = 1 - sum( ((lambda^H_norm)*exp(-lambda))/(factorial(H_norm)) ) # normalizing constant for poisson number of spillovers
        C = C / C1
      }
      return(1/sum(as.numeric( C * UD_prior$weights * (exp(-lambda*UD_prior$nodes)) )))
    }
    
    ## normalizing constant in gamma-poisson spillover past with 0 host jumps
    ## for H != 0, the solution will not be exact. use simulation only for this case
    norm_const3 <- function(par, UD_prior){ 
      k=as.numeric(par['k']); theta=as.numeric(par['theta']); H=as.numeric(par['H_crit']);
      if(min(k,theta)==0){return(1)}
      if(H != 0){print('Warning: this integral is not exact for non-zero values of H. ')}
      return(1/sum(as.numeric( UD_prior$weights * (1+theta*UD_prior$nodes)^(-k) ))) 
    }
    
    ##############################
    switch(past,
           
           #########################################
           #########################################
           {
             const = norm_const1(par=pars, UD_prior=UD_prior)
             switch(future,
                    #########################################
                    #########################################
                    {
                      N=as.numeric(pars['N']); M=as.numeric(pars['M']); H=as.numeric(pars['H_crit']);
                      vals = ((1-UD_prior$nodes)^M) * UD_prior$weights * choose(n=N, k=H)*(UD_prior$nodes^H)*(1-UD_prior$nodes)^(N-H)
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    {
                      N=as.numeric(pars['N']); lambda_f=as.numeric(pars['lambda_f']); H=as.numeric(pars['H_crit']);
                      vals = exp(-lambda_f*UD_prior$nodes) * UD_prior$weights * choose(n=N, k=H)*(UD_prior$nodes^H)*(1-UD_prior$nodes)^(N-H)
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    {
                      N=as.numeric(pars['N']); k_f=as.numeric(pars['k_f']); theta_f=as.numeric(pars['theta_f']); H=as.numeric(pars['H_crit']);
                      vals = ((1+theta_f*UD_prior$nodes)^(-k_f)) * UD_prior$weights * choose(n=N, k=H)*(UD_prior$nodes^H)*(1-UD_prior$nodes)^(N-H)
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    stop('error: invalid future value')
             )
             
           },
           #########################################
           #########################################
           {
             const = norm_const2(par=pars, UD_prior=UD_prior)
             switch(future,
                    #########################################
                    #########################################
                    {
                      lambda=as.numeric(pars['lambda']); M=as.numeric(pars['M']); H=as.numeric(pars['H_crit']);
                      vals = ((1-UD_prior$nodes)^M) * UD_prior$weights * exp(-lambda*UD_prior$nodes)*(lambda^H / factorial(H))
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    {
                      lambda=as.numeric(pars['lambda']); lambda_f=as.numeric(pars['lambda_f']); H=as.numeric(pars['H_crit']);
                      vals = exp(-lambda_f*UD_prior$nodes) * UD_prior$weights * exp(-lambda*UD_prior$nodes)*(lambda^H / factorial(H))
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    {
                      lambda=as.numeric(pars['lambda']); k_f=as.numeric(pars['k_f']); theta_f=as.numeric(pars['theta_f']); H=as.numeric(pars['H_crit']);
                      vals = ((1+theta_f*UD_prior$nodes)^(-k_f)) * UD_prior$weights * exp(-lambda*UD_prior$nodes)*(lambda^H / factorial(H))
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    stop('error: invalid future value')
                    
             )
             
           },
           #########################################
           #########################################
           {
             const = norm_const3(par=pars, UD_prior=UD_prior)
             switch(future,
                    #########################################
                    #########################################
                    {
                      k=as.numeric(pars['k']); theta=as.numeric(pars['theta']); M=as.numeric(pars['M']); H=as.numeric(pars['H_crit']);
                      if(H>0){print("Warning: solution is not exact for non-zero values of H")}
                      vals = ((1-UD_prior$nodes)^M) * UD_prior$weights * ((1+theta*UD_prior$nodes)^(-k))
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    {
                      k=as.numeric(pars['k']); theta=as.numeric(pars['theta']); lambda_f=as.numeric(pars['lambda_f']); H=as.numeric(pars['H_crit']);
                      if(H>0){print("Warning: solution is not exact for non-zero values of H")}
                      vals = exp(-lambda_f*UD_prior$nodes) * UD_prior$weights * ((1+theta*UD_prior$nodes)^(-k))
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    {
                      k=as.numeric(pars['k']); theta=as.numeric(pars['theta']); k_f=as.numeric(pars['k_f']); theta_f=as.numeric(pars['theta_f']); H=as.numeric(pars['H_crit']);
                      if(H>0){print("Warning: solution is not exact for non-zero values of H")}
                      vals = ((1+theta_f*UD_prior$nodes)^(-k_f)) * UD_prior$weights * ((1+theta*UD_prior$nodes)^(-k))
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    stop('error: invalid future value')
                    
             )
             
           },
           #########################################
           #########################################
           stop('error: invalid past value')
    )
    
  }else{ ## default to continuous / beta prior
    
    ##############################
    ##  functions to evaluate   ##
    ##  the re-scaling constant ##
    ## for the posterior of phi ##
    ##############################
    
    # use a, b vector inputs for all beta-type priors -- pull values from parameters object if no vectors or NA values supplied
    if(all(is.na(a)) & all(is.na(b))){ 
      print("using beta prior parameters from 'pars'")
      a = as.numeric(pars$a); b = as.numeric(pars$b); weights = 1 
    }
    
    ## integrand for the normalizing constant in N spillover past with H host jumps
    integrand1 <- function(x, par, a, b, weights){ 
      N=as.numeric(par['N']); H=as.numeric(par['H_crit']);
      return(as.numeric( betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * choose(n=N, k=H)*(x^H)*(1-x)^(N-H) )) 
    }
    
    ## integrand for the normalizing constant in poisson spillover past with H host jumps
    integrand2 <- function(x, par, a, b, weights){ 
      lambda=as.numeric(par['lambda']); H=as.numeric(par['H_crit']); C = 1;
      if(H>0){ 
        C = (lambda^H)/(factorial(H))
        H_norm = seq(0, H, by=1)[-(H+1)] # all values of N that must have been observed
        C1 = 1 - sum( ((lambda^H_norm)*exp(-lambda))/(factorial(H_norm)) ) # normalizing constant for poisson number of spillovers
        C = C / C1
      }
      return(as.numeric( C * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * (exp(-lambda*x)) )) 
    }
    
    ## integrand for the normalizing constant in gamma-poisson spillover past with 0 host jumps
    ## for H != 0, the solution will not be exact. use simulation only for this case
    integrand3 <- function(x, par, a, b, weights){ 
      k=as.numeric(par['k']); theta=as.numeric(par['theta']); H=as.numeric(par['H_crit']);
      if(H != 0){print('Warning: this integral is not exact for non-zero values of H. ')}
      return(as.numeric( betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * (1+theta*x)^(-k) )) 
    }
    
    ##############################
    switch(past,
           
           #########################################
           #########################################
           {
             const = 1/integrate(f = integrand1, par = pars, a = a, b = b, weights = weights, lower = 0, upper = 1, abs.tol = 0, stop.on.error = F)$value
             switch(future,
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par, a, b, weights){
                        N=as.numeric(par['N'])
                        M=as.numeric(par['M'])
                        return( ((1-x)^M) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * ((1-x)^N) )
                      }
                      value = (1-const*integrate(f = integrand, par = pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par, a, b, weights){
                        N=as.numeric(par['N']);
                        lambda_f=as.numeric(par['lambda_f'])
                        return( (exp(-lambda_f*x)) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * ((1-x)^N) )
                      }
                      value = (1-const*integrate(f = integrand, par = pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par, a, b, weights){
                        N=as.numeric(par['N']);
                        k_f=as.numeric(par['k_f']); theta_f=as.numeric(par['theta_f'])
                        return( ((1+theta_f*x)^(-k_f)) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * ((1-x)^N) )
                      }
                      value = (1-const*integrate(f = integrand, par = pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    stop('error: invalid future value')
                    
             )
             
           },
           #########################################
           #########################################
           {
             const = 1/integrate(f = integrand2, par = pars, a = a, b = b, weights = weights, lower = 0, upper = 1, abs.tol = 0, stop.on.error = F)$value
             switch(future,
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par, a = a, b = b, weights = weights){
                        lambda=as.numeric(par['lambda']);
                        M=as.numeric(par['M']);
                        return( ((1-x)^M) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * (exp(-lambda*x)) )
                      }
                      value = (1-const*integrate(f=integrand, par=pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par, a = a, b = b, weights = weights){
                        lambda=as.numeric(par['lambda']);
                        lambda_f=as.numeric(par['lambda_f']);
                        return( (exp(-lambda_f*x)) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * (exp(-lambda*x)) )
                      }
                      value = (1-const*integrate(f=integrand, par=pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par,a = a, b = b, weights = weights){
                        lambda=as.numeric(par['lambda']);
                        k_f=as.numeric(par['k_f']); theta_f=as.numeric(par['theta_f']);
                        return( ((1+theta_f*x)^(-k_f)) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * (exp(-lambda*x)) )
                      }
                      value = (1-const*integrate(f=integrand, par=pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    stop('error: invalid future value')
                    
             )
             
           },
           #########################################
           #########################################
           {
             const = 1/integrate(f = integrand3, par = pars, a = a, b = b, weights = weights, lower = 0, upper = 1, abs.tol = 0, stop.on.error = F)$value
             switch(future,
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par, a = a, b = b, weights = weights){
                        k=as.numeric(par['k']); theta=as.numeric(par['theta']);
                        M=as.numeric(par['M']);
                        return( const*((1-x)^M) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * ((1+theta*x)^(-k)) )
                      }
                      value = (1-const*integrate(f=integrand, par=pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par, a = a, b = b, weights = weights){
                        k=as.numeric(par['k']); theta=as.numeric(par['theta']);
                        lambda_f=as.numeric(par['lambda_f']);
                        return( const*(exp(-lambda_f*x)) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * ((1+theta*x)^(-k)) )
                      }
                      value = (1-const*integrate(f=integrand, par=pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par, a = a, b = b, weights = weights){
                        k=as.numeric(par['k']); theta=as.numeric(par['theta']);
                        k_f=as.numeric(par['k_f']); theta_f=as.numeric(par['theta_f']);
                        return( ((1+theta_f*x)^(-k_f)) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * ((1+theta*x)^(-k)) )
                      }
                      value = (1-const*integrate(f=integrand, par=pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    stop('error: invalid future value')
                    
             )
             
           },
           #########################################
           #########################################
           stop('error: invalid past value')
    )
    
  }
  
  return(value)
  
}


## exact solution to the model using the gamma function
## these results are only valid when the number of past and future spillover events is fixed
## this is true for both single and mixtures of beta distributions
## using lgamma exp trick should increase stability of solutions for large values of N,M
model_solution_exact <- function(past, future, pars, a=NA, b=NA, weights = c(), verbose = 0){
  
  
  stirling <- function(x){
    return( sqrt(2*pi*(x-1))*exp(-(x-1))*(x-1)^(x-1) )
  }
  lstirling <- function(x){
    return( log( sqrt(2*pi*(x-1))*((x-1)/exp(1))^(x-1)  ) )
  }
  
  #pars$N=5.47e14
  #pars$M=5.47e14
  
  if(pars$past_model != 0 && pars$future_model != 0){
    # print("invalid past or future model type")
    return(NA)
  }else{
    
    C_N = sum( weights * exp( (lgamma(a + b)+lgamma(a+pars$H_crit)+lgamma(b+pars$N-pars$H_crit)) - (lgamma(a)+lgamma(b)+lgamma(a+b+pars$N-pars$H_crit)) ) )
    
    T_k = weights * exp( (lgamma(a+b)+lgamma(b+pars$N+pars$M)+lgamma(a+pars$H_crit)) - (lgamma(a)+lgamma(b)+lgamma(a+b+pars$N+pars$M)) )
    
    # C_N_ap = sum( weights * exp( (lgamma(a + b)+lgamma(a+pars$H_crit)+lstirling(x=b+pars$N-pars$H_crit)) - (lgamma(a)+lgamma(b)+lstirling(a+b+pars$N-pars$H_crit)) ) )
    # T_k_ap = weights * exp( (lgamma(a+b)+lstirling(b+pars$N+pars$M)+lgamma(a+pars$H_crit)) - (lgamma(a)+lgamma(b)+lstirling(a+b+pars$N+pars$M)) )
    
    lim = sum(T_k/C_N)
    
    if(lim>1 || lim<0){lim=NA}

    
    ## ensure exact zero if no spillovers occur in future
    if(pars$M == 0){ return(0) }else{ return(1 - lim) }
    
  }
  
}


## R wrapper function for C++ simulation code
## n_reps == number of simulations to run
## parameters == a named data frame object containing all model parameters
## batch_name == naming tag for any files generated in the run
## output_mode == numeric value to modulate simulation data output. "0" returns all
#### simulation data while "1" returns only the simulated HJ probability
HJ_simulation <- function(n_reps, parameters, batch_name, output_mode = 1){
  
  if(parameters$H_crit > parameters$N){ parameters$N = parameters$H_crit } ## fail safe to ensure simulation of past data is actually possible
  
  parvec = c( format(n_reps, scientific = F),
              
              parameters$past_model, # 0, switch indicator 0: fixed N; 1: poisson N; 2: gamma-poisson N;
              parameters$N, # 1, fixed past spillovers -- used if past_model == 0;
              parameters$lambda, # 2, constant spillover rate 
              parameters$k, # 3, gamma shape parameter
              parameters$theta, # 4, gamma scale parameter 
              
              parameters$future_model, # 5, switch indicator 0: fixed M; 1: poisson M; 2: gamma-poisson M;
              parameters$M, # 6, fixed future spillovers -- used if future_model == 0;
              parameters$lambda_f, # 7, constant future spillover rate 
              parameters$k_f, # 8, gamma future shape parameter
              parameters$theta_f, # 9, gamma future scale parameter 
              
              parameters$redraw, # 10, redraw future spillover rate 0: no; 1: yes;
              parameters$verbose, # 11, verbose output level -- not currently used
              parameters$std_out, # 12, output type (0 - file; 1 - console) 
              output_mode, # 13, results output format (0 - full simulation, 1 - minimal)
              batch_name, # 14, naming tag for results file output
              
              parameters$prior_type, # 15, prior type (0 - beta, 1 - discrete, 2 - beta-mixture)
              parameters$a, # 16, beta shape1 prior parameter
              parameters$b, # 17, beta shape2 prior parameter
              parameters$H_crit # 18, number of permitted past host jumps. ideally 0
  )
  strvec = format(parvec, digits = 5)
  
  setwd("~/Desktop/Repos/HostJump_Model/src") ## call has to be from location of .exe file or parameter read-in fails???
  
  ## Run the model
  ## The path to the bTB cpp binary file must be set correctly in the sys call below:
  nm = paste0("./HJ_simulation.exe")
  r <- system2( nm, args = strvec, stdout = TRUE)
  
  ## capture model output of the simulation
  if(parameters$std_out == 0){
    ## load data from file
    out <- read.table(file = paste0("./simulation_files/", batch_name, "_fullSim.txt"), header = TRUE, sep = ';', check.names = FALSE) %>% mutate_all(as.numeric)
  }else{
    ## read output from console
    out <- read.table(text = r, header = TRUE, sep = ';', check.names = FALSE) %>% mutate_all(as.numeric)
  }
  
  setwd("..")
  
  return( out ) 
}


## derived solution to limit value
## provides exact limit when min(a) is unique, and provides an upper bound otherwise
## value is only exact when N,M are constant
lim_val <- function(pars, C=1, a=NA, b=NA, weights = c()){
  
  
  upr = 1
  lwr = 0  
  j_star = which(a %in% min(a))
  if(length(j_star) == 1){
    ## compute exact limit
    lwr = 1 - ( 1 / (1-C)^(min(a)) )
    upr = 1 - ( 1 / (length(a)*(1-C)^(min(a))) )
    
  }else{
    ## compute upper bound for the limit
    omega = which(b[j_star] == max(b[j_star]))
    
    lwr = 1 - ( 1 / (1-C)^(min(a)) )
    upr = 1 - ( 1 / (length(a)*(1-C)^(min(a))) )
    
  }
  
}

##########################
## Simulation functions ##
##   used to generate   ##
##  conceptual figures  ##
##########################

## vectorized funtion simulating time to host jump
T_HJ_sim <- function(par, reps) {
  N <- nrow(par)  # Number of parameter combinations
  results <- numeric(N)  # Initialize a vector to store results for each combination
  
  replicate_results <- function(phi, lambda) {
    HJ <- replicate(n = reps, expr = FALSE) 
    T_HJ <- replicate(n = reps, expr = 0) 
    rem <- reps
    
    if(phi == 0 | lambda == 0){ rem = 0; T_HJ = NA }
    
    while (rem > 0) {
      t <- rexp(rem, rate = lambda)
      T_HJ[HJ == FALSE] <- T_HJ[HJ == FALSE] + t
      HJ[HJ == FALSE] <- runif(n = rem) < phi
      rem <- sum(HJ == FALSE)
    }
    
    return(mean(T_HJ))
  }
  
  # Use mapply to apply replicate_results to each row of the data frame
  results <- mapply(replicate_results, par$phi, par$lambda)
  
  return(results)
}

## function to simulate wait times until the first host jump
T_HJ_dist <- function(par, reps){
  phi = par$phi
  lambda = par$lambda
  HJ = replicate(n = reps, expr = F)
  T_HJ = replicate(n = reps, expr = 0)
  rem = reps
  
  while(rem > 0){
    t = rexp(rem, rate = lambda) #time to next spillover
    T_HJ[HJ == F] = T_HJ[HJ == F] + t # update cumulative time to first host jump
    HJ[HJ == F] = runif(n = rem) < phi #indicator if HJ occured at this spillover
    rem = sum(HJ == F) # update number of remaining replicates
  }
  return(T_HJ)
}

##########################

######################################
######################################
#### Conceptual figure generation ####
######################################
######################################

lambda = seq(0,5, by = 0.1)[-1]
phi = seq(0,1, by = 0.01)[-1]
dat = expand.grid(lambda, phi)
names(dat) = c('lambda', 'phi')
dat$time = 1 / (dat$phi * dat$lambda)
dat$time_sim = T_HJ_sim(par = dat, reps = 10000)
dat$time[is.nan(dat$time)] = NA
# dat$time[dat$time < 10] = NA

## plot expected time to first host jump on log scale

p_log <- ggplot() +
  scale_fill_viridis(option = "B", discrete = F, direction = -1) +
  geom_tile(data = dat, aes(x = lambda, y = phi, fill = log(time, base = 10))) +
  geom_contour(data = dat, aes(x = lambda, y = phi, z = log(time_sim, base = 10)), color= '#FFFFFF', alpha = 0.55) +
  geom_contour(data = dat, aes(x = lambda, y = phi, z = log(time, base = 10)), color= '#041E42', alpha = 1, linetype = 'dashed') +
  #ggtitle("Expected Time to Host Jump") +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c('0','0.25','0.50','0.75','1')) +
  xlab("Spillover Rate") +
  #ylab(expression(phi)) +
  labs(fill = "log(Time)") + 
  coord_cartesian(expand = FALSE) +
  theme_bw() + 
  theme(plot.title = element_text(size = 23),
        #axis.title.x = element_text(size = 14, hjust = 0.12),
        #axis.title.y = element_text(size = 20, angle = 0, vjust = 0.8),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18, angle = 0, vjust = 0.54),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.key.width = unit(1.55, "cm"),
        aspect.ratio = 1.0)
p_log


## frequency/probability remaining plots
lambda = seq(0,5, by = 0.1)
phi = seq(0,1, by = 0.01)
T_lim = c(0,5,20,100) 
dat2 = expand.grid(lambda, phi, T_lim)
names(dat2) = c('lambda', 'phi', 'T_lim')
dat2$prob = 1-(1-exp(-dat2$lambda*dat2$phi*dat2$T_lim))

## fraction of pathogens remaining in zoonotic pool after T units of time
p2 <- ggplot() +
  scale_fill_viridis(option = "B", discrete = F, direction = 1) +
  geom_tile(data = dat2, aes(x = lambda, y = phi, fill = prob)) +
  geom_contour(data = dat2[dat2$T_lim > 0, ], aes(x = lambda, y = phi, z = prob), color = "#FFFFFF", breaks = seq(0.2,0.8,0.2), alpha = 0.55) +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c('0','0.25','0.50','0.75','1')) +
  xlab("Spillover Rate") +
  ylab(expression(phi)) +
  facet_wrap( ~ T_lim, nrow = 2 ) +
  coord_cartesian(expand = FALSE) +
  theme_bw() + 
  theme(
        legend.position = "bottom",
        legend.key.width = unit(1.55, "cm"),
        aspect.ratio = 1.0)
p2
# ggsave('p2.png', p2, bg = 'transparent')

ggarrange(p_log, p2)

######################################
######################################


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
## results for the Poisson version of the model can be obtained by setting
#### parameters "past_model" and "future_model" equal to 1. note that doing so
#### will make the "exact" model results invalid, and numerical integration must be used


####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 0, N = 1, lambda = 1, k = 1, theta = 1,
                  future_model = 0, M = 1, lambda_f = 1, k_f = 1, theta_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

spill_max = 20 # max mean spillover value
spill_step = 1 # step size -- must be an integer for N-M models
long_max = 1000 # maximum value for C:1 plots
long_step = 5 # step size for long C:1 line values
C = 1 # past:future spillover ratio -- slope of C:1 line

vals = expand.grid(seq(0, spill_max, by = spill_step), seq(0, spill_max, by = spill_step))
names(vals) = c('past_val', 'future_val')
vals_long = data.frame(seq(0, long_max, by = long_step), C*seq(0, long_max, by = long_step))
names(vals_long) = c('past_val', 'future_val')
vals = rbind(vals, vals_long); rm(vals_long)
## sum(duplicated(vals)) ## test count of duplicated rows
## vals[duplicated(vals),] ## test view duplicated rows
vals <- vals[!duplicated(vals),] ## removes rows that may be duplicated by merging vals_long

####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = 0.09,
                            b =  0.9,
                            weights = c(1) )

phi_vec = seq(from = 0, to = 1, by = 0.001)
prior = data.frame(phi = phi_vec, density = betaMix(x = phi_vec, weights = mixtureDistn$weights, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b )$density, 
                   prior = paste0("Scenario 1"))

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 1
run_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(run_out) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
run_out$scenario = "Scenario 1"

####################

####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(1/100),
                            b =  c(10),
                            weights = c(1) )

tmp = data.frame(phi = phi_vec, density = betaMix(x = phi_vec, weights = mixtureDistn$weights, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b )$density, 
                   prior = paste0("Scenario 2"))
prior = rbind(tmp, prior)

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 2 and merge results
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 2"
run_out <- rbind(run_out, tmp); rm(tmp)

####################


####################
## Define prior 3 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1, 100),
                            b =  c(10, 400),
                            weights = c(0.8, 0.2) )
tmp = data.frame(phi = phi_vec, density = betaMix(x = phi_vec, weights = mixtureDistn$weights, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b )$density, 
                 prior = paste0("Scenario 3"))
prior = rbind(tmp, prior)

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 3 and merge results
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 3"
run_out <- rbind(run_out, tmp); rm(tmp)

####################

sum(abs(run_out$sol-run_out$exact)) ## total error between numerical and exact solution over all scenarios
max(abs(run_out$sol-run_out$exact)) ## max point of error
which(abs(run_out$sol-run_out$exact) == max(abs(run_out$sol-run_out$exact)))

###############################
## plot prior  distributions ##
##  for the three scenarios  ##
###############################

prior_plots <- ggplot() +
  geom_line( data = prior, aes(x = phi, y = density), color = viridis(n=6)[3], linewidth = 1.3, alpha = 1 ) +
  #labs( x = paste0('Host Jump Probability (', expression(phi), ')'), y = 'density' ) + 
  labs( x = expression(phi), y = 'Density' ) + 
  facet_wrap( ~ prior) +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c('0','0.25','0.5','0.75','1')) +
  theme_classic() +
  theme(plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 10),
        aspect.ratio = 1.0)

prior_plots

###############################



###################
## Generate main ##
## results plots ##
###################
brk = seq(from = 0, to = max(max(run_out$exact), max(run_out$prob), max(run_out$sol)), by = 0.05)[-1] # contour value breakpoints
rng = 0 # range around the C:1 line as plot points
cond_l = C*run_out$future_val == run_out$past_val ## C:1 line -- i.e. future horizon is C times past horizon
cond = (C*run_out$future_val <= (run_out$past_val + rng)) & (C*run_out$future_val >= (run_out$past_val - rng)) ## range of values around C:1 line
cond_short = (C*run_out$future_val == run_out$past_val) & (run_out$past_val <= spill_max)  ## C:1 line values from heatmap

savePlots = FALSE ## option to save plots to ./plots in directory. Note that this will overwrite existing plots from other runs

## generate heatmap plots for exact solution and simulation results
## beyond the N-M model, the solution is obtained via numerical integration, 
## as the exact solution will not be valid in these cases
## 

{
  tmp = run_out[run_out$scenario %in% 'Scenario 1',] ## subset by scenario
  tmp = tmp[(tmp$past_val <= spill_max) & (tmp$future_val <= spill_max),] ## exclude long values
  cond_tmp = tmp$future_val == tmp$past_val ## C:1 line -- i.e. future horizon is C times past horizon
  
  int_map1 <- ggplot() +
    scale_fill_viridis(option = "plasma", discrete = F, breaks = c(0, 0.05, 0.1), labels = function(x) round(x, digits = 2)) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['sol']])) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('sol')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('exact')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_line(data = tmp[cond_tmp,], aes(x = past_val, y = future_val), color = "#FFFFFF", alpha = .6, linewidth = 1.25) +
    scale_x_continuous(breaks = seq(min(tmp$future_val), max(tmp$past_val), 5)) +
    scale_y_continuous(breaks = seq(from = min(tmp$future_val), to = max(tmp$future_val), by = 5)) +
    xlab(" ") +
    labs( fill = expression( paste(P(H[F]>0)) ) ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(title.vjust = 1.25)) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          plot.margin = unit(c(0.5,0,0,0), 'lines'),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 28),
          axis.text.x = element_text(size = 28),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 28),
          legend.title = element_text(size = 24),
          legend.position = "bottom",
          legend.title.align = 0.5,
          legend.text = element_text(size = 24),
          legend.key.width = unit(1.1, "cm"),
          aspect.ratio = 1.0)
  
  
  tmp = run_out[run_out$scenario %in% 'Scenario 2',] ## subset by scenario
  tmp = tmp[(tmp$past_val <= spill_max) & (tmp$future_val <= spill_max),] ## exclude long values
  cond_tmp = tmp$future_val == tmp$past_val ## C:1 line -- i.e. future horizon is C times past horizon
  int_map2 <- ggplot() +
    scale_fill_viridis(option = "plasma", discrete = F, breaks = c(0, 0.3, 0.6), labels = function(x) round(x, digits = 2)) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['sol']])) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('sol')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('exact')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_line(data = tmp[cond_tmp,], aes(x = past_val, y = future_val), color = "#FFFFFF", alpha = .6, linewidth = 1.25) +
    scale_x_continuous(breaks = seq(min(tmp$future_val), max(tmp$past_val), 5)) +
    scale_y_continuous(breaks = seq(from = min(tmp$future_val), to = max(tmp$future_val), by = 5)) +
    xlab(" ") +
    labs( fill = expression( paste(P(H[F]>0)) ) ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(title.vjust = 1.25)) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          plot.margin = unit(c(0.5,0,0,0), 'lines'),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 28),
          axis.text.x = element_text(size = 28),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 28, color = "#FFFFFF"),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 24),
          legend.position = "bottom",
          legend.title.align = 0.5,
          legend.text = element_text(size = 24),
          legend.key.width = unit(1.1, "cm"),
          aspect.ratio = 1.0)
  
  
  tmp = run_out[run_out$scenario %in% 'Scenario 3',] ## subset by scenario
  tmp = tmp[(tmp$past_val <= spill_max) & (tmp$future_val <= spill_max),] ## exclude long values
  cond_tmp = tmp$future_val == tmp$past_val*C ## C:1 line -- i.e. future horizon is C times past horizon
  int_map3 <- ggplot() +
    scale_fill_viridis(option = "plasma", discrete = F, breaks = c(0, 0.1, 0.2), labels = function(x) round(x, digits = 2)) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['sol']])) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('sol')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('exact')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_line(data = tmp[cond_tmp,], aes(x = past_val, y = future_val), color = "#FFFFFF", alpha = .6, linewidth = 1.25) +
    scale_x_continuous(breaks = seq(min(tmp$future_val), max(tmp$past_val), 5)) +
    scale_y_continuous(breaks = seq(from = min(tmp$future_val), to = max(tmp$future_val), by = 5)) +
    xlab(" ") +
    labs( fill = expression( paste(P(H[F]>0)) ) ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(title.vjust = 1.25)) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          plot.margin = unit(c(0.5,0,0,0), 'lines'),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 28),
          axis.text.x = element_text(size = 28),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 28, color = "#FFFFFF"),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 24),
          legend.position = "bottom",
          legend.title.align = 0.5,
          legend.text = element_text(size = 24),
          legend.key.width = unit(1.2, "cm"),
          aspect.ratio = 1.0)
}

int_map1

int_map2

int_map3

if(savePlots){
  ggsave(filename = paste("./figures/int_map1_", fname_tag, ".png"), plot = int_map1, width = 5.85, height = 5.85, units = "in")
  ggsave(filename = paste("./figures/int_map2_", fname_tag, ".png"), plot = int_map2, width = 5.85, height = 5.85, units = "in")
  ggsave(filename = paste("./figures/int_map3_", fname_tag, ".png"), plot = int_map3, width = 5.85, height = 5.85, units = "in")
}

## generate 1:1 plots 

if(pars$runSims){ fname_tagL = paste0(fname_tag, '_Sims') }else{ fname_tagL = paste0(fname_tag, '_noSims') } ## change filename tag for simulation values

{
  y_upr = max( c(run_out[cond_l,]$prob,run_out[cond_l,]$sol, run_out[cond_l,]$prob, 1 - (1/(C+1)^(pars$a))))*1.01 
  
  L_curve1 <- ggplot() +
    geom_line(data = run_out[cond_l,][run_out[cond_l,]$scenario == 'Scenario 1',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_l,][run_out[cond_l,]$scenario == 'Scenario 1',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, long_max, by=250) ) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +
    xlab("Spillover Rate") +
    ylab( expression( paste(P(H[F]>0)) ) ) +
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28),
      aspect.ratio = 1.0)
  
  curve1 <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == 'Scenario 1',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_short & run_out$scenario == 'Scenario 1',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, spill_max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    xlab("Spillover Rate") +
    ylab( expression( paste(P(H[F]>0)) ) ) +
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28),
      aspect.ratio = 1.0)
  
  
  L_curve2 <- ggplot() +
    geom_line(data = run_out[cond_l,][run_out[cond_l,]$scenario == 'Scenario 2',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_l,][run_out[cond_l,]$scenario == 'Scenario 2',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,long_max,by=250) ) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    xlab("Inherent rate of spillover (N = M)") +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      # axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.text.y = element_text(size = 28),
      axis.ticks.y = element_blank(),
      aspect.ratio = 1.0)
  
  
  curve2 <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == 'Scenario 2',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_short & run_out$scenario == 'Scenario 2',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, spill_max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    xlab("Inherent rate of spillover (N = M)") +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.ticks.y = element_blank(),
      aspect.ratio = 1.0)
 
  L_curve3 <- ggplot() +
    geom_line(data = run_out[cond_l,][run_out[cond_l,]$scenario == 'Scenario 3',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_l,][run_out[cond_l,]$scenario == 'Scenario 3',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,long_max,by=250) ) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      # axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.text.y = element_text(size = 28),
      axis.ticks.y = element_blank(),
      aspect.ratio = 1.0)
 
  
  curve3 <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == 'Scenario 3',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_short & run_out$scenario == 'Scenario 3',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, spill_max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.ticks.y = element_blank(),
      aspect.ratio = 1.0)
}

L_curve1
curve1

L_curve2
curve2

L_curve3
curve3

L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_l,], aes(x = past_val, y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
  geom_point(data = run_out[cond_l,], aes(x = past_val, y = .data[['sol']]), color = '#0C2340', size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Inherent rate of spillover') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)

curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
  geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['sol']]), color = '#0C2340', size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Inherent rate of spillover') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)

L_curves
curves

if(savePlots){
  ggsave(filename = paste("./figures/L_curve1_", fname_tagL, ".png"), plot = L_curve1, width = 3.85, height = 3.85, units = "in")
  ggsave(filename = paste("./figures/curve1_", fname_tagL, ".png"), plot = curve1, width = 5.85, height = 5.85, units = "in")
  
  ggsave(filename = paste("./figures/L_curve2_", fname_tagL, ".png"), plot = L_curve2, width = 3.85, height = 3.85, units = "in")
  ggsave(filename = paste("./figures/curve2_", fname_tagL, ".png"), plot = curve2, width = 5.85, height = 5.85, units = "in")
  
  ggsave(filename = paste("./figures/L_curve3_", fname_tagL, ".png"), plot = L_curve3, width = 3.85, height = 3.85, units = "in")
  ggsave(filename = paste("./figures/curve3_", fname_tagL, ".png"), plot = curve3, width = 5.85, height = 5.85, units = "in")
  
  ggsave(filename = paste("./figures/L_curves_", fname_tagL, ".png"), plot = L_curves, width = 5.85, height = 5.85, units = "in")
  ggsave(filename = paste("./figures/curves_", fname_tagL, ".png"), plot = curves, width = 10.00, height = 5.50, units = "in")
}

###################




##########################################################################################
##########################################################################################
## Simulation Runs #######################################################################
##########################################################################################
##########################################################################################

## results for the Poisson version of the model can be obtained by setting
#### parameters "past_model" and "future_model" equal to 1

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 0, N = 1, lambda = 1, k = 1, theta = 1,
                  future_model = 0, M = 1, lambda_f = 1, k_f = 1, theta_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = TRUE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

spill_max = 20 # max mean spillover value
spill_step = 1 # step size -- must be an integer for N-M models
long_max = 500 # maximum value for C:1 plots
long_step = 5 # step size for long C:1 line values
C = 1 # past:future spillover ratio -- slope of C:1 line

vals = expand.grid(seq(0, spill_max, by = spill_step), seq(0, spill_max, by = spill_step))
names(vals) = c('past_val', 'future_val')
vals_long = data.frame(seq(0, long_max, by = long_step), C*seq(0, long_max, by = long_step))
names(vals_long) = c('past_val', 'future_val')
vals = rbind(vals, vals_long); rm(vals_long)
## sum(duplicated(vals)) ## test count of duplicated rows
## vals[duplicated(vals),] ## test view duplicated rows
vals <- vals[!duplicated(vals),] ## removes rows that may be duplicated by merging vals_long

####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(10),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 1
run_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
run_out$scenario = "Scenario 1"

####################

####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(0.1),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 2 and merge results
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
tmp$scenario = "Scenario 2"
run_out <- rbind(run_out, tmp); rm(tmp)

####################


####################
## Define prior 3 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1, 100),
                            b =  c(10, 400),
                            weights = c(0.8, 0.2) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 3 and merge results
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
tmp$scenario = "Scenario 3"
run_out <- rbind(run_out, tmp); rm(tmp)

####################

sum(abs(run_out$sol-run_out$exact)) ## total error between numerical and exact solution over all scenarios
max(abs(run_out$sol-run_out$exact)) ## max point of error
which(abs(run_out$sol-run_out$exact) == max(abs(run_out$sol-run_out$exact)))

###################
## Generate main ##
## results plots ##
###################
brk = seq(from = 0, to = max(max(run_out$exact), max(run_out$prob), max(run_out$sol)), by = 0.05)[-1] # contour value breakpoints
rng = 0 # range around the C:1 line as plot points
cond_l = C*run_out$future_val == run_out$past_val ## C:1 line -- i.e. future horizon is C times past horizon
cond = (C*run_out$future_val <= (run_out$past_val + rng)) & (C*run_out$future_val >= (run_out$past_val - rng)) ## range of values around C:1 line
cond_short = (C*run_out$future_val == run_out$past_val) & (run_out$future_val <= spill_max) & (run_out$past_val <= spill_max)  ## C:1 line values from heatmap

savePlots = FALSE ## option to save plots to ./plots in directory. Note that this will overwrite existing plots from other runs

## generate heatmap plots for exact solution and simulation results
## beyond the N-M model, the solution is obtained via numerical integration, 
## as the exact solution will not be valid in these cases
## 

{
  tmp = run_out[run_out$scenario %in% 'Scenario 1',] ## subset by scenario
  tmp = tmp[(tmp$past_val <= spill_max) & (tmp$future_val <= spill_max),] ## exclude long values
  cond_tmp = tmp$future_val == tmp$past_val ## C:1 line -- i.e. future horizon is C times past horizon
  sim_map1 <- ggplot() +
    scale_fill_viridis(option = "plasma", discrete = F, breaks = c(0, 0.05, 0.1), labels = function(x) round(x, digits = 2)) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['prob']])) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('prob')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('sol')]]), color= '#FFFFFF', breaks = brk, alpha = 0.6, linewidth = 1.25, linetype = 'dashed') +
    geom_line(data = tmp[cond_tmp,], aes(x = past_val, y = future_val), color = "#FFFFFF", alpha = .6, linewidth = 1.25) +
    scale_x_continuous(breaks = seq(min(tmp$future_val), max(tmp$past_val), 5)) +
    scale_y_continuous(breaks = seq(from = min(tmp$future_val), to = max(tmp$future_val), by = 5)) +
    xlab(" ") +
    labs( fill = expression( paste(P(H[F]>0)) ) ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(title.vjust = 1.25)) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          plot.margin = unit(c(0.5,0,0,0), 'lines'),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 28),
          axis.text.x = element_text(size = 28),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 28),
          legend.title = element_text(size = 24),
          legend.position = "bottom",
          legend.title.align = 0.5,
          legend.text = element_text(size = 24),
          legend.key.width = unit(1.1, "cm"),
          aspect.ratio = 1.0)
  
  
  tmp = run_out[run_out$scenario %in% 'Scenario 2',] ## subset by scenario
  tmp = tmp[(tmp$past_val <= spill_max) & (tmp$future_val <= spill_max),] ## exclude long values
  cond_tmp = tmp$future_val == tmp$past_val ## C:1 line -- i.e. future horizon is C times past horizon
  sim_map2 <- ggplot() +
    scale_fill_viridis(option = "plasma", discrete = F, breaks = c(0, 0.3, 0.6), labels = function(x) round(x, digits = 2)) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['prob']])) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('prob')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('sol')]]), color= '#FFFFFF', breaks = brk, alpha = 0.6, linewidth = 1.25, linetype = 'dashed') +
    geom_line(data = tmp[cond_tmp,], aes(x = past_val, y = future_val), color = "#FFFFFF", alpha = .6, linewidth = 1.25) +
    scale_x_continuous(breaks = seq(min(tmp$future_val), max(tmp$past_val), 5)) +
    scale_y_continuous(breaks = seq(from = min(tmp$future_val), to = max(tmp$future_val), by = 5)) +
    xlab(" ") +
    labs( fill = expression( paste(P(H[F]>0)) ) ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(title.vjust = 1.25)) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          plot.margin = unit(c(0.5,0,0,0), 'lines'),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 28),
          axis.text.x = element_text(size = 28),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 28, color = "#FFFFFF"),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 24),
          legend.position = "bottom",
          legend.title.align = 0.5,
          legend.text = element_text(size = 24),
          legend.key.width = unit(1.1, "cm"),
          aspect.ratio = 1)
  
  
  tmp = run_out[run_out$scenario %in% 'Scenario 3',] ## subset by scenario
  tmp = tmp[(tmp$past_val <= spill_max) & (tmp$future_val <= spill_max),] ## exclude long values
  cond_tmp = tmp$future_val == tmp$past_val*C ## C:1 line -- i.e. future horizon is C times past horizon
  sim_map3 <- ggplot() +
    scale_fill_viridis(option = "plasma", discrete = F, breaks = c(0, 0.1, 0.2), labels = function(x) round(x, digits = 2)) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['prob']])) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('prob')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('sol')]]), color= '#FFFFFF', breaks = brk, alpha = 0.6, linewidth = 1.25, linetype = 'dashed') +
    geom_line(data = tmp[cond_tmp,], aes(x = past_val, y = future_val), color = "#FFFFFF", alpha = .6, linewidth = 1.25) +
    scale_x_continuous(breaks = seq(min(tmp$future_val), max(tmp$past_val), 5)) +
    scale_y_continuous(breaks = seq(from = min(tmp$future_val), to = max(tmp$future_val), by = 5)) +
    xlab(" ") +
    labs( fill = expression( paste(P(H[F]>0)) ) ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(title.vjust = 1.25)) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          plot.margin = unit(c(0.5,0,0,0), 'lines'),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 28),
          axis.text.x = element_text(size = 28),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 28, color = "#FFFFFF"),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 24),
          legend.position = "bottom",
          legend.title.align = 0.5,
          legend.text = element_text(size = 24),
          legend.key.width = unit(1.2, "cm"),
          aspect.ratio = 1.0)
}

sim_map1

sim_map2

sim_map3

if(savePlots){
  ggsave(filename = paste("./figures/int_map1_", fname_tag, ".png"), plot = int_map1, width = 5.85, height = 5.85, units = "in")
  ggsave(filename = paste("./figures/int_map2_", fname_tag, ".png"), plot = int_map2, width = 5.85, height = 5.85, units = "in")
  ggsave(filename = paste("./figures/int_map3_", fname_tag, ".png"), plot = int_map3, width = 5.85, height = 5.85, units = "in")
}

## generate 1:1 plots 

if(pars$runSims){ fname_tagL = paste0(fname_tag, '_Sims') }else{ fname_tagL = paste0(fname_tag, '_noSims') } ## change filename tag for simulation values

L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_l,], aes(x = past_val, y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Inherent rate of spillover') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)

if(pars$runSims){ L_curves <- L_curves + geom_point(data = run_out[cond_l,], aes(x = past_val, y = .data[['prob']]),  color = '#0C2340', size = 2.5)
}else{  L_curves <- L_curves + geom_point(data = run_out[cond_l,], aes(x = past_val, y = .data[['sol']]), color = '#0C2340', size = 2.5) }



curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Inherent rate of spillover') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)

if(pars$runSims){ curves <- curves + geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['prob']]),  color = '#0C2340', size = 2.5)
}else{  curves <- curves + geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['sol']]), color = '#0C2340', size = 2.5) }

L_curves
curves

if(savePlots){
  ggsave(filename = paste("./figures/L_curve1_", fname_tagL, ".png"), plot = L_curve1, width = 3.85, height = 3.85, units = "in")
  ggsave(filename = paste("./figures/curve1_", fname_tagL, ".png"), plot = curve1, width = 5.85, height = 5.85, units = "in")
  
  ggsave(filename = paste("./figures/L_curve2_", fname_tagL, ".png"), plot = L_curve2, width = 3.85, height = 3.85, units = "in")
  ggsave(filename = paste("./figures/curve2_", fname_tagL, ".png"), plot = curve2, width = 5.85, height = 5.85, units = "in")
  
  ggsave(filename = paste("./figures/L_curve3_", fname_tagL, ".png"), plot = L_curve3, width = 3.85, height = 3.85, units = "in")
  ggsave(filename = paste("./figures/curve3_", fname_tagL, ".png"), plot = curve3, width = 5.85, height = 5.85, units = "in")
  
  ggsave(filename = paste("./figures/L_curves_", fname_tagL, ".png"), plot = L_curves, width = 5.85, height = 5.85, units = "in")
  ggsave(filename = paste("./figures/curves_", fname_tagL, ".png"), plot = curves, width = 10.00, height = 5.50, units = "in")
}

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################




################################################################################
################################################################################
####     Runs and plots for different  past-future spillover slopes (c)     ####
################################################################################
################################################################################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 0, N = 1, lambda = 1, k = 1, theta = 1,
                  future_model = 0, M = 1, lambda_f = 1, k_f = 1, theta_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

long_max = 1000 # maximum value for 1:1 plots
long_step = 2 # step size for long C:1 line values
C = c(0.2, 0.5, 2, 5)
past_vals = c(0:19, seq(20,100,by=2), seq(100,long_max,by=5)[-1])
long_vals = data.frame(past = past_vals, future = past_vals, C = 1)

for(i in 1:length(C)){
  tmp = data.frame(past = past_vals, future = past_vals*C[i], C = C[i])
  long_vals = rbind(long_vals, tmp)
}
names(long_vals) = c('past_val', 'future_val', 'C')
####################


####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(10),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 1
long_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = long_vals)
names(long_out) = c('past_val', 'future_val', 'C', 'prob', 'sol', 'exact')
long_out$scenario = "Scenario 1"

####################

####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(0.1),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 2 and merge results
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = long_vals)
names(tmp) = c('past_val', 'future_val', 'C', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 2"
long_out <- rbind(long_out, tmp); rm(tmp)

####################


####################
## Define prior 3 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1, 100),
                            b =  c(10, 400),
                            weights = c(0.8, 0.2) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 3 and merge results
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = long_vals)
names(tmp) = c('past_val', 'future_val', 'C', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 3"
long_out <- rbind(long_out, tmp); rm(tmp)

####################

sum(abs(long_out$sol-long_out$exact)) ## total error between numerical and exact solution over all scenarios
max(abs(long_out$sol-long_out$exact)) ## max point of error
which(abs(long_out$sol-long_out$exact) == max(abs(long_out$sol-long_out$exact)))


###################
## Generate main ##
## results plots ##
###################

savePlots = FALSE
C = c(1, C) # add c=1 to vector of values for c to generate all limit lines
lim = 1-(1+long_out$C)^(-pars$a)
long_out$C = factor(long_out$C, levels = c(5,2,1,0.5,0.2))
y_upr = max( c(long_out$exact, long_out$sol, long_out$prob, 1 - (1/(C+1)^(pars$a))), na.rm=T)*1.01 
cond_mid = (long_out$past_val <= 100)

C_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = long_out[cond_mid,], aes(x = (past_val), y = .data[['sol']], color = C), linewidth = 1.25) +
  geom_point(data = long_out[cond_mid,], aes(x = past_val, y = .data[['sol']], color = C), size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0, 100, 25 )) +
  ylim(0, y_upr) +
  facet_wrap(.~scenario) + 
  theme_classic() +
  xlab('Number of past spillovers (N)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1) +
  guides(color = guide_legend(title="c")) +
  theme(# panel.grid = element_blank(),
    axis.line = element_line(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)


C_curves
if(savePlots){
  ggsave(filename = paste("./figures/C_curves_std_noSim.png"), plot = C_curves, width = 10.50, height = 5.00, units = "in")
}
###################


####################
## Calculate Rate ##
## of Convergence ##
##    Metrics     ##
####################

long_out$ROC = abs(long_out$exact-lim) / (lim)
long_out$ROC2 = NA
for(idx in which(long_out$past_val!=0)){
  long_out$ROC2[idx] = abs(long_out$exact[idx] - lim[idx]) / abs(long_out$exact[idx-1] - lim[idx])
}

conv_cond = (long_out$past_val > 100) 
cond_short = (long_out$past_val <= 20) 

L_ROC1 <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = long_out[conv_cond,], aes(x = past_val, y = .data[['ROC']], color = C), linewidth = 1.25) +
  geom_point(data = long_out[conv_cond,], aes(x = past_val, y = .data[['ROC']], color = C), size = 2) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  facet_wrap(~scenario) + 
  xlab('Number of past spillovers (N)') +
  ylab( 'Std_Dist' ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


L_ROC2 <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = long_out[conv_cond,], aes(x = past_val, y = .data[['ROC2']], color = C), linewidth = 1.25) +
  geom_point(data = long_out[conv_cond,], aes(x = past_val, y = .data[['ROC2']], color = C), size = 2) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  facet_wrap(~scenario) + 
  xlab('Number of past spillovers (N)') +
  ylab( 'ROC' ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1) +
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


L_ROC1
L_ROC2

if(savePlots){
  ggsave(filename = paste("./figures/C_ROC1.png"), plot = L_ROC1, width = 10.50, height = 5.50, units = "in")
  ggsave(filename = paste("./figures/C_ROC2.png"), plot = L_ROC2, width = 10.50, height = 5.50, units = "in")
}

####################


################################################################################
################################################################################
################################################################################
################################################################################
####    Runs and plots for different past-future spillover relationships    ####
################################################################################
################################################################################
################################################################################
################################################################################

###################
##   visualize   ##
## relationships ##
###################

## past-future spillover relationship parameters
long_max = 1000 # maximum value for 1:1 plots
long_step = 10 # step size for long C:1 line values
past_vals = c(0:10, seq(12,98,by=2), seq(100,long_max,by=long_step))
log_val = 1.1
pow_val = 2

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 0, N = 1, lambda = 1, k = 1, theta = 1,
                  future_model = 0, M = 1, lambda_f = 1, k_f = 1, theta_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

## data frame of input values (ex. N, M) uder different relationships
rel_vals = data.frame(past = past_vals, future = past_vals, fun = "linear")
rel_vals = rbind(rel_vals, data.frame(past = past_vals, future = past_vals^(1/pow_val), fun = "pow-") )
rel_vals = rbind(rel_vals, data.frame(past = past_vals, future = log(past_vals, base = log_val), fun = "log") )
rel_vals = rbind(rel_vals, data.frame(past = past_vals, future = log_val^past_vals, fun = "exp") )
rel_vals = rbind(rel_vals, data.frame(past = past_vals, future = past_vals^(pow_val), fun = "pow+") )
rel_vals$fun <- factor(rel_vals$fun, levels = c('exp', 'pow+', 'linear', 'pow-', 'log'))
rel_vals$future[rel_vals$future<0] = 0 # remove negative spillover values
if(pars$past_model == 0){ rel_vals$future <- round(rel_vals$future) } # round non-integer values for N-M model
names(rel_vals) = c('past_val', 'future_val', 'fun') # add column names

## visualization of past-future spillover relationships
rel_long <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = rel_vals, aes(x = past_val, y = future_val, color = fun), linewidth = 1.5) +
  theme_classic() +
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size = 24))

rel_short <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = rel_vals[rel_vals$past_val<=10,], aes(x = past_val, y = future_val, color = fun), linewidth = 1) +
  scale_x_continuous(breaks = c(0,5,10)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_text(size = 18))

rel_inset <- rel_long +
  inset_element(rel_short, left = 0.05, bottom = 0.5, right = .75, top = 1)
rel_inset

###################

####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(10),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 1
long_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = rel_vals)
names(long_out) = c('past_val', 'future_val', 'fun', 'prob', 'sol', 'exact')
long_out$scenario = "Scenario 1"

#####################


####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(0.1),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 2
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = rel_vals)
names(tmp) = c('past_val', 'future_val', 'fun', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 2"
long_out = rbind(long_out, tmp)

#####################


####################
## Define prior 3 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1, 100),
                            b =  c(10, 400),
                            weights = c(0.8, 0.2) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 3
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = rel_vals)
names(tmp) = c('past_val', 'future_val', 'fun', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 3"
long_out = rbind(long_out, tmp)
#####################


###################
## Generate main ##
## results plots ##
###################

long_outMAX = long_out
long_out = long_out[long_out$past_val<=100,]

y_uprMAX = max( c(long_outMAX$exact, long_outMAX$sol, long_outMAX$prob, 1 - (1/(C+1)^(pars$a))))*1.01 
y_upr = max( c(long_out$exact, long_out$sol, long_out$prob, 1 - (1/(C+1)^(pars$a))))*1.01 
y_upr_short = max( c(long_out$exact[long_out$past_val<=10], long_out$sol[long_out$past_val<=10], long_out$prob[long_out$past_val<=10], 1 - (1/(C+1)^(pars$a))))*1.01 

C = 1
savePlots = FALSE

## function legend labels
fun_label <- c(
  bquote(.(log_val) ^N),
  bquote(N^.(pow_val)),
  bquote(N),
  bquote( N^{"1/" ~ .(pow_val)} ),
  bquote( log[.(log_val)](N) )
)

L_curves <- ggplot() +
  geom_line(data = long_out, aes(x = past_val, y = .data[['sol']], color = fun), linewidth = 1.25) +
  scale_color_viridis(name = "f(N)", label = fun_label, discrete = T, option = 'H') +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = c(0,25,50,75,100)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4), limits = c(0,0.4)) +
  facet_wrap(.~scenario) + 
  xlab('Number of past spillovers (N)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)

if(pars$runSims){ L_curves <- L_curves + geom_point(data = long_out, aes(x = (past_val), y = .data[['prob']], color = fun),  size = 2.5)
}else{  L_curves <- L_curves + geom_point(data = long_out, aes(x = (past_val), y = .data[['sol']], color = fun), size = 2.5) }

L_curves


M_curves <- ggplot() +
  scale_color_viridis(name = "f(N)", label = fun_label, discrete = T, option = 'H') +
  geom_line(data = long_outMAX, aes(x = past_val, y = .data[['exact']], color = fun), linewidth = 1.25) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = c(0,250,500,750,1000)) +
  ylim(0, y_uprMAX) +
  facet_wrap(.~scenario) + 
  xlab('Number of past spillovers (N)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)

if(pars$runSims){ M_curves <- M_curves + geom_point(data = long_outMAX, aes(x = past_val, y = .data[['prob']], color = fun),  size = 2.5)
}else{  M_curves <- M_curves + geom_point(data = long_outMAX, aes(x = past_val, y = .data[['exact']], color = fun), size = 2.5) }

M_curves



M_lim0 <- ggplot() +
  scale_color_viridis(name = "f(N)", label = fun_label, discrete = T, option = 'H') +
  geom_line(data = long_outMAX, aes(x = past_val, y = .data[['sol']], color = fun), linewidth = 1.25) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = c(0,250,500,750,1000)) +
  scale_y_continuous(breaks = c(0,0.025, 0.05, 0.075, 0.1), limits = c(0,0.1)) +
  facet_wrap(.~scenario) + 
  xlab('Number of past spillovers (N)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)

if(pars$runSims){ M_lim0 <- M_lim0 + geom_point(data = long_outMAX, aes(x = past_val, y = .data[['prob']], color = fun),  size = 2.5)
}else{  M_lim0 <- M_lim0 + geom_point(data = long_outMAX, aes(x = past_val, y = .data[['sol']], color = fun), size = 2.5) }

M_lim0

if(savePlots){
  ggsave(filename = paste("./figures/R_curves_std_noSim.png"), plot = L_curves, width = 10.50, height = 5.00, units = "in")
  ggsave(filename = paste("./figures/M_curves_std_noSim.png"), plot = M_curves, width = 10.50, height = 5.00, units = "in")
  ggsave(filename = paste("./figures/M_lower_std_noSim.png"), plot = M_lim0, width = 10.50, height = 5.00, units = "in")
}
###################


####################################################################################
####################################################################################
####################################################################################
####################################################################################


############################
##  Evaluate  effects of  ##
## increasing b parameter ##
############################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 0, N = 1, lambda = 1, k = 1, theta = 1,
                  future_model = 0, M = 1, lambda_f = 1, k_f = 1, theta_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

spill_max = 20 # max mean spillover value
spill_step = 1 # step size -- must be an integer for N-M models
long_max = 1000 # maximum value for C:1 plots
long_step = 5 # step size for long C:1 line values
C = 1 # past:future spillover ratio -- slope of C:1 line

# vals = expand.grid(seq(0, spill_max, by = spill_step), seq(0, spill_max, by = spill_step))
# names(vals) = c('past_val', 'future_val')
# vals_long = data.frame(seq(0, long_max, by = long_step), C*seq(0, long_max, by = long_step))
# names(vals_long) = c('past_val', 'future_val')
# vals = rbind(vals, vals_long); rm(vals_long)
## sum(duplicated(vals)) ## test count of duplicated rows
## vals[duplicated(vals),] ## test view duplicated rows
# vals <- vals[!duplicated(vals),] ## removes rows that may be duplicated by merging vals_long
vec = c(seq(0,spill_max,by=spill_step), seq(spill_max, long_max, by = long_step)[-1])
vals = data.frame(vec, C*vec)
names(vals) = c('past_val', 'future_val')
####################

####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(pars$a),
                            b =  c(10),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 1
run_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
run_out$b = 10
names(run_out) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'b')
b_vec = c(1, 50, 100, 1000)
for(i in 1:length(b_vec)){
  mixtureDistn$b = b_vec[i]
  tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  tmp$b = b_vec[i]
  names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'b')
  run_out = rbind(run_out, tmp)
}
run_out$scenario = "Scenario 1"

####################

####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(pars$a),
                            b =  c(0.1),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 2 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
tmp2$b = 0.1
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'b')
b_vec = c(0.05, 0.01, 0.001)
for(i in 1:length(b_vec)){
  mixtureDistn$b = b_vec[i]
  tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  tmp$b = b_vec[i]
  names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'b')
  tmp2 = rbind(tmp2, tmp)
}
tmp2$scenario = "Scenario 2"
run_out = rbind(run_out, tmp2)

####################

####################
## Define prior 3 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(pars$a, 100),
                            b =  c(10, 400),
                            weights = c(0.8, 0.2) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 3 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
tmp2$b = 10
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'b')
b_vec = c(1, 50, 100, 1000)
for(i in 1:length(b_vec)){
  mixtureDistn$b = c(b_vec[i], 400)
  tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  tmp$b = b_vec[i]
  names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'b')
  tmp2 = rbind(tmp2, tmp)
}
tmp2$scenario = "Scenario 3"
run_out = rbind(run_out, tmp2)

####################

sum(abs(run_out$sol-run_out$exact)) ## total error between numerical and exact solution over all scenarios
max(abs(run_out$sol-run_out$exact)) ## max point of error
which(abs(run_out$sol-run_out$exact) == max(abs(run_out$sol-run_out$exact)))

###################
## Generate main ##
## results plots ##
###################
brk = seq(from = 0, to = max(max(run_out$exact), max(run_out$prob), max(run_out$sol)), by = 0.05)[-1] # contour value breakpoints
rng = 0 # range around the C:1 line as plot points
lim = 1-(C+1)^(-pars$a)
y_upr = max( c(run_out$prob, run_out$sol, run_out$prob, 1 - (1/(C+1)^(pars$a))))*1.01 
cond = (C*run_out$future_val <= (run_out$past_val + rng)) & (C*run_out$future_val >= (run_out$past_val - rng)) ## range of values around C:1 line
cond_mid = (run_out$past_val <= 100)
run_out$b = factor(run_out$b, levels = c(1000, 100, 50, 10, 1, 0.1, 0.05, 0.01, 0.001))
fname_tag = "Std"

savePlots = TRUE ## option to save plots to ./plots in directory. Note that this will overwrite existing plots from other runs

## generate 1:1 plots 

if(pars$runSims){ fname_tagL = paste0(fname_tag, '_Sims') }else{ fname_tagL = paste0(fname_tag, '_noSims') } ## change filename tag for simulation values

L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_mid,], aes(x = past_val, y = .data[['sol']], color = b), linewidth = 1.25) +
  geom_point(data = run_out[cond_mid,], aes(x = past_val, y = .data[['sol']], color = b), size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,100,by=25)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['sol']], color = b), linewidth = 1.25) +
  geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['sol']], color = b), size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


L_curves
curves

if(savePlots){
  ggsave(filename = paste("./figures/L_curvesb_", fname_tagL, ".png"), plot = L_curves, width = 10.50, height = 5.5, units = "in")
  ggsave(filename = paste("./figures/curvesb_", fname_tagL, ".png"), plot = curves, width = 10.50, height = 5.50, units = "in")
}
###################


####################
## Calculate Rate ##
## of Convergence ##
##    Metrics     ##
####################

run_out$ROC = abs(run_out$exact-lim) / (lim)
run_out$ROC2 = NA
for(idx in which(run_out$past_val!=0)){
  run_out$ROC2[idx] = abs(run_out$exact[idx] - lim) / abs(run_out$exact[idx-1] - lim)
}

conv_cond = (run_out$past_val > 100)

L_ROC1 <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[conv_cond,], aes(x = past_val, y = .data[['ROC']], color = b), linewidth = 1.25) +
  geom_point(data = run_out[conv_cond,], aes(x = past_val, y = .data[['ROC']], color = b), size = 2) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  facet_wrap(~scenario) + 
  xlab('Number of spillovers (N=M)') +
  ylab( 'Std_Dist' ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


L_ROC2 <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[conv_cond,], aes(x = past_val, y = .data[['ROC2']], color = b), linewidth = 1.25) +
  geom_point(data = run_out[conv_cond,], aes(x = past_val, y = .data[['ROC2']], color = b), size = 2) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  facet_wrap(~scenario) + 
  xlab('Number of spillovers (N=M)') +
  ylab( 'ROC' ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


L_ROC1
L_ROC2

if(savePlots){
  ggsave(filename = paste("./figures/L_bROC1_", fname_tagL, ".png"), plot = L_ROC1, width = 10.5, height = 5.5, units = "in")
  ggsave(filename = paste("./figures/L_bROC2", fname_tagL, ".png"), plot = L_ROC2, width = 10.50, height = 5.50, units = "in")
}


####################



##############################
##   Evaluate  effects of   ##
## increasing "a" parameter ##
##############################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 0, N = 1, lambda = 1, k = 1, theta = 1,
                  future_model = 0, M = 1, lambda_f = 1, k_f = 1, theta_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

spill_max = 20 # max mean spillover value
spill_step = 1 # step size -- must be an integer for N-M models
long_max = 1000 # maximum value for C:1 plots
long_step = 5 # step size for long C:1 line values
C = 1 # past:future spillover ratio -- slope of C:1 line

# vals = expand.grid(seq(0, spill_max, by = spill_step), seq(0, spill_max, by = spill_step))
# names(vals) = c('past_val', 'future_val')
# vals_long = data.frame(seq(0, long_max, by = long_step), C*seq(0, long_max, by = long_step))
# names(vals_long) = c('past_val', 'future_val')
# vals = rbind(vals, vals_long); rm(vals_long)
## sum(duplicated(vals)) ## test count of duplicated rows
## vals[duplicated(vals),] ## test view duplicated rows
# vals <- vals[!duplicated(vals),] ## removes rows that may be duplicated by merging vals_long
vec = c(seq(0,spill_max,by=spill_step), seq(spill_max, long_max, by = long_step)[-1])
vals = data.frame(vec, C*vec)
names(vals) = c('past_val', 'future_val')
####################

####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(10),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 1
run_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
run_out$a = 0.1
names(run_out) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'a')
a_vec = c(1, 0.5, 0.05, 0.01, 0.001)
for(i in 1:length(a_vec)){
  mixtureDistn$a = a_vec[i]
  tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  tmp$a = a_vec[i]
  names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'a')
  run_out = rbind(run_out, tmp)
}
run_out$scenario = "Scenario 1"

####################

####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(0.1),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 2 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
tmp2$a = 0.1
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'a')
for(i in 1:length(a_vec)){
  mixtureDistn$a = a_vec[i]
  tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  tmp$a = a_vec[i]
  names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'a')
  tmp2 = rbind(tmp2, tmp)
}
tmp2$scenario = "Scenario 2"
run_out = rbind(run_out, tmp2)

####################


####################
## Define prior 3 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1, 100),
                            b =  c(10, 400),
                            weights = c(0.8, 0.2) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 3 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
tmp2$a = 0.1
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'a')
for(i in 1:length(a_vec)){
  mixtureDistn$a = c(a_vec[i], 100)
  tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  tmp$a = a_vec[i]
  names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact', 'a')
  tmp2 = rbind(tmp2, tmp)
}
tmp2$scenario = "Scenario 3"
run_out = rbind(run_out, tmp2)

####################

sum(abs(run_out$sol-run_out$exact)) ## total error between numerical and exact solution over all scenarios
max(abs(run_out$sol-run_out$exact)) ## max point of error
which(abs(run_out$sol-run_out$exact) == max(abs(run_out$sol-run_out$exact)))


###################
## Generate main ##
## results plots ##
###################

brk = seq(from = 0, to = max(max(run_out$exact), max(run_out$prob), max(run_out$sol)), by = 0.05)[-1] # contour value breakpoints
rng = 0 # range around the C:1 line as plot points
lim = 1-(C+1)^(-run_out$a) ## vector of limit values used in convervence calculations
cond = (C*run_out$future_val <= (run_out$past_val + rng)) & (C*run_out$future_val >= (run_out$past_val - rng)) ## range of values around C:1 line
cond_short = (run_out$past_val <= spill_max)  ## C:1 line values from heatmap
cond_mid = (run_out$past_val <= 100)
y_upr = max( c(run_out$prob, run_out$sol, run_out$prob, 1 - (1/(C+1)^(pars$a))))*1.01 
run_out$a = factor(run_out$a, levels = c(1, 0.5, 0.1, 0.05, 0.01, 0.001))
a_vec = c(0.1, a_vec)

savePlots = TRUE ## option to save plots to ./plots in directory. Note that this will overwrite existing plots from other runs

## generate 1:1 plots 
fname_tag="std"
if(pars$runSims){ fname_tagL = paste0(fname_tag, '_Sims') }else{ fname_tagL = paste0(fname_tag, '_noSims') } ## change filename tag for simulation values

L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_mid,], aes(x = past_val, y = .data[['exact']], color = a), linewidth = 1.25) +
  geom_point(data = run_out[cond_mid,], aes(x = past_val, y = .data[['exact']], color = a), size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,100,by=25)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']], color = a), linewidth = 1.25) +
  geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']], color = a), size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


L_curves
curves

if(savePlots){
  ggsave(filename = paste("./figures/L_curvesa_", fname_tagL, ".png"), plot = L_curves, width = 10.50, height = 5.5, units = "in")
  ggsave(filename = paste("./figures/curvesa_", fname_tagL, ".png"), plot = curves, width = 10.50, height = 5.50, units = "in")
}

###################

####################
## Calculate Rate ##
## of Convergence ##
##    Metrics     ##
####################

run_out$ROC = abs(run_out$exact-lim) / (lim)
run_out$ROC2 = NA
for(idx in which(run_out$past_val!=0)){
  run_out$ROC2[idx] = abs(run_out$exact[idx] - lim[idx]) / abs(run_out$exact[idx-1] - lim[idx])
}

## test for logrithmic convergence
run_out$log_conv = NA
for(idx in which(run_out$past_val > 1)){
  run_out$log_conv[idx] = abs(run_out$exact[idx] - run_out$exact[idx-1] ) / abs(run_out$exact[idx-1] - lim[idx-2])
}

conv_cond = (run_out$past_val > 100)

L_ROC1 <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[conv_cond,], aes(x = past_val, y = .data[['ROC']], color = a), linewidth = 1.25) +
  geom_point(data = run_out[conv_cond,], aes(x = past_val, y = .data[['ROC']], color = a), size = 2) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  facet_wrap(~scenario) + 
  xlab('Number of spillovers (N=M)') +
  ylab( 'Std_Dist' ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


L_ROC2 <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[conv_cond,], aes(x = past_val, y = .data[['ROC2']], color = a), linewidth = 1.25) +
  geom_point(data = run_out[conv_cond,], aes(x = past_val, y = .data[['ROC2']], color = a), size = 2) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  facet_wrap(~scenario) + 
  xlab('Number of spillovers (N=M)') +
  ylab( 'ROC' ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1) +
  theme(
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)

L_ROC1
L_ROC2

if(savePlots){
  ggsave(filename = paste("./figures/L_aROC1_", fname_tagL, ".png"), plot = L_ROC1, width = 10.50, height = 5.5, units = "in")
  ggsave(filename = paste("./figures/L_aROC2_", fname_tagL, ".png"), plot = L_ROC2, width = 10.50, height = 5.50, units = "in")
}
####################

################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################



####################################################################################
##  Analyzing  different  ##########################################################
## beta mixture behaviors ##########################################################
####################################################################################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 0, N = 1, lambda = 1, k = 1, theta = 1,
                  future_model = 0, M = 1, lambda_f = 1, k_f = 1, theta_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

spill_max = 20 # max mean spillover value
spill_step = 1 # step size -- must be an integer for N-M models
long_max = 100000 # maximum value for C:1 plots
long_step = 100 # step size for long C:1 line values
C = 1 # past:future spillover ratio -- slope of C:1 line

vec = c(seq(0,spill_max,by=spill_step), seq(spill_max, long_max, by = long_step)[-1])
vals = data.frame(vec, C*vec)
names(vals) = c('past_val', 'future_val')

####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(10),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 1
run_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(run_out) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
run_out$scenario = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
scenario1 = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )

####################

####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.01),
                            b =  c(10),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 2 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp2$scenario = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
scenario2 = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
run_out = rbind(run_out, tmp2)

####################


####################
## Define prior 3 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1, 0.01),
                            b =  c(10, 10),
                            weights = c(0.1, 0.9) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}

## Run model with prior 3 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp2$scenario = 'Mixture'
scenario3 = 'Mixture'
run_out = rbind(run_out, tmp2)
####################


###################
## Generate main ##
## results plots ##
###################
brk = seq(from = 0, to = max(max(run_out$exact), max(run_out$prob), max(run_out$sol)), by = 0.05)[-1] # contour value breakpoints
rng = 0 # range around the C:1 line as plot points
cond = (C*run_out$future_val <= (run_out$past_val + rng)) & (C*run_out$future_val >= (run_out$past_val - rng)) ## range of values around C:1 line
cond_short = (run_out$past_val <= spill_max)  ## C:1 line values from heatmap
a_vec = c(0.1, 0.01)

savePlots = FALSE ## option to save plots to ./plots in directory. Note that this will overwrite existing plots from other runs

## generate 1:1 plots 

if(pars$runSims){ fname_tagL = paste0(fname_tag, '_Sims') }else{ fname_tagL = paste0(fname_tag, '_noSims') } ## change filename tag for simulation values

{
  y_upr = max( c(run_out$exact, 1 - (1/(C+1)^(pars$a))))*1.01 
  
  L_curve1 <- ggplot() +
    geom_line(data = run_out[run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,spill_max, y_upr/5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +
    xlab("Spillover Rate") +
    ylab( expression( paste(P(H[F]>0)) ) ) +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28),
      # axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      aspect.ratio = 1.0)
  
  if(pars$runSims){ L_curve1 <- L_curve1 + geom_point(data = run_out[run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5)
  }else{  L_curve1 <- L_curve1 + geom_point(data = run_out[run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) }
  
  
  
  curve1 <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, spill_max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    xlab("Spillover Rate") +
    ylab( expression( paste(P(H[F]>0)) ) ) +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28),
      # axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      aspect.ratio = 1.0)
  
  if(pars$runSims){ curve1 <- curve1 + geom_point(data = run_out[cond_short & run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5)
  }else{  curve1 <- curve1 + geom_point(data = run_out[cond_short & run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) }
  
  
  
  L_curve2 <- ggplot() +
    geom_line(data = run_out[run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,spill_max, y_upr/5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    xlab("Inherent rate of spillover (N = M)") +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      aspect.ratio = 1.0)
  
  if(pars$runSims){ L_curve2 <- L_curve2 + geom_point(data = run_out[run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5)
  }else{  L_curve2 <- L_curve2 + geom_point(data = run_out[run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) }
  
  
  curve2 <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, spill_max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    xlab("Inherent rate of spillover (N = M)") +
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      aspect.ratio = 1.0)
  
  if(pars$runSims){ curve2 <- curve2 + geom_point(data = run_out[cond_short & run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5)
  }else{  curve2 <- curve2 + geom_point(data = run_out[cond_short & run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) }
  
  L_curve3 <- ggplot() +
    geom_line(data = run_out[run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,spill_max, y_upr/5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      aspect.ratio = 1.0)
  
  if(pars$runSims){ L_curve3 <- L_curve3 + geom_point(data = run_out[run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5)
  }else{  L_curve3 <- L_curve3 + geom_point(data = run_out[run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) }
  
  
  curve3 <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, spill_max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    theme_classic() +
    theme(
      plot.title = element_blank(),
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 28),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 28, color = "#FFFFFF"),
      axis.ticks.y = element_blank(),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      aspect.ratio = 1.0)
  
  if(pars$runSims){ curve3 <- curve3 + geom_point(data = run_out[cond_short & run_out$scenario == 'Mixture'], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5)
  }else{  curve3 <- curve3 + geom_point(data = run_out[cond_short & run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) }
  
  ## Merge Scenarios 1 and 3 plots
  
  L_curve_merge <- ggplot() +
    geom_line(data = run_out[run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_line(data = run_out[run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_line(data = run_out[run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,spill_max, y_upr/5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    theme_classic() +
    theme(# panel.grid = element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.line = element_line(),
      strip.text = element_text(size = 20),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      # axis.ticks.y = element_blank(),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      aspect.ratio = 1.0)
  
  if(pars$runSims){ curve_merge <- curve_merge + 
    geom_point(data = run_out[run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5) +
    geom_point(data = run_out[run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5) +
    geom_point(data = run_out[run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5)
  }else{  curve_merge <- curve_merge + 
    geom_point(data = run_out[run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) +
    geom_point(data = run_out[run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) +
    geom_point(data = run_out[run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) }
  
  
  curve_merge <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_line(data = run_out[cond_short & run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_line(data = run_out[cond_short & run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, spill_max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    theme_classic() +
    theme(# panel.grid = element_blank(),
      panel.spacing = unit(1, "lines"),
      axis.line = element_line(),
      strip.text = element_text(size = 20),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      # axis.ticks.y = element_blank(),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      aspect.ratio = 1.0)
  
  if(pars$runSims){ curve_merge <- curve_merge + 
    geom_point(data = run_out[cond_short & run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5) +
    geom_point(data = run_out[cond_short & run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5) +
    geom_point(data = run_out[cond_short & run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['prob']]), color = '#0C2340', size = 2.5)
  }else{  curve_merge <- curve_merge + 
    geom_point(data = run_out[cond_short & run_out$scenario == scenario1,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) +
    geom_point(data = run_out[cond_short & run_out$scenario == scenario2,], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) +
    geom_point(data = run_out[cond_short & run_out$scenario == 'Mixture',], aes(x = (past_val), y = .data[['exact']]), color = '#0C2340', size = 2.5) }
}

L_curve1
curve1

L_curve2
curve2

L_curve3
curve3

L_curve_merge
curve_merge

L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out, aes(x = past_val, y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max, y_upr/5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Inherent rate of spillover') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)

if(pars$runSims){ L_curves <- L_curves + geom_point(data = run_out, aes(x = past_val, y = .data[['prob']]), color = '#0C2340', size = 2.5)
}else{  L_curves <- L_curves + geom_point(data = run_out, aes(x = past_val, y = .data[['exact']]), color = '#0C2340', size = 2.5) }



curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color = '#0C2340', linewidth = 1.25) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(a_vec)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max, y_upr/5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab('Inherent rate of spillover') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    aspect.ratio = 1.0)

if(pars$runSims){ curves <- curves + geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['prob']]), color = '#0C2340', size = 2.5)
}else{  curves <- curves + geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color = '#0C2340', size = 2.5) }

L_curves
curves

if(savePlots){
  ggsave(filename = paste("./figures/L_curve1a_", fname_tagL, ".png"), plot = L_curve1, width = 3.85, height = 3.85, units = "in")
  ggsave(filename = paste("./figures/curve1a_", fname_tagL, ".png"), plot = curve1, width = 5.85, height = 5.85, units = "in")
  
  ggsave(filename = paste("./figures/L_curve2a_", fname_tagL, ".png"), plot = L_curve2, width = 3.85, height = 3.85, units = "in")
  ggsave(filename = paste("./figures/curve2a_", fname_tagL, ".png"), plot = curve2, width = 5.85, height = 5.85, units = "in")
  
  ggsave(filename = paste("./figures/L_curve3a_", fname_tagL, ".png"), plot = L_curve3, width = 3.85, height = 3.85, units = "in")
  ggsave(filename = paste("./figures/curve3a_", fname_tagL, ".png"), plot = curve3, width = 5.85, height = 5.85, units = "in")
  
  ggsave(filename = paste("./figures/L_curvesa_", fname_tagL, ".png"), plot = L_curves, width = 5.85, height = 5.85, units = "in")
  ggsave(filename = paste("./figures/curvesa_", fname_tagL, ".png"), plot = curves, width = 10.00, height = 5.50, units = "in")
}

###################


#########################
## Prior Mean-Variance ##
##       Analyses      ##
#########################


####################
## Specify  Model ##
##   Parameters   ##
####################
## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 0, N = 1, lambda = 1, k = 1, theta = 1,
                  future_model = 0, M = 1, lambda_f = 1, k_f = 1, theta_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

spill_max = 20 # max mean spillover value
spill_step = 1 # step size -- must be an integer for N-M models
long_max = 100 # maximum value for C:1 plots
long_step = 2 # step size for long C:1 line values
C = 1 # past:future spillover ratio -- slope of C:1 line

vec = c(seq(0,spill_max,by=spill_step), seq(spill_max, long_max, by = long_step)[-1])
vals = data.frame(vec, C*vec)
names(vals) = c('past_val', 'future_val')
####################

phi_bar = 1/100
a_vec = seq(0,1,by =0.0001)[-1] # exclude a=0
a_vec = a_vec[-length(a_vec)] # exclude a=1
b_vec = ((1/phi_bar) - 1)*a_vec # all corresponding values of b that give a mean of phi_bar
all( round(a_vec/(a_vec+b_vec), digits = 15) == 1/100) ## ensure all means are "numerically" equal 

prior_vals = data.frame(a=a_vec, b=b_vec)
prior_vals$mean = (prior_vals$a) / (prior_vals$a + prior_vals$b)
prior_vals$var = (prior_vals$a*prior_vals$b) / ( (prior_vals$a + prior_vals$b)^2 * (prior_vals$a + prior_vals$b + 1) )
prior_vals$CV = sqrt(prior_vals$var) / prior_vals$mean

if( which(prior_vals$b == min(prior_vals$b)) == which(prior_vals$a == min(prior_vals$a)) ){
  minTrue = prior_vals[which(prior_vals$b == min(prior_vals$b)),]
  prior_vals$mean[which(prior_vals$b == min(prior_vals$b))]
  prior_vals$var[which(prior_vals$b == min(prior_vals$b))]
}

prior_vals = prior_vals[!(b_vec<1),] ## remove b<1
apply(prior_vals, MARGIN = 2, FUN = min)
apply(prior_vals, MARGIN = 2, FUN = max)

## bool test for matching extreme values of a and b
if( which(prior_vals$b == min(prior_vals$b)) == which(prior_vals$a == min(prior_vals$a)) ){
  minVals = prior_vals[which(prior_vals$b == min(prior_vals$b)),]
  prior_vals$mean[which(prior_vals$b == min(prior_vals$b))]
  prior_vals$var[which(prior_vals$b == min(prior_vals$b))]
}
if( which(prior_vals$b == max(prior_vals$b)) == which(prior_vals$a == max(prior_vals$a)) ){
  maxVals = prior_vals[which(prior_vals$b == max(prior_vals$b)),]
  prior_vals$mean[which(prior_vals$b == max(prior_vals$b))]
  prior_vals$var[which(prior_vals$b == max(prior_vals$b))]
}


####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################
mixtureDistn <- data.frame( a = minVals['a'],
                            b =  minVals['b'],
                            weights = c(1) )
phi_vec = seq(from = 0, to = 1, by = 0.0001)
prior = data.frame(phi = phi_vec, density = betaMix(x = phi_vec, weights = mixtureDistn$weights, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b )$density, 
                   prior = paste0("Scenario 1"))
prior = data.frame(phi = phi_vec, density = dbeta(x = phi_vec, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b ), 
                   prior = paste0("Scenario 1"))

## Run model with prior 1
run_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(run_out) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
run_out$prior = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
run_out$scenario = 'Scenario 1'
run_out$vals = paste0('a=', as.numeric(minVals['a']),'; b=', as.numeric(minVals['b']))
run_out$lim = as.numeric(1 - (C+1)^(-minVals['a']))
scenario1 = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )

####################

 
####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################
mixtureDistn <- data.frame( a = maxVals['a'],
                            b =  maxVals['b'],
                            weights = c(1) )
tmp = data.frame(phi = phi_vec, density = betaMix(x = phi_vec, weights = mixtureDistn$weights, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b )$density, 
                 prior = paste0("Scenario 2"))
tmp = data.frame(phi = phi_vec, density = dbeta(x = phi_vec, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b ), 
                 prior = paste0("Scenario 2"))
prior = rbind(tmp, prior)

## Run model with prior 2 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp2$prior = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
tmp2$scenario = 'Scenario 2'
tmp2$vals = paste0('a=', as.numeric(maxVals['a']),'; b=', as.numeric(maxVals['b']))
tmp2$lim = as.numeric(1 - (C+1)^(-maxVals['a']))
scenario2 = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
run_out = rbind(run_out, tmp2)

####################

####################
## Define prior 3 ##
## parameters and ##
## run the model  ##
####################
a_mid = 0.3
mixtureDistn <- data.frame( a = prior_vals$a[prior_vals$a == a_mid],
                            b =  prior_vals$b[prior_vals$a == a_mid],
                            weights = c(1) )
tmp = data.frame(phi = phi_vec, density = betaMix(x = phi_vec, weights = mixtureDistn$weights, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b )$density, 
                 prior = paste0("Scenario 3"))
tmp = data.frame(phi = phi_vec, density = dbeta(x = phi_vec, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b ), 
                 prior = paste0("Scenario 3"))
prior = rbind(tmp, prior)

## Run model with prior 3 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp2$prior = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
tmp2$scenario = 'Scenario 3'
tmp2$vals = paste0('a=', as.numeric(prior_vals$a[prior_vals$a == a_mid]),'; b=', as.numeric(prior_vals$b[prior_vals$a == a_mid]))
tmp2$lim = 1 - (C+1)^(-a_mid)
scenario3 = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
run_out = rbind(run_out, tmp2)

####################

prior_plots <- ggplot() +
  geom_line( data = prior, aes(x = phi, y = density), color = viridis(n=6)[3], linewidth = 1.3, alpha = 1 ) +
  geom_vline(xintercept = phi_bar, color = 'red', linewidth = 1.2) +
  labs( x = expression(phi), y = 'Density' ) + 
  facet_wrap( ~ prior) +
  scale_x_continuous(breaks=seq(0,0.1,by=0.05), limits = c(0,0.1)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme_classic() +
  theme(panel.spacing = unit(2, "lines"),
        plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 10),
        aspect.ratio = 1.0)
prior_plots

a_vals = data.frame(a = c(minVals$a, maxVals$a, a_mid), scenario = c("Scenario 1", "Scenario 2", "Scenario 3"))
cond_short = (run_out$past_val <= spill_max)  ## C:1 line values from heatmap
y_upr = max( c(run_out$prob, run_out$sol, run_out$prob, c(1 - (C+1)^(-a_vals$a)) ) )*1.01 


L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out, aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out, aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out, aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,100,by=25)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~vals) + 
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~vals) + 
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)

y_upr1 = max( c(run_out$prob[run_out$scenario == "Scenario 1"], run_out$sol[run_out$scenario == "Scenario 1"], run_out$prob[run_out$scenario == "Scenario 1"], as.numeric(1 - (C+1)^(-minVals['a'])) ) )*1.01 

curve1 <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short & run_out$scenario == "Scenario 1",], aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out[cond_short & run_out$scenario == "Scenario 1",], aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out[cond_short & run_out$scenario == "Scenario 1",], aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(limits = c(0, y_upr1)) +     
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    axis.line = element_line(),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 28),
    aspect.ratio = 1.0)

L_curves
curves
curve1

if(savePlots){
  ggsave(filename = paste("./figures/curve1_abvals.png"), plot = curve1, width = 10.5, height = 5.85, units = "in")
  ggsave(filename = paste("./figures/L_curves_abvals.png"), plot = L_curves, width = 10.5, height = 5.85, units = "in")
}

#########################

########################################################################
########################################################################
########################################################################
########################################################################

####################
## Define mixture ##
##  prior 1 and   ##
## run the model  ##
####################
mixtureDistn <- data.frame( a = c(as.numeric(minVals['a']), 100),
                            b =  c(as.numeric(minVals['b']), 400),
                            weights = c(0.8, 0.2) )
phi_vec = seq(from = 0, to = 1, by = 0.0001)
prior = data.frame(phi = phi_vec, density = betaMix(x = phi_vec, weights = mixtureDistn$weights, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b )$density, 
                   prior = paste0("Scenario 1"))

## Run model with prior 1
run_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(run_out) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
run_out$scenario = 'Scenario 1'
run_out$vals = paste0('a=', as.numeric(minVals['a']),'; b=', as.numeric(minVals['b']))
run_out$lim = as.numeric(1 - (C+1)^(-minVals['a']))
scenario1 = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )

####################


####################
## Define mixture ##
##  prior 2 and   ##
## run the model  ##
####################
mixtureDistn <- data.frame( a = c(as.numeric(maxVals['a']), 100),
                            b =  c(as.numeric(maxVals['b']), 400),
                            weights = c(0.8, 0.2) )
tmp = data.frame(phi = phi_vec, density = betaMix(x = phi_vec, weights = mixtureDistn$weights, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b )$density, 
                 prior = paste0("Scenario 2"))
prior = rbind(tmp, prior)

## Run model with prior 2 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp2$scenario = 'Scenario 2'
tmp2$vals = paste0('a=', as.numeric(maxVals['a']),'; b=', as.numeric(maxVals['b']))
tmp2$lim = as.numeric(1 - (C+1)^(-maxVals['a']))
scenario2 = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
run_out = rbind(run_out, tmp2)

####################

####################
## Define mixture ##
##  prior 3 and   ##
## run the model  ##
####################
a_mid = 0.3
mixtureDistn <- data.frame( a = c(as.numeric(prior_vals$a[prior_vals$a == a_mid]), 100),
                            b = c(as.numeric(prior_vals$b[prior_vals$a == a_mid]), 400),
                            weights = c(0.8, 0.2) )
tmp = data.frame(phi = phi_vec, density = betaMix(x = phi_vec, weights = mixtureDistn$weights, shape1 = mixtureDistn$a, shape2 = mixtureDistn$b )$density, 
                 prior = paste0("Scenario 3"))
prior = rbind(tmp, prior)

## Run model with prior 3 and merge results
tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp2$scenario = 'Scenario 3'
tmp2$vals = paste0('a=', as.numeric(prior_vals$a[prior_vals$a == a_mid]),'; b=', as.numeric(prior_vals$b[prior_vals$a == a_mid]))
tmp2$lim = 1 - (C+1)^(-a_mid)
scenario3 = paste0('beta(', mixtureDistn$a, '; ', mixtureDistn$b, ')' )
run_out = rbind(run_out, tmp2)

####################


prior_plots <- ggplot() +
  geom_line( data = prior, aes(x = phi, y = density), color = viridis(n=6)[3], linewidth = 1.3, alpha = 1 ) +
  geom_vline(xintercept = phi_bar*.8 + (0.1*.2), color = 'red', linewidth = 1.2) +
  labs( x = expression(phi), y = 'Density' ) + 
  facet_wrap( ~ prior) +
  scale_x_continuous(breaks=seq(0,0.3,by=0.1), limits = c(0,0.3)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme_classic() +
  theme(panel.spacing = unit(2, "lines"),
        plot.title = element_text(size = 28),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 10),
        aspect.ratio = 1.0)
prior_plots

a_vals = data.frame(a = c(minVals$a, maxVals$a, a_mid), scenario = c("Scenario 1", "Scenario 2", "Scenario 3"))
cond_short = (run_out$past_val <= spill_max)  ## C:1 line values from heatmap
y_upr = max( c(run_out$prob, run_out$sol, run_out$prob, c(1 - (C+1)^(-a_vals$a)) ) )*1.01 


L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out, aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out, aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out, aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,100,by=25)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~vals) + 
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~vals) + 
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_text(size = 20),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    # axis.ticks.y = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)

y_upr1 = max( c(run_out$prob[run_out$scenario == "Scenario 1"], run_out$sol[run_out$scenario == "Scenario 1"], run_out$prob[run_out$scenario == "Scenario 1"], as.numeric(1 - (C+1)^(-minVals['a'])) ) )*1.01 

curve1 <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short & run_out$scenario == "Scenario 1",], aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out[cond_short & run_out$scenario == "Scenario 1",], aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out[cond_short & run_out$scenario == "Scenario 1",], aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(limits = c(0, y_upr1)) +     
  xlab('Number of spillovers (N=M)') +
  ylab( expression( paste(P(H[F]>0)) ) ) +
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1)+
  theme(# panel.grid = element_blank(),
    axis.line = element_line(),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 28),
    aspect.ratio = 1.0)

L_curves
curves
curve1

if(savePlots){
  ggsave(filename = paste("./figures/curve1_abvals_mixture.png"), plot = curve1, width = 10.5, height = 5.85, units = "in")
  ggsave(filename = paste("./figures/L_curves_abvals_mixture.png"), plot = L_curves, width = 10.5, height = 5.85, units = "in")
}





