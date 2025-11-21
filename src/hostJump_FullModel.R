#### UPDATED FOR POISSON TIME-EXPLICIT MODEL ####

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
library(colorspace)
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
  mixtureDistn$cum_weights <- cumsum(mixtureDistn$weights)
  
  ## ensure path matches the location in the cpp code
  write.table(x = mixtureDistn, sep = ',', row.names = F, col.names = F, file = "./src/simulation_files/betaMix_pars.csv")
  
  return(0)
  
}


## generate full output file for model runs over all 3 priors used in the main text
## includes functionality to plot different past-future spillover relationships
## pars should be a data frame containing all model parameters
## mixtureDistn should be a data frame containing the parameters for the prior distribution
## vals should be a data frame of values for the past and future parameters to iterate over (i.e., N, M, lambda)
fullModel_eval <- function(pars, mixtureDistn, vals){
  
  run_out <- vals
  run_out$prob = 0 # create simulation probability column
  run_out$sol = 0 # create analytic probability column
  run_out$exact = 0 # create exact probability column
  
  for(index in 1:dim(vals)[1]){
    
    ## set values being tested
    pars[names(vals)] <- vals[index, names(vals)]
    
    ## run simulation model
    if(pars$runSims == TRUE){
      out <- HJ_simulation(n_reps = pars$reps, parameters = pars, batch_name = 'test', output_mode = 0)
      run_out$prob[index] = out$total_prob[1] # store results -- the [1] prevents errors when long output is requested 
    }else{
      run_out$prob[index] = 0
    }
    
    ## evaluate model using numerical integration
    val <- model_solution(past = pars$past_model+1, future = pars$future_model+1, pars = pars, a = mixtureDistn$a, b = mixtureDistn$b, weights = mixtureDistn$weights)
    run_out$sol[index] = val
    
    ## exact analytic solution via beta binomial or confluent hypergeometric function
    ex_val <- model_solution_exact(past = pars$past_model, future = pars$future_model, pars = pars, a = mixtureDistn$a, b = mixtureDistn$b, weights = mixtureDistn$weights)
    run_out$exact[index] = ex_val
    
  }
  
  return(run_out)
  
}


## Numerical model solution (integration)
# UD_prior == data frame with columns "nodes" and "weights" to be used as the prior distribution
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
      lambda=as.numeric(par['t_p'])*as.numeric(par['lambda']); H=as.numeric(par['H_crit']); C = 1;
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
                      N=as.numeric(pars['N']); lambda_f=as.numeric(par['t_f'])*as.numeric(pars['lambda_f']); H=as.numeric(pars['H_crit']);
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
                      lambda=as.numeric(par['t_f'])*as.numeric(pars['lambda']); M=as.numeric(pars['M']); H=as.numeric(pars['H_crit']);
                      vals = ((1-UD_prior$nodes)^M) * UD_prior$weights * exp(-lambda*UD_prior$nodes)*(lambda^H / factorial(H))
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    {
                      lambda=as.numeric(par['t_f'])*as.numeric(pars['lambda']); lambda_f=as.numeric(pars['lambda_f']); H=as.numeric(pars['H_crit']);
                      vals = exp(-lambda_f*UD_prior$nodes) * UD_prior$weights * exp(-lambda*UD_prior$nodes)*(lambda^H / factorial(H))
                      value = ( 1-const*sum( vals ) )
                    },
                    #########################################
                    #########################################
                    {
                      lambda=as.numeric(par['t_f'])*as.numeric(pars['lambda']); k_f=as.numeric(pars['k_f']); theta_f=as.numeric(pars['theta_f']); H=as.numeric(pars['H_crit']);
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
      lambda=as.numeric(par['lambda'])*as.numeric(par['t_p']); H=as.numeric(par['H_crit']); C = 1;
      if(H>0){ 
        C = (lambda^H)/(factorial(H))
        H_norm = seq(0, H, by=1)[-(H+1)] # all values of N that must have been observed
        C1 = 1 - sum( ((lambda^H_norm)*exp(-lambda))/(factorial(H_norm)) ) # normalizing constant for poisson number of spillovers
        C = C / C1
      }
      return(as.numeric( C * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * (exp(-lambda*x)) )) 
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
                        lambda_f=as.numeric(par['lambda_f'])*as.numeric(par['t_f'])
                        return( (exp(-lambda_f*x)) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * ((1-x)^N) )
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
                        lambda=as.numeric(par['lambda'])*as.numeric(par['t_p']);
                        M=as.numeric(par['M']);
                        return( ((1-x)^M) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * (exp(-lambda*x)) )
                      }
                      value = (1-const*integrate(f=integrand, par=pars, a = a, b = b, weights = weights, lower=0, upper=1, abs.tol=0, stop.on.error=F)$value)
                    },
                    #########################################
                    #########################################
                    {
                      integrand <- function(x, par, a = a, b = b, weights = weights){
                        lambda=as.numeric(par['lambda'])*as.numeric(par['t_p']);
                        lambda_f=as.numeric(par['lambda_f'])*as.numeric(par['t_f']);
                        return( (exp(-lambda_f*x)) * betaMix(x = x, weights = weights, shape1 = a, shape2 = b)$density * (exp(-lambda*x)) )
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
  
  ## confluent hypergeometric solution
  if(pars$past_model == 1 && pars$future_model == 1){
    
    strvec = format( c( t=1, length(a), pars$lambda*pars$t_p, pars$lambda_f*pars$t_f, a, b, weights), digits = 5)
    
    setwd("~/Desktop/Repos/HostJump_Model/src") ## call has to be from location of .exe file or parameter read-in fails???
    
    ## Run the model
    ## The path to the bTB cpp binary file must be set correctly in the sys call below:
    nm = paste0("./poissonExact.exe")
    r <- system2( nm, args = strvec, stdout = TRUE)
    out <- read.table(text = r, header = TRUE, sep = ';', check.names = FALSE) %>% mutate_all(as.numeric)
    
    setwd("..")
    
    return( as.numeric(out) ) 
    
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
              parameters$lambda*parameters$t_p, # 2, constant spillover rate 
              0, # 3, gamma shape parameter -- removed
              0, # 4, gamma scale parameter -- removed
              
              parameters$future_model, # 5, switch indicator 0: fixed M; 1: poisson M; 2: gamma-poisson M;
              parameters$M, # 6, fixed future spillovers -- used if future_model == 0;
              parameters$lambda_f*parameters$t_f, # 7, constant future spillover rate 
              0, # 8, gamma future shape parameter -- removed
              0, # 9, gamma future scale parameter -- removed 
              
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

## frequency/probability remaining plots
lambda = seq(0,5, by = 0.1)
phi = seq(0,1, by = 0.01)
T_lim = c(0,5,20,100) 
dat2 = expand.grid(lambda, phi, T_lim)
names(dat2) = c('lambda', 'phi', 'T_lim')
dat2$prob = exp(-dat2$lambda*dat2$phi*dat2$T_lim)


# New facet label names for T
T.labs <- c("T[P]==0", 
            "T[P]==5", 
            "T[P]==20", 
            "T[P]==100")

names(T.labs) = c("0","5","20","100")

## fraction of pathogens remaining in zoonotic pool after T units of time
annotations_df <- data.frame(
  x = rep(c(4.3, 0.7, 4.3, 0.7),times=4),   # x positions of the annotations
  y = rep(c(0.8, 0.8, 0.1, 0.1),times=4),       # y positions of the annotations
  label = rep(c("A", "C", "B", "D"),times=4),    # Labels
  T_lim = rep(c("0","5","20","100"), each = 4),
  color = c(rep("black",times=4), 
          c(rep("white",times=3),'black'),
          #rep("white",times=4),
          c(rep("white",times=4)),
            c(rep("white",times=4)))  # Different colors for each facet
)

annotations_df$T_lim <- factor(annotations_df$T_lim, levels = c("0", "5", "20", "100"))
dat2$T_lim <- factor(dat2$T_lim, levels = c("0", "5", "20", "100"))

## conceptual figure 1E
p2 <- ggplot() +
  scale_fill_viridis(option = "B", discrete = F, direction = 1) +
  geom_tile(data = dat2, aes(x = lambda, y = phi, fill = prob)) +
  geom_contour(data = dat2[dat2$T_lim != "0", ], aes(x = lambda, y = phi, z = prob), color = "#FFFFFF", breaks = seq(0.2,0.8,0.2), alpha = 0.55) +
  # geom_vline(xintercept = 1, color = '#4490FEFF', linetype = 'dashed', linewidth = 0.75) +
  # geom_vline(xintercept = 4, color = '#EA4F0DFF', linetype = 'dashed', linewidth = 0.75) +
  geom_text(data = annotations_df, aes(x = x, y = y, label = label, color = color), size = 5, fontface = "bold") +
  scale_color_identity() +  # Use colors directly from the `color` column
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels=c('0','0.25','0.50','0.75','1')) +
  # labs(x = expression( "Spillover rate "*(lambda)), 
  #      y = expression(phi), 
  #      fill = expression("P("*H[P]*"=0)")) + 
  labs(x = expression( "spillover rate "*(lambda)), 
       y = expression("host jump prob.\n per spillover"*(phi)), 
       fill = expression("fraction of remaining\nzoonotic pathogens")) + 
  facet_wrap( ~ T_lim, nrow = 1, labeller = as_labeller(T.labs, default = label_parsed)) +
  coord_cartesian(expand = FALSE) +
  theme_bw() + 
  theme(
    legend.position = "bottom",
    axis.title.y = element_text( vjust = 0.5, size = 16),
    axis.title.x = element_text(size = 14),
    legend.key.width = unit(1.55, "cm"),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(vjust = 0.8, size = 12),
    legend.title.align = 0.1,
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.spacing = unit(1, "lines"),
    aspect.ratio = 1.0)
p2

ggsave(filename = paste("./figures/conFig_E_dist.png"), plot = p2, width = 10.00, height = 4, units = "in")


############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################


####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## past/future _model: 0 = count model, 1 = poisson model
## prior_type = 2 uses the mixture model which effectively also captures the single beta prior
LAMBDA = 1
pars = data.frame(past_model = 1, N = 1, lambda = LAMBDA, t_p = 1,
                  future_model = 1, M = 1, lambda_f = LAMBDA, t_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

max = 20 # max mean spillover value
step = 0.1 # step size -- must be an integer for N-M models
long_max = 100 # maximum value for C:1 plots
long_step = 5 # step size for long C:1 line values
C = 1; # past:future spillover ratio -- slope of C:1 line

vals <- rbind(
  expand.grid(seq(0, max, by = step), seq(0, max, by = step)),
  data.frame(Var1 = seq(0, long_max, long_step), 
             Var2 = C * seq(0, long_max, long_step))
) %>%
  dplyr::distinct()  # removes duplicates

names(vals) <- c('lambda', 'lambda_f')

####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

## define beta mixture distribution
mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(10),
                            weights = c(1) )

## generate dataframe for prior plots (fig. 2A-C)
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

## define beta mixture distribution
mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(0.1),
                            weights = c(1) )

## bind to dataframe for prior plots (fig. 2A-C)
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

## define beta mixture distribution
mixtureDistn <- data.frame( a = c(0.1, 100),
                            b =  c(10, 400),
                            weights = c(0.8, 0.2) )

## bind to dataframe for prior plots (fig. 2A-C)
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

###############################


#######################
##  Generate Prior   ##
## Plots (Fig. 2A-C) ##
#######################

y_upr = max(prior$density[!is.infinite(prior$density)])
{
  ## SCENARIO 1 ----
  tmp = prior[prior$prior %in% "Scenario 1",]
  pr1 <- ggplot() + 
    geom_line(data = tmp, aes(x = phi, y = density), color = viridis(n = 6)[3], linewidth = 1.5) +
    labs(
      x = " ",
      y = expression(atop(phantom("Density"), "Density"))
    ) + 
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c('0','0.25','0.5','0.75','1')) +
    scale_y_continuous(breaks = c(0,20,40,60), labels = c('0','20','40','60'),
                       limits = c(0, y_upr)) +
    annotate(
      "text", x = 0.4, y = 40,
      label = "paste(pi, '(', phi, ') = Beta(0.1, 10)')",
      parse = TRUE, color = "black", alpha = 1, size = 18 / 2.845
    ) +
    theme_classic() + 
    theme(
      panel.border = element_blank(),
      plot.margin = unit(c(0.5,0,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 24),
      axis.text.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      axis.text.y = element_text(size = 24),
      aspect.ratio = 1.0
    )
  
  ## SCENARIO 2 ----
  tmp = prior[prior$prior %in% "Scenario 2",]
  pr2 <- ggplot() + 
    geom_line(data = tmp, aes(x = phi, y = density), color = viridis(n = 6)[3], linewidth = 1.5) +
    labs(
      x = expression(paste("Per-spillover host jump probability (", phi, ")")),
      y = " "
    ) + 
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c('0','0.25','0.5','0.75','1')) +
    scale_y_continuous(breaks = c(0,20,40,60), labels = c('0','20','40','60'),
                       limits = c(0, y_upr)) +
    annotate(
      "text", x = 0.4, y = 40,
      label = "paste(pi, '(', phi, ') = Beta(0.1, 0.1)')",
      parse = TRUE, color = "black", alpha = 1, size = 18 / 2.845
    ) +
    theme_classic() + 
    theme(
      panel.border = element_blank(),
      plot.margin = unit(c(0.5,0,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 24),
      axis.text.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      #axis.text.y = element_text(size = 24),
      axis.text.y = element_blank(),
      aspect.ratio = 1.0
    )
  
  ## SCENARIO 3 ----
  tmp = prior[prior$prior %in% "Scenario 3",]
  pr3 <- ggplot() + 
    geom_line(data = tmp, aes(x = phi, y = density), color = viridis(n = 6)[3], linewidth = 1.5) +
    labs(
      x = " ",
      y = " "
    ) + 
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c('0','0.25','0.5','0.75','1')) +
    scale_y_continuous(breaks = c(0,20,40,60), labels = c('0','20','40','60'),
                       limits = c(0, y_upr)) +
    annotate(
      "text", x = 0.5, y = 40,
      label = "atop(pi(phi) == 0.8 %.% 'Beta(0.1, 10)', '+' ~ 0.2 %.% 'Beta(100, 400)')",
      parse = TRUE, color = "black", alpha = 1, size = 18 / 2.845
    ) +
    theme_classic() + 
    theme(
      panel.border = element_blank(),
      plot.margin = unit(c(0.5,0,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 24),
      axis.text.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      #axis.text.y = element_text(size = 24),
      axis.text.y = element_blank(),
      aspect.ratio = 1.0
    )
  
}

wrap_plots(list(pr1, pr2, pr3), ncol = 3)

(pr1 + pr2 + pr3) +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = '',
    tag_suffix = '',
    tag_position = 'topleft'  
  ) &
  theme(
    plot.tag = element_text(size = 22, face = "bold", hjust = 0, vjust = 1),
    plot.tag.position = c(0.9, 0.98)  # fine-tune if needed
  )
####################


###################
## Generate main ##
## results plots ##
###################
brk = seq(from = 0, to = max(max(run_out$exact), max(run_out$prob), max(run_out$sol)), by = 0.05)[-1] # contour value breakpoints
rng = 0 # range around the C:1 line as plot points
cond_l = C*run_out$future_val == run_out$past_val ## C:1 line -- i.e. future horizon is C times past horizon
cond_long = paste(run_out$past_val, run_out$future_val) %in% 
  paste(seq(0, long_max, long_step), C * seq(0, long_max, long_step))
cond = (C*run_out$future_val <= (run_out$past_val + rng)) & (C*run_out$future_val >= (run_out$past_val - rng)) ## range of values around C:1 line
cond_short = (C*run_out$future_val == run_out$past_val) & (run_out$past_val <= max)  ## C:1 line values from heatmap
cond_short_P = paste(run_out$past_val, run_out$future_val) %in% 
  paste(seq(0, max, 1), C * seq(0, max, 1))

savePlots = FALSE ## option to save plots to ./plots in directory. Note that this will overwrite existing plots from other runs

## generate heatmap plots (Fig. 2D-F)
## using .data[['sol']] gives the integral solution
## using .data[['exact']] gives the analytic solution via the confluent hypergeometric function or the gamma function (count model)
## using .data[['prob']] gives the simulation values (default 0 if pars$runSims == FALSE)
{
  tmp = run_out[run_out$scenario %in% 'Scenario 1',] ## subset by scenario
  tmp = tmp[(tmp$past_val <= max) & (tmp$future_val <= max),] ## exclude long values
  cond_tmp = tmp$future_val == tmp$past_val ## C:1 line -- i.e. future horizon is C times past horizon
  
  int_map1 <- ggplot() +
    scale_fill_viridis(option = "plasma", discrete = F, breaks = c(0, 0.05, 0.1), labels = function(x) round(x, digits = 2)) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['sol']])) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('sol')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('exact')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_line(data = tmp[cond_tmp,], aes(x = past_val, y = future_val), color = "#FFFFFF", alpha = .6, linewidth = 1.25) +
    scale_x_continuous(breaks = seq(min(tmp$future_val), max(tmp$past_val), 5)) +
    scale_y_continuous(breaks = seq(from = min(tmp$future_val), to = max(tmp$future_val), by = 5)) +
    annotate(
      "text", x = 9, y = 10.2, 
      label = "(lambda*T[P] == c*lambda*T[F])", 
      parse = TRUE, color = "white", alpha = 0.6, size = 14 / 2.845, angle = 45
    ) +
    labs( 
      x = " ",
      y = expression(atop("Expected number of", paste("future spillovers (", lambda, T[F], ")"))),
      fill = expression( paste(P(H[F]>0)) ) 
    ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(
      title.vjust = 1.2,
      barwidth = unit(4, "cm"),
      barheight = unit(0.6, "cm")
      )
    ) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          plot.margin = unit(c(0.5,0,0,0), 'lines'),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 24),
          axis.text.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          axis.text.y = element_text(size = 24),
          legend.title = element_text(size = 20),
          legend.position = "bottom",
          legend.title.align = 0.5,
          legend.text = element_text(size = 20),
          legend.key.width = unit(1.1, "cm"),
          aspect.ratio = 1.0)
  
  
  tmp = run_out[run_out$scenario %in% 'Scenario 2',] ## subset by scenario
  tmp = tmp[(tmp$past_val <= max) & (tmp$future_val <= max),] ## exclude long values
  cond_tmp = tmp$future_val == tmp$past_val ## C:1 line -- i.e. future horizon is C times past horizon
  int_map2 <- ggplot() +
    scale_fill_viridis(option = "plasma", discrete = F, breaks = c(0, 0.3, 0.6), labels = function(x) round(x, digits = 2)) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['sol']])) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('sol')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('exact')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_line(data = tmp[cond_tmp,], aes(x = past_val, y = future_val), color = "#FFFFFF", alpha = .6, linewidth = 1.25) +
    scale_x_continuous(breaks = seq(min(tmp$future_val), max(tmp$past_val), 5)) +
    scale_y_continuous(breaks = seq(from = min(tmp$future_val), to = max(tmp$future_val), by = 5)) +
    annotate(
      "text", x = 9, y = 10.2, 
      label = "(lambda*T[P] == c*lambda*T[F])", 
      parse = TRUE, color = "white", alpha = 0.6, size = 14 / 2.845, angle = 45
    ) + 
    labs( 
      x = expression(paste("Expected number of past spillovers (", lambda, T[P], ")")),
      fill = expression( paste(P(H[F]>0)) ) 
      ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(
      title.vjust = 1.2,
      barwidth = unit(4, "cm"),
      barheight = unit(0.6, "cm")
      )
    ) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          plot.margin = unit(c(0.5,0,0,0), 'lines'),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 24),
          axis.text.x = element_text(size = 24),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 24, color = "#FFFFFF"),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 20),
          legend.position = "bottom",
          legend.title.align = 0.5,
          legend.text = element_text(size = 20),
          legend.key.width = unit(1.1, "cm"),
          aspect.ratio = 1.0)
  
  
  tmp = run_out[run_out$scenario %in% 'Scenario 3',] ## subset by scenario
  tmp = tmp[(tmp$past_val <= max) & (tmp$future_val <= max),] ## exclude long values
  cond_tmp = tmp$future_val == tmp$past_val*C ## C:1 line -- i.e. future horizon is C times past horizon
  int_map3 <- ggplot() +
    scale_fill_viridis(option = "plasma", discrete = F, breaks = c(0, 0.1, 0.2), labels = function(x) round(x, digits = 2)) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['sol']])) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('sol')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_contour(data = tmp, aes(x = past_val, y = future_val, z = .data[[c('exact')]]), color= '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.25) +
    geom_line(data = tmp[cond_tmp,], aes(x = past_val, y = future_val), color = "#FFFFFF", alpha = .6, linewidth = 1.25) +
    scale_x_continuous(breaks = seq(min(tmp$future_val), max(tmp$past_val), 5)) +
    scale_y_continuous(breaks = seq(from = min(tmp$future_val), to = max(tmp$future_val), by = 5)) +
    annotate(
      "text", x = 9, y = 10.2, 
      label = "(lambda*T[P] == c*lambda*T[F])", 
      parse = TRUE, color = "white", alpha = 0.6, size = 14 / 2.845, angle = 45
    ) + 
    labs( 
      fill = expression( paste(P(H[F]>0)) ), 
      x = " "
      ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(
      title.vjust = 1.2,
      barwidth = unit(4, "cm"),
      barheight = unit(0.6, "cm")
      )
    ) +
    theme(panel.grid = element_blank(), panel.border = element_blank(), 
          plot.margin = unit(c(0.5,0,0,0), 'lines'),
          plot.title = element_blank(),
          axis.title.x = element_text(size = 24),
          axis.text.x = element_text(size = 24),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 24, color = "#FFFFFF"),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 20),
          legend.position = "bottom",
          legend.title.align = 0.5,
          legend.text = element_text(size = 20),
          legend.key.width = unit(1.2, "cm"),
          aspect.ratio = 1.0)
  
  
  ## facet plot using a single color scale
  tmp = run_out[(run_out$past_val <= max) & (run_out$future_val <= max),]
  brk = log10(brk)
  int_maps <- ggplot() +
    scale_fill_viridis(
      option = "plasma", discrete = FALSE, 
      #breaks = c(0, 0.3, 0.6), 
      labels = function(x) round(x, digits = 2)
    ) +
    geom_tile(data = tmp, aes(x = past_val, y = future_val, fill = .data[['sol']])) +
    geom_contour(
      data = tmp, 
      aes(x = past_val, y = future_val, z = .data[['sol']]), 
      color = '#041E42', breaks = brk, alpha = 0.6, linewidth = 1.0
    ) +
    geom_contour(
      data = tmp,
      aes(x = past_val, y = future_val, z = .data[['exact']]),
      color = '#FFFFFF', breaks = brk, alpha = 0.6, linewidth = 1.0, linetype = 'dashed'
    ) +
    geom_line(
      data = tmp[cond_tmp,], 
      aes(x = past_val, y = future_val), 
      color = "#FFFFFF", alpha = .6, linewidth = 1.25
    ) +
    facet_wrap(~scenario, scales = "fixed") + 
    scale_x_continuous(
      breaks = seq(min(tmp$past_val), max, 5)
    ) +
    scale_y_continuous(
      breaks = seq(from = min(tmp$future_val), to = max, by = 5)
    ) +
    labs(
      x = expression(paste("Expected number of past spillovers (", lambda, T[P], ")")),
      y = expression(atop("Expected number of", paste("future spillovers (", lambda, T[F], ")"))),
      fill = expression(P(H[F] > 0))
    ) +
    annotate(
      "text", x = 9, y = 10.2, 
      label = "(lambda*T[P] == c*lambda*T[F])", 
      parse = TRUE, color = "white", alpha = 0.6, size = 10 / 2.845, angle = 45
    ) + 
    coord_cartesian(expand = FALSE) +
    theme_bw() + 
    guides(fill = guide_colorbar(
      title.vjust = 1.2,
      barwidth = unit(4, "cm"),
      barheight = unit(0.6, "cm")
      )
    ) +
    theme(
      panel.grid = element_blank(), 
      panel.border = element_blank(), 
      plot.title = element_blank(),
      axis.title.x = element_text(size = 14),
      axis.text.x  = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.y  = element_text(size = 14),
      legend.title = element_text(size = 14, vjust = -0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.2, "cm"),
      aspect.ratio = 1.0,
      panel.spacing = unit(1, "lines")
    )
  rm(tmp)
  int_maps
}



int_map1

int_map2

int_map3

wrap_plots(list(int_map1, int_map2, int_map3), ncol = 3)

## generate 1:1 plots (Fig. 2G-I) 
fname_tag = ""
if(pars$runSims){ fname_tagL = paste0(fname_tag, '_Sims') }else{ fname_tagL = paste0(fname_tag, '_noSims') } ## change filename tag for simulation values

{
  y_upr = max( c(run_out[cond_l,]$prob,run_out[cond_l,]$sol, run_out[cond_l,]$prob, 1 - (1/(C+1)^(pars$a))))*1.01 
  
  L_curve1 <- ggplot() +
    geom_line(data = run_out[cond_l,][run_out[cond_l,]$scenario == 'Scenario 1',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_long,][run_out[cond_long,]$scenario == 'Scenario 1',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, long_max, by=long_max) ) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +
    labs( 
      x = expression( paste("spillover rate (", lambda, ")") ),
      y = expression( paste(P(H[F]>0)) ) 
    ) + 
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 20),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 20),
      aspect.ratio = 1.0)
  
  curve1 <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == 'Scenario 1',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_short_P & run_out$scenario == 'Scenario 1',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    labs( 
      x = expression( paste("spillover rate (", lambda, ")") ),
      y = expression( atop( "Probability of a host", paste("jump  ", P(H[F]>0)) )) 
    ) + 
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 24, color = "#FFFFFF"),
      axis.text.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      axis.text.y = element_text(size = 24),
      aspect.ratio = 1.0)
  
  
  L_curve2 <- ggplot() +
    geom_line(data = run_out[cond_l,][run_out[cond_l,]$scenario == 'Scenario 2',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_long,][run_out[cond_long,]$scenario == 'Scenario 2',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,long_max,by=long_max) ) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    labs( 
      x = expression( paste("spillover rate (", lambda, ")") ),
      y = expression( paste(P(H[F]>0)) ) 
    ) + 
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 20),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 20),
      aspect.ratio = 1.0)

  
  curve2 <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == 'Scenario 2',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_short_P & run_out$scenario == 'Scenario 2',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    labs( 
      x = expression( paste("Rate of spillover (", lambda, ")") ),
      y = " " 
    ) + 
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 24),
      axis.text.x = element_text(size = 24),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      aspect.ratio = 1.0)
  
  L_curve3 <- ggplot() +
    geom_line(data = run_out[cond_l,][run_out[cond_l,]$scenario == 'Scenario 3',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_long,][run_out[cond_long,]$scenario == 'Scenario 3',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,long_max,by=long_max) ) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    theme_classic() +
    labs( 
      x = expression( paste("spillover rate (", lambda, ")") ),
      y = expression( paste(P(H[F]>0)) ) 
    ) + 
    theme(
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 20),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 20),
      aspect.ratio = 1.0)
  
  
  curve3 <- ggplot() +
    geom_line(data = run_out[cond_short & run_out$scenario == 'Scenario 3',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', linewidth = 1.25) +
    geom_point(data = run_out[cond_short_P & run_out$scenario == 'Scenario 3',], aes(x = (past_val), y = .data[['sol']]), color = '#0C2340', size = 2.5) +
    geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0, max, by = 5)) +
    scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
    labs( 
      x = " ",
      y = " " 
    ) + 
    theme_classic() +
    theme(
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 24, color = "#FFFFFF"),
      axis.text.x = element_text(size = 24),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      aspect.ratio = 1.0)
}

L_curve1
curve1

L_curve2
curve2

L_curve3
curve3

wrap_plots(list(curve1, curve2, curve3), ncol = 3)
wrap_plots(list(L_curve1, L_curve2, L_curve3), ncol = 3)

library(ggplot2)
library(cowplot)
library(patchwork)

#--- Step 1: Combine main plots (no insets yet)
base_plot <- wrap_plots(list(curve1, curve2, curve3), ncol = 3)

#--- Step 2: Draw that combined layout into a drawable canvas
# Using cowplot::ggdraw allows absolute positioning over the *entire combined figure*
final_plot <- ggdraw(base_plot)

#--- Step 3: Define inset size and positions (in relative figure coordinates)
# Values between 0 and 1 refer to the full figure area (not per-panel)
# Adjust inset_x and inset_y to fine-tune placement
inset_width  <- 0.3   # 18% of full figure width
inset_height <- 0.3    # 30% of full figure height
y_offset     <- 0.55   # vertical offset from bottom
x_positions  <- c(0.32, 0.61, .898)  # roughly the horizontal centers of 3 panels

#--- Step 4: Overlay each inset at a fixed position
final_plot <- final_plot +
  draw_plot(L_curve1, x = x_positions[1] - inset_width / 2, y = y_offset,
            width = inset_width, height = inset_height) +
  draw_plot(L_curve2, x = x_positions[2] - inset_width / 2, y = y_offset,
            width = inset_width, height = inset_height) +
  draw_plot(L_curve3, x = x_positions[3] - inset_width / 2, y = y_offset,
            width = inset_width, height = inset_height)

final_plot


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
  scale_x_continuous(breaks = seq(0,max,by=5)) +
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


##########################
## Generate figure 2J-L ##
##  P(HJ) vs.  T_P for  ##
##  multiple  rates of  ##
##       spillover      ##
##########################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## past/future _model: 0 = count model, 1 = poisson model
## prior_type = 2 uses the mixture model which effectively also captures the single beta prior
LAMBDA = 1
lambda_vec = c( 5, 10, 100)
pars = data.frame(past_model = 1, N = 1, lambda = LAMBDA, t_p = 1,
                  future_model = 1, M = 1, lambda_f = LAMBDA, t_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

max = 20 # max mean spillover value
step = 0.1 # step size -- must be an integer for N-M models
C = 1; pars$lambda_f = C*pars$lambda_f; # past:future spillover ratio -- slope of C:1 line

##vals = data.frame(c(seq(0,20,0.5), seq(20,100,2)[-1]), c(seq(0,20,0.5), seq(20,100,2)[-1]))
vals = expand.grid(seq(0,max,step), c(5))
names(vals) = c('t_p', 't_f') ## names must be adjusted to match parameter name in "pars" data frame

####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = 0.1,
                            b =  10.0,
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
run_out$lambda = pars$lambda
for(i in 1:length(lambda_vec)){
  pars$lambda = lambda_vec[i]
  pars$lambda_f = C*lambda_vec[i]
  tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
  tmp2$scenario = "Scenario 1"
  tmp2$lambda = pars$lambda
  run_out <- rbind(run_out, tmp2)
}
rm(tmp2)
pars$lambda = LAMBDA ## reset lambda to initial value
pars$lambda_f = C*LAMBDA ## reset lambda to initial value
####################

####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(0.1),
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
tmp$lambda = pars$lambda
for(i in 1:length(lambda_vec)){
  pars$lambda = lambda_vec[i]
  pars$lambda_f = C*lambda_vec[i]
  tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
  tmp2$scenario = "Scenario 2"
  tmp2$lambda = pars$lambda
  tmp <- rbind(tmp, tmp2)
}
run_out <- rbind(run_out, tmp); rm(tmp); rm(tmp2)
pars$lambda = LAMBDA ## reset lambda to initial value
pars$lambda_f = C*LAMBDA ## reset lambda to initial value
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
tmp$lambda = pars$lambda
for(i in 1:length(lambda_vec)){
  pars$lambda = lambda_vec[i]
  pars$lambda_f = C*lambda_vec[i]
  tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
  tmp2$scenario = "Scenario 3"
  tmp2$lambda = pars$lambda
  tmp <- rbind(tmp, tmp2)
}
run_out <- rbind(run_out, tmp); rm(tmp); rm(tmp2)
pars$lambda = LAMBDA ## reset lambda to initial value
pars$lambda_f = C*LAMBDA ## reset lambda to initial value
####################

sum(abs(run_out$sol-run_out$exact)) ## total error between numerical and exact solution over all scenarios
max(abs(run_out$sol-run_out$exact)) ## max point of error

###############################

###################
## Generate main ##
## results plots ##
###################
brk = seq(from = 0, to = max(max(run_out$exact), max(run_out$prob), max(run_out$sol)), by = 0.05)[-1] # contour value breakpoints
rng = 0 # range around the C:1 line as plot points
lim = 1-(C+1)^(-pars$a)
y_upr = max( c(run_out$prob, run_out$sol, run_out$prob, 1 - (1/(C+1)^(pars$a))))*1.01 
cond = (C*run_out$future_val <= (run_out$past_val + rng)) & (C*run_out$future_val >= (run_out$past_val - rng)) ## range of values around C:1 line
cond_short = (run_out$past_val <= 5)
cond_mid = (run_out$past_val <= 100)
run_out$lambda = factor(run_out$lambda, levels = c(1,5,10,100))
run_out$future_val = factor(run_out$future_val, levels = c(20))

fname_tag = "Std"

savePlots = TRUE ## option to save plots to ./plots in directory. Note that this will overwrite existing plots from other runs

## generate figures 2J-L
{
  brks = seq(0,by = round((y_upr+.05)/2, digits = 1), length = 3)
  labamba_curve1 <- ggplot() +
    scale_color_manual(values = rev(sequential_hcl(7, palette = "Heat")[1:(length(lambda_vec)+1)])) +
    geom_line(data = run_out[run_out$scenario == "Scenario 1", ], aes(x = past_val, y = .data[['exact']], color = lambda), linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,20,by=5)) +
    scale_y_continuous(breaks = brks, limits = c(0, max(brks))) +     
    theme_classic() +
    labs(x = expression( " "), 
         y = expression( atop( "Probability of a host", paste("jump  ", P(H[F]>0)) )), 
         color = expression(lambda)) + 
    theme(
      legend.position = "none",
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 24, color = "#FFFFFF"),
      axis.text.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      axis.text.y = element_text(size = 24),
      aspect.ratio = 1.0)
  
  labamba_curve2 <- ggplot() +
    scale_color_manual(values = rev(sequential_hcl(7, palette = "Heat")[1:(length(lambda_vec)+1)])) +
    geom_line(data = run_out[run_out$scenario == "Scenario 2", ], aes(x = past_val, y = .data[['exact']], color = lambda), linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,20,by=5)) +
    scale_y_continuous(breaks = brks, limits = c(0, max(brks))) +     
    theme_classic() +
    labs(x = expression( "Past spillover window "*(T[P])), 
         y = expression( atop( "Probability of a host", paste("jump  ", P(H[F]>0)) )), 
         color = expression(lambda)) + 
    theme(
      legend.position = "none",
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 24),
      axis.text.x = element_text(size = 24),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      aspect.ratio = 1.0)
  
  labamba_curve3 <- ggplot() +
    scale_color_manual(values = rev(sequential_hcl(7, palette = "Heat")[1:(length(lambda_vec)+1)])) +
    geom_line(data = run_out[run_out$scenario == "Scenario 3", ], aes(x = past_val, y = .data[['exact']], color = lambda), linewidth = 1.25) +
    scale_x_continuous(breaks = seq(0,20,by=5)) +
    scale_y_continuous(breaks = brks, limits = c(0, max(brks))) +     
    theme_classic() +
    labs(x = expression( " "), 
         y = expression( atop( "Probability of a host", paste("jump  ", P(H[F]>0)) )), 
         color = expression(lambda)) + 
    theme(
      plot.margin = unit(c(0.2,1,0,0), 'lines'),
      plot.title = element_blank(),
      axis.title.x = element_text(size = 24, color = "#FFFFFF"),
      axis.text.x = element_text(size = 24),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none",
      # legend.title = element_text(size = 20),
      # legend.text = element_text(size = 18),
      aspect.ratio = 1.0)
  
}

wrap_plots(list(labamba_curve1, labamba_curve2, labamba_curve3), ncol = 3)

labamba_curve1
labamba_curve2
labamba_curve3



if(pars$runSims){ fname_tagL = paste0(fname_tag, '_Sims') }else{ fname_tagL = paste0(fname_tag, '_noSims') } ## change filename tag for simulation values

## facet version of fig 2J-L on x scale T_P = 0:20
L_curves <- ggplot() +
  scale_color_manual(values = rev(sequential_hcl(7, palette = "Heat")[1:(length(lambda_vec)+1)])) +
  # scale_color_manual(values = viridis(5, option = "H")[1:(length(lambda_vec)+1)] ) +
  geom_line(data = run_out, aes(x = past_val, y = .data[['exact']], color = lambda), linewidth = 1.25) +
  #geom_point(data = run_out, aes(x = past_val, y = .data[['exact']], color = lambda), size = 2.5) +
  scale_x_continuous(breaks = seq(0,20,by=5)) +
  scale_y_continuous(breaks = seq(0,by = round((y_upr+.05)/2, digits = 1), length = 3)) +     
  facet_wrap(~scenario) + 
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1) +
  labs(x = expression( "Past spillover window "*(T[P])), 
       y = expression( atop( "Probability of a host", paste("jump  ", P(H[F]>0)) )), 
       color = expression(lambda)) + 
  theme(
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)

## facet version of fig 2J-L on x scale T_P = 0:5
curves <- ggplot() +
  scale_color_manual(values = rev(sequential_hcl(7, palette = "Heat")[1:(length(lambda_vec)+1)])) +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']], color = lambda), linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,5,by=1)) +
  scale_y_continuous(breaks = seq(0,by = round((y_upr+.05)/2, digits = 1), length = 3)) +     
  facet_wrap(~scenario) + 
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1) +
  labs(x = expression( "Past spillover window "*(T[P])), 
       y = expression( atop( "Probability of a host", paste("jump  ", P(H[F]>0)) )), 
       color = expression(lambda)) + 
  theme(
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)

L_curves
curves

if(savePlots){
  ggsave(filename = paste("./figures/TP_figs.png"), plot = curves, width = 10.00, height = 5.5, units = "in")
  ggsave(filename = paste("./figures/TP_figs_long.png"), plot = L_curves, width = 10.00, height = 5.5, units = "in")
}



###
### Extract values for figure R1
###

x = as.data.frame(matrix(nrow=0,ncol=4)); 
names(x) <- c('lambda','t_spill','T_P','y_pos')

T_F = 10
lambda_lo = 0.8
lambda_hi = 2
TP_vec = c(2,10)
Lsp = 0.001

# Assign y-positions for each T_P and lambda combo
y_map = rbind(
  data.frame(T_P = 10, lambda = lambda_lo, y_pos = 1.0+Lsp*3),
  data.frame(T_P = 10, lambda = lambda_hi, y_pos = 1.0+Lsp*2),
  data.frame(T_P = 2, lambda = lambda_lo, y_pos = 1.0+Lsp),
  data.frame(T_P = 2, lambda = lambda_hi, y_pos = 1.0)
)

for(i in 1:nrow(y_map)){
  T_P = y_map$T_P[i]
  lambda_val = y_map$lambda[i]
  y_pos = y_map$y_pos[i]
  
  M = round(lambda_val * (T_P + T_F))
  if (M > 0){
    t_vals = seq(-1*T_P, T_F, length.out = M+2)[-1] 
    t_vals = t_vals[-length(t_vals)]
    x_F = data.frame(lambda = lambda_val, t_spill = t_vals, T_P = T_P, y_pos = y_pos)
  } else {
    x_F = data.frame(lambda = lambda_val, t_spill = (T_F - T_P)*0.5, T_P = T_P, y_pos = y_pos)
  }
  x <- rbind(x, x_F)
}

x$lambda_val <- as.factor(x$lambda)

df <- y_map
df$negLim = -1 * df$T_P
top = 1+Lsp*4

figR1 <- ggplot() +
  # geom_rect(data = x, aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf), fill = "gray60", alpha = 0.5) +
  # geom_rect(data = x, aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf), fill = "gray90", alpha = 0.5) +
  # geom_point(data = x, aes(x = t_spill, y = y_pos, color = lambda_val), size = 4, shape = 1, stroke = 1.2)  +
  # scale_color_manual(values = c('#4490FEFF', '#EA4F0DFF'), breaks = c("0.8","2")) +
  geom_point(data = x, aes(x = t_spill, y = y_pos), size = 4, shape = 1, stroke = 1.2)  +
  annotate("segment", x = 0, xend = 0, y = 1.0-Lsp, yend = top, size = 1, linetype = 'dashed', alpha = 0.95, color = "gray0") +
  geom_point(data = df, aes(x = negLim, y = y_pos), color = "#000000", size = 13, shape = "*") +
  geom_segment(data = df, aes(x = -1 * T_P, xend = T_F, y = y_pos, yend = y_pos), size = 1) +
  scale_x_continuous(limits = c(-1.15 * max(TP_vec), T_F)) +
  scale_y_continuous(limits = c(1.0-4*Lsp, top)) +
  # annotate("text", x = -2.5, y = 1+Lsp*4, label = "past", size = 6, hjust = 0.5) +
  # annotate("text", x = 2.5, y = 1+Lsp*4, label = "future", size = 6, hjust = 0.5) +
  # annotate("text", x = 0, y = top, label = "present", size = 6, hjust = 0.5) +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        panel.spacing = unit(1, "lines")
  )

figR1

dat = run_out[run_out$past_val %in% c(2, 15),]
dat = dat[dat$lambda %in% c(1,10),c(1,5,6,7)]
dat$exact = format(dat$exact, digits = 2) 
write.csv(dat, '~/Desktop/R1_tableRAW.csv')
dat = split(dat, dat$scenario)




################################################################################
################################################################################
################################################################################
################################################################################

##########################
## Test alternate plot  ##
## types for Fig. 2 J-L ##
##########################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## past/future _model: 0 = count model, 1 = poisson model
## prior_type = 2 uses the mixture model which effectively also captures the single beta prior
## model parameters -- not all are used depending on the values of past/future _model
## past/future _model: 0 = count model, 1 = poisson model
## prior_type = 2 uses the mixture model which effectively also captures the single beta prior
LAMBDA = 1
pars = data.frame(past_model = 1, N = 1, lambda = LAMBDA, t_p = 1,
                  future_model = 1, M = 1, lambda_f = LAMBDA, t_f = 20,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

max = 20 # max mean spillover value
step = 0.5 # step size -- must be an integer for N-M models
long_max = 100 # maximum value for C:1 plots
long_step = 5 # step size for long C:1 line values
C = 1; # past:future spillover ratio -- slope of C:1 line
tp_vals = c(1,2,5,20)

vals = expand.grid(c(seq(0, max, by = step), seq(max, long_max, by = long_step)[-1]), tp_vals)
names(vals) = c('lambda', 't_p')
vals$lambda_f = vals$lambda*C ## assume 1:C correlation between rates of past and future spillover 
names(vals) = c('lambda', 't_p', 'lambda_f') ## names must be adjusted to match parameter name in "pars" data frame


####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = 0.1,
                            b =  10.0,
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
names(run_out) = c('lambda', 't_p', 'lambda_f', 'prob', 'sol', 'exact')
run_out$scenario = "Scenario 1"
####################

####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(1/10),
                            b =  c(1/10),
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
names(tmp) = c('lambda', 't_p', 'lambda_f', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 2"

run_out <- rbind(run_out, tmp); rm(tmp);
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
names(tmp) = c('lambda', 't_p', 'lambda_f', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 3"

run_out <- rbind(run_out, tmp); rm(tmp)
####################

sum(abs(run_out$sol-run_out$exact)) ## total error between numerical and exact solution over all scenarios
max(abs(run_out$sol-run_out$exact)) ## max point of error
which(abs(run_out$sol-run_out$exact) == max(abs(run_out$sol-run_out$exact)))

###############################

###################
## Generate main ##
## results plots ##
###################
brk = seq(from = 0, to = max(max(run_out$exact), max(run_out$prob), max(run_out$sol)), by = 0.05)[-1] # contour value breakpoints
rng = 0 # range around the C:1 line as plot points
lim = 1-(C+1)^(-pars$a)
y_upr = max( c(run_out$prob, run_out$sol, run_out$prob, 1 - (1/(C+1)^(pars$a))))*1.01 
cond_short = (run_out$lambda <= 5)
cond_mid = (run_out$lambda <= 20)
run_out$t_p = factor(run_out$t_p, levels = tp_vals)
# run_out$lambda = factor(run_out$lambda, levels = c(1,5,10,100))
# run_out$future_val = factor(run_out$future_val, levels = c(20))

fname_tag = "Std"

savePlots = TRUE ## option to save plots to ./plots in directory. Note that this will overwrite existing plots from other runs

## generate 1:1 plots 

if(pars$runSims){ fname_tagL = paste0(fname_tag, '_Sims') }else{ fname_tagL = paste0(fname_tag, '_noSims') } ## change filename tag for simulation values

L_curves <- ggplot() +
  scale_color_manual(values = rev(sequential_hcl(8, palette = "Heat")[2:(length(tp_vals)+1)])) +
  # scale_color_manual(values = viridis(5, option = "H")[1:(length(lambda_vec)+1)] ) +
  geom_line(data = run_out, aes(x = lambda_f, y = .data[['exact']], color = t_p), linewidth = 1.25) +
  # geom_point(data = run_out, aes(x = lambda_f, y = .data[['exact']], color = t_p), size = 2.5) +
  scale_x_continuous(breaks = seq(0, long_max, by=25)) +
  #scale_y_continuous(breaks = seq(0, by = round((y_upr+.05)/2, digits = 1), length = 3)) +     
  scale_y_continuous(breaks = seq(0,by = round((y_upr+.05)/2, digits = 1), length = 3), limits = c(0, 0.8)) +     
  facet_wrap(~scenario) + 
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1) +
  labs(x = expression("Rate of Spillover ("*lambda[p]*"="*lambda[f]*")"), 
       y = expression( paste(P(H[F]>0)) ), 
       color = expression(t[P])) + 
  theme(
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


curves <- ggplot() +
  scale_color_manual(values = rev(sequential_hcl(8, palette = "Heat")[2:(length(tp_vals)+1)])) +
  geom_line(data = run_out[cond_mid,], aes(x = lambda_f, y = .data[['exact']], color = t_p), linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,20,by=5)) +
  #scale_y_continuous(breaks = seq(0,by = round((y_upr+.05)/2, digits = 1), length = 3)) +     
  scale_y_continuous(breaks = seq(0,by = round((y_upr+.05)/2, digits = 1), length = 3), limits = c(0, 0.8)) +     
  facet_wrap(~scenario) + 
  theme_classic() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = 1) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth = 1) +
  labs(x = expression("Rate of Spillover ("*lambda[p]*"="*lambda[f]*")"), 
       y = expression( paste(P(H[F]>0)) ), 
       color = expression(t[P])) + 
  theme(
    panel.spacing = unit(1, "lines"),
    axis.line = element_line(),
    strip.text = element_blank(),
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)

L_curves
curves

if(savePlots){
  ggsave(filename = paste("./figures/TP_figs_v2.png"), plot = curves, width = 10.00, height = 5.5, units = "in")
  ggsave(filename = paste("./figures/TP_figs_long_v2.png"), plot = L_curves, width = 10.00, height = 5.5, units = "in")
}

################################################################################################################################################################################################################################
################################################################################################################################################################################################################################
################################################################################################################################################################################################################################
################################################################################################################################################################################################################################

## Changing C vs changing t_p prior to asymptotic behavior

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## past/future _model: 0 = count model, 1 = poisson model
## prior_type = 2 uses the mixture model which effectively also captures the single beta prior
C_vec = c(0.01, 0.1, 2, 5, 10, 100)
pars = data.frame(past_model = 1, N = 1, lambda = 1, t_p = 1,
                  future_model = 1, M = 1, lambda_f = 1, t_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

max = 10 # max mean spillover value
step = 0.2 # step size -- must be an integer for N-M models
long_max = 10 # maximum value for C:1 plots
long_step = 5 # step size for long C:1 line values
C = 1; pars$lambda_f = C*pars$lambda_f; # past:future spillover ratio -- slope of C:1 line

## vals = data.frame(c(seq(0,20,0.5), seq(20,100,2)[-1]), c(seq(0,20,0.5), seq(20,100,2)[-1]))
## vals = expand.grid(c(0, 5,20), c(seq(0,20,0.1), seq(20,100,10)[-1], seq(100,1000,100)[-1]))
vals = expand.grid(c(0, 5, 20, 100), c(1, 5, 10, 100))
names(vals) = c('t_p', 't_f') ## names must be adjusted to match parameter name in "pars" data frame

####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = 0.1,
                            b =  10.0,
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
run_out$lambda = pars$lambda
for(i in 1:length(lambda_vec)){
  pars$lambda = lambda_vec[i]
  pars$lambda_f = C*lambda_vec[i]
  tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
  tmp2$scenario = "Scenario 1"
  tmp2$lambda = pars$lambda
  run_out <- rbind(run_out, tmp2)
}
rm(tmp2)
pars$lambda = LAMBDA ## reset lambda to initial value
pars$lambda_f = 1*LAMBDA ## reset lambda to initial value
####################

####################
## Define prior 2 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(1/10),
                            b =  c(1/10),
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
tmp$lambda = pars$lambda
for(i in 1:length(lambda_vec)){
  pars$lambda = lambda_vec[i]
  pars$lambda_f = C*lambda_vec[i]
  tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
  tmp2$scenario = "Scenario 2"
  tmp2$lambda = pars$lambda
  tmp <- rbind(tmp, tmp2)
}
run_out <- rbind(run_out, tmp); rm(tmp); rm(tmp2)
pars$lambda = LAMBDA ## reset lambda to initial value
pars$lambda_f = C*LAMBDA ## reset lambda to initial value
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
tmp$lambda = pars$lambda
for(i in 1:length(lambda_vec)){
  pars$lambda = lambda_vec[i]
  pars$lambda_f = C*lambda_vec[i]
  tmp2 <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
  names(tmp2) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
  tmp2$scenario = "Scenario 3"
  tmp2$lambda = pars$lambda
  tmp <- rbind(tmp, tmp2)
}
run_out <- rbind(run_out, tmp); rm(tmp); rm(tmp2)
pars$lambda = LAMBDA ## reset lambda to initial value
pars$lambda_f = C*LAMBDA ## reset lambda to initial value
####################

sum(abs(run_out$sol-run_out$exact)) ## total error between numerical and exact solution over all scenarios
max(abs(run_out$sol-run_out$exact)) ## max point of error
which(abs(run_out$sol-run_out$exact) == max(abs(run_out$sol-run_out$exact)))

###############################
run_out$lambda = factor(run_out$lambda, levels = c(0.01, 0.1, 1, 2, 5, 10, 100, 500))
run_out$past_val = factor(run_out$past_val, levels = c(0, 5, 20, 100))
run_out$future_val = factor(run_out$future_val, levels = c(1, 5, 10, 100))
run_out_0 = run_out[run_out$past_val != 0, ]

## grouped by t_past
ggplot(data = run_out, aes(x = future_val, y = exact, color = past_val)) +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  facet_wrap(scenario~., dir = "v") 

ggplot(data = run_out_0, aes(x = future_val, y = exact, color = past_val)) +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  facet_wrap(scenario~., dir = "v") 

## grouped by lambda
ggplot(data = run_out, aes(x = future_val, y = exact, color = lambda)) +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  facet_wrap(scenario~., dir = "v") 

ggplot(data = run_out_0, aes(x = future_val, y = exact, color = lambda)) +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_boxplot() +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  facet_wrap(scenario~., dir = "v") 



####################################################################################
####################################################################################
####################################################################################
####################################################################################


############################
##  Evaluate  effects of  ##
## increasing b parameter ##
##    for SI Fig. S4 B    ##
############################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 1, N = 1, lambda = 1, t_p = 1,
                  future_model = 1, M = 1, lambda_f = 1, t_f = 1,
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
names(vals) = c('lambda', 'lambda_f')
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
cond_short = (run_out$past_val <= 20)
cond_mid = (run_out$past_val <= 100)
run_out$b = factor(run_out$b, levels = c(1000, 100, 50, 10, 1, 0.1, 0.05, 0.01, 0.001))
fname_tag = "Std"

savePlots = TRUE ## option to save plots to ./plots in directory. Note that this will overwrite existing plots from other runs

## generate 1:1 plots 

if(pars$runSims){ fname_tagL = paste0(fname_tag, '_Sims') }else{ fname_tagL = paste0(fname_tag, '_noSims') } ## change filename tag for simulation values

L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_mid,], aes(x = past_val, y = .data[['sol']], color = b), linewidth = 1.25) +
  geom_point(data = run_out[cond_mid,], aes(x = past_val, y = .data[['exact']], color = b), size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,100,by=25)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    aspect.ratio = 1.0)


curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['sol']], color = b), linewidth = 1.25) +
  geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']], color = b), size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~scenario) + 
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
##  (Fig. S5 B)   ##
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
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
##     for SI Fig. S4 A     ##
##############################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 1, N = 1, lambda = 1, t_p = 1,
                  future_model = 1, M = 1, lambda_f = 1, t_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

spill_max = 20 # max mean spillover value
spill_step = 0.5 # step size -- must be an integer for N-M models
long_max = 1000 # maximum value for C:1 plots
long_step = 5 # step size for long C:1 line values
C = 1 # past:future spillover ratio -- slope of C:1 line

vec = c(seq(0,spill_max,by=spill_step), seq(spill_max, long_max, by = long_step)[-1])
vals = data.frame(vec, C*vec)
names(vals) = c('lambda', 'lambda_f')
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
##   Fig. S4 A   ##
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
  #xlab('Number of spillovers (N=M)') +
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
##  (Fig. S5 A)   ##
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
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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


#############################################
#############################################
####     Runs and plots for different    ####
####   past-future spillover slopes (c)  ####
####           (SI Fig. S5 C)            ####
#############################################
#############################################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## note that past/future _model are designated using C++ indexing
pars = data.frame(past_model = 1, N = 1, lambda = 1, t_p = 1,
                  future_model = 1, M = 1, lambda_f = 1, t_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

long_max = 1000 # maximum value for 1:1 plots
long_step = 2 # step size for long C:1 line values
C = c(0.2, 0.5, 2, 5)
past_vals = c(seq(0,20,0.5)[-41], seq(20,100,by=2), seq(100,long_max,by=5)[-1])
long_vals = data.frame(past = past_vals, future = past_vals, c = 1)

for(i in 1:length(C)){
  tmp = data.frame(past = past_vals, future = past_vals*C[i], c = C[i])
  long_vals = rbind(long_vals, tmp)
}
names(long_vals) = c('lambda', 'lambda_f', 'c')
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
names(long_out) = c('past_val', 'future_val', 'c', 'prob', 'sol', 'exact')
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
names(tmp) = c('past_val', 'future_val', 'c', 'prob', 'sol', 'exact')
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
names(tmp) = c('past_val', 'future_val', 'c', 'prob', 'sol', 'exact')
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

savePlots = TRUE
C = c(1, C) # add c=1 to vector of values for c to generate all limit lines
lim = 1-(1+as.numeric(long_out$c))^(-pars$a)
long_out$c = factor(long_out$c, levels = c(5,2,1,0.5,0.2))
y_upr = max( c(long_out$exact, long_out$sol, long_out$prob, 1 - (1/(C+1)^(pars$a))), na.rm=T)*1.01 
cond_mid = (long_out$past_val <= 100)

C_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = long_out[cond_mid,], aes(x = (past_val), y = .data[['sol']], color = c), linewidth = 1.25) +
  geom_point(data = long_out[cond_mid,], aes(x = past_val, y = .data[['sol']], color = c), size = 2.5) +
  geom_abline(slope = 0, intercept = 1 - (1/(C+1)^(pars$a)), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0, 100, 25 )) +
  ylim(0, y_upr) +
  facet_wrap(.~scenario) + 
  theme_classic() +
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
##  (Fig. S5 C)   ##
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
  geom_line(data = long_out[conv_cond,], aes(x = past_val, y = .data[['ROC']], color = c), linewidth = 1.25) +
  geom_point(data = long_out[conv_cond,], aes(x = past_val, y = .data[['ROC']], color = c), size = 2) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  facet_wrap(~scenario) + 
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
  geom_line(data = long_out[conv_cond,], aes(x = past_val, y = .data[['ROC2']], color = c), linewidth = 1.25) +
  geom_point(data = long_out[conv_cond,], aes(x = past_val, y = .data[['ROC2']], color = c), size = 2) +
  scale_x_continuous(breaks = seq(0,long_max,by=250)) +
  facet_wrap(~scenario) + 
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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


################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################



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
pars = data.frame(past_model = 1, N = 1, lambda = 1, t_p = 1,
                  future_model = 1, M = 1, lambda_f = 1, t_f = 1,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

spill_max = 20 # max mean spillover value
spill_step = 0.5 # step size -- must be an integer for N-M models
long_max = 100 # maximum value for C:1 plots
long_step = 2 # step size for long C:1 line values
C = 1 # past:future spillover ratio -- slope of C:1 line

vec = c(seq(0,spill_max,by=spill_step), seq(spill_max, long_max, by = long_step)[-1])
vals = data.frame(vec, C*vec)
names(vals) = c('lambda', 'lambda_f')
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

## visualize priors
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

## SI Fig. S6 -- extended x-axis
L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out, aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out, aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out, aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,100,by=25)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~vals) + 
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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

## SI Fig. S6
curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~vals) + 
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
  ggsave(filename = paste("./figures/L_curves_abvals.png"), plot = L_curves, width = 10.5, height = 5.85, units = "in")
}

#########################

######################################
##    Repeat using these priors     ##
## as the right-skewes component of ##
##     the mixture distribution     ##
##          (SI figure S7)          ##
######################################

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

## visualize prior
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

## fig. S7 -- extended x-axis
L_curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out, aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out, aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out, aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,100,by=25)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~vals) + 
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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

## fig. S7
curves <- ggplot() +
  scale_color_viridis(discrete = T, option = 'H') +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color= '#041E42', linewidth = 1.25) +
  geom_point(data = run_out[cond_short,], aes(x = past_val, y = .data[['exact']]), color= '#041E42', size = 2.5) +
  geom_line(data = run_out[cond_short,], aes(x = past_val, y = .data[['lim']]), color = '#C99700', linetype = 'dashed', linewidth = 1.25) +
  scale_x_continuous(breaks = seq(0,spill_max,by=5)) +
  scale_y_continuous(breaks = seq(0, y_upr , by = round(y_upr/2, digits = 2)), limits = c(0, y_upr)) +     
  facet_wrap(~vals) + 
  xlab( expression( "Spillover rate "*(lambda) ) ) +
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
  ggsave(filename = paste("./figures/L_curves_abvals_mixture.png"), plot = L_curves, width = 10.5, height = 5.85, units = "in")
}


####################################
### Review Response: Application ###
###    to Lassa Virus - table    ###
####################################

####################
## Specify  Model ##
##   Parameters   ##
####################

## model parameters -- not all are used depending on the values of past/future _model
## past/future _model: 0 = count model, 1 = poisson model
## prior_type = 2 uses the mixture model which effectively also captures the single beta prior
LAMBDA = 100000
t_pvec = c(75, 1000, 10)
pars = data.frame(past_model = 1, N = 1, lambda = LAMBDA, t_p = 1,
                  future_model = 1, M = 1, lambda_f = LAMBDA, t_f = 100,
                  redraw = 1, verbose = 0, std_out = 1, runSims = FALSE, reps = 50000,
                  prior_type = 2, a = .1, b = 10, H_crit = 0
)

1-((pars$t_f/t_pvec)+1)^(-0.1)

vals <- data.frame(LAMBDA, LAMBDA)

names(vals) <- c('lambda', 'lambda_f')
####################
## Define prior 1 ##
## parameters and ##
## run the model  ##
####################

mixtureDistn <- data.frame( a = c(0.1),
                            b =  c(10),
                            weights = c(1) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){betaMix_data(mixtureDistn)}

## Run model with prior 1 -- T_P1
pars$t_p = t_pvec[1]
run_out <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(run_out) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
run_out$scenario = "Scenario 1"
run_out$t_p = pars$t_p

## Run model with prior 1 -- T_P2
pars$t_p = t_pvec[2]
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 1"
tmp$t_p = pars$t_p
run_out <- rbind(run_out, tmp); rm(tmp)

## Run model with prior 1 -- T_P3
pars$t_p = t_pvec[3]
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 1"
tmp$t_p = pars$t_p
run_out <- rbind(run_out, tmp); rm(tmp)

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

## Run model with prior 2 -- T_P1
pars$t_p = t_pvec[1]
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 2"
tmp$t_p = pars$t_p
run_out <- rbind(run_out, tmp); rm(tmp)

## Run model with prior 2 -- T_P2
pars$t_p = t_pvec[2]
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 2"
tmp$t_p = pars$t_p
run_out <- rbind(run_out, tmp); rm(tmp)

## Run model with prior 2 -- T_P3
pars$t_p = t_pvec[3]
tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 2"
tmp$t_p = pars$t_p
run_out <- rbind(run_out, tmp); rm(tmp)

####################


####################
## Define prior 3 ##
## parameters and ##
## run the model  ##
####################

## note that for these particular parameter values, 
## numerical underflow issues in c++ model evaluation 
## necessitated the use of a simulation-only approach

mixtureDistn <- data.frame( a = c(0.1, 100),
                            b =  c(10, 400),
                            weights = c(0.8, 0.2) )

## generate parameter file for beta mixture for simulation use
if(pars$prior_type == 2 && pars$runSims){
  betaMix_data(mixtureDistn)
}


## Run model with prior 3 -- T_P1
pars$reps = 500000
pars$t_p = t_pvec[1]

out <- HJ_simulation(n_reps = pars$reps, parameters = pars, batch_name = 'test', output_mode = 0)
out$total_prob[1] # store results -- the [1] prevents errors when long output is requested 

tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 3"
tmp$t_p = pars$t_p
run_out <- rbind(run_out, tmp); rm(tmp)

## Run model with prior 3 -- T_P2
pars$t_p = t_pvec[2]

out <- HJ_simulation(n_reps = pars$reps, parameters = pars, batch_name = 'test', output_mode = 0)
out$total_prob[1] # store results -- the [1] prevents errors when long output is requested 

tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 3"
tmp$t_p = pars$t_p
run_out <- rbind(run_out, tmp); rm(tmp)

## Run model with prior 3 -- T_P3
pars$t_p = t_pvec[3]

out <- HJ_simulation(n_reps = pars$reps, parameters = pars, batch_name = 'test', output_mode = 0)
out$total_prob[1] # store results -- the [1] prevents errors when long output is requested 

tmp <- fullModel_eval(pars = pars, mixtureDistn = mixtureDistn, vals = vals)
names(tmp) = c('past_val', 'future_val', 'prob', 'sol', 'exact')
tmp$scenario = "Scenario 3"
tmp$t_p = pars$t_p
run_out <- rbind(run_out, tmp); rm(tmp)

####################
