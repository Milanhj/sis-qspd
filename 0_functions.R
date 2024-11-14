
# All Functions for SIS models:

library(tidyverse)
library(deSolve)


# Deterministic SIS ------------------------------------------------------------

# (S(t), I(t))

sis_det <- function(t, x, params){
  with(as.list(c(params, x)), {
    dS <- -(beta/N)*S*I + (gamma)*I
    dI <- (beta/N)*S*I - (gamma)*I
    res <- c(dS, dI)
    list(res)
  })
} # end function



# CTMC SIS ---------------------------------------------------------------------


sis_ctmc <- function(N, tmax, I, beta, gamma, nsims, prop_i = NULL){
  
  # Initial Infected as a proportion of N
  if (is.null(prop_i)) { 
    I0 <- I
  } else { 
    # Initial infected fixed for all N
    I0 <- floor(N*prop_i)
  }
  
  # initialize a list to store state paths
  sis_paths <- vector("list", length = nsims)
  
  for (k in 1:nsims) {
    # define variables
    t <- 0
    i <- I0
    s <- N-I0
    
    # counter
    j <- 1
    while (i[j] > 0 & t[j] < tmax){
      # generate two variables from uniform distribution
      U <- runif(2)
      
      # rate of an event occurring: b(i)t + d(i)t
      rate <- (beta/N)*i[j]*s[j] + (gamma)*i[j]
      # probability of infection
      probi <- (beta*s[j]/N) / (beta*s[j]/N + gamma)
      
      # exponentially distributed wait time
      # add time step dependent on the rate of an event occurring
      t[j+1] <- t[j] - log(U[1])/rate
      
      # rate competition
      if (U[2] <= probi){
        # add new infection
        i[j+1] <- i[j]+1
        s[j+1] <- s[j]-1
        
      } else{
        # remove infection
        i[j+1] <- i[j]-1
        s[j+1] <- s[j]+1
        
      } # end else
      
      # add to counter
      j <- j+1
      
    } # end while
    
    # store state path
    sis_paths[[k]] <- i
    
  } # end for
  
  return(sis_paths) 
  
} # end function



# SIS Cross-section -----------------------------------------------------

# Get a cross-sectional distribution of I from simulations
crsect_distrib <- function(paths, threshold = length(paths)/2){ # x_scale
  
  # find longest path to choose cross-section  
  lengths <- vector()
  for (k in 1:length(paths)){
    lengths[[k]] <- length(paths[[k]])
  }
  # longest path
  lmax <- max(lengths)
  
  # counter
  x <- 1
  while (x < lmax) {
    # create a vector of I values, NA represents simulations where I=0
    Ix <- rep(NA, length(paths))
    for (p in 1:length(paths)){
      # extract each path individually
      path <- paths[[p]]
      # take the value of I at the x selected above and store
      Ix[p] <- path[x]
    }
    # Remove NA
    Ix <- Ix[!is.na(Ix)]
    # Stop loop when the number of live infections drops below threshold
    if (length(Ix) <= threshold) break
    # Add 1 to counter
    x <- x+1
    
  } # end while
  
  return(Ix)
  
} # end function




# Approximations ---------------------------------------------------------------



# First approximation of quasi-stationary distribution
approx_p1 <- function(N, beta, gamma) {
  
  # State space vector
  v <- seq(1, N, 1)
  # New infection: 
  bi <- beta*v*(N-v)/N
  # Death or recovery:
  di <- (gamma)*v
  # Calculate R0
  R0 <- beta/(gamma)
  
  # Calculate initial value of p1
  vals <- vector()
  for (k in v){
    vals[k] <- (factorial(N-1) / (k*factorial(N-k))) * ((R0/N)^(k-1))
  }  
  # sum from k=1 to N
  p1_1 <- sum(vals)^-1
  
  # Store p1_1 as first value in p1
  p1 <- p1_1
  
  # Calculate and store p for all i
  for (i in 1:N) {
    p1[i+1] <- p1[i]*(bi[i]/di[i+1])
  }
  # Remove the value at N+1
  p1 <- p1[-c(N+1)]
  # Transpose p1
  p1 <- t(p1)
  
  # Store output as a list and return list
  out1 <- list(R0 = R0, p1 = p1)
  return(out1)
  
} # end function



# Second approximation of quasi-stationary distribution
approx_p2 <- function(N, beta, gamma) {
  
  # State space vector
  v <- seq(1, N, 1)
  # New infection: 
  bi <- beta*v*(N-v)/N
  # Recovery:
  di <- (gamma)*v
  # Calculate R0
  R0 <- beta/(gamma)
  
  # Calculate initial value of p2
  vals <- vector()
  for (k in v){
    vals[k] <- (factorial(N-1) / (factorial(N-k))) * ((R0/N)^(k-1))
  } 
  # sum from k = 1 to N
  p2_1 <- sum(vals)^-1
  
  # Store p2_1 as first value in p2
  p2 <- p2_1
  
  # Calculate and store p for all i
  for (i in 2:N) {
    p2[i] <- p2_1 * (factorial(N-1)/factorial(N-i)) * (R0/N)^(i-1)
  }
  
  # Transpose p2
  p2 <- t(p2)
  
  # Store output as a list and return list
  out2 <- list(R0 = R0, p2 = p2)
  return(out2)
  
} # end function


# Generator matrix  for approximation 1
Q1 <- function(N, beta, gamma){
  
  # Transition matrix
  Q1 <- matrix(0, nrow = N, ncol = N)
  # State space vector
  v <- seq(1, N, 1)
  # New infection: 
  bi <- beta*v*(N-v)/N
  # Death or recovery:
  di <- (gamma)*v
  
  # Fill matrix
  # starts from position [2,2], ends at  fill in other values manually
  for (i in 2:(N-1)) {
    Q1[i,i] <- -(bi[i] + di[i])   # diagonal entries
    Q1[i-1,i] <- di[i]            # superdiagonal entries
    Q1[i+1,i] <- bi[i]            # subdiagonal entries
  }
  
  # manually set first and last row values 
  # assume d(1) = 0 (no chance of letting system go to 0)
  Q1[1,1] <- -bi[1] # state 0 to 0
  Q1[2,1] <- bi[1] # 
  Q1[1,2] <- di[2] # 
  Q1[N-1,N] <- di[N]
  Q1[N,N] <- -di[N] # state N to N
  
  return(Q1)
}




# Generator matrix approximation 2
Q2 <- function(N, beta, gamma){
  
  # Transition matrix
  Q2 <- matrix(0, nrow = N, ncol = N)
  # State space vector
  v <- seq(1, N, 1)
  # New infection: 
  bi <- beta*v*(N-v)/N
  # Death or recovery:
  di <- (gamma)*v
  
  
  # Fill matrix
  # starts from position [2,2], ends at  fill in other values manually
  for (i in 2:(N-1)) {
    Q2[i,i] <- -(bi[i] + di[i-1])  # diagonal entries
    Q2[i-1,i] <- di[i-1]           # superdiagonal entries
    Q2[i+1,i] <- bi[i]             # subdiagonal entries
  }
  
  # manually set first and last row values 
  # assume d(1) = 0 (no chance of letting system go to 0)
  Q2[1,1] <- -bi[1] # state 0 to 0
  Q2[2,1] <- bi[1] # 
  Q2[1,2] <- di[1] # 
  Q2[N-1,N] <- di[N-1]
  Q2[N,N] <- -di[N] # state N to N
  
  return(Q2)
}




# Plotting functions -----------------------------------------------------------


## for simulation ------------------------------------

# Plot state paths from each simulation for a single value of R
plot_sims <- function(paths, title, N, tmax){ 
  
  # Find longest path for plotting xlim
  lengths <- vector()
  for (k in 1:length(paths)){
    lengths[[k]] <- length(paths[[k]])
  }
  # longest path 
  lmax <- max(lengths)
  
  # First simulation
  y <- paths[[1]]
  # Initialize plot
  plot(y, main = title_text, ylim = c(0,N), xlim = c(0, lmax),
       ylab = "Infected individuals", xlab = "time", 
       type = "l", lwd = 0.7)
  # Loop through all simulations
  for (i in 2:length(paths)) {
    yi <- paths[[i]]
    lines(yi, col = alpha("black", 0.2), lwd = 0.7)
  }
  
} # end function




## for approximation --------------------------------


# Exponential Maximum-likelihood
neg_ll_exp <- function(lambda, x) {
  if (lambda <= 0) return(Inf) # validate the params
  ll <- sum(log(dexp(x, rate = lambda)))
  return(-ll)
}

# Normal Maximum-likelihood
neg_ll_normal <- function(theta, x) {
  mu <- theta[1] # mean
  sd <- theta[2] # standard deviation
  # Validate
  if (sd < 0 | mu >= max(x) | sd >= max(x)) return(Inf) 
  ll <- sum(log(dnorm(x, mean = mu, sd = sd)))
  return(-ll)
}



# Plot the approximations of the PDF at quasi-equilibrium sequentially
plot_qspd <- function(approx, N, plot_dist = FALSE){
  
  if (plot_dist == TRUE){
    # add the lines for exp/normal density functions
    
    for(i in length(out_p1):1){
      
      print(i)
      # objects for plotting
      dat <- approx[[i]][[2]]
      title <- names(approx)[i]
      xlab <- "Infected"
      states <- seq(1,N,1)
      
      # value of R0
      r <- approx[[i]][[1]]
      
      # Plot Distributions
      if (r < 1){ # plot when R0 < 1
        
        # Generate values using pdf to fit normal
        exp_vals <- sample(states, replace = TRUE, prob = dat[1,])
        # initial param values from poisson distribution using mean&sd of generated values
        lambda <- 1/rpois(1, mean(exp_vals))
        # estimate and store exponential parameter
        theta_exp <- optim(lambda, neg_ll_exp, method = "BFGS", 
                           hessian = TRUE, x = exp_vals)$par 
        
        # plot values when R0 < 1
        plot(states, dat, type = "h", lwd = 3, lend = 2,
             main = title, xlab = xlab, ylab = "Density"
        )
        # exponential distribution
        lines(states, dexp(states, rate = theta_exp), col = "red", lwd = 2)
        
      } else{
        
        if (r > 1){ # Plot when R0 > 1
          
          # Generate values using pdf to fit normal
          norm_vals <- sample(states, replace = TRUE, prob = dat[1,])
          # initial param values from poisson distribution using mean&sd of generated values
          theta <- rpois(1, mean(norm_vals))
          theta[2] <- rpois(1, sd(norm_vals))
          
          # store mean and standard deviation
          theta_norm <- optim(theta, neg_ll_normal, method = "BFGS", 
                              hessian = TRUE, x = norm_vals)$par 
          
          # plot values when R0 > 1
          plot(states, dat, type = "h", lwd = 3, lend = 2,
               main = title, xlab = xlab, ylab = "Density"
          )
          # normal distribution
          lines(states, dnorm(states, mean = theta_norm[1],
                              sd = theta_norm[2]),
                col = "blue", lwd = 2)
          
        } else {
          
          # Generate values using pdf to fit normal
          exp_vals <- sample(states, replace = TRUE, prob = dat[1,])
          # initial param values from poisson distribution using mean+sd of generated values
          lambda <- 1/rpois(1, mean(exp_vals))
          # store rate
          theta_exp <- optim(lambda, neg_ll_exp, method = "BFGS", 
                             hessian = TRUE, x = exp_vals)$par #lambda
          
          # Fit normal
          norm_vals <- sample(states, replace = TRUE, prob = dat[1,])
          # initial parameter values from a poisson distribution
          theta <- rpois(1, mean(norm_vals))
          theta[2] <- rpois(1, sd(norm_vals))
          
          # store mean and standard deviation
          theta_norm <- optim(theta, neg_ll_normal, method = "BFGS",  
                              hessian = TRUE, x = norm_vals)$par 
          
          # plot values when R0 = 1
          plot(states, dat, type = "h", lwd = 3, lend = 2,
               main = title, xlab = xlab, ylab = "Density"
          )
          # exponential distribution
          lines(states, dexp(states, rate = theta_exp),
                col = "red", lwd = 2)
          # normal distribution
          lines(states, dnorm(states, mean = theta_norm[1], sd = theta_norm[2]),
                col = "blue", lwd = 2)
          
        } # end if else
        
      } # end else 1
      
    } # end for
    
    # Just qspd plots
  } else { 
    
    for (i in length(out_p1):1){
      
      # objects for plotting
      dat <- approx[[i]][[2]]
      title <- names(approx)[i]
      xlab <- "Infected"
      states <- seq(1,N,1)
      
      # plot distribution infected
      plot(states, dat, type = "h", lwd = 3, lend = 2,
           main = title, xlab = xlab, ylab = "Density"
      )
      
    } # end for
    
  } # end else
  
} # end function





# Scale R with N --------------------------------------------------------------


# Scale R with N
scale_r <- function(N = c(25, 50, 100, 150), beta0, gamma0, scale_param = "beta"){
  
  # Calculate R0 for initial value of N
  R0 <- beta0/gamma0
  # Vector of beta for each N
  betas <- c(beta0)
  # Vector of gamma for each N
  gammas <- c(gamma0)
  # Vector of R for each N
  Rs <- c(R0)
  # Data frame to store new epidemic parameters
  params <- data.frame(
    N = N[1],
    beta = beta0,
    gamma = gamma0,
    R = R0
  )
  
  # Scale R0 proportional to N
  for (n in 1:length(N)){
    # Calculate the percent increase/decrease between each N
    inc_n <- (N[n+1] - N[n])/N[n]
    
    # Return an error if neither beta nor gamma are provided
    if (scale_param != "beta" & scale_param != "gamma"){ 
      stop("Invalid parameter input. Cannot be scaled")
    }
    
    # Scale beta
    if (scale_param == "beta"){
      beta_n <- betas[n]
      betas[n+1] <- beta_n + (beta_n*inc_n)    # next beta reflects change in N
      gammas[n+1] <- gammas[n]                 # no change in gamma
      Rs[n+1] <- Rs[n] + (inc_n*Rs[n])         # next R reflects change in beta and N
      
      # Or scale gamma  
    } else{
      gamma_n <- gammas[n]
      gammas[n+1] <- gamma_n + (gamma_n*inc_n)  # next beta reflects change in N
      betas[n+1] <- betas[n]                    # no change in beta
      Rs[n+1] <- Rs[n] + (inc_n*Rs[n])          # next R reflects change in gamma and N
      
    } # end else  
    
    # Store new beta, gamma, and R for next position of N
    params[n+1,] <-  c(N[n+1], betas[n+1], gammas[n+1], Rs[n+1])
    
  } # end for loop
  
  # Remove extraneous final row of df (comes from storing at n+1) 
  params <- params[-(length(N)+1),]
  
  return(params)
  
} # end function



# Uses the parameter outputs from scale_r() to run simulations for each N and R combo
multi_sis <- function(params, I = NULL, prop_i = 0.25, tmax = 200, sims = 5){
  
  # Initial Infected as a proportion of N
  if (is.null(I)) {
    I0 <- N*0.25 } else {
      # Initial infected fixed for all N
      I0 <- I}
  
  # Store each N
  N <- params$N 
  
  # List with each element as the simulations from each N
  all_sims <- vector("list", length = length(N))
  
  # Simulations
  for (k in 1:length(all_sims)) {
    
    # Get the correct parameters
    N_k <- N[k]
    # N_params <- list_params[[k]]
    N_params <- params[k,]
    B <- N_params[,2]
    g <- N_params[,3]
    
    # Store state paths for single N
    sis_paths <- vector("list")
    for (sim in 1:sims) {
      
      # define variables
      t <- 0
      i <- I0
      s <- N_k-I0
      
      # counter
      j <- 1
      while (i[j] > 0 & t[j] < tmax){
        # generate two variables from uniform distribution
        U <- runif(2)
        # rate of an event occurring: b(i)t + d(i)t
        rate <- (B/N_k)*i[j]*s[j] + g*i[j]
        # probability of infection
        probi <- (B*s[j]/N_k) / (B*s[j]/N_k + g)
        # exponentially distributed wait time
        # add time step dependent on the rate of an event occurring
        t[j+1] <- t[j] - log(U[1])/rate
        # rate competition
        if (U[2] <= probi){
          # add new infection
          i[j+1] <- i[j]+1
          s[j+1] <- s[j]-1
          
        } else{
          # remove infection
          i[j+1] <- i[j]-1
          s[j+1] <- s[j]+1
        }
        # add to counter
        j <- j+1
        
      } # end while loop
      
      # store state path
      sis_paths[[sim]] <- i
      
    } # end simulation for loop
    
    all_sims[[k]] <- sis_paths 
    
  } # end N for loop 
  
  return(all_sims)
  
} # end function



# Zeros overtime ----------------------------------------------------

num_zeros <- function(state_paths, by = 1){
  
  # find longest path to choose sequence of x  
  lengths <- vector()
  for (n in 1:length(state_paths)){
    lengths[[n]] <- length(state_paths[[n]])
  }
  # maximum value for x
  lmax <- max(lengths)
  # x values to calculate at
  seq_indices <- seq(0, lmax, by = by)
  
  # Rows = paths
  # Columns = values at each x in seq_indices
  # 0 and NA reflects an extinct infection
  point_vals <- matrix(0, nrow = length(state_paths), ncol = length(seq_indices))
  
  for (p in 1:length(state_paths)) {
    # extract path p
    path <- state_paths[[p]]
    # Fill in each column for corresponding value in the sequence of cross-sections 
    # Each row contains a unique path 
    for (i in 1:length(seq_indices)) {
      point_vals[p,i] <- path[i]
    }
  }
  
  # vector to hold the number of paths that are extinct at each x position
  num_extinct <- vector()
  # count the number of zeros at each x (in each column)
  for (j in 1:ncol(point_vals)){
    v <- point_vals[,j]
    v <- v[v == 0]
    v <- v[is.na(v)]
    num_extinct[j] <- length(v)
  }
  
  out <- list(num_extinct, seq_indices)
  return(out)
  
}



# Trachoma Prevalence  -----------------------------------------------------------


## Maximum-Likelihood for Discrete ------------------

# Geometric
neg_ll_geom <- function(pp, a, b) {
  # Validate parameters
  if (pp < 0 | pp > 1) return(Inf) 
  # Calculate negative log-likelihood
  nll <- -sum(log(dgeom(a, prob = 1/(1 + b*pp))))
  return(nll)
}


# Negative Binomial
neg_ll_nbinom <- function(theta, a, b) {
  # Parameter values
  n <- theta[1]
  pp <- theta[2]
  # need finite values for parameters
  if (pp < 0 | pp > 1 | n < 0) return(Inf)
  # Total kids surveyed
  n_total <- sum(b)
  # Density calculations
  full_dens <- dnbinom(a, size = n, prob = n/(n + b*pp))
  dens_n <- dnbinom(seq(0, n_total, 1), size = n, prob = n/(n + b*pp))
  # Truncate observed density with proportion of density below n
  truc_dens <- full_dens/sum(dens_n)
  # Negative log-likelihood
  nll <- -sum(log(truc_dens))
  
  return(nll)
}



# Poisson
neg_ll_pois <- function(pp, a, b) {
  # Validate Parameters
  if (pp < 0 | pp > 1) return(Inf) 
  # Total kids surveyed
  n_total <- sum(b)
  # Density Calculations
  full_dens <- dpois(a, pp*b)
  dens_n <- dpois(seq(0, n_total, 1), pp*b)
  # Truncate observed density with proportion of density below n
  truc_dens <- full_dens/sum(dens_n)
  # Negative log-likelihood
  nll <- -sum(log(truc_dens))
  
  return(nll)
}



# Beta-Geometric
neg_ll_betageom <- function(theta, x, n) {
  # Parameters
  alpha <- theta[1]
  mu <- theta[2]
  # Calculate beta
  beta <- mu * n * (alpha - 1)
  # Total kids
  n_total <- sum(n)
  # Validate alpha and beta (must be positive)
  if (alpha <= 0 | any(beta <= 0) |
      mu <= 0 | mu >= 1) return(Inf)
  full_dens <- dbetageom(x, shape1 = alpha, shape2 = beta)
  dens_n <- dbetageom(seq(0, n_total, 1), shape1 = alpha, shape2 = beta)
  truc_dens <- full_dens/sum(dens_n)
  # Negative log-likelihood
  nll <- -sum(log(truc_dens))
  return(nll)
}



## Cramér-von Mises ----------------------------------

# All functions combined in cvm() to return p-values in a single function call

# Data with identifying variable in [,1]
# Parameter variable in [,2]
# Numerator (TF) [,3]
# Denominator (village size) [,4]

# Returns List (all districts) of lists (district) of lists of p+matrices (one sim)
sim <- function(data, n_sims, pp0 = 0.5){
  
  # Create a vector of districts to sample
  names <- as.vector(unlist(unique(data[,1])))
  
  # Name of variable to filter with
  id_var <- names(data[,1])  
  # Parameter variable
  param_var <- names(data[,2])
  # Name of numerator variable
  num_var <- names(data[,3])
  # Name of denominator variable
  denom_var <- names(data[,4])
  
  # Empty list to hold simulations for all districts
  id_sims <- list()
  
  for (i in 1:length(names)){ # loop through all districts
    # District i
    # List for holding the sims for single district
    sims <- vector("list", length = n_sims+1)
    # District Name
    d_i <- names[i]
    # New data set with values for single ID in bootstrap
    filter_data <- data %>% 
      filter(!!sym(id_var) == d_i) %>% 
      mutate(prop_pos = !!sym(num_var)/!!sym(denom_var))
    
    # Observed Parameter p 
    pp <- unique(filter_data[[2]])
    # Store observed p and observed data together
    sim_vals_1 <- list(pp, filter_data)
    # Add observed data at first position in district's simulations list
    sims[[1]] <- sim_vals_1
    # Vector of the denominators... want to keep these the same as actual data
    denoms <- filter_data[[4]]
    
    for (j in 2:length(sims)){ # 
      
      # Matrix to hold single simulation
      sim_matrix <- matrix(nrow = length(denoms), ncol = 3,
                           dimnames = list(NULL, c("num", "denom", "prop")))
      # Fill first col with simulated tf data
      # Using shape of 1 and the pp from fitting that district
      numerators <- rgeom(n = nrow(sim_matrix),
                          prob = 1/(1 + denoms*pp))
      sim_matrix[,1] <- numerators
      # Store the actual denominators in second col of matrix
      sim_matrix[,2] <- denoms
      # Store proportion positive in third col
      sim_matrix[,3] <- sim_matrix[,1]/sim_matrix[,2]
      
      # Estimate p for that simulation
      tryCatch(
        { # function call
          out_i <- optim(pp0, neg_ll_geom, 
                         method = "BFGS", hessian = TRUE, 
                         a = numerators,  b = denoms
          )
        },  # end function call
        error = function(cond){
          # Error message
          message(conditionMessage(cond))
          # Choose a return value in case of error
          NA
        } # end error
        
      ) # end tryCatch
      
      # Take p from optim output
      p_i <- as.numeric(out_i$par)
      # Store observed p and observed data together in a list
      sim_vals_i <- list(p_i, sim_matrix)
      # Store one p and sim in list for that district
      sims[[j]] <- sim_vals_i
      
    } # end for loop 2 (individual sim)
    id_sims[[i]] <- sims
  } # end for loop 1
  
  # return list of lists
  return(id_sims)
  
} # end function




# Function to take cross-sectional value at each x based on bin size
# Stores the theoretical and observed at each point together as rows in matrix
# Calculates a single value for one district


# Calculates the distance from theoretical for a single curve
crssct_cprob <- function(dat, p, step = 0.001){ # needs observed data as proportion
  
  # Vector of positions to take cross-section between 0 and 1
  x_step <- seq(0, 1, by = step)
  
  # Matrix to hold simulated values
  m_values <- matrix(
    nrow = length(x_step), ncol = 4, 
    dimnames = list(NULL, c("x", "theoretical", "observed", "sqr_dist")))
  
  # Observed: ECDF function
  ecdf_obs <- ecdf(dat)
  
  for (i in 1:length(x_step)){
    # Value at step i
    step_val <- x_step[i] # i
    step_val_scaled <- step_val*100
    
    # Theoretical cdf value at x
    y_theor <- pgeom(step_val_scaled, 1 / (1 + 100*p)) # calculate with pgeom()
    # Observed cdf value at x
    y_obs <- ecdf_obs(step_val)
    # Summed squared distance between theoretical and observed
    sqr_dist <- (y_theor - y_obs)^2
    
    # Store Values
    m_values[i,1] <- step_val
    m_values[i,2] <- y_theor
    m_values[i,3] <- y_obs
    m_values[i,4] <- sqr_dist
    
  } # end for loop (matrix)
  
  # Calculate the sum of the squared distances to store at second position of output
  sum_sqr_dist <- sum(m_values[,4])
  crssct_vals_sum_diff <- list(m_values, sum_sqr_dist)
  return(crssct_vals_sum_diff)
  
} # end function



# For a single *district*: spits out a vector of distances for each simulation 
sum_sqr_dists <- function(sims, step){
  
  district_vals <- list()
  # Observed p
  pp <- sims[[1]][[1]]
  # Observed proportions
  prop_obs <- sims[[1]][[2]]$prop_pos
  # Call function to calculate distance
  calc_obs <- crssct_cprob(prop_obs, pp , step = 0.001)
  
  # Vector to hold summed distances for each of a district's sims
  sim_dist <- vector()
  for (i in 2:length(sims)){
    # Simulated estimation for p
    p_sim <- sims[[i]][[1]]
    # Simulated proportions to calculate the ecdf
    props <- sims[[i]][[2]][,3]
    # Call function to calculate distance
    calc <- crssct_cprob(props, p_sim , step = 0.001)
    # Store the summed distance
    sim_dist[i-1] <- calc[[2]]
    
  } # end for loop
  
  # Store observed value at 1st potistion in output list
  district_vals[[1]] <- calc_obs[[2]]
  # Store vector of simmed values at 1st potistion in output list
  district_vals[[2]] <- sim_dist
  
  names(district_vals) <- c("obs_diff", "sim_diffs") 
  return(district_vals) 
  
} # end function


get_p_val <- function(sims, name){
  
  # Vector to store proportions together
  prop_higher <- list()
  # loop through each district
  for (i in 1:length(sims)){
    obs_diff <- sims[[i]][[1]]
    sim_diff <- sims[[i]][[2]]
    n <- length(sim_diff)
    # Proportion of simulated differences greater than observed value
    prop <- length(sim_diff[sim_diff > obs_diff])/n
    prop_higher[[i]] <- prop
  } # end loop
  # Add district names
  names(prop_higher) <- name 
  return(prop_higher)
  
} # end function


### final combined function that calls all together
# Simulate, then calculate sum of squares distances, then calculate and return p-values
cvm <- function(data, pp0 = 0.5, n_sims, stp = 0.001){
  
  all_sims <- sim_districts(data, n_sims, pp0)
  # all_sims <- boot_then_sim(data, b, n_sims, pp0)
  # Matrix to hold p values
  cvm_p <- vector()
  # For loop to fill list of country bootstrapping
  for (i in 1:length(all_sims)){
    sim_group <- all_sims[[i]]
    
    print("calculate cvm stat")
    # Calculate distances from theoretical CDF
    sim_diffs <- district_sum_sqr_dists(sim_group, stp)
    
    print("get p")
    # Store p values for a district
    # The proportion of simulated distances higher than bootstrap for a district
    obs_diff <- sim_diffs[[1]]
    sim_diff <- sim_diffs[[2]]
    n <- length(sim_diff)
    # Proportion of simulated differences greater than observed value
    cvm_p[i] <- length(sim_diff[sim_diff > obs_diff])/n
    
  } # end loop
  
  # Return vector of p-values
  return(cvm_p)

} # end function
