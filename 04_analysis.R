
# Look at simulation and approximation results

library(tidyverse)

set.seed(301)


# Load Objects -----------------------------------------------------------------


# Approximations
load("outputs/approximations/approx_p1_p2.rda")

# Simulations 
load("outputs/simulations/sim20_n100_tmax150.rda")         # n = 20
# load("outputs/simulations/sim50_n100_tmax200.rda")       # n = 50

# Cross-sections
load("outputs/simulations/crsct_5000_n100_tmax200.rda")    # cross-sections



# Functions --------------------------------------------------------------------


## Simulation -------------------------------------


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
  # Initialize plot-- first index in solid black
  plot(y, main = title, ylim = c(0,N), xlim = c(0, lmax),
       ylab = "I", xlab = "time", 
       type = "l", lwd = 0.7)
  # Loop through all simulations
  for (i in 2:length(paths)) {
    yi <- paths[[i]]
    lines(yi, col = alpha("black", 0.2), lwd = 0.7)
  }
  
} # end function





## Approximation -----------------------------------------------


### Maximum-likelihood ---------------------


# Fit Expontential 
neg_ll_exp <- function(lambda, x) {
  if (lambda <= 0) return(Inf) # validate the params
  ll <- sum(log(dexp(x, rate = lambda)))
  return(-ll)
}

# Fit Normal 
neg_ll_normal <- function(theta, x) {
  mu <- theta[1] # mean
  sd <- theta[2] # standard deviation
  # Validate
  if (sd < 0 | mu >= max(x) | sd >= max(x)) return(Inf) 
  ll <- sum(log(dnorm(x, mean = mu, sd = sd)))
  return(-ll)
}


### Plotting -------------------------------


# Plot the approximations of the PDF at quasi-equilibrium sequentially
plot_qspd <- function(approx, N, plot_dist = FALSE){
  
  if (plot_dist == TRUE){
    
    for(i in length(approx):1){
      
      print(i)
      # value of R0
      r <- approx[[i]][[1]]
      # objects for plotting
      dat <- approx[[i]][[2]]
      title <- names(approx)[i]
      xlab <- "Infected"
      states <- seq(1,N,1)
      
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



# Plot Simulations -------------------------------------------------------------


names(sims_fixed_i)

for (i in 1:length(sims_fixed_i)){
  sim <- sims_fixed_i[[i]]
  title <- str_c(names(sims_fixed_i)[i], " Fixed I0")
  plot_sims(paths = sim, title, N = 100, tmax = 300)
}


  

# names(sims_prop_i)
# 
# for (i in 1:length(sims_prop_i)){
#   sim <- sims_prop_i[[i]]
#   title <- str_c(names(sims_prop_i)[i], " I0 as 0.2 of N")
#   plot_sims(paths = sim, title, N = 100, tmax = 300)
# }



# Plot Approximations ----------------------------------------------------------


set.seed(301) # again for reproducibility 


# Approximation 1
plot_qspd(approx = out_p1, N = 100, plot_dist = FALSE)
plot_qspd(approx = out_p1, N = 100, plot_dist = TRUE)



# Approximation 2
plot_qspd(approx = out_p2, N = 100, plot_dist = FALSE)
plot_qspd(approx = out_p2, N = 100, plot_dist = TRUE)







# Cross-sections ---------------------------------------------------------------


# Plot histogram of cross-section from SIS simulations
for (i in 1:length(crsct_5000)){
  
  # Cross-section of simulations where there are still 3000 live infections
  Ix <- crsct_5000[[i]]
  
  # Plot the distribution of infections in cross-section
  title <- names(crsct_5000)[i]
  hist(Ix, breaks = 100, xlim = c(0,100), freq = FALSE,
       main = title)
  
}



