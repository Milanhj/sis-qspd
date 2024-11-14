

# SIS Continuous Time Markov Chain  


# Load Packages ----------------------------------------------------------------

library(tidyverse)
library(doMC)


# Set-up parallel processing to run in background
registerDoMC(cores = 12)

# Set seed
set.seed(3012)



# Function ---------------------------------------------------------------------


# `sis_ctmc()`: generates nsims simulations from a stochastic SIS model
  #> Option to scale initial conditions by N 
  #> Arguments:
    #> beta = effective contact rate
    #> gamma = recovery rate
    #> N = population size
    #> I = fixed number for initially infected OR..
    #> prop_i = proportion of N to be initially infected
    #> nsims = number of state paths to simulate
  #> Returns a list of simulates state paths



sis_ctmc <- function(N, tmax, I, beta, gamma, nsims, prop_i = NULL){
  
  # Initial Infected as a proportion of N
  if (is.null(prop_i)) { 
    I0 <- I
    } else { 
      # Initial infected fixed for all N
      I0 <- floor(N*prop_i)
      }
  
  # Initialize a list to store state paths
  sis_paths <- vector("list", length = nsims)
  
  for (k in 1:nsims) {
    # Define variables
    t <- 0
    i <- I0
    s <- N-I0
    
    # Counter
    j <- 1
    while (i[j] > 0 & t[j] < tmax){
      
      # Generate two variables from uniform distribution
      U <- runif(2)
      
      # Rate of an event occurring: b(i)t + d(i)t
      rate <- (beta/N)*i[j]*s[j] + (gamma)*i[j]
      # Probability of infection
      probi <- (beta*s[j]/N) / (beta*s[j]/N + gamma) # i cancels out
      
      # Exponentially distributed wait time -log()
      # Add time step dependent on the rate of an event occurring
      t[j+1] <- t[j] - log(U[1])/rate
      
      # Rate competition
      if (U[2] <= probi){
        # Add new infection
        i[j+1] <- i[j]+1
        s[j+1] <- s[j]-1
        
      } else{
        # Remove infection
        i[j+1] <- i[j]-1
        s[j+1] <- s[j]+1
        
      } # end else
      
      # add to counter
      j <- j+1
    
    } # end while
    
    # Store state path
    sis_paths[[k]] <- i
    
  } # end for
  
  # Return a list with 1:nsims simulations
  return(sis_paths) 
  
} # end function



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
    } # end for
    # Remove NA for correct length of live infections
    Ix <- Ix[!is.na(Ix)]
    # Stop loop when the number of live infections drops below threshold
    if (length(Ix) <= threshold) break
    # Add 1 to counter
    x <- x+1
  } # end while
  
  return(Ix)
  
} # end function



# Simulations ------------------------------------------------------------------


# Parameters
params <- tibble(
  beta = c(0.7, 0.8, 0.9, 1, 0.75, 1, 1.5),
  gamma = c(1, 1, 1, 1, 0.5, 0.5, 0.5)
  ) %>% 
  mutate(R0 = round(beta/(gamma), digits = 2)) %>% 
  relocate(R0)


# Population size
N <- 100

# Number of simulations and time
nsims <- 20
tmax <- 150




## Fixed I0 --------------------------------------


# Initial conditions
I0 <- 2

# List for holding simulations for all values of R0
sims_fixed_i <- list()
for (i in 1:nrow(params)){

  # parameter values
  beta <- params[i,"beta"]
  gamma <- params[i,"gamma"]

  # function call
  sims_fixed_i[[i]] <- sis_ctmc(N = N, tmax = tmax, I = I0,
                                beta = beta, gamma = gamma,
                                nsims = nsims
                                )

} # end loop

# Name each index in main list by the value of R0
names(sims_fixed_i) <- c(str_c("R0 = ", as.character(params$R0)))



## Proportional I0 ----------------------------------


# Initial conditions
# prop_i <- 0.2
# 
# 
# # List for holding simulations for all values of R0
# sims_prop_i <- list()
# for (i in 1:nrow(params)){
#   
#   # parameter values
#   beta <- params[i,"beta"]
#   gamma <- params[i,"gamma"]
#   
#   # function call
#   sims_prop_i[[i]] <- sis_ctmc(N = N, tmax = tmax, I = NULL,
#                                 beta = beta, gamma = gamma, 
#                                 nsims = nsims, prop_i = prop_i
#   )
#   
# } # end loop
# 
# # Name each index in main list by the value of R0
# names(sims_prop_i) <- c(str_c("R0 = ", as.character(params$R0)))



# Cross-section prevalence -----------------------------------------------------


## Simulate -----------------------------------------


# List for holding simulations for all values of R0
# sims_fixed_i_5000 <- list()
# for (i in 1:nrow(params)){
#   
#   # parameter values
#   beta <- params[i,"beta"]
#   gamma <- params[i,"gamma"]
#   
#   # function call
#   sims_fixed_i_5000[[i]] <- sis_ctmc(N = N, tmax = tmax, I = I0,
#                                beta = beta, gamma = gamma, 
#                                nsims = 5000, prop_i = NULL
#   )
#   
# } # end loop
# 
# # Name each index in main list by the value of R0
# names(sims_fixed_i_5000) <- c(str_c("R0 = ", as.character(params$R0)))


# simulations where it goes away
# equilibrium of all at 0


## Grab values ------------------------------------


# crsct_5000 <- list()
# for (i in 1:length(sims_fixed_i_5000)){
#   # Cross-section of simulations where there are still 3000 live infections
#   crsct_5000[[i]] <- crsect_distrib(sims_fixed_i_5000[[i]], threshold = 3000)
# 
# }
# 
# names(crsct_5000) <- names(sims_fixed_i_5000)




# Save -------------------------------------------------------------------------


# Save simulations to analyze in a separate script (04_analysis.R)

# 10 simulations
# save(sims_fixed_i, # sims_prop_i,
#      file = "outputs/simulations/sim10_n100_tmax150.rda"
#      )


# 20 simulations
save(sims_fixed_i, # sims_prop_i,
     file = "outputs/simulations/sim20_n100_tmax150.rda"
)


# 50 simulations
# save(sims_fixed_i, sims_prop_i,
#      file = "outputs/simulations/sim50_n100_tmax200.rda")

# Save cross-section I (5000 sims)
# save(crsct_5000, file = "outputs/simulations/crsct_5000_n100_tmax200.rda")



