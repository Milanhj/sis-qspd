
# Quasi-stationary Probability Distribution


# Approximate kolmogorov equations to solve for quasi-stationary distribution indirectly
  #> q∗ cannot be solved directly from the system of equations 
  #> use iterative scheme to solve indirectly 
  #> two different methods for approximation p1 and p2





# Load Packages ----------------------------------------------------------------

library(tidyverse)


# Functions --------------------------------------------------------------------


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






# Parameter table --------------------------------------------------------------


# Parameters
params <- tibble(
  beta = c(0.5, 0.7, 0.8, 0.9, 0.95, 1, 0.75, 1, 1.5, 2),
  gamma = c(1, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5)
) %>% 
  mutate(R0 = round(beta/(gamma), digits = 2)) %>% 
  relocate(R0)


# Population size
N <- 100
# Vector of states
states <- seq(1, N, 1)



# Approximate ------------------------------------------------------------------


## first method ----------------------------


# List for holding simulations for all values of R0
out_p1 <- list()
for (i in 1:nrow(params)){
  
  # parameter values
  beta <- params[[i,"beta"]]
  gamma <- params[[i,"gamma"]]
  
  # function call
  out_p1[[i]] <- approx_p1(N = N, beta, gamma)
  
} # end loop

# Name each index in main list by the value of R0
names(out_p1) <- c(str_c("R0 = ", as.character(params$R0)))


# Check that p1 sums to 1
all(
  for (i in 1:length(out_p1)){
    # round to deal with precision issues
    round(sum(out_p1[[i]][[2]]), 2) == 1
  }
)




## second method --------------------------------------

# List for holding simulations for all values of R0
out_p2 <- list()
for (i in 1:nrow(params)){
  
  # parameter values
  beta <- params[[i,"beta"]]
  gamma <- params[[i,"gamma"]]
  
  # function call
  out_p2[[i]] <- approx_p2(N = N, beta, gamma)
  
} # end loop

# Name each index in main list by the value of R0
names(out_p2) <- c(str_c("R0 = ", as.character(params$R0)))


# Check that p1 sums to 1
all(
  for (i in 1:length(out_p2)){
    # round to deal with precision issues
    round(sum(out_p2[[i]][[2]]), 2) == 1
  }
)



# Save -------------------------------------------------------------------------


save(out_p1, out_p2, params, file = "outputs/approximations/approx_p1_p2.rda")

