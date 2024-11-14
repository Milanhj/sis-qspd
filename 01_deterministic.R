
# SIS deterministic model


# Load Packages ----------------------------------------------------------------

library(tidyverse)
library(deSolve)


# DE Function ------------------------------------------------------------------

# (S(t), I(t))

sis_det <- function(t, x, params){
  with(as.list(c(params, x)), {
    dS <- -(beta/N)*S*I + (gamma)*I
    dI <- (beta/N)*S*I - (gamma)*I
    res <- c(dS, dI)
    list(res)
  })
} # end function


# Solve R0 >1 ------------------------------------------------------------------


## Params and Conditions -------------------------------


# Parameters
N <- 100                        # population size
beta <- 1                       # contact rate
gamma <- 0.5                    # recovery rate

# Initial conditions
S0 <- 90                        # initial susceptible 
I0 <- 10                        # initial infected


# Check R0
round(beta/(gamma), digits = 2)


# Define parameters for function call
params <- c(N=100, beta=beta, gamma=gamma)

# Define Initial conditions
xstart <- c(S=S0, I=I0)


# Define Times
times <- seq(0, 100)



##  Solve --------------------------------------------

sis_out <- as.data.frame(lsoda(xstart, times, sis_det, params))



## Plot ----------------------------------------------

# I over time
with(sis_out, plot(time, I, type = "l", ylim = c(0,100), main = "R0 > 1"))     
# S over time
with(sis_out, lines(time, S, col = "red"))    


# Legend
legend("topright",
       legend = c("S", "I"),
       lty = 1,
       col =  c("red", "black"),
       #"darkgreen"),
       #colors,
       border = FALSE,
       xjust = 1,
       yjust = 1,
       lwd = 1.5,
       # x = 0.4,
       x.intersp = 0.8,
       y.intersp = 0.97,
       # pch = c(15, 15, 15),
       bty = "n",
       #pt.cex = 1.4,
       cex = 0.85, #1.12, #1.15,
       text.col = "black",
       horiz = F ,
       inset = c(0.04, 0.04)
)




# Solve R0 < 1 -----------------------------------------------------------------


## Params and Conditions -------------------------------


# Parameters
N <- 100                        # population size
beta <- 0.95                     # contact rate
gamma <- 1                      # recovery rate

# Initial conditions
S0 <- 80                        # initial susceptible 
I0 <- 20                        # initial infected


# Check R0
round(beta/(gamma), digits = 2)


# Define parameters for function call
params <- c(N=100, beta=beta, gamma=gamma)

# Define Initial conditions
xstart <- c(S=S0, I=I0)


# Define Times
times <- seq(0, 100)



##  Solve ----------------------------

sis_out <- as.data.frame(lsoda(xstart, times, sis_det, params))



## Plot ------------------------------

# I over time
with(sis_out, plot(time, I, type = "l", ylim = c(0,100), main = "R0 < 1"))     
# S over time
with(sis_out, lines(time, S, col = "red"))    
# Legend
# legend("topright",
#        horiz = F,
#        legend = c("S", "I"),
#        lty = 1,
#        col =  c("red", "black"),
#        #"darkgreen"),
#        #colors,
#        border = FALSE,
#        xjust = 1,
#        yjust = 1,
#        lwd = 1.5,
#        # x = 0.4,
#        x.intersp = 0.8,
#        y.intersp = 0.97,
#        # pch = c(15, 15, 15),
#        bty = "n",
#        #pt.cex = 1.4,
#        #cex = 0.5, #1.12, #1.15,
#        text.col = "black",
#        inset = c(0.03, 0.2)
# )


par(mfrow = c(1,2))
