
# Quasi-Stationary Distribution of an SIS Model

Simulating and approximating the quasi-stationary distribution of a stochastic $SIS$ model.

<br>

## Project Structure

-   ***0_functions.R*** : Collected SIS functions. Also includes functions used to asses trachoma prevalence distributions with a Cramér-von Mises statistic and maximum-likelihood estimations for geometric, negative binomial, Poisson, beta-geometric, exponential, and normal distributions.

-   ***01_deterministic.R*** : Coding a deterministic SIS model.

-   ***02_ctmc.R*** : Simulating from a CTMC SIS model.

-   ***03_approximation.R*** : Approximating the quasi-stationary distribution.

-   ***04_analysis.R*** : Plotting simulated sample paths. Comparing a cross-sectional distribution of infections from simulations to the approximated probability density functions at quasi-equilibrium.

<br>

## Deterministic SIS Model

Here I use a Susceptible-infectious-susceptible ($SIS$) process to model the transmission dynamics of an infectious disease in a closed population, assuming that infected individuals return to a fully susceptible state upon recovery.

Suppose a population of $N$ individuals, with an initial number, $I_0$, of infected individuals, and $S_0$ of initially susceptible individuals. Infected individuals, $I$, come in contact with, and successfully infect susceptible individuals, $S$, at rate $\beta$, and recover at rate $\gamma$.

The dynamics of a deterministic $SIS$ model is given by a system of differential equations:

$$
\begin{align*}
\frac{dS}{dt} &= -\frac{\beta}{N}SI + {\gamma}I \\
\frac{dI}{dt} &= \frac{\beta}{N}SI - {\gamma}I
\end{align*}
$$

Where $\beta$ is the effective contact rate, $\gamma$ is the recovery rate, and $N$ is the population size. Because the model assumes a closed population, $N = S(t) + I(t)$.

Together, $\beta$ and $\gamma$ define the basic reproduction number, $R_0$, where $R_0 = \frac{\beta}{\gamma}$. $R_0$ is a critical parameter reflecting the number of secondary cases expected from a single index case in a completely susceptible population. There is a high probability that an epidemic will occur when $R_0$ is greater than $1$. No epidemic is possible when $R_0$ is less than $1$.


<br>

## Stochastic SIS Model

Now, suppose we introduce probabilities that a new infection or recovery will occur. A stochastic $SIS$ model can be represented as a continuous time markov chain (CTMC), in which the state of the system at any time, $t$, depends only on the state at time $t-1$.

Take the number of susceptible, $S(t)$, and infected individuals, $I(t)$, at time $t$, with $I$ in $\{0, 1, 2, ... N\}$. The infinitesimal transition probabilities for $I(t+1)$ in a CTMC can be defined as,

$$
p_{ji}({\Delta}t) = 
\begin{cases} {b(i){\Delta}t + o{\Delta}t}, & \mbox{j = i + 1} \\ 
{d(i){\Delta}t + o{\Delta}t}, & \mbox{j = i - 1} \\
{1 - [b(i) + d(i)]{\Delta}t + o{\Delta}t}, & \mbox{j = i} \\
o{\Delta}t, & \mbox{otherwise} \end{cases}
$$

Where $b(i) = \frac{{\beta}i(N-i)}{N}{\Delta}t$ denotes a new infection, and $d(i) = {\gamma}(i){\Delta}t$ a recovery.


Here, zero is an absorbing state, meaning that all simulations will eventually converge to a zero-infection state regardless of $R_0$. However, the model can reach a stable distribution when conditioned on non-extinction.

<br>

## Quasi-Stationary Probability Distribution

True equilibrium of a stochastic $SIS$ model is $0$ for all values of $R_0$, but the system can come to a separate, quasi-equilibrium if prevented from reaching the disease-free state. This quasi-stationary distribution approaches a geometric distribution when $R_0 <1$, and a normal distribution when $R_0 >1$. If a disease is disappearing ($R_0 <1$), the distribution of infectious cases should be approximately geometric.

I use two methods from Näsell 1996 to approximate the quasi-stationary distribution for different $R_0$. Both methods use an iterative scheme to approximate the Kolmogorov forward equations and solve for quasi-stationary distribution indirectly. We can then use these approximations to plot the probability density at quasi-equilibrium for different values of $R$.

The first, $p^1$, and second, $p^2$, approximations are given by,

$$
\begin{align*}
p_i^1 &= p_1^1 \frac{ (N-1)! }{ i(N-i)! } \left( \frac{R_0}{N} \right)^{i-1}{ ,} \qquad{i = 2, 3,... N} \\
p_i^2 &= p_1^2 \frac{ (N-1)! }{ (N-i)! } \left( \frac{R_0}{N} \right)^{i-1}{ ,} \qquad{i = 2, 3,... N}
\end{align*}
$$

where $p_1^1$,

$$
p_1^1 = \left[\sum^N_{k=1} \frac{ (N-1)! }{ k(N-k)! } \left( \frac{R_0}{N} \right)^{k-1} \right]^{-1}
$$

and $p_1^2$

$$
p_1^2 = \left[\sum^N_{k=1} \frac{ (N-1)! }{ (N-k)! } \left(\frac{R_0}{N} \right)^{k-1} \right]^{-1}
$$

<br> <br>


# Trachoma Surveillance

The trachoma modeling team at Proctor has been using TF prevalence surveys to test this theory with data. We expect that, in areas where trachoma is disappearing, the distribution of cases should be approximately geometric (discrete data) or exponential (continuous data).

There are 3 survey types for trachoma surveillance: baseline, impact (TIS), and trachoma surveillance surveys (or TSS), which are administered two years after a TIS measurement below 5% TF. TSS provide a unique opportunity to explore these distributions because they detail trachoma prevalence in areas that have already met the criteria for control.

We found strong evidence that TF distributions are monotonically decreasing in a majority of the districts surveyed. This result is consistent with past findings on leprosy and onchocerciasis, as well as trachoma done at a smaller scale. We were able to fit a geometric distribution to a much larger database, including TSS from around 10,000 villages across about 300 districts in 6 countries. 


<br>

## Implications for Control

A disease is likely disappearing in areas where prevalence fits a geometric distribution. If prevalence is normally distributed, it is unlikely that the disease is disappearing, and suggests areas of super-critical transmission $(R>1)$.

Using distributions to assess progress towards elimination is encouraging for three reasons:

1.  **A geometric distribution is defined by a single parameter,** $p$. Here, $p$ represents the probability of screening positive for TF. $p$ is a sufficient statistic, meaning that it contains all the information in the data for the geometric model. If we have an estimate for $p$, we are able to predict the entire distribution of cases without having to screen every village.

2.  **Outliers don't necessarily indicate program failure or hotspots for transmission.** A geometric distribution is relatively heavy-tailed, which means higher prevalence villages are expected by random chance, even if a disease is in decline.

3. **A geometric distribution of cases may imply control.** An infectious disease may be considered controlled when sub-critical transmission has been achieved $(R<1)$. An observed geometric distribution provides evidence that $R$ is less than $1$, given that the theoretical distribution under sub-critical transmission is geometric.


<br> <br> 
