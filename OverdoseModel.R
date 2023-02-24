---
      title: "Thesis Model"
author: "Liza Eubanks"

output:
      pdf_document: default
urlcolor: blue
---
      
      
library(tidyverse)
library(ggrepel)
library(deSolve)


theme_set(theme_minimal())
options(pillar.min_character_chars = 15)


##Let's write down our assumptions

## Everyone in the simulation is either Susceptible (S), Using (U_1), Active Recovory (A), Passive Recovery (P), or Removed (R), so we will have 5 dependent variables now, and will need an initial condition for each.
## Each group of people is homogeneous, we are not going to be concerned with age or co-morbidities or anything else
## The population size is constant, $N$
## We arent going to consider birth or death out of the susceptible population
## Transmission happens at rate $\lambda_1$
      
##Let's set the parameters

mu = .00728
lambda_1 = .36
alpha_1 = .006467
alpha_2 = .01681
beta_1 = .00455
lambda_2 = .65
beta_2 = .00151667

##Let's write down our IVP:

##\begin{align}
##\displaystyle\frac{dS}{dt} &= \Lambda - \muS \lambda_1SU_1 \\
##\displaystyle\frac{dU_1}{dt} &=   \lambda_1SU_1 - \muU_1 - \alpha_1U_1 - \alpha_2U_1 + 
                 ## \lambda_2U_1P + \gammaR \\
##\displaystyle\frac{dA}{dt} &= \alpha_1U_1 - \muA - \beta_1A \\
##\displaystyle\frac{dP}{dt} &= \alpha_2U_1 - \lambda_2U_1P - \beta_2P - \muP \\
##\displaystyle\frac{dR}{dt} &= \beta_1A + \beta_2P - \muR - \gammaR
##\end{align}

##where $N = S+U_1+A+P+R$. We will use initial conditions of $S=999,999$, $U_1=1$, $A=0$, $P=0$, $R=0$.

##So first, let's set up our parameters vector using the `c` function and naming the parameters into the vector `parameters` again.
  
parameters <- c(mu = mu, lambda_1 = lambda_1, 
                alpha_1 = alpha_1, alpha_2 = alpha_2, 
                beta_1 = beta_1, lambda_2 = lambda_2, beta_2 = beta_2)

##Now let's set our initial values into the vector `state.`

state <- c(S = .9, 
           U_1 = .1,
           A = 0,
           P = 0,
           R = 0)

##Let's use time from 0 to 250 with a timestep of $\frac{1}{2^4}$ for our output.

times <- seq(0, 250, by = 1/2^4)

##Now we need to setup our function `SIRModel` that contains our differential equations. Put each right hand side equation of the differential equations system on a new line. Make sure to return all three differential equations, $dS$, $dI$, and $dR$.
  
SIRModel <- function(t, state, parameters){
      with(as.list(c(state, parameters)), {
            ## put the DEs here
            dS <- -mu*S - lambda_1*S*U_1
            dU_1 <- lambda_1*S*U_1 - mu*U_1 - alpha_1*U_1 - alpha_2*U_1 + 
                  lambda_2*U_1*P 
            dA <- alpha_1*U_1 - mu*A - beta_1*A
            dP <- alpha_2*U_1 - lambda_2*U_1*P - beta_2*P - mu*P
            dR <- beta_1*A + beta_2*P - mu*R 
            
            ## return the differential equation. This gets called by ode for step of the simulation 
            list(c(dS, dU_1, dA, dP, dR))
      }) ## end the with statement
}

##Plot all 5 trajectories $S(t)$, $U_1(t)$, $A(t)$, $P(t)$, and $R(t)$ for times 0 to 250 days on the same plot. Color them all differently.

SIROut <- ode(y=state, times = times, func = SIRModel, parms = parameters)
q<-as.data.frame(SIROut) %>% ggplot() + geom_line(aes(x=time, y=U_1), color="red", size=1) + geom_line(aes(x=time, y=S), color="blue", size=1) + geom_line(aes(x=time, y=A), color="green", size=1)+ geom_line(aes(x=time, y=P), color="orange", size=1)+ geom_line(aes(x=time, y=R), color="purple", size=1)
print(q)
