# Load Library
library(here)
library(tidyverse)
library(kableExtra)
library(deSolve)
library(reshape2)

#Load data
rff <- readRDS(here("data/rff_2024.rds"))

?here

# Define parameters -- THESE ARE THE MODEL PARAMETERS
parms <- c(alpha= 2.123,        # alpha = daily contacts
           beta=0.581,        # beta = probability of infection on contact
           sigma=0.609,       # sigma = rate of recovery per day
           mu = 0.0,        # mu =  per capita birth and death rate
           omega = 0.0)     # omega = rate of immune loss per day

# Initial conditions --  THESE ARE THE CONDITIONS AT THE START OF THE SIMULATION
init <- c(S=90,           # number initially susceptible
          I=1,            # number initially infectious
          R=0)            # initially immune or "recovered"

# Define model equations -- do not change -- or change with care!
# These are the model equations.  They are written as a function called sir_ode.  
# "parms" and "init" input the parameters and inital conditions into the equations 
sir_ode <- function(times,init,parms){
    with(as.list(c(parms,init)), {
        # ODEs
        dS <- mu*(S+I+R) + omega*R -alpha*beta*I*S/(S+I+R) - mu*S 
        dI <- alpha*beta*I*S/(S+I+R)-sigma*I - mu*I
        dR <- sigma*I  - omega*R - mu*R
        list(c(dS,dI,dR))
    })
}

# This creates the output from model equations.  
#If you want to run the model for longer, change the second term eg: seq(0,200,...)
times <- seq(0,100,length.out=100)
sir_out <- lsoda(init,times,sir_ode,parms)
sir_out_long <- melt(as.data.frame(sir_out),"time")

#Plotting the model output
SIR <-ggplot(sir_out_long,aes(x=time,y=value,colour=variable,group=variable))+
    geom_line(lwd=2)+             #Add line
    xlab("Time")+ylab("Number")   #Add labels

print(SIR)

parms2 <- c(alpha= 2.123,        # alpha = daily contacts
           beta=0.581,        # beta = probability of infection on contact
           sigma=0.609,       # sigma = rate of recovery per day
           mu = 0.0,        # mu =  per capita birth and death rate
           omega = 1/60)     # omega = rate of immune loss per day


sirs_ode <- function(times,init,parms2){
    with(as.list(c(parms2,init)), {
        # ODEs
        dS <- mu*(S+I+R) + omega*R -alpha*beta*I*S/(S+I+R) - mu*S 
        dI <- alpha*beta*I*S/(S+I+R)-sigma*I - mu*I
        dR <- sigma*I  - omega*R - mu*R
        list(c(dS,dI,dR))
    })
}

# This creates the output from model equations.  
#If you want to run the model for longer, change the second term eg: seq(0,200,...)
times <- seq(0,100,length.out=100)
sirs_out <- lsoda(init,times,sirs_ode,parms2)
sirs_out_long <- melt(as.data.frame(sirs_out),"time")

SIRS <-ggplot(sirs_out_long,aes(x=time,y=value,colour=variable,group=variable))+
    geom_line(lwd=2)+             #Add line
    xlab("Time")+ylab("Number") #Add labels

print(SIRS)

parms3 <- c(alpha= 2.123,        # alpha = daily contacts
            beta=0.581,        # beta = probability of infection on contact
            sigma=0.609,       # sigma = rate of recovery per day
            mu = 0.0,        # mu =  per capita birth and death rate
            omega = 1/15)     # omega = rate of immune loss per day


sirs_ode <- function(times,init,parms3){
    with(as.list(c(parms3,init)), {
        # ODEs
        dS <- mu*(S+I+R) + omega*R -alpha*beta*I*S/(S+I+R) - mu*S 
        dI <- alpha*beta*I*S/(S+I+R)-sigma*I - mu*I
        dR <- sigma*I  - omega*R - mu*R
        list(c(dS,dI,dR))
    })
}

# This creates the output from model equations.  
#If you want to run the model for longer, change the second term eg: seq(0,200,...)
times <- seq(0,100,length.out=100)
sirs2_out <- lsoda(init,times,sirs_ode,parms3)
sirs2_out_long <- melt(as.data.frame(sirs2_out),"time")

SIRS2 <-ggplot(sirs2_out_long,aes(x=time,y=value,colour=variable,group=variable))+
    geom_line(lwd=2)+             #Add line
    xlab("Time")+ylab("Number") #Add labels

print(SIRS2)

# Define model parameters
parms <- c(beta = 1.444,   # beta = daily effective contacts (alpha*beta)
           sigma = 0.973,  # sigma = rate of exposed to infectious per day
           gamma = 0.918,  # gamma = rate of recovery from infectious to recovered per day
           mu = 0.0,       # mu = per capita birth and death rate
           omega = 0.0)     # omega = rate of immune loss per day

# Initial conditions
init <- c(S = 90,           # number initially susceptible
          E = 2,            # number initially exposed
          I = 1,            # number initially infectious
          R = 0)            # initially immune or "recovered"

# Define SEIR model equations
seir_ode <- function(times, init, parms) {
    with(as.list(c(parms, init)), {
        dS <- mu * (S + E + I + R) + omega * R - beta * I * S / (S + E + I + R) - mu * S
        dE <- beta * I * S / (S + E + I + R) - sigma * E - mu * E
        dI <- sigma * E - gamma * I - mu * I
        dR <- gamma * I - omega * R - mu * R
        list(c(dS, dE, dI, dR))
    })
}

# Model time and solution
times <- seq(0, 100, length.out = 100)
seir_out <- lsoda(init, times, seir_ode, parms)
seir_out_long <- melt(as.data.frame(seir_out), "time")

# Plot the SEIR model output
SEIR <- ggplot(seir_out_long, aes(x = time, y = value, colour = variable, group = variable)) +
    geom_line(lwd = 2) +
    xlab("Time") + ylab("Number")

# Save the SEIR plot as a PNG image
ggsave(filename = "SEIR.png", plot = SEIR, width = 8, height = 6, dpi = 300)
