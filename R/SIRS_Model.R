# Load Library
library(here)
library(tidyverse)
library(kableExtra)
library(deSolve)
library(reshape2)

#Load data
rff <- readRDS(here("/Users/ravikineu/Desktop/EPI 569/FINAL/EPI-569-Final-Project/data/rff_2024.rds"))


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
ggplot(sir_out_long,aes(x=time,y=value,colour=variable,group=variable))+
    geom_line(lwd=2)+             #Add line
    xlab("Time")+ylab("Number")   #Add labels
