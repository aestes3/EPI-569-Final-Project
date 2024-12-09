# Aidan Estes
# R0 and Rt 

# Load Library
library(here)
library(tidyverse)
library(kableExtra)
library(deSolve)
library(reshape2)
library(dplyr)
library(ggplot2)
library(EpiEstim)
#install.packages("EpiNow2")
library(EpiNow2)

#Load data
rff <- readRDS(here("data/rff_2024.rds"))


##### CREATION OF R0 and SIR MODEL #####

# Calculate mean daily contacts
#mean(rff$`How many people did you expose?`, na.rm = TRUE)


# METHOD 1: CALC PROB OF INFECTION ON CONTACT
#tot_exposed <- sum(rff$`How many people did you expose?`, na.rm = TRUE) #155
#tot_infected <- sum(rff$`How many people did you infect?`, na.rm = TRUE) #90
# 90 / 155 = 0.581

##### Proper way ######
# CALC ALPHA
alpha1 <- mean(rff$`How many people did you expose?`, na.rm = TRUE)

# METHOD 2: Calc prob of infection on contact (BETA)
rff2 <- rff
# Create column for infected / exposed for each person
# How many ppl did you infect / how many did you expose
# Take mean of that column
rff2$prob_of_infect <- (rff$`How many people did you infect?`) / (rff$`How many people did you expose?`)
beta1 <- mean(rff2$prob_of_infect, na.rm = TRUE)

# Calc rate of recovery

rff2$delta <- (rff2$Date_of_Recovery - rff2$Date_of_Onset)
delta1 <- mean(rff2$delta, na.rm= TRUE)
#Avg Duration of illness 1.643 days
# 1/d = 0.609

#R0 = alpha * beta * delta = 2.075
R0 <- 2.075261


# Define parameters -- THESE ARE THE MODEL PARAMETERS
parms <- c(alpha= as.numeric(alpha1),        # alpha = daily contacts
           beta= as.numeric(beta1),        # beta = probability of infection on contact
           sigma= 1/ as.numeric(delta1),       # sigma = rate of recovery per day
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

################################
# Start Rt                     #
################################
rff3 <- rff

# Plotting the epidemic curve
epicurve <- as.data.frame(rff3 %>% group_by(Date_of_Onset) %>% summarize(`How many people did you infect?`=n()))

plot1 <- ggplot(data=epicurve, aes(x=Date_of_Onset, y=`How many people did you infect?`, labs=FALSE)) + 
    geom_bar(stat="identity", fill="gray45") +
    scale_x_date(date_breaks = "2 days", date_labels = "%m/%d") +
    scale_y_continuous(breaks = seq(0,15, by = 2)) +
    theme(axis.text.x = element_text(angle = 45)) +
    labs(x = "Symptom Onset Date", y = "Number of Cases", title = "Epi Curve for Rollins Fall Fever")
plot1

# First, convert the data to a format that can be used in the EpiEstim wallinga_teunis function 
# Data must be in the following format: 
### 1 column for symptom onset dates in ascending order, including dates on which 0 cases were reported, titled "dates"
### 1 column for case counts (incidence) titled "I"
### Note: to calculate an Rt estimate for day day 1 of the outbreak, we must start our epi curve 2 days prior the first symptom onset date

epicurve2 <- epicurve %>% arrange(Date_of_Onset) %>% rename(dates = Date_of_Onset, I = `How many people did you infect?`)

#Original code used 2 days before min of date onset; min is 2024-09-09

all.dates <- as.data.frame(seq(as.Date("2024-09-07"), by = "day", length.out = 45))
names(all.dates) <- "dates"

epicurve.epiestim <- merge(x=epicurve2, y=all.dates, by="dates", all="TRUE")
epicurve.epiestim <- epicurve.epiestim %>% mutate(I = ifelse(is.na(I), 0, I)) 


# Next, run the code below to estimate Rt, along with 95% confidence intervals for Rt estimates
# This requires that we specify the mean and standard deviation of the serial interval  
# An offset gamma distribution will be used for the serial interval (by default)

#################################################
#Create Serial interval using transmission pairs#
#################################################

transmission_pairs_2024 <- readRDS(here("data/transmission_pairs_2024.rds"))
transmission <- na.omit(transmission_pairs_2024)
#Create serial interval column
transmission$case_1_sick <- as.POSIXct.Date(rff3$Date_of_Onset[1])
transmission$diff_in_days<- difftime(transmission$Date_Time_Onset, transmission$case_1_sick, units = c("days"))

#Find SI Mean
mean_si <- as.numeric(mean(transmission$diff_in_days))
std_si <- as.numeric(sd(transmission$diff_in_days))

#Find duration of outbreak
max(rff$Date_of_Onset, na.rm=TRUE) - min(rff$Date_of_Onset, na.rm=TRUE) #33 days

estimates <- wallinga_teunis(epicurve.epiestim$I, 
                             method="parametric_si",
                             config = list(t_start = seq(3, 37), #Change to 33 days 
                                           t_end = seq(3, 37), #Needed to change to fix issue with row length
                                           mean_si = mean_si, 
                                           std_si = std_si, 
                                           n_sim = 1000))


# You can examine the serial interval distribution using the code below

plot(estimates$si_distr, xlab="Serial Interval (Days)", ylab="Proportion")


#Then, use the code below to plot the R(t) estimates and 95% CIs over the epi curve to examine trends:
    
epicurve_no_na <- na.omit(epicurve)
plot2.data <- cbind(epicurve_no_na, estimates$R$`Mean(R)`,
                    estimates$R$`Quantile.0.025(R)`, estimates$R$`Quantile.0.975(R)`)
names(plot2.data) <- c("dates", "I", "R", "lowerCI", "upperCI")

plot2 <- ggplot(data=plot2.data, aes(x=dates, y=I, labs=FALSE)) + 
    geom_bar(stat="identity", fill="gray45") +
    scale_x_date(date_breaks = "2 days", date_labels = "%m/%d") +
    scale_y_continuous(breaks = seq(0,15, by = 2)) +
    theme(axis.text.x = element_text(angle = 45)) +
    labs(x = "Symptom Onset Date", 
         y = "Number of Cases (bars) and Rt (line;95% CI)",
         title = "Epi Curve for Rollins Fall Fever with Rt") +
    geom_hline(aes(yintercept=1), colour="red", linetype="dashed", size=0.5) +
    geom_errorbar(data=plot2.data, aes(ymax=upperCI, ymin=lowerCI, width=0.6),stat="identity", size=0.8, show.legend=FALSE) +
    geom_line(data=plot2.data[!is.na(plot2.data$R),],aes(x=dates, y=R), color='blue', size=0.5) +
    geom_point(data = plot2.data, aes(x=dates, y=R), size=1.2, show.legend=FALSE) 
plot2
