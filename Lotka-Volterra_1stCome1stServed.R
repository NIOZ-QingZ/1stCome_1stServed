# Load required libraries
library(deSolve)
library(ggplot2)
library(reshape2)  # For reshaping the data for ggplot

# Define the phytoplankton competition model function
phytoplankton_model <- function(t, state, parameters) {
  # Unpack state variables (phytoplankton species populations)
  P1 <- state[1]  # Population of phytoplankton species 1
  P2 <- state[2]  # Population of phytoplankton species 2
  
  # Unpack parameters
  r1 <- parameters["r1"]  # Growth rate of species 1
  r2 <- parameters["r2"]  # Growth rate of species 2
  alpha12 <- parameters["alpha12"]  # Effect of species 2 on species 1
  alpha21 <- parameters["alpha21"]  # Effect of species 1 on species 2
  K1 <- parameters["K1"]  # Carrying capacity for species 1 (nutrient/light limitation)
  K2 <- parameters["K2"]  # Carrying capacity for species 2 (nutrient/light limitation)
  t_intro2 <- parameters["t_intro2"]  # Time when species 2 is introduced
  
  # Apply condition for the introduction of species 2
  if (t < t_intro2) {
    P2 <- 0  # Before introduction time, species 2 is absent
  }
  
  # Define phytoplankton-specific growth equations
  dP1_dt <- r1 * P1 * (1 - (P1 + alpha12 * P2) / K1)  # Logistic growth with competition for species 1
  dP2_dt <- r2 * P2 * (1 - (P2 + alpha21 * P1) / K2)  # Logistic growth with competition for species 2
  
  return(list(c(dP1_dt, dP2_dt)))
}

####################
## dynamics of species 1 and 2
###################

# Define initial conditions and parameters
initial_state <- c(P1 = 10, P2 = 10)  # Start with species 1 present, species 2 absent

t_intro2 = 0  
parameters <- c(
    r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
    r2 = 0.8,      # Growth rate of phytoplankton species 2
    alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
    alpha21 = 1, # Competition effect of species 1 on species 2
    K1 = 100,      # Carrying capacity for species 1 (related to nutrient/light availability)
    K2 = 90,      # Carrying capacity for species 2 (related to nutrient/light availability)
    t_intro2 = t_intro2  # Introduce species 2 at time = 10
  )
  
# Define the time points for the simulation
times <- seq(0, 1000, by = 1)

# Run the simulation using ode solver
output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)

# Convert the output to a data frame
output_df <- as.data.frame(output)

plot(output_df$time, output_df$P1, col="red", ylab="species abundance",ylim=c(0,100))
lines(output_df$time, output_df$P2, col="green")
legend("right",legend=c("species 1 (K1=100)","species 2 (K2=90)"), col=c("red","green"), lty=c(1,1),text.col = c("red","green"))

################ 
### senstivity analysis of individual parameter
##############

#### 
## for introducing days of t
sen_t_intro <- NULL

t_seq = seq(0,10, by=0.5)

for (t_intro2 in t_seq) {
  
  # Define initial conditions and parameters
  initial_state <- c(P1 = 10, P2 = 10)  # Start with species 1 present, species 2 absent
  
  parameters <- c(
    r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
    r2 = 0.8,      # Growth rate of phytoplankton species 2
    alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
    alpha21 = 1, # Competition effect of species 1 on species 2
    K1 = 90,      # Carrying capacity for species 1 (related to nutrient/light availability)
    K2 = 90,      # Carrying capacity for species 2 (related to nutrient/light availability)
    t_intro2 = t_intro2  # Introduce species 2 at time = 10
  )
  
  # Define the time points for the simulation
  times <- seq(0, 1000, by = 1)
  
  # Run the simulation using ode solver
  output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
  
  # Convert the output to a data frame
  output_df <- as.data.frame(output)
  
  # Melt the data frame for ggplot
  output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
  
  P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
  P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
  
  diff_P1_to_P2 = P1_end_biomass-P2_end_biomass
  P1_perc = P1_end_biomass
  
  sen_t_intro = rbind(sen_t_intro, data.frame(t_intro2=t_intro2, 
                                              diff_P1_to_P2=diff_P1_to_P2,
                                              P1_end_biomass=P1_end_biomass,
                                              P2_end_biomass=P2_end_biomass,
                                              sum_biomass=P1_end_biomass+P2_end_biomass))
}

# plot(sen_t_intro$t_intro2, sen_t_intro$diff_P1_to_P2)

plot(sen_t_intro$t_intro2, sen_t_intro$P1_end_biomass/sen_t_intro$sum_biomass, ylab="percentage of species 1", xlab="introducing day of species 2 (x days after species 1)")

######
## for carry capacity K
######

sen_K <- NULL

K2_to_K1_seq = seq(0.98,1.05, by=0.0005)

for (K2_to_K1 in K2_to_K1_seq) {
  
  # Define initial conditions and parameters
  initial_state <- c(P1 = 10, P2 = 10)  # Start with species 1 present, species 2 absent
  
  parameters <- c(
    r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
    r2 = 0.8,      # Growth rate of phytoplankton species 2
    alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
    alpha21 = 1, # Competition effect of species 1 on species 2
    K1 = 90,      # Carrying capacity for species 1 (related to nutrient/light availability)
    K2 = K2_to_K1*K1,      # Carrying capacity for species 2 (related to nutrient/light availability)
    t_intro2 = 0  # Introduce species 2 at time x
  )
  
  # Define the time points for the simulation
  times <- seq(0, 1000, by = 1)
  
  # Run the simulation using ode solver
  output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
  
  # Convert the output to a data frame
  output_df <- as.data.frame(output)
  
  # Melt the data frame for ggplot
  output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
  
  P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
  P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
  
  diff_P1_to_P2 = P1_end_biomass-P2_end_biomass
  P1_perc = P1_end_biomass
  
  sen_K = rbind(sen_K, data.frame(K2_to_K1=K2_to_K1, 
                                              diff_P1_to_P2=diff_P1_to_P2,
                                              P1_end_biomass=P1_end_biomass,
                                              P2_end_biomass=P2_end_biomass,
                                              sum_biomass=P1_end_biomass+P2_end_biomass))
}


plot(sen_K$K2_to_K1, sen_K$P1_end_biomass/sen_K$sum_biomass, ylab="percentage of species 1", xlab="carry capacity K2/K1")
abline(h=0.5, col="red")
abline(v=1, col="red")

######
## for initial states P
######

sen_P <- NULL

P2_to_P1_seq = seq(0.01,10, by=0.15)

for (P2_to_P1 in P2_to_P1_seq) {
  
  # Define initial conditions and parameters
  P1 = 10
  P2 = P1*P2_to_P1
  initial_state <- c(P1 = P1, P2 = P2)  # Start with species 1 present, species 2 absent
  
  parameters <- c(
    r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
    r2 = 0.8,      # Growth rate of phytoplankton species 2
    alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
    alpha21 = 1, # Competition effect of species 1 on species 2
    K1 = 90,      # Carrying capacity for species 1 (related to nutrient/light availability)
    K2 = 90,      # Carrying capacity for species 2 (related to nutrient/light availability)
    t_intro2 = 0  # Introduce species 2 at time x
  )
  
  # Define the time points for the simulation
  times <- seq(0, 1000, by = 1)
  
  # Run the simulation using ode solver
  output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
  
  # Convert the output to a data frame
  output_df <- as.data.frame(output)
  
  # Melt the data frame for ggplot
  output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
  
  P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
  P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
  
  sen_P = rbind(sen_P, data.frame(P2_to_P1=P2_to_P1, 
                                  P1_end_biomass=P1_end_biomass,
                                  P2_end_biomass=P2_end_biomass,
                                  sum_biomass=P1_end_biomass+P2_end_biomass))
}


plot(sen_P$P2_to_P1, sen_P$P1_end_biomass/sen_P$sum_biomass, ylab="percentage of species 1", xlab="initial biomass P2/P1")
abline(h=0.5, col="red")
abline(v=1, col="red")

#####################
###### sensitivity analysis of all parameters:
#####################
sen_tIntro_K_P <- NULL

t_seq = seq(0,10, by=0.5)
K2_to_K1_seq = seq(0.98,1.03, by=0.001)
P2_to_P1_seq = seq(0.01,10, by=0.15)

for (t_intro2 in t_seq) {
  for (K2_to_K1 in K2_to_K1_seq) {
    for (P2_to_P1 in P2_to_P1_seq) {
    
      # Define initial conditions and parameters
      P1 = 10
      P2 = P1*P2_to_P1
      initial_state <- c(P1 = P1, P2 = P2)  # Start with species 1 present, species 2 absent
      
      parameters <- c(
      r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
      r2 = 0.8,      # Growth rate of phytoplankton species 2
      alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
      alpha21 = 1, # Competition effect of species 1 on species 2
      K1 = 90,      # Carrying capacity for species 1 (related to nutrient/light availability)
      K2 = K2_to_K1*K1,      # Carrying capacity for species 2 (related to nutrient/light availability)
      t_intro2 = t_intro2  # Introduce species 2 at time x
    )
    
    # Define the time points for the simulation
    times <- seq(0, 1000, by = 1)
    
    # Run the simulation using ode solver
    output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
    
    # Convert the output to a data frame
    output_df <- as.data.frame(output)
    
    # Melt the data frame for ggplot
    output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
    
    P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
    P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
    
    sen_tIntro_K_P = rbind(sen_tIntro_K_P, 
                         data.frame(t_intro2=t_intro2,
                                    K2_to_K1=K2_to_K1,
                                    P2_to_P1=P2_to_P1,
                                    P1_end_biomass=P1_end_biomass,
                                    P2_end_biomass=P2_end_biomass,
                                    sum_biomass=P1_end_biomass+P2_end_biomass))
  }
}
}

head(sen_tIntro_K_P)
sen_tIntro_K_P$P1_Perc = sen_tIntro_K_P$P1_end_biomass/sen_tIntro_K_P$sum_biomass

######################
#### visualization
######################
require("plot3D")
## Loading required package: plot3D
## Warning: package 'plot3D' was built under R version 3.6.3
par(cex=2)
optimization.df <- read.csv("…\\all.regime.csv",header = T) # the sensitivity analysis results.
optimization.df$W_T <- optimization.df$W_T*149.6 # million m3
optimization.df$DOC_C <- 100-optimization.df$DOC_C*100 # percentage

x <- sen_tIntro_K_P$t_intro2
y <- sen_tIntro_K_P$K2_to_K1
z <- sen_tIntro_K_P$P2_to_P1
P1_Perc <- sen_tIntro_K_P$P1_Perc

#png("…\\Optimization.png",width = 16,height = 12,units = "cm",res = 300)
scatter3D(x,y,z,
          colvar = P1_Perc,
          xlab="Spec 2 introducing d",
          ylab="Carrying cap K2/K1",
          zlab="Init biomass P2/P1",
          clab = c("Spec 1 percentage"),
          pch=20,phi=0,theta = 60,bty="g",main="Optimization",
          ticktype="simple",cex=1)

# scatter3D(x=df.regime2$Q_thres,y=df.regime2$DOC_thres,z=df.regime2$W_T*149.6,colvar=df.regime2$DOC_C*100,add = T,type="l",colkey = F,lwd=10)

write.csv(sen_tIntro_K_P,".../1stCome1stServed/sen_tIntro_K_P.csv")
