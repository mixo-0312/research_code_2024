


rm(list=ls())
# Load required libraries


library(tidyverse)# plot and tidying

library(deSolve)  # For solving differential equations

library(readr) # read csv


# Read the CSV file
South_Africa_data <- read_csv("South_Africa_data.csv")

# Rename the columns (assuming 4 columns in total)
colnames(South_Africa_data) <- c("Index", "Date", "Infections", "Deaths")

# Convert all columns to character, add a new row
South_Africa_data <- South_Africa_data %>%
  mutate(across(everything(), as.character)) %>%  # Convert all columns to character to prevent data loss
  add_row(Index = "0", Date = "05/03/2020", Infections = "1", Deaths = "0", .before = 1)

# Convert columns back to numeric types where applicable
South_Africa_data <- South_Africa_data %>%
  mutate(
    Index = as.numeric(Index),
    Infections = as.numeric(Infections),
    Deaths = as.numeric(Deaths),
    Date = as.Date(Date, format = "%d/%m/%Y")  # Convert Date to Date type with the correct format
  )

# Filter the data to include only the date range from 05/03/2020 to 12/03/2020
week_one <- South_Africa_data %>%
  filter(Date >= as.Date("2020-03-05") & Date <= as.Date("2020-03-12"))

# View the result
week_one







# Filter the data to include only the date range from 05/03/2020 to 17/03/2020
week_two <- South_Africa_data %>%
  filter(Date >= as.Date("2020-03-12") & Date <= as.Date("2020-03-19"))

# View the result
week_two


# Filter the data to include only the date range from 05/03/2020 to 17/03/2020
week_three <- South_Africa_data %>%
  filter(Date >= as.Date("2020-03-19") & Date <= as.Date("2020-03-26"))

# View the result
week_three



# Filter the data to include only the date range from 05/03/2020 to 17/03/2020
#week_four <- South_Africa_data %>% mutate(Infections = as.numeric(Infections), Deaths=as.numeric(Deaths)) %>% filter(Date >= as.Date("2020-03-26") & Date <= as.Date("2020-04-02"))
week_four <- South_Africa_data %>%  filter(Date >= as.Date("2020-03-26") & Date <= as.Date("2020-04-02"))


# View the result
week_four




# View the results
print(week_one)
print(week_two)



# Define the SIRD model
sird_model <- function(time, state, parameters) {
  S <- state[1]
  I <- state[2]
  R <- state[3]
  D <- state[4]
  N <- S + I + R + D  # Total population size
  
  beta <- parameters["beta"]
  gamma <- parameters["gamma"]
  alpha <- parameters["alpha"]
  
  dS <- -beta * S * I / N
  dI <- beta * S * I / N - alpha * I - gamma * I
  dR <- alpha * I
  dD <- gamma * I
  
  list(c(dS, dI, dR, dD))
}

# Initial conditions for week 1
initial_state <- c(S = 59309999, I = 1, R = 0, D = 0)
abc_rejection <- function(real_data, tolerance, alpha, num_accepted_params) {
  accepted_parameters <- list()
  accepted_count <- 0
  
  for (i in 1:(max_num_of_iterations)) {  # Control the number of total iterations, increase if acceptance is slow
    # Generate random parameters for beta and gamma
    beta <- runif(1, 0, 1)  # Transmission rate
    gamma <- runif(1, 0, 1)  # Death rate
    
    # Parameters for this simulation
    parameters <- c(beta = beta, gamma = gamma, alpha = alpha)
    
    # Time points (7 days for week 1)
    times <- seq(0, 7, by = 1)
    
    # Solve the ODE system for this parameter set
    out <- ode(y = initial_state, times = times, func = sird_model, parms = parameters)
    
    # Extract the infections and deaths from the simulation
    simulated_infections <- out[, "I"]
    simulated_deaths <- out[, "D"]
    
    # Compare the real data to the simulated data (use Euclidean distance)
    distance <- sqrt(sum((simulated_infections - week_one$Infections)^2) +
                       sum((simulated_deaths - week_one$Deaths)^2))
    
    # Accept or reject the parameters based on the tolerance
    if (distance < tolerance) {
      accepted_parameters[[length(accepted_parameters) + 1]] <- c(beta, gamma)
      accepted_count <- accepted_count + 1
    }
    
    
    
    # Break the loop when we reach the desired number of accepted parameters
    if (accepted_count >= num_accepted_params) {
      break
    }
  }
  
  return(do.call(rbind, accepted_parameters))  # Return accepted parameters
}

# Set the tolerance, alpha, and number of accepted parameters
tolerance <-4
alpha <- 0.07  # Recovery rate
max_num_of_iterations<-10000000
num_accepted_params <- 1000  # Number of accepted parameter sets

# Using week_one for the real data
set.seed(12345L)
accepted_params_w1 <- abc_rejection(real_data = week_one, tolerance = tolerance, alpha = alpha, num_accepted_params = num_accepted_params)


#saveRDS(accepted_params_w1, "accepted_params_week1.rds")
## Reading the RDS file
#accepted_params_w1 <- readRDS("accepted_params_week1.rds")


#saveRDS(accepted_params, "accepted_params_week2.rds")
#saveRDS(accepted_params, "accepted_params_week3.rds")
# Display the accepted parameters

print(accepted_params_w1)





# Calculate the mean of the accepted beta and gamma
mean_beta <- mean(accepted_params_w1[,1])
mean_gamma <- mean(accepted_params_w1[,2])

cat("Estimated beta: ", mean_beta, "\n")
cat("Estimated gamma: ", mean_gamma, "\n")
# Generate simulated data using the inferred parameters
parameters <- c(beta = mean_beta, gamma = mean_gamma, alpha = alpha)
times <- seq(0, 7, by = 1)
simulated_output <- ode(y = initial_state, times = times, func = sird_model, parms = parameters)

# Extract results
simulated_infections <- simulated_output[, "I"]
simulated_deaths <- simulated_output[, "D"]




ggplot(week_one, aes(x = Date)) +
  geom_point(aes(y = Infections, color = "Observed Infections"), size = 4, shape = 16) +  # Filled circles
  geom_line(aes(y = simulated_infections, color = "Simulated Infections"), linetype = "dashed", size=1) +
  
  labs(title = "Observed Infections vs Model predicted Infections for Week 1", y = "Number of Individuals") +
  scale_color_manual(values = c("Observed Infections" = "skyblue", "Simulated Infections" = "red" )) +
  
  theme_minimal() +
  theme(panel.grid = element_blank())
ggplot(week_one, aes(x = Date)) +
  geom_point(aes(y = Deaths, color = "Observed Deaths"), size = 4,shape = 16) +
  geom_line(aes(y = simulated_deaths, color = "Simulated Deaths"), linetype = "dashed", size=1) +
  labs(title = "Observed Deaths  vs Simulated Deaths for Week 1", y = "Number of Individuals") +
  scale_color_manual(values = c("Observed Deaths" = "skyblue", "Simulated Deaths" = "red"))+
  
  theme_minimal()+
  theme(panel.grid = element_blank())  


# Plotting Observed Infections with filled circles and custom axis names
plot(week_one$Date, week_one$Infections, 
     pch = 21,         # Circle with fill
     bg = "lightblue",      # Fill color
     col = "black",    # Border color
     cex = 1.5,        # Point size
     xlab = "Date",    # X-axis label
     ylab = "Number of Individuals",  # Y-axis label
     main = "Observed Infections vs simulated Infections Week 1")  # Title

# Adding the simulated infections line
lines(week_one$Date, simulated_infections, col = "red", lty = "dashed", lwd = 2)
# Add a legend to differentiate between observed and simulated deaths
legend("topright", legend = c("Observed Infections", "Simulated Infections"), 
       col = c("skyblue", "red"), lty = c(1, 2), lwd = 2)

# Optionally, you can also adjust the axis limits or other graphical parameters as needed.
# Plot Observed Deaths
plot(week_one$Date, week_one$Deaths, 
     pch = 21,         # Circle with fill
     bg = "lightblue",      # Fill color
     col = "black",    # Border color
     cex = 1.5,        # Point size
     xlab = "Date",    # X-axis label
     ylab = "Number of Individuals",  # Y-axis label
     main = "Observed Deaths vs Simulated Deaths for Week 1",  # Title
     ylim = range(c(week_one$Deaths, simulated_deaths))  # Ensure both lines fit within the y-axis range
)

# Add Simulated Deaths Line (Dashed)
lines(week_one$Date, simulated_deaths, col = "red", lty = "dashed", lwd = 2)

# Add a legend to differentiate between observed and simulated deaths
legend("topright", legend = c("Observed Deaths", "Simulated Deaths"), 
       col = c("skyblue", "red"), lty = c(1, 2), lwd = 2)


# Convert the matrix to a data frame for easier manipulation
accepted_params_w1 <- as.data.frame(accepted_params_w1)

# Rename the columns (for better readability)
colnames(accepted_params_w1) <- c("beta", "gamma")


# Summary statistics
summary_stats <- accepted_params_w1 %>%
  summarise(
    mean_beta = mean(beta),
    median_beta = median(beta),
    sd_beta = sd(beta),
    mean_gamma = mean(gamma),
    median_gamma = median(gamma),
    sd_gamma = sd(gamma)
  )

# Calculate 95% credible intervals using the 2.5% and 97.5% quantiles
credible_intervals <- accepted_params_w1 %>%
  summarise(
    beta_lower = quantile(beta, 0.025),
    beta_upper = quantile(beta, 0.975),
    gamma_lower = quantile(gamma, 0.025),
    gamma_upper = quantile(gamma, 0.975)
  )

# Display the results
summary_stats
credible_intervals

# Convert the accepted_params to a data frame for easier plotting
accepted_params_df <- data.frame(beta = accepted_params_w1[,1], gamma = accepted_params_w1[,2])

# Plot histogram for beta
ggplot(accepted_params_df, aes(x=beta)) +
  geom_histogram(binwidth=0.002, fill="lightcoral", color="white") +
  labs(title="Beta Distribution for week 1 estimates", x="Beta", y="Frequency") +
  theme_minimal()+
  theme(panel.grid = element_blank())  


# Plot histogram for gamma
ggplot(accepted_params_df, aes(x=gamma)) +
  geom_histogram(binwidth=0.002, fill="darkred", color="white") +
  labs(title="Gamma Distribution for week 1 estimates", x="Gamma", y="Frequency") +
  theme_minimal()+
  theme(panel.grid = element_blank())  


#########################################

'WEEK TWO'

# Calculate the mean of the accepted beta and gamma
mean_beta <- mean(accepted_params_w1[,1])

mean_gamma <- mean(accepted_params_w1[,2])

cat("Estimated beta: ", mean_beta, "\n")
cat("Estimated gamma: ", mean_gamma, "\n")
# Generate simulated data using the inferred parameters
parameters <- c(beta = mean_beta, gamma = mean_gamma, alpha = alpha)
times <- seq(0, 7, by = 1)
simulated_output <- ode(y = initial_state, times = times, func = sird_model, parms = parameters)

end_of_week1<-as.data.frame(tail(simulated_output,1))

initial_state_week2 <- c(S = end_of_week1[,2], I = end_of_week1[,3], R = end_of_week1[,4], D = end_of_week1[,5])



abc_rejection_week2 <- function(real_data, tolerance, alpha, num_accepted_params) {
  accepted_parameters <- list()
  accepted_count <- 0
  
  for (i in 1:(max_num_of_iterations)) {  # Control the number of total iterations, increase if acceptance is slow
    # Generate random parameters for beta and gamma
    beta <- runif(1, 0, 1)  # Transmission rate
    gamma <- runif(1, 0, 1)  # Death rate
    
    # Parameters for this simulation
    parameters <- c(beta = beta, gamma = gamma, alpha = alpha)
    
    # Time points (7 days for week 1)
    times <- seq(0, 7, by = 1)
    
    # Solve the ODE system for this parameter set
    out <- ode(y = initial_state_week2, times = times, func = sird_model, parms = parameters)
    
    # Extract the infections and deaths from the simulation
    simulated_infections <- out[, "I"]
    simulated_deaths <- out[, "D"]
    
    # Compare the real data to the simulated data (use Euclidean distance)
    distance <- sqrt(sum((simulated_infections - week_two$Infections)^2) +
                       sum((simulated_deaths - week_two$Deaths)^2))
    
    # Accept or reject the parameters based on the tolerance
    if (distance < tolerance) {
      accepted_parameters[[length(accepted_parameters) + 1]] <- c(beta, gamma)
      print(accepted_count)
      accepted_count <- accepted_count + 1
      print(accepted_count)
    }
    
    
    
    # Break the loop when we reach the desired number of accepted parameters
    if (accepted_count >= num_accepted_params) {
      break
    }
  }
  
  return(do.call(rbind, accepted_parameters))  # Return accepted parameters
}

# Set the tolerance, alpha, and number of accepted parameters
tolerance <- 21
alpha <- 0.07  # Recovery rate
max_num_of_iterations<-10000000
num_accepted_params <- 500  # Number of accepted parameter sets

# Using week_two for the real data
set.seed(12345L)

accepted_params_w2 <- abc_rejection_week2(real_data = week_two, tolerance = tolerance, alpha = alpha, num_accepted_params = num_accepted_params)


#saveRDS(accepted_params_w2, "accepted_params_week2.rds")
## Reading the RDS file
#accepted_params_w2 <- readRDS("accepted_params_week2.rds")


# Display the accepted parameters
print(accepted_params_w2)

# Calculate the mean of the accepted beta and gamma
mean_beta <- mean(accepted_params_w2[,1])
mean_gamma <- mean(accepted_params_w2[,2])


cat("Estimated beta: ", mean_beta, "\n")
cat("Estimated gamma: ", mean_gamma, "\n")
# Generate simulated data using the inferred parameters
parameters <- c(beta = mean_beta, gamma = mean_gamma, alpha = alpha)
times <- seq(0, 7, by = 1)
simulated_output <- ode(y = initial_state_week2, times = times, func = sird_model, parms = parameters)

# Extract results
simulated_infections <- simulated_output[, "I"]
simulated_deaths <- simulated_output[, "D"]




# Plot real vs simulated data for week two
ggplot(week_two, aes(x = Date)) +
  geom_line(aes(y = Infections, color = "Observed Infections"), size = 1) +
  geom_line(aes(y = simulated_infections, color = "Simulated Infections"), linetype = "dashed", size=1) +
  
  labs(title = "Observed Infections vs Model predicted Infections for Week 2", y = "Number of Individuals") +
  scale_color_manual(values = c("Observed Infections" = "blue", "Simulated Infections" = "red" ))+
  
  theme_minimal()+
  theme(panel.grid = element_blank())  

ggplot(week_two, aes(x = Date)) +
  geom_line(aes(y = Deaths, color = "Observed Deaths"), size = 1) +
  geom_line(aes(y = simulated_deaths, color = "Simulated Deaths"), linetype = "dashed", size=1) +
  labs(title = "Observed Deaths  vs Simulated Deaths for Week 2", y = "Number of Individuals") +
  scale_color_manual(values = c("Observed Deaths" = "black", "Simulated Deaths" = "red"))+
  
  theme_minimal()+
  theme(panel.grid = element_blank())  

# Plotting Observed Infections with filled circles and custom axis names
plot(week_two$Date, week_two$Infections, 
     pch = 21,         # Circle with fill
     bg = "lightblue",      # Fill color
     col = "black",    # Border color
     cex = 1.5,        # Point size
     xlab = "Date",    # X-axis label
     ylab = "Number of Individuals",  # Y-axis label
     main = "Observed Infections vs simulated Infections for Week 2")  # Title

# Adding the simulated infections line
lines(week_two$Date, simulated_infections, col = "red", lty = "dashed", lwd = 2)
# Add a legend to differentiate between observed and simulated deaths
legend("topright", legend = c("Observed Infections", "Simulated Infections"), 
       col = c("skyblue", "red"), lty = c(1, 2), lwd = 2)

# Optionally, you can also adjust the axis limits or other graphical parameters as needed.
# Plot Observed Deaths
plot(week_two$Date, week_two$Deaths, 
     pch = 21,         # Circle with fill
     bg = "lightblue",      # Fill color
     col = "black",    # Border color
     cex = 1.5,        # Point size
     xlab = "Date",    # X-axis label
     ylab = "Number of Individuals",  # Y-axis label
     main = "Observed Deaths vs Simulated Deaths for Week 2",  # Title
     ylim = range(c(week_two$Deaths, simulated_deaths))  # Ensure both lines fit within the y-axis range
)

# Add Simulated Deaths Line (Dashed)
lines(week_two$Date, simulated_deaths, col = "red", lty = "dashed", lwd = 2)

# Add a legend to differentiate between observed and simulated deaths
legend("topright", legend = c("Observed Deaths", "Simulated Deaths"), 
       col = c("skyblue", "red"), lty = c(1, 2), lwd = 2)


# Convert the matrix to a data frame for easier manipulation
accepted_params_w2 <- as.data.frame(accepted_params_w2)

# Rename the columns (for better readability)
colnames(accepted_params_w2) <- c("beta", "gamma")

# Summary statistics
summary_stats <- accepted_params_w2 %>%
  summarise(
    mean_beta = mean(beta),
    median_beta = median(beta),
    sd_beta = sd(beta),
    mean_gamma = mean(gamma),
    median_gamma = median(gamma),
    sd_gamma = sd(gamma)
  )

# Calculate 95% credible intervals using the 2.5% and 97.5% quantiles
credible_intervals <- accepted_params_w2 %>%
  summarise(
    beta_lower = quantile(beta, 0.025),
    beta_upper = quantile(beta, 0.975),
    gamma_lower = quantile(gamma, 0.025),
    gamma_upper = quantile(gamma, 0.975)
  )
# Display the results

summary_stats
credible_intervals




# Convert the accepted_params to a data frame for easier plotting
accepted_params_df <- data.frame(beta = accepted_params_w2[,1], gamma = accepted_params_w1[,2])

# Plot histogram for beta
ggplot(accepted_params_df, aes(x=beta)) +
  geom_histogram(binwidth=0.002, fill="lightcoral", color="white") +
  labs(title="Beta Distribution for week 2 estimates", x="Beta", y="Frequency") +
  theme_minimal()+
  theme(panel.grid = element_blank())  


# Plot histogram for gamma
ggplot(accepted_params_df, aes(x=gamma)) +
  geom_histogram(binwidth=0.002, fill="darkred", color="white") +
  labs(title="Gamma Distribution for week 2 estimates", x="Gamma", y="Frequency") +
  theme_minimal()+
  theme(panel.grid = element_blank())

##################################################################################
'Week three'



# Calculate the mean of the accepted beta and gamma
mean_beta <- mean(accepted_params_w2[,1])
mean_gamma <- mean(accepted_params_w2[,2])

cat("Estimated beta: ", mean_beta, "\n")
cat("Estimated gamma: ", mean_gamma, "\n")
# Generate simulated data using the inferred parameters
parameters <- c(beta = mean_beta, gamma = mean_gamma, alpha = alpha)
times <- seq(0, 7, by = 1)
simulated_output <- ode(y = initial_state_week2, times = times, func = sird_model, parms = parameters)

end_of_week2<-as.data.frame(tail(simulated_output,1))

initial_state_week3 <- c(S = end_of_week2[,2], I = end_of_week2[,3], R = end_of_week2[,4], D = end_of_week2[,5])



abc_rejection_week3 <- function(real_data, tolerance, alpha, num_accepted_params) {
  accepted_parameters <- list()
  accepted_count <- 0
  
  for (i in 1:(max_num_of_iterations)) {  # Control the number of total iterations, increase if acceptance is slow
    # Generate random parameters for beta and gamma
    beta <- runif(1, 0, 1)  # Transmission rate
    gamma <- runif(1, 0, 1)  # Death rate
    
    # Parameters for this simulation
    parameters <- c(beta = beta, gamma = gamma, alpha = alpha)
    
    # Time points (7 days for week 1)
    times <- seq(0, 7, by = 1)
    
    # Solve the ODE system for this parameter set
    out <- ode(y = initial_state_week3, times = times, func = sird_model, parms = parameters)
    
    # Extract the infections and deaths from the simulation
    simulated_infections <- out[, "I"]
    simulated_deaths <- out[, "D"]
    
    # Compare the real data to the simulated data (use Euclidean distance)
    distance <- sqrt(sum((simulated_infections - week_three$Infections)^2) +
                       sum((simulated_deaths - week_three$Deaths)^2))
    
    # Accept or reject the parameters based on the tolerance
    if (distance < tolerance) {
      accepted_parameters[[length(accepted_parameters) + 1]] <- c(beta, gamma)
      accepted_count <- accepted_count + 1
      print(accepted_count)
    }
    
    
    
    # Break the loop when we reach the desired number of accepted parameters
    if (accepted_count >= num_accepted_params) {
      break
    }
  }
  
  return(do.call(rbind, accepted_parameters))  # Return accepted parameters
}

# Set the tolerance, alpha, and number of accepted parameters
tolerance <- 80
alpha <- 0.07  # Recovery rate
max_num_of_iterations<-10000000
num_accepted_params <- 1000  # Number of accepted parameter sets, or 200 or 500

# Using week_two for the real data
accepted_params_w3 <- abc_rejection_week3(real_data = week_three, tolerance = tolerance, alpha = alpha, num_accepted_params = num_accepted_params)


#saveRDS(accepted_params_w3, "accepted_params_week3.rds")
## Reading the RDS file
#accepted_params_w3 <- readRDS("accepted_params_week3.rds")


# Display the accepted parameters
print(accepted_params_w3)

# Calculate the mean of the accepted beta and gamma
mean_beta <- mean(accepted_params_w3[,1])
mean_gamma <- mean(accepted_params_w3[,2])

cat("Estimated beta: ", mean_beta, "\n")
cat("Estimated gamma: ", mean_gamma, "\n")
# Generate simulated data using the inferred parameters
parameters <- c(beta = mean_beta, gamma = mean_gamma, alpha = alpha)
times <- seq(0, 7, by = 1)
simulated_output <- ode(y = initial_state_week3, times = times, func = sird_model, parms = parameters)

# Extract results
simulated_infections <- simulated_output[, "I"]
simulated_deaths <- simulated_output[, "D"]




# Plot real vs simulated data for week three
ggplot(week_three, aes(x = Date)) +
  geom_line(aes(y = Infections, color = "Observed Infections"), size = 1) +
  geom_line(aes(y = simulated_infections, color = "Simulated Infections"), linetype = "dashed", size=1) +
  
  labs(title = "Observed Infections vs simulated Infections for Week 3", y = "Number of Individuals") +
  scale_color_manual(values = c("Observed Infections" = "blue", "Simulated Infections" = "red" ))+
  
  theme_minimal()+
  theme(panel.grid = element_blank())  

ggplot(week_three, aes(x = Date)) +
  geom_line(aes(y = Deaths, color = "Observed Deaths"), size = 1) +
  geom_line(aes(y = simulated_deaths, color = "Simulated Deaths"), linetype = "dashed", size=1) +
  labs(title = "Observed Deaths  vs Simulated Deaths for Week 3", y = "Number of Individuals") +
  scale_color_manual(values = c("Observed Deaths" = "black", "Simulated Deaths" = "red"))+
  
  theme_minimal()+
  theme(panel.grid = element_blank()) 


# Plotting Observed Infections with filled circles and custom axis names
plot(week_three$Date, week_three$Infections, 
     pch = 21,         # Circle with fill
     bg = "lightblue",      # Fill color
     col = "black",    # Border color
     cex = 1.5,        # Point size
     xlab = "Date",    # X-axis label
     ylab = "Number of Individuals",  # Y-axis label
     main = "Observed Infections vs Model Predicted Infections for Week 3")  # Title

# Adding the simulated infections line
lines(week_three$Date, simulated_infections, col = "red", lty = "dashed", lwd = 2)
# Add a legend to differentiate between observed and simulated deaths
legend("topright", legend = c("Observed Infections", "Simulated Infections"), 
       col = c("skyblue", "red"), lty = c(1, 2), lwd = 2)

# Optionally, you can also adjust the axis limits or other graphical parameters as needed.
# Plot Observed Deaths
plot(week_three$Date, week_three$Deaths, 
     pch = 21,         # Circle with fill
     bg = "lightblue",      # Fill color
     col = "black",    # Border color
     cex = 1.5,        # Point size
     xlab = "Date",    # X-axis label
     ylab = "Number of Individuals",  # Y-axis label
     main = "Observed Deaths vs Simulated Deaths for Week 3",  # Title
     ylim = range(c(week_three$Deaths, simulated_deaths))  # Ensure both lines fit within the y-axis range
)

# Add Simulated Deaths Line (Dashed)
lines(week_three$Date, simulated_deaths, col = "red", lty = "dashed", lwd = 2)

# Add a legend to differentiate between observed and simulated deaths
legend("topright", legend = c("Observed Deaths", "Simulated Deaths"), 
       col = c("skyblue", "red"), lty = c(1, 2), lwd = 2)


# Convert the matrix to a data frame for easier manipulation
accepted_params_w3 <- as.data.frame(accepted_params_w3)

# Rename the columns (for better readability)
colnames(accepted_params_w3) <- c("beta", "gamma")

# Summary statistics
summary_stats <- accepted_params_w3 %>%
  summarise(
    mean_beta = mean(beta),
    median_beta = median(beta),
    sd_beta = sd(beta),
    mean_gamma = mean(gamma),
    median_gamma = median(gamma),
    sd_gamma = sd(gamma)
  )

# Calculate 95% credible intervals using the 2.5% and 97.5% quantiles
credible_intervals <- accepted_params_w3 %>%
  summarise(
    beta_lower = quantile(beta, 0.025),
    beta_upper = quantile(beta, 0.975),
    gamma_lower = quantile(gamma, 0.025),
    gamma_upper = quantile(gamma, 0.975)
  )

# Display the results
summary_stats
credible_intervals



# Convert the accepted_params to a data frame for easier plotting
accepted_params_df <- data.frame(beta = accepted_params_w3[,1], gamma = accepted_params_w3[,2])

# Plot histogram for beta
ggplot(accepted_params_df, aes(x=beta)) +
  geom_histogram(binwidth=0.0003, fill="lightcoral", color="white") +
  labs(title="Beta Distribution for week 3", x="Beta", y="Frequency") +
  theme_minimal()+
  theme(panel.grid = element_blank())  


# Plot histogram for gamma
ggplot(accepted_params_df, aes(x=gamma)) +
  geom_histogram(binwidth=0.0003, fill="darkred", color="white") +
  labs(title="Gamma Distribution for week 3", x="Gamma", y="Frequency") +
  theme_minimal()+
  theme(panel.grid = element_blank())

##########################################################################################################
'week four'

# Calculate the mean of the accepted beta and gamma
mean_beta <- mean(accepted_params_w3[,1])
mean_gamma <- mean(accepted_params_w3[,2])

cat("Estimated beta: ", mean_beta, "\n")
cat("Estimated gamma: ", mean_gamma, "\n")
# Generate simulated data using the inferred parameters
parameters <- c(beta = mean_beta, gamma = mean_gamma, alpha = alpha)
times <- seq(0, 7, by = 1)
simulated_output <- ode(y = initial_state_week3, times = times, func = sird_model, parms = parameters)

end_of_week3<-as.data.frame(tail(simulated_output,1))

initial_state_week4 <- c(S = end_of_week3[,2], I = end_of_week3[,3], R = end_of_week3[,4], D = end_of_week3[,5])



abc_rejection_week4 <- function(real_data, tolerance, alpha, num_accepted_params) {
  accepted_parameters <- list()
  accepted_count <- 0
  
  for (i in 1:(max_num_of_iterations)) {  # Control the number of total iterations, increase if acceptance is slow
    # Generate random parameters for beta and gamma
    beta <- runif(1, 0, 1)  # Transmission rate
    gamma <- runif(1, 0, 1)  # Death rate
    
    # Parameters for this simulation
    parameters <- c(beta = beta, gamma = gamma, alpha = alpha)
    
    # Time points (7 days for week 1)
    times <- seq(0, 7, by = 1)
    
    # Solve the ODE system for this parameter set
    out <- ode(y = initial_state_week4, times = times, func = sird_model, parms = parameters)
    
    # Extract the infections and deaths from the simulation
    simulated_infections <- out[, "I"]
    simulated_deaths <- out[, "D"]
    
    # Compare the real data to the simulated data (use Euclidean distance)
    distance <- sqrt(sum((simulated_infections - week_four$Infections)^2) +
                       sum((simulated_deaths - week_four$Deaths)^2))
    
    # Accept or reject the parameters based on the tolerance
    if (distance < tolerance) {
      accepted_parameters[[length(accepted_parameters) + 1]] <- c(beta, gamma)
      accepted_count <- accepted_count + 1
      print(accepted_count)
    }
    
    
    
    # Break the loop when we reach the desired number of accepted parameters
    if (accepted_count >= num_accepted_params) {
      break
    }
  }
  
  return(do.call(rbind, accepted_parameters))  # Return accepted parameters
}

# Set the tolerance, alpha, and number of accepted parameters
tolerance <- 150
alpha <- 0.07  # Recovery rate
max_num_of_iterations<-10000000
num_accepted_params <- 1000  # Number of accepted parameter sets, or 200 or 500

# Using week_two for the real data
accepted_params_w4 <- abc_rejection_week4(real_data = week_four, tolerance = tolerance, alpha = alpha, num_accepted_params = num_accepted_params)


#saveRDS(accepted_params_w4, "accepted_params_week4.rds")
## Reading the RDS file
#Week_four1000params <- readRDS("accepted_params_week4.rds")
#accepted_params_w4 <- readRDS("accepted_params_week4.rds")

# Display the accepted parameters
print(accepted_params_w4)

# Calculate the mean of the accepted beta and gamma
mean_beta <- mean(accepted_params_w4[,1])
mean_gamma <- mean(accepted_params_w4[,2])

cat("Estimated beta: ", mean_beta, "\n")

cat("Estimated gamma: ", mean_gamma, "\n")
# Generate simulated data using the inferred parameters
parameters <- c(beta = mean_beta, gamma = mean_gamma, alpha = alpha)
times <- seq(0, 7, by = 1)
simulated_output <- ode(y = initial_state_week4, times = times, func = sird_model, parms = parameters)

# Extract results
simulated_infections <- simulated_output[, "I"]
simulated_deaths <- simulated_output[, "D"]




# Plot real vs simulated data for week four
ggplot(week_four, aes(x = Date)) +
  geom_line(aes(y = Infections, color = "Observed Infections"), size = 1) +
  geom_line(aes(y = simulated_infections, color = "Simulated Infections"), linetype = "dashed", size=1) +
  
  labs(title = "Observed Infections vs simulated Infections for Week 4", y = "Number of Individuals") +
  scale_color_manual(values = c("Observed Infections" = "blue", "Simulated Infections" = "red" ))+
  
  theme_minimal()+
  theme(panel.grid = element_blank())  

ggplot(week_four, aes(x = Date)) +
  geom_line(aes(y = Deaths, color = "Observed Deaths"), size = 1) +
  geom_line(aes(y = simulated_deaths, color = "Simulated Deaths"), linetype = "dashed", size=1) +
  labs(title = "Observed Deaths  vs Simulated Deaths for Week 4", y = "Number of Individuals") +
  scale_color_manual(values = c("Observed Deaths" = "black", "Simulated Deaths" = "red"))+
  
  theme_minimal()+
  theme(panel.grid = element_blank()) 



# Plotting Observed Infections with filled circles and custom axis names
plot(week_four$Date, week_four$Infections, 
     pch = 21,         # Circle with fill
     bg = "lightblue",      # Fill color
     col = "black",    # Border color
     cex = 1.5,        # Point size
     xlab = "Date",    # X-axis label
     ylab = "Number of Individuals",  # Y-axis label
     main = "Observed Infections vs Model Predicted Infections for Week 4")  # Title

# Adding the simulated infections line
lines(week_four$Date, simulated_infections, col = "red", lty = "dashed", lwd = 2)
# Add a legend to differentiate between observed and simulated deaths
legend("topright", legend = c("Observed Infections", "Simulated Infections"), 
       col = c("skyblue", "red"), lty = c(1, 2), lwd = 2)

# Optionally, you can also adjust the axis limits or other graphical parameters as needed.
# Plot Observed Deaths
plot(week_four$Date, week_four$Deaths, 
     pch = 21,         # Circle with fill
     bg = "lightblue",      # Fill color
     col = "black",    # Border color
     cex = 1.5,        # Point size
     xlab = "Date",    # X-axis label
     ylab = "Number of Individuals",  # Y-axis label
     main = "Observed Deaths vs Simulated Deaths for Week 4",  # Title
     ylim = range(c(week_four$Deaths, simulated_deaths))  # Ensure both lines fit within the y-axis range
)

# Add Simulated Deaths Line (Dashed)
lines(week_four$Date, simulated_deaths, col = "red", lty = "dashed", lwd = 2)

# Add a legend to differentiate between observed and simulated deaths
legend("topright", legend = c("Observed Deaths", "Simulated Deaths"), 
       col = c("skyblue", "red"), lty = c(1, 2), lwd = 2)


# Convert the matrix to a data frame for easier manipulation
accepted_params_w4 <- as.data.frame(accepted_params_w4)

# Rename the columns (for better readability)
colnames(accepted_params_w4) <- c("beta", "gamma")

# Summary statistics
summary_stats <- accepted_params_w4 %>%
  summarise(
    mean_beta = mean(beta),
    median_beta = median(beta),
    sd_beta = sd(beta),
    mean_gamma = mean(gamma),
    median_gamma = median(gamma),
    sd_gamma = sd(gamma)
  )

# Calculate 95% credible intervals using the 2.5% and 97.5% quantiles
credible_intervals <- accepted_params_w4 %>%
  summarise(
    beta_lower = quantile(beta, 0.025),
    beta_upper = quantile(beta, 0.975),
    gamma_lower = quantile(gamma, 0.025),
    gamma_upper = quantile(gamma, 0.975)
  )

# Display the results
summary_stats
credible_intervals


# Convert the accepted_params to a data frame for easier plotting
accepted_params_df <- data.frame(beta = accepted_params_w4[,1], gamma = accepted_params_w4[,2])

# Plot histogram for beta
ggplot(accepted_params_df, aes(x=beta)) +
  geom_histogram(binwidth=0.0001, fill="lightcoral", color="white") +
  labs(title="Beta Distribution for week 4", x="Beta", y="Frequency") +
  theme_minimal()+
  theme(panel.grid = element_blank())  


# Plot histogram for gamma
ggplot(accepted_params_df, aes(x=gamma)) +
  geom_histogram(binwidth=0.0001, fill="darkred", color="white") +
  labs(title="Gamma Distribution for week 4", x="Gamma", y="Frequency") +
  theme_minimal()+
  theme(panel.grid = element_blank())


