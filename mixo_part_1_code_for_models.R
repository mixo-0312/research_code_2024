rm(list=ls())
install.packages("reshape2")

library(tidyverse)  # General data manipulation and plotting tools (includes dplyr, ggplot2, etc.)
library(deSolve)    # For solving differential equations (e.g., SIR/SIRD models)
library(reshape2)   # For reshaping data (e.g., melt for long-format)

sird_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -Beta * S * I / N
    dI <- Beta * S * I / N - Alpha * I-gamma*I
    dR <- Alpha * I 
    dD <- gamma * I
    
    list(c(dS, dI, dR,dD))
  })
}


# parameters

parameters <- c(Beta =0.75,Alpha=0.25, gamma =0.1,N = 10000)
initial_state <- c(S = 9999, I = 1, R = 0,D=0)  # initial conditions
times <- seq(0,100,1)

modeloutputsird<- as.data.frame( ode(y=initial_state, times=times, func=sird_model, parms = parameters))

ggplot(data = modeloutputsird, mapping = aes(x = time)) +
  geom_line(aes(y = R), color = 'blue') +
  geom_line(aes(y = I), color = 'red') +
  geom_line(aes(y = S), color = 'green') +
  geom_line(aes(y = D), color = 'black') +
  labs(x = "Time", y = "Population", title =  paste0("Dynamics of SIRD Model"))+
  theme_minimal() +
  theme(legend.position = "top")



sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -Beta * S * I / N
    dI <- Beta * S * I / N - Alpha * I
    dR <- Alpha * I 
    
    
    list(c(dS, dI, dR))
  })
}


# parameters

parameters <- c(Beta =0.75,Alpha=0.25,N = 10000)
initial_state <- c(S = 9999, I = 1, R = 0)  # initial conditions
times <- seq(0,100,1)

modeloutputsir<- as.data.frame( ode(y=initial_state, times=times, func=sir_model, parms = parameters))

ggplot(data = modeloutputsir, mapping = aes(x = time)) +
  geom_line(aes(y = S), color = 'blue') +
  geom_line(aes(y = I), color = 'red') +
  geom_line(aes(y = R), color = 'green') +
  labs(x = "Time", y = "Population", title =  paste0("Dynamics of SIR Model"))+
  theme_minimal() +
  theme(legend.position = "top")



# Function to perform sensitivity analysis and return a combined data frame
sensitivity_analysis <- function(model, params_list, initial_state, times) {
  results <- list()
  
  for (params in params_list) {
    model_output <- as.data.frame(ode(y = initial_state, times = times, func = model, parms = params))
    model_output$Label <- paste0("Beta=", params["Beta"], "_Alpha=", params["Alpha"],
                                 if ("gamma" %in% names(params)) paste0("_Gamma=", params["gamma"]) else "")
    
    # Rename the 'time' column to 'Time' for consistency
    colnames(model_output)[colnames(model_output) == "time"] <- "Time"
    
    results[[model_output$Label[1]]] <- model_output
  }
  
  # Combine all data frames into one
  combined_results <- do.call(rbind, results)
  
  return(combined_results)
}

'SIR MODEL results'
# Define parameter sets for SIR and SIRD models
sir_params_list_1 <- list(
  c(Beta = 0.75, Alpha = 0.25, N = 10000),
  c(Beta = 0.50, Alpha = 0.25, N = 10000),
  c(Beta = 1.00, Alpha = 0.25, N = 10000)
)


sir_params_list_2 <- list(
  c(Beta = 0.75, Alpha = 0.15, N = 10000),
  c(Beta = 0.75, Alpha = 0.25, N = 10000),
  c(Beta = 0.75, Alpha = 0.35, N = 10000)
)


# Adjusted plotting function to exclude 'Time' as a compartment
plot_sensitivity_results <- function(results, title) {
  # Melt the data, but exclude 'Time' from being treated as a compartment
  results_melt <- melt(results, id.vars = c("Time", "Label"), variable.name = "Compartment", value.name = "Value")
  
  # Plot
  ggplot(results_melt, aes(x = Time, y = Value, color = Label)) +
    geom_line() +
    facet_wrap(~ Compartment, scales = "free_y") +
    labs(title = title, x = "Time", y = "Population") +
    theme_minimal() +
    theme(panel.grid = element_blank())  # removes grid
}

# Now, you can use this function as before to generate the plots
sir_results_1 <- sensitivity_analysis(sir_model, sir_params_list_1, initial_state = c(S = 9999, I = 1, R = 0), times = times)
plot_sensitivity_results(sir_results_1, "SIR Model Parameter Sensitivity (Beta Variations)")

sir_results_2 <- sensitivity_analysis(sir_model, sir_params_list_2, initial_state = c(S = 9999, I = 1, R = 0), times = times)
plot_sensitivity_results(sir_results_2, "SIR Model Parameter Sensitivity (Alpha Variations)")


'SIRD MODEL'


###### SIRD MODEL RESULTS
# Define three sets of parameter variations for the SIRD model
sird_params_list_1 <- list(
  c(Beta = 0.75, Alpha = 0.25, gamma = 0.1, N = 10000),  # Variation 1
  c(Beta = 0.50, Alpha = 0.25, gamma = 0.05, N = 10000), # Variation 2
  c(Beta = 1.00, Alpha = 0.25, gamma = 0.2, N = 10000)   # Variation 3
)

sird_params_list_2 <- list(
  c(Beta = 0.75, Alpha = 0.15, gamma = 0.1, N = 10000),  # Variation 4
  c(Beta = 0.75, Alpha = 0.25, gamma = 0.05, N = 10000), # Variation 5
  c(Beta = 0.75, Alpha = 0.35, gamma = 0.2, N = 10000)   # Variation 6
)

sird_params_list_3 <- list(
  c(Beta = 0.85, Alpha = 0.2, gamma = 0.08, N = 10000),  # Variation 7
  c(Beta = 0.65, Alpha = 0.3, gamma = 0.1, N = 10000),   # Variation 8
  c(Beta = 0.9, Alpha = 0.4, gamma = 0.15, N = 10000)    # Variation 9
)

# Run sensitivity analysis and plot for SIRD model
sird_results_1 <- sensitivity_analysis(sird_model, sird_params_list_1, initial_state = c(S = 9999, I = 1, R = 0, D = 0), times = times)
plot_sensitivity_results(sird_results_1, "SIRD Model Parameter Sensitivity (Beta & Gamma Variations)")

sird_results_2 <- sensitivity_analysis(sird_model, sird_params_list_2, initial_state = c(S = 9999, I = 1, R = 0, D = 0), times = times)
plot_sensitivity_results(sird_results_2, "SIRD Model Parameter Sensitivity (Alpha & Gamma Variations)")



# Run sensitivity analysis for each parameter list (three separate analyses)
sird_results_1 <- sensitivity_analysis(sird_model, sird_params_list_1, initial_state = c(S = 9999, I = 1, R = 0, D = 0), times = times)
sird_results_2 <- sensitivity_analysis(sird_model, sird_params_list_2, initial_state = c(S = 9999, I = 1, R = 0, D = 0), times = times)
sird_results_3 <- sensitivity_analysis(sird_model, sird_params_list_3, initial_state = c(S = 9999, I = 1, R = 0, D = 0), times = times)

# Plot the results for each sensitivity analysis
plot_sensitivity_results(sird_results_1, "SIRD Model Sensitivity Analysis (Variation 1)")
plot_sensitivity_results(sird_results_2, "SIRD Model Sensitivity Analysis (Variation 2)")
plot_sensitivity_results(sird_results_3, "SIRD Model Sensitivity Analysis (Variation 3)")