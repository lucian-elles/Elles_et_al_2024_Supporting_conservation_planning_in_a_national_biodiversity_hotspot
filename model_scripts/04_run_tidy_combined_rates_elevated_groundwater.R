# Elles et al. (2024). Supporting conservation planning in a national 
# biodiversity hotspot -Projecting species composition across a groundwater level
# gradient using a demographic forest model

# This Script is used to run the PPA Model "01_PPA-model_tidy.R" for the 
# "Combined Rates" to project the species composition in the LFF with and 
# without an increase in Groundwater level. To simulate no increase in 
# groundwater level, the demographic rates for  moist, intermediate and dry 
# sites are used, together with the respective initial state from 2020. When 
# using this script make sure to have the initial state and demographic rates 
# for each scenario (dry, intermediate...etc) in a separate folder.



library(tidyverse)
library(ggpubr)
library(patchwork)


############ Run combined Initial state Moist

#set the Source for the mainfolder where the demographic rates are saved
mainfolders <-list("demographic_rates/rates_combined_census_moist","demographic_rates/rates_combined_census_wet","demographic_rates/rates_combined_census_v_wet")
#set the Source for the respective demographic rate data file (must be located within the mainfolder)
spdata_paths <- list("rates_combined_census_moist.txt","rates_combined_census_wet.txt","rates_combined_census_v_wet.txt")
#set the Source for the respect initial state data file (must be located within the mainfolder)
initdata_paths <- list("PPA_initial_state_moist20.txt","PPA_initial_state_moist20.txt","PPA_initial_state_moist20.txt")


all_plots <- list()
for (x in 1:3) {
  id <- ((x - 1) %/% 5) + 1
  print(id)
  source("model_scripts/01_PPA_model_tidy.R")
  results <- main_cohorts_tidy(mainfolder = mainfolders[[x]],
                               spdata_path = spdata_paths[[x]],
                               initdata_path = initdata_paths[[x]],
                               growth = "rates",
                               mort = "rates")
  
  # Extract relevant information from results
  spdata_name <- sub(".txt", "", spdata_paths[[x]])  # Remove the file extension
  
  # Create a plot title with the spdata name
  plot_title <- paste("Plot for", spdata_name)
  
  # Define the path for saving the plot in the model_plot directory
  plot_filename <- file.path(paste0(spdata_name, "_init_moist.png"))
  
  png(file = plot_filename, width = 900, height = 600)  # Adjust width and height as needed
  plot(results[[4]], main = plot_title)
  dev.off()
  
  # Store the plot in the list
  all_plots[[x]] <- results[[4]]
}

# Combine all plots within one figure
my_plots <- wrap_plots(all_plots, ncol = 3)

plot(my_plots)

############ Run combined Initial state Middle

#set the Source for the mainfolder where the demographic rates are saved
mainfolders <-list("demographic_rates/rates_combined_census_middle","demographic_rates/rates_combined_census_moist","demographic_rates/rates_combined_census_wet","demographic_rates/rates_combined_census_v_wet")
#set the Source for the respective demographic rate data file (must be located within the mainfolder)
spdata_paths <- list("rates_combined_census_middle.txt","rates_combined_census_moist.txt","rates_combined_census_wet.txt","rates_combined_census_v_wet.txt")
#set the Source for the respect initial state data file (must be located within the mainfolder)
initdata_paths <- list("PPA_initial_state_middle20.txt","PPA_initial_state_middle20.txt","PPA_initial_state_middle20.txt","PPA_initial_state_middle20.txt")

all_plots <- list()
for (x in 1:4) {
  id <- ((x - 1) %/% 5) + 1
  print(id)
  source("model_scripts/01_PPA_model_tidy.R")
  results <- main_cohorts_tidy(mainfolder = mainfolders[[x]],
                               spdata_path = spdata_paths[[x]],
                               initdata_path = initdata_paths[[x]],
                               growth = "rates",
                               mort = "rates")
  
  # Extract relevant information from results
  spdata_name <- sub(".txt", "", spdata_paths[[x]])  # Remove the file extension
  
  # Create a plot title with the spdata name
  plot_title <- paste("Plot for", spdata_name)
  
  # Define the path for saving the plot in the model_plot directory
  plot_filename <- file.path(paste0(spdata_name, "_init_intermediate.png"))
  
  png(file = plot_filename, width = 900, height = 600)  # Adjust width and height as needed
  plot(results[[4]], main = plot_title)
  dev.off()
  
  # Store the plot in the list
  all_plots[[x]] <- results[[4]]
}
# Combine all plots within one figure
my_plots <- wrap_plots(all_plots, ncol = 4)

plot(my_plots)

############ Run combined Initial state Moist

#set the Source for the mainfolder where the demographic rates are saved
mainfolders <-list("demographic_rates/rates_combined_census_dry","demographic_rates/rates_combined_census_middle","demographic_rates/rates_combined_census_moist","demographic_rates/rates_combined_census_wet","demographic_rates/rates_combined_census_v_wet")
#set the Source for the respective demographic rate data file (must be located within the mainfolder)
spdata_paths <- list("rates_combined_census_dry.txt","rates_combined_census_middle.txt","rates_combined_census_moist.txt","rates_combined_census_wet.txt","rates_combined_census_v_wet.txt")
#set the Source for the respect initial state data file (must be located within the mainfolder)
initdata_paths <- list("PPA_initial_state_dry20.txt","PPA_initial_state_dry20.txt","PPA_initial_state_dry20.txt","PPA_initial_state_dry20.txt","PPA_initial_state_dry20.txt")

all_plots <- list()
for (x in 1:5) {
  id <- ((x - 1) %/% 5) + 1
  print(id)
  source("model_scripts/01_PPA_model_tidy.R")
  results <- main_cohorts_tidy(mainfolder = mainfolders[[x]],
                               spdata_path = spdata_paths[[x]],
                               initdata_path = initdata_paths[[x]],
                               growth = "rates",
                               mort = "rates")
  
  # Extract relevant information from results
  spdata_name <- sub(".txt", "", spdata_paths[[x]])  # Remove the file extension
  
  # Create a plot title with the spdata name
  plot_title <- paste("Plot for", spdata_name)
  
  # Define the path for saving the plot in the model_plot directory
  plot_filename <- file.path(paste0(spdata_name, "_init_dry.png"))
  
  png(file = plot_filename, width = 900, height = 600)  # Adjust width and height as needed
  plot(results[[4]], main = plot_title)
  dev.off()
  
  # Store the plot in the list
  all_plots[[x]] <- results[[4]]
}
# Combine all plots within one figure
my_plots <- wrap_plots(all_plots, ncol = 5)

plot(my_plots)