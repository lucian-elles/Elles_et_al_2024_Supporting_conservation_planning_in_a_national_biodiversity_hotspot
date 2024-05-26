# Elles et al. (2024). Supporting conservation planning in a national 
# biodiversity hotspot -Projecting species composition across a groundwater level
# gradient using a demographic forest model

# This Script is used to run the PPA Models "02_PPA-model_tidy_no_FAH.R" and 
# the "03_PPA-model_tidy_no_SAH.R".

# The ""02_PPA-model_tidy_no_FAH.R" uses demographic rates from 2013 to 2020 for 
# moist, intermediate and dry sites in the LFF. The initial State from 2013 of 
# the respective site (dry, intermediate..etc) is used. The model excludes 
# Acer campestre (FAH). For more information see the method section in the thesis. 


# The ""03_PPA-model_tidy_no_SAH.R" uses demographic rates from 2016 to 2020 for 
# moist, intermediate and dry sites in the LFF. The initial State from 2016 of 
# the respective site (dry, intermediate..etc) is used. The model excludes 
# Acer platanoides (SAH). For more information see the method section in the thesis.

# For validation the basal area of the sites (dry, intermdeiate and moist) in 2020
# is displayed in the generated plots
 
# When using this script make sure to have the initial state, demographic rates
# and basal area in 2020 for each scenario (dry, intermediate and moist) togehter 
# in a folder seperated from other scenarios.


library(tidyverse)
library(ggpubr)
library(patchwork)

# Run the Model 6 Plots moist, intermediate (middle) and dry for 2013 - 2020 
# (Initial State 2013) and 2016 - 2020 (Initial State 2016) for 100 years  

#set the Source for the scripts for the PPA Model are located 
sources <- list("model_scripts/02_PPA_model_tidy_no_FAH.R","model_scripts/02_PPA_model_tidy_no_FAH.R","model_scripts/02_PPA_model_tidy_no_FAH.R",
                "model_scripts/03_PPA_model_tidy_no_SAH.R","model_scripts/03_PPA_model_tidy_no_SAH.R","model_scripts/03_PPA_model_tidy_no_SAH.R")
#set the Source for the mainfolder where the demographic rates are saved
mainfolders <-list("demographic_rates/rates_13_moist_no_FAH", "demographic_rates/rates_13_middle_no_FAH","demographic_rates/rates_13_dry_no_FAH",
                  "demographic_rates/rates_16_moist_no_SAH", "demographic_rates/rates_16_middle_no_SAH","demographic_rates/rates_16_dry_no_SAH")
#set the Source for the respective demographic rate data file (must be located within the mainfolder)
spdata_paths <- list("rates_13_moist_no_FAH.txt", "rates_13_middle_no_FAH.txt","rates_13_dry_no_FAH.txt",
                     "rates_16_moist_no_SAH.txt", "rates_16_middle_no_SAH.txt","rates_16_dry_no_SAH.txt")
#set the Source for the respect initial state data file (must be located within the mainfolder)
initdata_paths <- list("PPA_initial_state_moist13_no_FAH.txt", "PPA_initial_state_middle13_no_FAH.txt","PPA_initial_state_dry13_no_FAH.txt",
                       "PPA_initial_state_moist16_no_SAH.txt", "PPA_initial_state_middle16_no_SAH.txt","PPA_initial_state_dry16_no_SAH.txt")
#basal area in 2020 for model validation
badata_paths <- list("ba_dat_moist_13.txt", "ba_dat_middle_13.txt","ba_dat_dry_13.txt",
                     "ba_dat_moist_16.txt", "ba_dat_middle_16.txt","ba_dat_dry_16.txt")



all_plots <-list()
for (x  in 1:6 ){
  id = ((x-1) %/% 5)+1
  print(id)
  source(sources[[x]])
  results=main_cohorts_tidy(mainfolder = mainfolders[[x]],
                            spdata_path = spdata_paths[[x]] ,
                            initdata_path = initdata_paths[[x]],
                            growth = "rates",
                            mort = "rates")
  a <- list(results[[1]], results[[2]], results[[3]])
  plot_temp <- results[[4]]
  badata = read_table(paste(mainfolders[[x]], badata_paths[[x]], sep="/"))
  plot_temp <- plot_temp +
    geom_point(data = badata, aes(x = time, y = ba, color = sp), size = 7) +
    scale_color_manual(values = c("SAH" = "#228B22",
                                  "BAH" = "#00FF00",
                                  "HBU" = "#FFA500",
                                  "GES" = "blue",
                                  "SEI" = "red",
                                  "UL" = "#17becf",
                                  "WLI" = "grey",
                                  "FAH" = "#808000")) 
  plot(plot_temp)
  # Extract relevant information from results
  spdata_name <- sub(".txt", "", spdata_paths[[x]])  # Remove the file extension
  
  # Create a plot title with the spdata name
  plot_title <- paste("Plot for", spdata_name)
  
  # Define the path for saving the plot in the model_plot directory
  plot_filename <- file.path( paste0(spdata_name, "_Validierung.png"))
  
  png(file = plot_filename, width = 900, height = 600)  # Adjust width and height as needed
  plot(plot_temp, main = plot_title)
  dev.off()
  all_plots[[x]] <- plot_temp
}
# Combine all plots within one figure
my_plots <- wrap_plots(all_plots, ncol = 3)

plot(my_plots)


