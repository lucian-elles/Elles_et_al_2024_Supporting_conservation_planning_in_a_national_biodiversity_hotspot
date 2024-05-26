# Elles et al. (2024). Supporting conservation planning in a national 
# biodiversity hotspot -Projecting species composition across a groundwater level
# gradient using a demographic forest model
# 
# This script computes the demographic rates specific to groundwater-dependent 
# species, focusing on eight key species within the Leipzig Floodplain Forest 
# (LFF). Drawing from forest inventory data gathered as part of the 
# "Lebendige Luppe" project, encompassing 8491 trees surveyed between 2013 and 
# 2020, and 2016 to 2020, it incorporates spatial groundwater modeling to 
# determine the distance to the groundwater table for each tree, linked by their
# coordinates. The calculation of tree height and crown area relies on 
# allometric equations adapted from Holzwarth et al. (2015). The necessary 
# functions are encapsulated in the R script functions_lelu, available in the 
# electronic supplements. Note that the utilized data is not publicly accessible;
# for inquiries, please reach out to the authors.


library(tidyverse)
library(ggplot2)
library(openxlsx)

# library(readxl)
source("functions_lelu.R")

theme_set(theme_bw()) # for ggplot

# prepare the data
dat = readxl::read_excel("data/LL_Data_subplot_02_05_2023.xlsx", na = "NA") %>% # Add Forest inventory data here
  select(treeID = BAUM_CODE, 
         sp = BAUMART17, 
         sp20 = BAUMART20, 
         dbh = BHD17, 
         dbh20 = BHD20, 
         height = hG17, 
         height20 = hG20, 
         gx = Rechtswert, 
         gy = Hochwert,
         plot = PLOTID,
         splot= SubplotNr,
         subplot = SubplotNr,
         plot_area = Area,
         period = veg_period,
         status = Status,
         nutzbar_wachstum = Nutzbar_Wachstum,
         nutzbar_mort = Nutzbar_Mort,
         nutzbar_rec = Nutzbar_Recruit) %>% 
  
# convert datatypes to numeric (some are wrong)
  mutate(dbh = as.numeric(dbh),
       dbh20 = as.numeric(dbh20),
       height = as.numeric(height),
       height20 = as.numeric(height20)) %>% 

  # use info from 2020 if initial entry was NA: 
  mutate(sp = replace(sp, is.na(sp), sp20[is.na(sp)]), 
         height = replace(height, is.na(height), height20[is.na(height)])) %>%
  # create separate plot, tree and stem columns; replace DBH=0 by NA
  mutate(plotID = str_sub(treeID, 2, 3), 
         tree_no = str_sub(unlist(map(str_split(treeID, fixed(".")), 2)), 2, -1), 
         stem_no = str_sub(treeID, -1,-1), 
         dbh = replace(dbh, dbh == 0, NA), 
         dbh20 = replace(dbh20, dbh20 == 0, NA))
  #Combine Plotid and Subplots
dat = dat %>% 
 unite("plotID", c(plotID,splot), sep = "_")

#create a column for the original species 
dat = dat %>% 
  mutate(sp_orig = sp)
#replace Species in dat$sp by alternative species with more accurate allometrie data
dat = dat %>% 
  mutate(sp = case_when(
    sp == "Aes_hip" ~ "Til_sp.",
    sp == "Aln_glu" ~ "Til_sp.",
    sp == "Bet_pen" ~ "Til_sp.",
    sp == "Cor_ave" ~ "Til_sp.",
    sp == "Cor_san" ~ "Til_sp.",
    sp == "Cra_spe" ~ "Til_sp.",
    sp == "Euo_eur" ~ "Til_sp.",
    sp == "Fra_pen" ~ "Fra_exc",
    sp == "Lar_dec" ~ "Til_sp.",
    sp == "Mal_syl" ~ "Til_sp.",
    sp == "Pop_spe" ~ "Til_sp.",
    sp == "Pru_pad" ~ "Pru_avi",
    sp == "Que_rub" ~ "Que_sp.",
    sp == "Rha_car" ~ "Til_sp.",
    sp == "Rob_pse" ~ "Til_sp.",
    sp == "Sam_nig" ~ "Til_sp.",
    sp == "Ulm_spe" ~ "Ulm_sp.",
    TRUE ~ sp
  ))

# read in the allometry data (adapted from Holzwarth et al 2015, DOI 10.1098/rsos.140541, for detailed inforamtion see appendix)
allom = read_csv("data/allometries_sigmoidal_MES.csv", na = "NA") %>% 
  mutate(sp = paste(str_sub(sp_full, 1,3), str_sub(map(str_split(sp_full, fixed(" ")),2),1,3), sep="_"))



# bind crown area and height columns to dat: 
dat = dat %>% left_join(allom, by = "sp") %>% 
  mutate(ca = param1 / (1 + exp(-(inflection - dbh)/steepness)) + param2 / (1 + exp((inflection - dbh)/steepness))) %>% 
        mutate( height = replace(height, is.na(height), (h_a[is.na(height)]+h_b[is.na(height)])*(dbh[is.na(height)]/100)^h_c[is.na(height)]))

# assign trees to canopy layers: 
dat = dat %>% 
  left_join(assign.crown.layer.in.single.census(dat, nr.of.layers = 2), by = "treeID")

# calculate individual growth and add, if tree survived or is a recruit: 

dat = dat %>% 
  mutate(growth = dbh20-dbh, 
         growth_annual = ifelse(period == 2013, growth/7, growth/4), 
         mort = ifelse(status == "tot", 1, 0), 
         recr = ifelse(status == "NEU", 1, 0))



#calculate rates 

#assign original tree species to the trees, remove trees that are not relevant (NOt_rel) for the scenarios
#For the case of a scenario where we want to have a look at Pru_avi (Vogelkirsche) change last row
dat = dat %>% 
  mutate( sp = sp_orig) %>% 
  mutate(sp = case_when(
    sp == "Ace_cam" ~ "FAH",
    sp == "Ace_pla" ~ "SAH",
    sp == "Ace_pse" ~ "BAH",
    sp == "Car_bet" ~ "HBU",
    sp == "Fra_exc" ~ "GES",
    sp == "Que_rob" ~ "SEI",
    sp == "Aes_hip" ~ "Not_rel",
    sp == "Aln_glu" ~ "Not_rel",
    sp == "Bet_pen" ~ "Not_rel",
    sp == "Cor_ave" ~ "Not_rel",
    sp == "Cor_san" ~ "Not_rel",
    sp == "Cra_spe" ~ "Not_rel",
    sp == "Euo_eur" ~ "Not_rel",
    sp == "Fra_pen" ~ "GES",
    sp == "Lar_dec" ~ "Not_rel",
    sp == "Mal_syl" ~ "Not_rel",
    sp == "Pop_spe" ~ "Not_rel",
    sp == "Pru_pad" ~ "Not_rel",
    sp == "Que_rub" ~ "Not_rel",
    sp == "Que_sp." ~ "SEI",
    sp == "Rha_car" ~ "Not_rel",
    sp == "Rob_pse" ~ "Not_rel",
    sp == "Sam_nig" ~ "Not_rel",
    sp == "Ulm_spe" ~ "UL",
    sp == "Ulm_min" ~ "UL",
    sp == "Ulm_gla" ~ "UL",
    sp == "Til_cor" ~ "WLI",
    sp == "Til_pla" ~ "WLI",
    sp == "Fag_syl" ~ "Not_rel",
    sp == "Pru_avi" ~ "Not_rel",
    
    TRUE ~ sp
  ))

table(dat$sp)
#load "Distance to Groundwater table" ( in th following "gwfa") dataset and join it to the table 
gwfa = readxl::read_excel("data/groundwater_to_surface_distances_2013_2021.xlsx", sheet = 2 ,na = "NA") %>% 
  mutate(mean13_20 = ((GWFA_mean_2013+GWFA_mean_2014+GWFA_mean_2015+GWFA_mean_2016+GWFA_mean_2017+GWFA_mean_2018+GWFA_mean_2019+GWFA_mean_2020)/8),
         mean16_20 = ((GWFA_mean_2016+GWFA_mean_2017+GWFA_mean_2018+GWFA_mean_2019+GWFA_mean_2020)/5)) %>% 
  select(treeID = BAUM_CODE,
         mean13_20,
         mean16_20)
dat = dat %>% 
  left_join(gwfa, by = "treeID") %>%
  
  #Give trees the dustance to groundwater table value (gwfa) for there specific vegetation period
  mutate(gwfa_veg_mean = ifelse(period == 2013, mean13_20, mean16_20))

#Get the median distance to groundwater table on a plot level 
gwfa_plot <- readxl::read_excel("data/Feuchteklassen_stand_2018.xlsx", na = "NA") %>% 
  select(plot_median = median, plot = Plot_Nr_1, gwfa_blocks_plot = ...15) %>% 
  mutate(gwfa_blocks_plot = case_when(
    gwfa_blocks_plot == "feucht" ~ "moist",
    gwfa_blocks_plot == "mittel" ~ "middle",
    gwfa_blocks_plot == "trocken" ~ "dry",
    
    TRUE ~ gwfa_blocks_plot
  ))

#Bind the distance to groundwater table on plot level to the data 
dat = dat %>% 
  left_join(gwfa_plot, by = "plot") 

# Calculate the area of the inspected groundwater class (block) 
# (moist, intermediate and dry)
blocks_area_combined = dat %>% 
  group_by(plotID, gwfa_blocks_plot) %>% 
  summarise( area_subplots = mean(plot_area)) %>% 
  group_by(gwfa_blocks_plot) %>% 
  summarise(area_gwfa_block = sum(area_subplots))

blocks_area_2013 = dat %>% 
  filter(period == 2013) %>% 
  group_by(plotID, gwfa_blocks_plot, period) %>% 
  summarise( area_subplots = mean(plot_area)) %>% 
  group_by(gwfa_blocks_plot) %>% 
  summarise(area_gwfa_block = sum(area_subplots))

blocks_area_2017 = dat %>%
  filter(period == 2016) %>% 
  group_by(plotID, gwfa_blocks_plot, period) %>% 
  summarise( area_subplots = mean(plot_area)) %>% 
  group_by(gwfa_blocks_plot) %>% 
  summarise(area_gwfa_block = sum(area_subplots))

blocks_area_average = dat %>% 
  group_by(plotID) %>% 
  summarise( area_subplots = mean(plot_area)) %>% 
  summarise(area = sum(area_subplots))

blocks_area_average_13 = dat %>%
  filter(period == 2013) %>%
  group_by(plotID) %>% 
  summarise( area_subplots = mean(plot_area)) %>% 
  summarise(area = sum(area_subplots))

blocks_area_average_16 = dat %>%
  filter(period == 2016) %>%
  group_by(plotID) %>% 
  summarise( area_subplots = mean(plot_area)) %>% 
  summarise(area = sum(area_subplots))


gwfa_plot1 = gwfa_plot %>% 
  group_by(gwfa_blocks_plot) %>% 
  summarize(mean = mean(plot_median), max = max(plot_median), min = min(plot_median))

# prepare a species table:
sp_table = dat %>% 
  select(sp,inflection,steepness,param1, param2) %>% 
  unique() %>% 
  arrange(sp) %>% 
  filter(sp!="NA" & sp!="Not_rel") %>%
  drop_na()


#######################Calculate the recruitment rates per species and hectar in the periods###############################
r_13_moist = dat %>%
  filter(is.na(nutzbar_rec) & gwfa_blocks_plot == "moist"& period!=2016) %>% 
  mutate(cens_int = ifelse(period == 2013, 7, 4), 
         recr_year = recr/cens_int) %>%
  group_by(sp, period, gwfa_blocks_plot) %>%
  summarise(rec_year = sum(recr_year), 
            rec_year_ha = rec_year/(2.7425))#Value can be taken from data blocks_area_2013

r_13_middle = dat %>%
  filter(is.na(nutzbar_rec) & gwfa_blocks_plot == "middle"& period!=2016) %>% 
  mutate(cens_int = ifelse(period == 2013, 7, 4), 
         recr_year = recr/cens_int) %>%
  group_by(sp, period, gwfa_blocks_plot) %>%
  summarise(rec_year = sum(recr_year), 
            rec_year_ha = rec_year/(2.7661))#Value can be taken from data blocks_area_2013

r_13_dry = dat %>%
  filter(is.na(nutzbar_rec) & gwfa_blocks_plot == "dry"& period!=2016) %>% 
  mutate(cens_int = ifelse(period == 2013, 7, 4), 
         recr_year = recr/cens_int) %>%
  group_by(sp, period, gwfa_blocks_plot) %>%
  summarise(rec_year = sum(recr_year), 
            rec_year_ha = rec_year/(2.29154))#Value can be taken from data blocks_area_2013

r_16_moist = dat %>%
  filter(is.na(nutzbar_rec) & gwfa_blocks_plot == "moist"& period!=2013) %>% 
  mutate(cens_int = ifelse(period == 2013, 7, 4), 
         recr_year = recr/cens_int) %>%
  group_by(sp, period, gwfa_blocks_plot) %>%
  summarise(rec_year = sum(recr_year), 
            rec_year_ha = rec_year/(2.2516))#Value can be taken from data blocks_area_2017

r_16_middle = dat %>%
  filter(is.na(nutzbar_rec) & gwfa_blocks_plot == "middle"& period!=2013) %>% 
  mutate(cens_int = ifelse(period == 2013, 7, 4), 
         recr_year = recr/cens_int) %>%
  group_by(sp, period, gwfa_blocks_plot) %>%
  summarise(rec_year = sum(recr_year), 
            rec_year_ha = rec_year/(2.2535))#Value can be taken from data blocks_area_2017

r_16_dry = dat %>%
  filter(is.na(nutzbar_rec) & gwfa_blocks_plot == "dry"& period!=2013) %>% 
  mutate(cens_int = ifelse(period == 2013, 7, 4), 
         recr_year = recr/cens_int) %>%
  group_by(sp, period, gwfa_blocks_plot) %>%
  summarise(rec_year = sum(recr_year), 
            rec_year_ha = rec_year/(2.765))#Value can be taken from data blocks_area_2017



########################## Use a linear model (regression) to predict the annual growth #######################

predictions <- data.frame(
  sp = character(),
  cl = character(),
  period = character(),
  gwfa_veg_mean = numeric(),
  prediction = numeric()
)

dat_growth <-  dat %>% 
  filter(sp != "Not_rel" & !is.na(cl) & is.na(nutzbar_wachstum) & !is.na(growth_annual) & growth_annual >= -1 & growth_annual <= 2.5)

groups <- dat_growth %>%select(sp, cl, period) %>%  group_by(sp, cl, period) %>% slice(1) %>% ungroup()

for (i in 1:nrow(groups)) {
  group <- groups[i, ]
  group_data <- dat_growth %>% filter(sp == group$sp, cl == group$cl, period == group$period)
  group_model <- lm(growth_annual ~ gwfa_veg_mean, data = group_data)
  new_data <- data.frame(
    sp = group$sp,
    cl = group$cl,
    period = group$period,
    gwfa_veg_mean = c(0.3,0.8,1.3, 1.8, 2.3)
  )
  new_data$prediction <- predict(group_model, newdata = new_data, type = "response")
  predictions <- rbind(predictions, new_data)
}


##################Use the predicted data to calculate the rates############################

growth_rates = predictions %>% 
  mutate(growth_annual = prediction) %>% 
  mutate(gwl = gwfa_veg_mean) %>%
  mutate(gwl = case_when(
    gwl == 0.3 ~ "very wet",
    gwl == 0.8 ~ "wet",
    gwl == 1.3 ~ "moist",
    gwl == 1.8 ~ "middle",
    gwl == 2.3 ~ "dry",
   
    
    TRUE ~ as.character(gwl)
  )) %>%
  select(-prediction, -gwfa_veg_mean )


################### Predict Mortality Values for different GWL ############################
predictions <- data.frame(
  sp = character(),
  cl = character(),
  period = character(),
  gwfa_veg_mean = numeric(),
  prediction = numeric()
)
dat_mort <- dat %>%
  filter(
    !is.na(cl) &
      is.na(nutzbar_mort) &
      sp != "Not_rel")


groups <- dat_mort %>%select(sp, cl, period) %>%  group_by(sp, cl, period) %>% slice(1) %>% ungroup()

for (i in 1:nrow(groups)) {
  group <- groups[i, ]
  group_data <- dat_mort %>% filter(sp == group$sp, cl == group$cl, period == group$period)
  group_model <- glm(mort ~ gwfa_veg_mean, family = binomial, data = group_data)
  new_data <- data.frame(
    sp = group$sp,
    cl = group$cl,
    period = group$period,
    gwfa_veg_mean = c(0.3,0.8,1.3, 1.8, 2.3)
  )
  new_data$prediction <- predict(group_model, newdata = new_data, type = "response")
  predictions <- rbind(predictions, new_data)
}

##################Use the predicted data to calculate the rates############################

mort_rates = predictions %>% 
  mutate(mort = prediction) %>% 
  mutate(mort_annual = round(-log(1-mort)/ ifelse(period == 2013, 7, 4), digits=4)) %>%
  mutate(gwl = gwfa_veg_mean) %>% 
  mutate(gwl = case_when(
    gwl == 0.3 ~ "very wet",
    gwl == 0.8 ~ "wet",
    gwl == 1.3 ~ "moist",
    gwl == 1.8 ~ "middle",
    gwl == 2.3 ~ "dry",
    TRUE ~ as.character(gwl)
  )) %>% 
  select(-mort, -gwfa_veg_mean, -prediction)


##################### create the rates for 2013 - 2020 ######################################

g_13_wet = growth_rates %>% filter(period == 2013 & gwl == "wet") %>% select(-gwl)
g_13_dry = growth_rates %>% filter(period == 2013 & gwl == "dry") %>% select(-gwl)
g_13_middle = growth_rates %>%  filter(period == 2013 & gwl == "middle")%>% select(-gwl)
g_13_moist = growth_rates %>%  filter(period == 2013 & gwl == "moist")%>% select(-gwl)
g_13_v_wet = growth_rates %>% filter(period == 2013 & gwl == "very wet") %>% select(-gwl)

m_13_wet = mort_rates %>%  filter(period == 2013 & gwl == "wet")%>% select(-gwl)
m_13_dry = mort_rates %>%  filter(period == 2013 & gwl == "dry")%>% select(-gwl)
m_13_middle = mort_rates %>%  filter(period == 2013 & gwl == "middle")%>% select(-gwl)
m_13_moist = mort_rates %>%  filter(period == 2013 & gwl == "moist")%>% select(-gwl)
m_13_v_wet = mort_rates %>%  filter(period == 2013 & gwl == "very wet")%>% select(-gwl)

#####################   create the rates for 2016 - 2020 ######################################

g_16_wet = growth_rates %>% filter(period == 2016 & gwl == "wet") %>% select(-gwl)
g_16_dry = growth_rates %>% filter(period == 2016 & gwl == "dry") %>% select(-gwl)
g_16_middle = growth_rates %>%  filter(period == 2016 & gwl == "middle")%>% select(-gwl)
g_16_moist = growth_rates %>%  filter(period == 2016 & gwl == "moist")%>% select(-gwl)
g_16_v_wet = growth_rates %>% filter(period == 2016 & gwl == "very wet") %>% select(-gwl)

m_16_wet = mort_rates %>%  filter(period == 2016 & gwl == "wet")%>% select(-gwl)
m_16_dry = mort_rates %>%  filter(period == 2016 & gwl == "dry")%>% select(-gwl)
m_16_middle = mort_rates %>%  filter(period == 2016 & gwl == "middle")%>% select(-gwl)
m_16_moist = mort_rates %>%  filter(period == 2016 & gwl == "moist")%>% select(-gwl)
m_16_v_wet = mort_rates %>%  filter(period == 2016 & gwl == "very wet")%>% select(-gwl)

##### Assumption regarding Quercus robur #####

# There are only 9 recorded oak individuals in the understory that survived 
# until the second forest inventory in 2020. In 2013, there were only 3 individuals,
# and in 2016, there were 6 individuals, respectively. Therefore, it's not advisable
# to use a linear model with such few data points. The results indicate that for 
# very wet areas, the growth rate is negative. Accordingly, assumptions need to be made. 
# Therefore we use the mean value of all oak individuals for all groundwater level classifications.

# mean_oak_growth <- dat %>% 
#             filter( sp == "SEI" & cl == 2 & status == "2017 und 2020 vorhanden") %>% 
#   summarise(mean = mean(growth_annual))

g_13_dry = g_13_dry %>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))
g_13_middle = g_13_middle%>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))
g_13_moist = g_13_moist%>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))
g_13_wet = g_13_wet%>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))
g_13_v_wet = g_13_v_wet%>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))

g_16_dry = g_16_dry %>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))
g_16_middle = g_16_middle%>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))
g_16_moist = g_16_moist%>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))
g_16_wet = g_16_wet%>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))
g_16_v_wet = g_16_v_wet%>%
  mutate(growth_annual = if_else(sp == "SEI" & cl == 2 ,0.2416, growth_annual))

##### Assumption 2 regarding Quercus robur #####

# Due to low abundance of Quercus robur in the understory especially in 2013-2020
# the mean mortality rate of Quercus robur in 2016 - 2020 is applied as an substitute 
# independent on plots or distance to groundwater.


m_13_dry = m_13_dry %>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))
m_13_middle = m_13_middle%>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))
m_13_moist = m_13_moist%>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))
m_13_wet = m_13_wet%>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))
m_13_v_wet = m_13_v_wet%>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))

m_16_dry = m_16_dry %>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))
m_16_middle = m_16_middle%>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))
m_16_moist = m_16_moist%>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))
m_16_wet = m_16_wet%>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))
m_16_v_wet = m_16_v_wet%>%
  mutate(mort_annual = if_else(sp == "SEI" & cl == 2 ,0.0828, mort_annual))


#######################Calculate the rates for 2013 Census##############
rates_13_dry = g_13_dry %>% 
  full_join(m_13_dry) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_13_dry %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_13_middle = g_13_middle %>% 
  full_join(m_13_middle) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_13_middle %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_13_moist = g_13_moist %>% 
  full_join(m_13_moist) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_13_moist %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2)

# make the table for the extreme scenarios use the recruitment of moist for wet and dry for very dry
rates_13_v_wet = g_13_v_wet %>% 
  full_join(m_13_v_wet) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_13_moist %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_13_wet = g_13_wet %>% 
  full_join(m_13_wet) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_13_moist %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 
#######################Calculate the rates for 2016 Census##############

rates_16_dry = g_16_dry %>% 
  full_join(m_16_dry) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_16_dry %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_16_middle = g_16_middle %>% 
  full_join(m_16_middle) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_16_middle %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_16_moist = g_16_moist %>% 
  full_join(m_16_moist) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_16_moist %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2)

# make the table for the extreme scenarios use the recruitment of moist for wet and dry for very dry
rates_16_v_wet = g_16_v_wet %>% 
  full_join(m_16_v_wet) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_16_moist %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_16_wet = g_16_wet %>% 
  full_join(m_16_wet) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_16_moist %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2)



############################## save the dataframes ############################

######Assumptions regarding Acer campestre and Acer platanoides##########

# In 2013, only very few individuals of the Field Maple (FAH) Acer campstre 
# were mapped. Therefore, the data foundation for calculating demographic data 
# is inadequate.
# Solution approach: Therefore Acer campestre is filtered out when creating the 
# demographic rates for 2013. 
# For calculating the "Combined rates" the value from 2016 is taken as an substitute



# In 2016, only very few individuals of the Field Maple (SAH) Acer platanoides 
# were mapped. Therefore, the data foundation for calculating demographic data 
# is inadequate.
# Solution approach: Therefore Acer platanoides is filtered out when creating the 
# demographic rates for 2016. 
# For calculating the "Combined rates" the value from 2013 is taken as an substitute


rates_13_dry_no_FAH <- rates_13_dry %>% 
  filter(sp != "FAH")
rates_13_middle_no_FAH<- rates_13_middle %>% 
  filter(sp != "FAH")
rates_13_moist_no_FAH<- rates_13_moist %>% 
  filter(sp != "FAH")
rates_13_v_wet_no_FAH<- rates_13_v_wet %>% 
  filter(sp != "FAH")
rates_13_wet_no_FAH<- rates_13_wet %>% 
  filter(sp != "FAH")
rates_16_dry_no_SAH  <- rates_16_dry %>% 
  filter(sp != "SAH")
rates_16_middle_no_SAH <- rates_16_middle %>% 
  filter(sp != "SAH")
rates_16_moist_no_SAH <- rates_16_moist %>% 
  filter(sp != "SAH")
rates_16_v_wet_no_SAH <- rates_16_v_wet %>% 
  filter(sp != "SAH")
rates_16_wet_no_SAH <- rates_16_wet %>% 
  filter(sp != "SAH")

write.table(rates_13_dry_no_FAH, "demographic_rates/rates_13_dry_no_FAH/rates_13_dry_no_FAH.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(rates_13_middle_no_FAH, "demographic_rates/rates_13_middle_no_FAH/rates_13_middle_no_FAH.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(rates_13_moist_no_FAH, "demographic_rates/rates_13_moist_no_FAH/rates_13_moist_no_FAH.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(rates_16_dry_no_SAH, "demographic_rates/rates_16_dry_no_SAH/rates_16_dry_no_SAH.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(rates_16_middle_no_SAH, "demographic_rates/rates_16_middle_no_SAH/rates_16_middle_no_SAH.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(rates_16_moist_no_SAH, "demographic_rates/rates_16_moist_no_SAH/rates_16_moist_no_SAH.txt", sep="\t", row.names=FALSE, quote = FALSE)
########################### calculate the average rates for both census intervals ########################

#Take the mean value of annual growth between the two periods
g_c_dry <- g_13_dry %>% 
  left_join(g_16_dry, by = c("sp","cl")) %>% 
  mutate(growth_annual = rowMeans(.[,c("growth_annual.x","growth_annual.y")]),
         growth_annual = ifelse(sp == "FAH", growth_annual.y, growth_annual),
         growth_annual = ifelse(sp == "SAH", growth_annual.x, growth_annual))

g_c_middle <- g_13_middle %>% 
  left_join(g_16_middle, by = c("sp","cl")) %>% 
  mutate(growth_annual = rowMeans(.[,c("growth_annual.x","growth_annual.y")]),
         growth_annual = ifelse(sp == "FAH", growth_annual.y, growth_annual),
         growth_annual = ifelse(sp == "SAH", growth_annual.x, growth_annual))

g_c_moist <- g_13_moist %>% 
  left_join(g_16_moist, by = c("sp","cl")) %>% 
  mutate(growth_annual = rowMeans(.[,c("growth_annual.x","growth_annual.y")]),
         growth_annual = ifelse(sp == "FAH", growth_annual.y, growth_annual),
         growth_annual = ifelse(sp == "SAH", growth_annual.x, growth_annual))

g_c_v_wet <- g_13_v_wet %>% 
  left_join(g_16_v_wet, by = c("sp","cl")) %>% 
  mutate(growth_annual = rowMeans(.[,c("growth_annual.x","growth_annual.y")]),
         growth_annual = ifelse(sp == "FAH", growth_annual.y, growth_annual),
         growth_annual = ifelse(sp == "SAH", growth_annual.x, growth_annual))

g_c_wet <- g_13_wet %>% 
  left_join(g_16_wet, by = c("sp","cl")) %>% 
  mutate(growth_annual = rowMeans(.[,c("growth_annual.x","growth_annual.y")]),
         growth_annual = ifelse(sp == "FAH", growth_annual.y, growth_annual),
         growth_annual = ifelse(sp == "SAH", growth_annual.x, growth_annual))

#Take the mean value of annual mortality between the two periods
m_c_dry <- m_13_dry %>% 
  left_join(m_16_dry, by = c("sp","cl")) %>% 
  mutate(mort_annual = rowMeans(.[,c("mort_annual.x","mort_annual.y")]),
         mort_annual = ifelse(sp == "FAH", mort_annual.y, mort_annual),
         mort_annual = ifelse(sp == "SAH", mort_annual.x, mort_annual))

m_c_middle <- m_13_middle %>% 
  left_join(m_16_middle, by = c("sp","cl")) %>% 
  mutate(mort_annual = rowMeans(.[,c("mort_annual.x","mort_annual.y")]),
         mort_annual = ifelse(sp == "FAH", mort_annual.y, mort_annual),
         mort_annual = ifelse(sp == "SAH", mort_annual.x, mort_annual))


m_c_moist <- m_13_moist %>% 
  left_join(m_16_moist, by = c("sp","cl")) %>%
  mutate(mort_annual = rowMeans(.[,c("mort_annual.x","mort_annual.y")]),
         mort_annual = ifelse(sp == "FAH", mort_annual.y, mort_annual),
         mort_annual = ifelse(sp == "SAH", mort_annual.x, mort_annual))

m_c_v_wet <- m_13_v_wet %>% 
  left_join(m_16_v_wet, by = c("sp","cl")) %>% 
  mutate(mort_annual = rowMeans(.[,c("mort_annual.x","mort_annual.y")]),
         mort_annual = ifelse(sp == "FAH", mort_annual.y, mort_annual),
         mort_annual = ifelse(sp == "SAH", mort_annual.x, mort_annual))

m_c_wet <- m_13_wet %>% 
  left_join(m_16_wet, by = c("sp","cl")) %>% 
  mutate(mort_annual = rowMeans(.[,c("mort_annual.x","mort_annual.y")]),
         mort_annual = ifelse(sp == "FAH", mort_annual.y, mort_annual),
         mort_annual = ifelse(sp == "SAH", mort_annual.x, mort_annual))

#Calculate the recruitment rates for both periods mean between rates from both periods

r_c_moist <- r_13_moist %>%
  drop_na(sp) %>% 
  filter(sp != "Not_rel") %>% 
  left_join(r_16_moist, by = c("sp")) %>% 
  mutate(rec_year_ha = (rec_year_ha.x + rec_year_ha.y) / 2,
         rec_year_ha = ifelse(sp == "FAH", rec_year_ha.y, rec_year_ha),
         rec_year_ha = ifelse(sp == "SAH", rec_year_ha.x, rec_year_ha))

r_c_middle <- r_13_middle %>% 
  drop_na(sp) %>% 
  filter(sp != "Not_rel") %>%
  left_join(r_16_middle, by = c("sp")) %>% 
  mutate(rec_year_ha = (rec_year_ha.x + rec_year_ha.y) / 2,
         rec_year_ha = ifelse(sp == "FAH", rec_year_ha.y, rec_year_ha),
         rec_year_ha = ifelse(sp == "SAH", rec_year_ha.x, rec_year_ha))

r_c_dry <- r_13_dry %>% 
  drop_na(sp) %>% 
  filter(sp != "Not_rel") %>%
  left_join(r_16_dry, by = c("sp")) %>% 
  mutate(rec_year_ha = (rec_year_ha.x + rec_year_ha.y) / 2,
         rec_year_ha = ifelse(sp == "FAH", rec_year_ha.y, rec_year_ha),
         rec_year_ha = ifelse(sp == "SAH", rec_year_ha.x, rec_year_ha))




#create the combined rates
rates_c_dry = g_c_dry %>% 
  full_join(m_c_dry) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_c_dry %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_c_middle = g_c_middle %>% 
  full_join(m_c_middle) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_c_middle %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_c_moist = g_c_moist %>% 
  full_join(m_c_moist) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_c_moist %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_c_v_wet = g_c_v_wet %>% 
  full_join(m_c_v_wet) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_c_moist %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2) 

rates_c_wet = g_c_wet %>% 
  full_join(m_c_wet) %>%
  pivot_wider(id_cols = "sp", names_from = cl, values_from = c("growth_annual","mort_annual")) %>% 
  full_join(r_c_moist %>% select(sp, rec_year_ha)) %>% filter(sp!="NA"& sp!="Not_rel") %>%
  left_join(sp_table) %>% 
  arrange(sp) %>% 
  select(sp, G1 = growth_annual_1, G2 = growth_annual_2,
         mu1 = mort_annual_1, mu2 = mort_annual_2, rec_ha = rec_year_ha, inflection,
         steepness,param1, param2)

######################## Save the combined rates ###############################

write.table(rates_c_dry, "demographic_rates/rates_combined_census_dry/rates_combined_census_dry.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(rates_c_middle, "demographic_rates/rates_combined_census_middle/rates_combined_census_middle.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(rates_c_moist, "demographic_rates/rates_combined_census_moist/rates_combined_census_moist.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(rates_c_v_wet, "demographic_rates/rates_combined_census_v_wet/rates_combined_census_v_wet.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(rates_c_wet, "demographic_rates/rates_combined_census_wet/rates_combined_census_wet.txt", sep="\t", row.names=FALSE, quote = FALSE)


