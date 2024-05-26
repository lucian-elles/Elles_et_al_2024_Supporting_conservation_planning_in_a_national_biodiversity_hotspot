# Elles et al. (2024). Supporting conservation planning in a national 
# biodiversity hotspot -Projecting species composition across a groundwater level
# gradient using a demographic forest model

# This Script incorporates the PPA Model arranged to be parameterized by the 
# demographic rates from 2016 - 2020 and the Inital state from 2016. This script
# excludes Acer platanoides (SAH). To run this script use the the R scripts 
# "05_run_tidy_current_conditions13_16_and_validation" located in the folder 
# "model_scripts".

main_cohorts_tidy = function(
    mainfolder,       # main folder in which files for this simulation are stored
    spdata_path,      # path to the file with vital rates and allometry
    initdata_path,    # path to the file with initial data, or NA if simulation is to start from bare ground
    growth = "rates", 
    mort = "rates"    # mortality based on demographic rates ("rates") or based on Holzwarth_2012 ("holzwarth")? 
    ) {
 # define basic settings: 
  deltaT = 1         # model timestep
  maxT = 100        # maximum simulation time
  dnot = 0.1         # dbh with which recruits are "born"
  PA = 10000         # simulation area in m² 
  cutT = 1           # timesteps, at which stats shall be recorded
  mincohortN = 0.01  # minimum number of trees in a cohort, before cohort gets removed
  
  spdata = read_table(paste(mainfolder, spdata_path, sep="/"), show_col_types = F) %>% 
    mutate(sp = as.integer(factor(sp)))
  
  if(is.na(initdata_path)) {
    initdata = NA
  } else {
    initdata = as_tibble(read.table(paste(mainfolder, initdata_path, sep="/"), header = F))
  }
  
  # prepare initial data: 
  data = prepare.initdata(initdata, spdata, dnot, deltaT, mincohortN, PA)
  
  # record data to save: 
  all_cohorts_out = data
  stats_out = calc.stats(data,time=0)
  dbh_dist = derive.dbh.distribution(data, time=0)
  
  # MAIN LOOP: .............................................................................
  for(time.tmp in seq(deltaT,maxT,by=deltaT)) {
    
    # Step 1: Mortality
    if(mort == "rates") {
      data = data %>% 
        mutate(n = ifelse(cl == 1, 
                          n * (1 - spdata$mu1[sp]) ^ deltaT, 
                          n * (1 - spdata$mu2[sp]) ^ deltaT)) %>% 
        filter(n > mincohortN)}
    
    if(mort == "holzwarth") {
      data = data %>% 
        mutate(n = n * (1-mortality.holzwarth(dbh))^deltaT) %>% 
        filter(n > mincohortN)
    }
    
    if(mort == "francis.lutz.farrior") {
      data = data %>% 
        mutate(n = ifelse(cl == 2, 
                          n * (1 - spdata$mu2[sp]) ^ deltaT,
                          ifelse(dbh > 70, 
                                 n * (1 - spdata$mu1b[sp]) ^ deltaT, 
                                 n * (1 - spdata$mu1s[sp]) ^ deltaT))) %>% 
        filter(n > mincohortN)
    }
    
    # Step 2: Growth
    if(growth == "rates") {
      data = data %>% 
        mutate(dbh = ifelse(cl == 1, 
                            dbh + spdata$G1[sp] * deltaT, 
                            dbh + spdata$G2[sp] * deltaT))
    }
    if(growth == "modelled") {
      data = data %>% 
        mutate(dbh = dbh + growth.modelled(dbh) * deltaT)
    }
    
    # Step 3: Recruitment
    data = data %>% 
      bind_rows(recruitment(spdata, dnot, deltaT, PA, time.tmp))
    
    # Step 4: Assign canopy layer
    data = assign.layers(data, spdata, PA) %>% 
      filter(cl < 3,         # kill everything higher than layer 3
             n > mincohortN) # filter for mincohortN again
    
    # Step 5: Record stats
    all_cohorts_out = bind_rows(all_cohorts_out, data %>% mutate(time = time.tmp))
    if(time.tmp %% cutT == 0) {
      stats_out = bind_rows(stats_out, calc.stats(data, time.tmp))
      dbh_dist = bind_rows(dbh_dist, derive.dbh.distribution(data, time.tmp))
    }
    
  } #end main loop...........................................................................
  # name trees 
  all_cohorts_out<- all_cohorts_out %>%  mutate(sp_name = case_when(sp == 2 ~ "Acer campestre",
                                                                    sp == 1 ~ "Acer pseudoplatanus",
                                                                    sp == 4 ~ "Carpinus betulus",
                                                                    sp == 3 ~ "Fraxinus excelsior",
                                                                    sp == 5 ~ "Quercus robur",
                                                                    sp == 6 ~ "Ulmus spp.",
                                                                    sp == 7 ~ "Tilia spp.",
                                                                    
                                                                    TRUE ~ as.character(sp)
  ))
  #safe output
  write_csv(all_cohorts_out, paste0(mainfolder, "/all_cohorts_out.csv"))
  write_csv(stats_out, paste0(mainfolder, "/stats_out.csv"))
  write_csv(dbh_dist, paste0(mainfolder, "/dbh_dist_out.csv"))
  
  plot_temp = plot.output(all_cohorts_out, stats_out, dbh_dist, spdata, mainfolder, maxT)
  
  return(list(all_cohorts_out, stats_out, dbh_dist, plot_temp))
}


prepare.initdata = function(initdata, spdata, dnot, deltaT, mincohortN, PA) {
  # create initial data (if NA) or get it in shape: 
  if(is_tibble(initdata)) {
    data = initdata %>% mutate(time=0)
    names(data) = c("dbh", "n", "cl", "sp", "time")
    data = data %>% mutate(sp = as.integer(factor(sp, 
                                                  levels = c("BAH", "FAH","GES", "HBU", "SEI", "UL", "WLI"))))
  } else if(is.na(initdata)) {
    data = recruitment(spdata, dnot, deltaT, PA, time=0)
  }
  
  # add crown area columns and filter for mincohortN: 
  data = data %>% 
    mutate(ca_ind = get.crown.area(dbh, sp, spdata), 
           ca_cohort = ca_ind * n)
  
  # assign canopy layers: 
  data = assign.layers(data, spdata, PA)
  
  # group cohorts, that are now in the same canopy layer: 
  data = data %>% 
    group_by(dbh, cl, sp, time) %>% 
    summarise(n = sum(n), 
              ca_ind = unique(ca_ind), 
              ca_cohort = sum(ca_cohort)) %>% 
    select(dbh, n, cl, sp, time, ca_ind, ca_cohort) %>% 
    ungroup() %>% 
    filter(n >= mincohortN) %>% 
    assign.layers(spdata, PA) #assign canopy layers anew, after filtering for mincohortN
  
  return(data)
}

get.crown.area = function(dbh, sp, spdata) {
  spdata$param1[sp] / (1 + exp(-(spdata$inflection[sp]-dbh)/spdata$steepness[sp])) + spdata$param2[sp] / (1 + exp((spdata$inflection[sp]-dbh)/spdata$steepness[sp]))
}

assign.layers = function(data, spdata, PA) {
  
  # first, update crown area
  data = data %>% 
    mutate(ca_ind = get.crown.area(dbh, sp, spdata), 
           ca_cohort = ca_ind * n) %>% 
    arrange(desc(ca_ind), cl) %>% 
    mutate(cumca = cumsum(ca_cohort), 
           cl = ceiling(cumca/PA)) # get rough cl
  
  layers = max(data$cl)
  
  if(layers > 1) {
    # split data in the rough layers
    dat.tmp = data %>% group_by(cl) %>% group_split()
    
    # split the first cohort in the bottom layer to fill the leftover open canopy space
    for(i in seq(1, layers-1)){
      tosplit = dat.tmp[[i+1]][1,]
      opencan = i*PA - last(dat.tmp[[i]]$cumca)
      split_intop = tosplit %>% mutate(n = (n * opencan) / ca_cohort)
      dat.tmp[[i]] = bind_rows(dat.tmp[[i]], split_intop) %>% 
        mutate(cl = i)
      dat.tmp[[i+1]][1,]$n = dat.tmp[[i+1]][1,]$n - split_intop$n
    }
    
    # bind back together and update crown area
    data = bind_rows(dat.tmp) %>% 
      mutate(ca_ind = get.crown.area(dbh, sp, spdata), 
             ca_cohort = ca_ind * n, 
             cumca = cumsum(ca_cohort))
  }
  
  return(data)
}

mortality.holzwarth = function(dbh) {
  beta_e0 = 1.8
  beta_e1 = -2.1
  beta_e2 = -1.4
  
  beta_l0 = -8.9
  beta_l1 = 0.052
  
  y1 = -3.4
  y2 = 2.1
  y3 = -0.00035
  y4 = 2.5
  
  dinc = exp(y1 + y2*(1-exp((-y3)*(dbh^y4))))
  
  logit = function(x){ 1/(1+exp((-x))) }
  
  mort_e = logit(beta_e0 + beta_e1 * log(dbh+8) + beta_e2*dinc)
  mort_l = logit(beta_l0 + beta_l1 * dbh)
  
  mort = mort_e + mort_l
  final.mort = pmin(mort, 0.02)
  
  return(final.mort)
}

growth.modelled = function(dbh) {
  growth = 0.0824 + 0.004210263*dbh
  return(growth)
}

recruitment = function(spdata, dnot, deltaT, PA, time) {
  # function to generate a tibble with newly recruited trees in one model timestep
  tibble(dbh = dnot, 
         n = spdata$rec_ha * PA / 10000 * deltaT, 
         cl = NA, 
         sp = spdata$sp, 
         time = time)
}

calc.stats = function(data, time) {
  out = data %>% 
    summarise(numtrees = sum(n), 
              ba = sum((dbh/200)^2 * pi * n), 
              timb_vol = NA, 
              maxdbh = max(dbh), 
              densVLT = sum(n[dbh > 80]), 
              totcarea  = sum(ca_cohort), 
              Dstar = min(dbh[cl == 1]), 
              numcohorts = n()) %>% 
    mutate(time = time)
}

derive.dbh.distribution = function(data, time, binwidth=5) {
  dbh_dist = data %>% 
    mutate(dbhcut = cut(dbh, 
                        breaks = seq(0,300,binwidth), 
                        labels = seq(0+binwidth/2, 300-binwidth/2, binwidth))) %>% 
    group_by(dbhcut) %>% 
    summarise(n = sum(n, na.rm = T)) %>% 
    transmute(dbh = as.numeric(as.character(dbhcut)), 
              n, 
              rel_freq = n/sum(n), 
              rel_freq = replace(rel_freq, rel_freq == 0, NA), 
              time = time)
}
# create the Output-plot of the basal area m²/ha for 7 species (excluding Acer platanoides)  for 100 years
plot.output = function(all_cohorts_out, stats_out, dbh_dist, spdata, mainfolder, maxT) {
  
  theme_set(theme_bw())
  theme_update(panel.grid = element_blank(), 
               legend.position = "right")
  
all_cohorts_out$sp <- as.character(all_cohorts_out$sp)
  
  p4 = all_cohorts_out %>%
    filter(cl==1) %>% 
    group_by(time,sp) %>% 
    summarise(ba = sum((dbh/200)^2 * pi * n)) %>%
    group_by(time) %>% 
    mutate(sum_ba = sum(ba)) %>% 
    mutate(sp = case_when(sp == 2 ~ "FAH",
                          sp == 1 ~ "BAH",
                          sp == 4 ~ "HBU",
                          sp == 3 ~ "GES",
                          sp == 5 ~ "SEI",
                          sp == 6 ~ "UL",
                          sp == 7 ~ "WLI",
                          
                          TRUE ~ as.character(sp)
    )) %>% 
    ggplot() + 
    geom_line(aes(x = time, y = ba, color = sp), size = 2) + 
    scale_color_manual(values = c("SAH" = "#228B22",
                                  "BAH" = "#00FF00",
                                  "HBU" = "#FFA500",
                                  "GES" = "blue",
                                  "SEI" = "red",
                                  "UL" = "#17becf",
                                  "WLI" = "grey",
                                  "FAH" = "#808000"))+
    theme(legend.position = "none") +
    ylab("none") + 
    xlab("none") +
    guides(color = guide_legend(override.aes = list(size = 8))) +
    theme(legend.text = element_text(size = 14))+
    theme(
      axis.text = element_text(size = 46, color = "black"),  # Adjust size and style
      axis.title = element_text(size = 1, color = "black", face = "bold"), # Adjust size and style
      axis.ticks.length=unit(.25, "cm"),
      axis.ticks = element_line(size = 2), # Add ticks to the axes  # Adjust size and style
      axis.text.x = element_text(margin = margin(t = 2)),  # Adjust margin for x-axis text
      axis.text.y = element_text(margin = margin(r = 2))  # Adjust margin for y-axis text# Add ticks to the axes  # Adjust size and style
    )+
    labs(color = "Species") +  # Change the legend title here
    ylim(0, 17)+  # Set the y-axis limits from 0 to 15
    xlim(0, 100)  # Set the y-axis limits from 0 to 15

  
  return(p4)
  
}
