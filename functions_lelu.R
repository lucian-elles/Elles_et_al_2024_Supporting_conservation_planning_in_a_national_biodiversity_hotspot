assign.crown.layer.in.single.census = function(df, nr.of.layers){
  
  df = df %>% 
    # remove all dead, prior and non-located stems
    filter(!is.na(ca) & 
             !is.na(gx) & 
             !is.na(gy) &
             !is.na(plotID) &
             !is.na(dbh)) %>% 
    # correct crown area, individuals with crown area > subplot area get the value of subplot area
    mutate(ca_corrected = replace(ca, ca > 625, 625)) %>% 
    # order trees by height
    arrange(plotID, desc(height)) %>% 
    # assign the individual trees to crown layers
    group_by(plotID) %>%
    mutate(cl = helper.function.assign.trees.to.crown.layers(crown.cover.corrected = ca_corrected,
                                                             subplot.area = plot_area[1],
                                                             allometry = "ppa")) %>% 
    ungroup()
  
  df.sum = df %>% 
    group_by(plotID) %>% 
    summarise(max_cl = max(cl), 
              l1 = sum(cl == 1), 
              l2 = sum(cl == 2), 
              l3plus = sum(cl >2))
  
  # change the maximum number of crown levels (i.e. change all higher values to the maximum)  
  df <- df %>% 
    mutate(cl_true = cl, 
           cl = replace(cl, cl > nr.of.layers, nr.of.layers))
  
  # print number of crown layers
  print(paste("Median number of possible crown layers =", median(df.sum$max_cl)))
  print(paste(c("Range of possible crown layers =", range(df.sum$max_cl)), collapse = " "))
  print(paste("Nr. of live trees total:", nrow(df)))
  print(paste("Nr. of trees in layers 3 and higher:", sum(df.sum$l3plus)))
  
  crownlayers <- df %>% 
    select(treeID, ca_corrected, cl)
  
  # return to global environment
  return(crownlayers)
}


helper.function.assign.trees.to.crown.layers = function(crown.cover.corrected, 
                                                        subplot.area, 
                                                        allometry){
  # arguments:
  # crown.cover.corrected must be a numeric vector that is ordered from high to low values
  # subplot area = area of the subplot used for PPA layer assignment
  # allometry = on which basis should trees be assigned to crown layers:
  # "ppa" = based on crown - dbh allometry
  # "even" = evenly assignment the same number of individuals to every group
  # nr.of.layers = how many layers should be assigned, will only be used for allometry = "even"
  
  crown.layer.assigned = rep(NA, length(crown.cover.corrected))
  
  temp.cover.level = 0
  
  if(allometry == "ppa"){
    while(length(which(is.na(crown.layer.assigned))) > 0){
      
      temp.cover.level = temp.cover.level + 1
      
      # calculate the cummulative sum of all individuals that have not yet been assigned to a crown layer  
      cumsum.temp = rep(NA, length(crown.cover.corrected))
      cumsum.temp[which(is.na(crown.layer.assigned))] = cumsum(crown.cover.corrected[which(is.na(crown.layer.assigned))])
      #cat("which(is.na(crown.layer.assigned)) : ", which(is.na(crown.layer.assigned)))
      
      # assign those individuals whose cumsum is < plot area to the same crown layer
      if(length(which(cumsum.temp < subplot.area)) > 0){
        crown.layer.assigned[which(cumsum.temp < subplot.area)] = temp.cover.level
      }else{
        crown.layer.assigned[min(which(is.na(crown.layer.assigned)))] = temp.cover.level}
      
      # stop the loop if all idividuals are already assigned to a crown layer
      if(length(which(is.na(crown.layer.assigned))) == 0){break}
      
      # if more than 50% of the next individual would fit within the current crown layer assign it to this layer
      plot.area.minus.area.covered.by.temp.layer = subplot.area - cumsum.temp[max(which(crown.layer.assigned == temp.cover.level))]
      if((plot.area.minus.area.covered.by.temp.layer / crown.cover.corrected[min(which(is.na(crown.layer.assigned)))]) >= 0.5){
        crown.layer.assigned[min(which(is.na(crown.layer.assigned)))] = temp.cover.level}
    }}
  
  if(allometry == "even"){
    crown.layer.assigned = rep(1:nr.of.layers, each =  round(length(crown.cover.corrected) / nr.of.layers))
    crown.layer.assigned = crown.layer.assigned[1:length(crown.cover.corrected)]
  }
  return(crown.layer.assigned)
}

adjust.recruitment.rate = function(r, n_years, mu){
  ### function to adjust the recruitment rate including unobserved recruits due to mortality ###
  # r = annual recruitment rate
  # n_years = timespan between censusses in years
  # mu = mortality rate
  
  div = 0
  i = 0
  while(n_years-i != 0){
    div = div + (1-mu)^(n_years-i)
    i = i+1
  }
  Fadj = (r*n_years) / div
  return(Fadj)
}
