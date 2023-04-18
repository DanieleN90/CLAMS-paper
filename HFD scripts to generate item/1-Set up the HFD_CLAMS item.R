setwd("~")


library(readxl)
library(tidyverse)
library(stringr)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)

# so the structure is a little different from the cold challenge, so I will have to change a couple of things
# first, let's get the names of the mice and their relative group 

mice_names <- read_xlsx("GTT ITT.xlsx")[,1:2]
mice_names$Mouse

# so apparently this file has inverted names/groups, so the mice under 22C in this files are actually the 30C and viceversa

born_at_30 <- mice_names$Mouse[mice_names$...1 == "22C born"]
born_at_22 <- mice_names$Mouse[mice_names$...1 == "30C born"]
#it seems like the 30C group is missing one mouse from that file
# I will manually add the name of this mouse as I find in 201904-06 CLAMS HFD challenge record.xlsx
# it is F980, and that's because her data is missing after the third week
born_at_30 <- c(born_at_30, "F980")




Weeks <- dir(pattern = "2019.0")
for(i in 1:8){
  setwd(paste(getwd(), Weeks[i], sep = "/"))
  week_data <- read.table(dir(pattern = "Oxymax"), header = T, fill = T, sep = ",")
  week_data <- select(week_data, Subject, Interval, Date.Time, Light.Dark, Volume.O2, RER, Heat, Feed.Weight.1, Feed.Acc..1, X.Total, X.Ambulatory)
  this_week <- list()
  
  #format the data
  for(v in colnames(week_data)){
    if(v %in% c("Interval", "Volume.O2", "RER", "Heat", "Feed.Weight.1", "Feed.Acc..1", "X.Total", "X.Ambulatory")){
      week_data[,v] <- as.numeric(week_data[,v])
    }
    if(v == "Date.Time"){
      week_data$Date.Time <- strptime(week_data$Date.Time, "%m/%d/%Y %I:%M:%S %p")
    }
  }
  for (j in 1:length(unique(week_data$Subject))) {
    my_data <- week_data[week_data$Subject == unique(week_data$Subject)[j],]
    

    ## check how often is the cage sampled and if it is consistent through the file
    difference <- my_data$Date.Time[2] - my_data$Date.Time[1]
    for(b in 3:nrow(my_data)){
      difference <- c(difference, my_data$Date.Time[b] - my_data$Date.Time[b-1])
    }
    var(difference)
    unique(difference)
    # if the sampling is consistent there will be no variance (all values are the same), and unique() will return the only interval value in minutes
    
    # define light time and dark time
    # create check points for each specific exp
    days_in_series <- unique(as.Date(my_data$Date.Time))
    lights_off <- c()
    
    
    for (g in 1:nrow(my_data)){
      current_day <- as.Date(my_data$Date.Time[g])
      if(my_data$Date.Time[g] > strptime(paste(current_day, "07:00:00 EDT"), "%Y-%m-%d %H:%M:%S") & my_data$Date.Time[g] < strptime(paste(current_day, "19:00:00 EDT"), "%Y-%m-%d %H:%M:%S")){
        lights_off <- c(lights_off, FALSE)
      }else{
        lights_off <- c(lights_off, TRUE)
      }
    }
    
    my_datalights_off <- lights_off

    # add extra days at beginning and end to trim later
    putative_dark_ranges <- data.frame(
      from = strptime(c(paste(days_in_series, "19:00:00 EDT")), "%Y-%m-%d %H:%M:%S"),
      to = strptime(c(paste(days_in_series+1, "07:00:00 EDT")), "%Y-%m-%d %H:%M:%S")
    )
    
    
    dark_ranges <- data.frame(from=rep(0,length(days_in_series)), to=rep(0,length(days_in_series)))
    
    k=1
    m=1
    n=1
    while (k <= nrow(my_data)) {
      if(my_data$Date.Time[k] > putative_dark_ranges[m,n]){
        if(n%%2 == 0){
          dark_ranges[m,n] <- my_data$Interval[k-1]
          k <- k + 1
          n <- n - 1
          m <- m + 1
        }else{
          dark_ranges[m,n] <- my_data$Interval[k]
          n <- n + 1
          k <- k + 1
        }
      }else{k <- k + 1}
    }
    
    #trim the last row if the difference is 0
    if(dark_ranges[nrow(dark_ranges),2] - dark_ranges[nrow(dark_ranges),1] == 0){
      dark_ranges <- dark_ranges[-nrow(dark_ranges),]
    }
    
    dark_ranges$ymin <- rep(-Inf, nrow(dark_ranges))
    dark_ranges$ymax<- rep(Inf, nrow(dark_ranges))
    
    
    
    
    mouse_id <- unique(my_data$Subject)
    #take these two parameters from 201904-06 CLAMS HFD challenge record.xlsx
    # I will just put 0 as a place holder, will come back to weight and fat mass later
    mouse_weight <- 0
    mouse_fat_mass <- 0
    mouse_sex <- substring(mouse_id, 1, 1)
    if(mouse_id %in% born_at_22){
      mouse_group <- "22C"
    }else{
      mouse_group <- "30C"
    }
    
    
    this_mouse <- list(mouse_id, mouse_weight, mouse_sex, mouse_group, my_data)
    names(this_mouse) <- c("ID", "Weight", "Sex", "Group", "Data")
    assign(mouse_id, this_mouse)
    if(mouse_sex == "F"){
      this_week[[j]] <- get(ls()[startsWith(ls(), "F")])
    }else{
      this_week[[j]] <- get(ls()[startsWith(ls(), "M")])
    }
    names(this_week[j]) <- mouse_id
    rm(list = c(mouse_id))
    
    
    
  }
  
  assign(paste("Week", i, sep = "_"), this_week)
  rm(this_week)
  
  setwd("..")
}



CLAMS_HFD <- list(Week_1, Week_2, Week_3, Week_4, Week_5, Week_6, Week_7, Week_8)

### let's set the names right

for(i in 1:8){
  setwd(paste(getwd(), Weeks[i], sep = "/"))
  mouse_names <- c()
  for(j in 1:length(CLAMS_HFD[[i]])){
    mouse_id <- CLAMS_HFD[[i]][[j]]$ID
    mouse_names <- c(mouse_names, mouse_id)
  }
  names(CLAMS_HFD[[i]]) <- mouse_names
  setwd("..")
}

names(CLAMS_HFD) <- c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7", "Week8")



save(CLAMS_HFD, file = "CLAMS_HFD.RData")

rm(list= ls()[!(ls() %in% c("CLAMS_HFD", "born_at_22", "born_at_30"))])



