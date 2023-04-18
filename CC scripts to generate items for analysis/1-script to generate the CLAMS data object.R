#load up this script to use the function
source("/compensate_negative_values.R")


library(readxl)
library(tidyverse)
library(stringr)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)

raised_at_22 <- c("M744", "F745", "F746", "F747", "F749", "F750", "M751", "M754", "M755", "M756", "M757", "F758")
raised_at_30 <- c("M773", "M775", "F781", "M783", "M791", "M793", "M795", "F786", "F787", "F788", "F792", "F794")

### this line depends on how you organize the files in your folders
### in our case, each folder correspond to one week
folders <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)[-1]   

#this value of threshold has been estimated as the 97.5th percentile of all the feeding events (thresholded at 5g per feeding event)
threshold <- 1.32

for(i in 1:7){
  setwd(folders[i]) #go inside the folder/week
  this_week <- list() #this list will contain each mouse for any given week
  #let's count how many samples we have for the current week
  if(dir.exists("./Females")){ #this loop will execute only if the mice hav been subdivided in males and females
    setwd("./Females")
    females_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    setwd("..")
    setwd("./Males")
    males_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    setwd("..")
    mice_number <- females_number + males_number
    #create the list to append all the plots of the week to (we are going to plot 8 parameters per mouse)
    #now let's loop first through the females folder
    setwd("./Females")
    females <- list.files(getwd(),pattern='2018')
    for(j in 1:length(females)){
      #read in the csv file and keep the relevant data
      my_data <- read.table(females[j], header = T, skip = 22, fill = T, sep = ",")
      my_data <- my_data[-1,]
      my_data <- my_data[-(dim(my_data)[1]),]
      my_data <- select(my_data, INTERVAL, CHAN, DATE.TIME, FEED1, FEED1.ACC, XTOT, RER, HEAT, VO2)
      
      #format the data
      for(o in colnames(my_data)){
        if(o %in% c("INTERVAL", "FEED1", "FEED1.ACC", "XTOT", "RER", "HEAT", "VO2")){
          my_data[,o] <- as.numeric(my_data[,o])
        }
        if(o == "DATE.TIME"){
          my_data$DATE.TIME <- strptime(my_data$DATE.TIME, "%m/%d/%Y %I:%M:%S %p")
        }
      }
      
      #reset rownames
      rownames(my_data) <- 1:nrow(my_data)
      
      #remove NAs
      while (is.na(my_data$DATE.TIME[nrow(my_data)])) {
        my_data <- my_data[-nrow(my_data),]
      }
      
      
      ## check how often is the cage sampled and if it is consistent through the file
      difference <- my_data$DATE.TIME[2] - my_data$DATE.TIME[1]
      for(p in 3:nrow(my_data)){
        difference <- c(difference, my_data$DATE.TIME[p] - my_data$DATE.TIME[p-1])
      }
      var(difference)
      unique(difference)
      # if the sampling is consistent there will be no variance (all values are the same), and unique() will return the only interval value in minutes
      
      # define light time and dark time
      # create check points for each specific exp
      days_in_series <- unique(as.Date(my_data$DATE.TIME))
      lights_off <- c()
      
      
      for (q in 1:nrow(my_data)){
        current_day <- as.Date(my_data$DATE.TIME[q])
        if(my_data$DATE.TIME[q] > strptime(paste(current_day, "07:00:00 EDT"), "%Y-%m-%d %H:%M:%S") & my_data$DATE.TIME[q] < strptime(paste(current_day, "19:00:00 EDT"), "%Y-%m-%d %H:%M:%S")){
          lights_off <- c(lights_off, FALSE)
        }else{
          lights_off <- c(lights_off, TRUE)
        }
      }
      
      my_data$LIGHTS.OFF <- lights_off
      
      
      
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
        if(my_data$DATE.TIME[k] > putative_dark_ranges[m,n]){
          if(n%%2 == 0){
            dark_ranges[m,n] <- my_data$INTERVAL[k-1]
            k <- k + 1
            n <- n - 1
            m <- m + 1
          }else{
            dark_ranges[m,n] <- my_data$INTERVAL[k]
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
      
      
      # fix the feeding values (check the relative script of Material and Methods section for description)
      my_data$NEW.FEED1 <- compensate_negative_values(my_data$FEED1)
      my_data$NEW.FEED1[my_data$NEW.FEED1 >= threshold] <- threshold
      # recalculate cumulative food intake accordingly
      my_data$NEW.FEED.ACC <- cumsum(my_data$NEW.FEED1)
      
      
      #make plots for each variable for visualization
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = FEED1), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Food intake (g)") +
        theme_classic()
      tot_plots[[1 + (4*(j-1))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = FEED1.ACC), size = 0.1) +
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Cumulative food intake (g)") +
        theme_classic()
      tot_plots[[2 + (4*(j-1))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = NEW.FEED1), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Food intake fixed (g)") +
        theme_classic()
      tot_plots[[3 + (4*(j-1))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = NEW.FEED.ACC), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Cumulative food intake fixed (g)") +
        theme_classic()
      tot_plots[[4 + (4*(j-1))]] <- temp_plot
      
      
      
      
      
      
      
      #store the number of females in a new variable for later
      bb <- j
      mouse_id <- read.table(females[j], skip = 6, fill = T, sep = ",")[4,2]
      mouse_weight <- as.numeric(read.table(females[j], skip = 6, fill = T, sep = ",")[5,2])
      mouse_cage <- as.numeric(read.table(females[j], skip = 6, fill = T, sep = ",")[1,2])
      mouse_sex <- substring(mouse_id, 1, 1)
      if(mouse_id %in% raised_at_22){
        mouse_group <- "22C"
      }else{
        mouse_group <- "30C"
      }
      
      #collect all the variable for this mouse in one list
      this_mouse <- list(mouse_id, mouse_weight, mouse_sex, mouse_group, mouse_cage, my_data)
      names(this_mouse) <- c("ID", "Weight", "Sex", "Group", "Cage", "Data")
      #name the list with this mouse name
      assign(mouse_id, this_mouse)
      # add the mouse to the relative week
      this_week[[j]] <- get(ls()[startsWith(ls(), "F")])
      names(this_week[j]) <- mouse_id
      rm(list = c(mouse_id))
      
      
    }
    
    # now let's go back to the main folder of the week and loop through the males
    setwd("..")
    setwd("./Males")
    males <- list.files(getwd(),pattern='2018')
    for(j in 1:length(males)){
      my_data <- read.table(males[j], header = T, skip = 22, fill = T, sep = ",")
      my_data <- my_data[-1,]
      my_data <- my_data[-(dim(my_data)[1]),]
      my_data <- select(my_data, INTERVAL, CHAN, DATE.TIME, FEED1, FEED1.ACC, XTOT, RER, HEAT, VO2)
      
      #format the data
      for(o in colnames(my_data)){
        if(o %in% c("INTERVAL", "FEED1", "FEED1.ACC", "XTOT", "RER", "HEAT", "VO2")){
          my_data[,o] <- as.numeric(my_data[,o])
        }
        if(o == "DATE.TIME"){
          my_data$DATE.TIME <- strptime(my_data$DATE.TIME, "%m/%d/%Y %I:%M:%S %p")
        }
      }
      
      #reset rownames
      rownames(my_data) <- 1:nrow(my_data)
      
      
      
      while (is.na(my_data$DATE.TIME[nrow(my_data)])) {
        my_data <- my_data[-nrow(my_data),]
      }
      
      
      ## check how often is the cage sampled and if it is consistent through the file
      difference <- my_data$DATE.TIME[2] - my_data$DATE.TIME[1]
      for(p in 3:nrow(my_data)){
        difference <- c(difference, my_data$DATE.TIME[p] - my_data$DATE.TIME[p-1])
      }
      var(difference)
      unique(difference)
      # if the sampling is consistent there will be no variance (all values are the same), and unique() will return the only interval value in minutes
      
      # define light time and dark time
      # create check points for each specific exp
      days_in_series <- unique(as.Date(my_data$DATE.TIME))
      lights_off <- c()
      
      
      for (q in 1:nrow(my_data)){
        current_day <- as.Date(my_data$DATE.TIME[q])
        if(my_data$DATE.TIME[q] > strptime(paste(current_day, "07:00:00 EDT"), "%Y-%m-%d %H:%M:%S") & my_data$DATE.TIME[q] < strptime(paste(current_day, "19:00:00 EDT"), "%Y-%m-%d %H:%M:%S")){
          lights_off <- c(lights_off, FALSE)
        }else{
          lights_off <- c(lights_off, TRUE)
        }
      }
      
      my_data$LIGHTS.OFF <- lights_off
      
      # 
      
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
        if(my_data$DATE.TIME[k] > putative_dark_ranges[m,n]){
          if(n%%2 == 0){
            dark_ranges[m,n] <- my_data$INTERVAL[k-1]
            k <- k + 1
            n <- n - 1
            m <- m + 1
          }else{
            dark_ranges[m,n] <- my_data$INTERVAL[k]
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
      
      
      my_data$NEW.FEED1 <- compensate_negative_values(my_data$FEED1)
      my_data$NEW.FEED1[my_data$NEW.FEED1 >= threshold] <- threshold
      my_data$NEW.FEED.ACC <- cumsum(my_data$NEW.FEED1)
      
      
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = FEED1), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Food intake (g)") +
        theme_classic()
      tot_plots[[1 + (4*(bb + (j-1)))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = FEED1.ACC), size = 0.1) +
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Cumulative food intake (g)") +
        theme_classic()
      tot_plots[[2 + (4*(bb + (j-1)))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = NEW.FEED1), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Food intake fixed (g)") +
        theme_classic()
      tot_plots[[3 + (4*(bb + (j-1)))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = NEW.FEED.ACC), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Cumulative food intake fixed (g)") +
        theme_classic()
      tot_plots[[4 + (4*(bb + (j-1)))]] <- temp_plot
      
      # all_feeding_events <-c(all_feeding_events, my_data$NEW.FEED1)
      
      mouse_id <- read.table(males[j], skip = 6, fill = T, sep = ",")[4,2]
      mouse_weight <- as.numeric(read.table(males[j], skip = 6, fill = T, sep = ",")[5,2])
      mouse_cage <- as.numeric(read.table(males[j], skip = 6, fill = T, sep = ",")[1,2])
      mouse_sex <- substring(mouse_id, 1, 1)
      if(mouse_id %in% raised_at_22){
        mouse_group <- "22C"
      }else{
        mouse_group <- "30C"
      }
      
      
      this_mouse <- list(mouse_id, mouse_weight, mouse_sex, mouse_group, mouse_cage, my_data)
      names(this_mouse) <- c("ID", "Weight", "Sex", "Group", "Cage", "Data")
      assign(mouse_id, this_mouse)
      this_week[[bb+j]] <- get(ls()[startsWith(ls(), "M")])
      names(this_week[bb+j]) <- mouse_id
      rm(list = c(mouse_id))
      
      
      
    }
    # now we get back to the week folder
    setwd("..")
    #print the plots
    # ggexport(plotlist = tot_plots, filename = paste("Week", i, "plots food fixed final.pdf"),
    #          nrow = 2, ncol = 2)
    #and now back to the main folder 
    
    
    
    setwd("..")
    
  }else{
    
    
    # now the normal loop for folders where males and females are not split
    mice_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    #create the list to append all the plots of the week to (we are going to plot 8 parameters per mouse)
    tot_plots <- vector(mode = "list", length = mice_number*4)
    
    mice <- list.files(getwd(),pattern='2018')
    for(j in 1:length(mice)){
      my_data <- read.table(mice[j], header = T, skip = 22, fill = T, sep = ",")
      my_data <- my_data[-1,]
      my_data <- my_data[-(dim(my_data)[1]),]
      my_data <- select(my_data, INTERVAL, CHAN, DATE.TIME, FEED1, FEED1.ACC, XTOT, RER, HEAT, VO2)
      
      #format the data
      for(o in colnames(my_data)){
        if(o %in% c("INTERVAL", "FEED1", "FEED1.ACC", "XTOT", "RER", "HEAT", "VO2")){
          my_data[,o] <- as.numeric(my_data[,o])
        }
        if(o == "DATE.TIME"){
          my_data$DATE.TIME <- strptime(my_data$DATE.TIME, "%m/%d/%Y %I:%M:%S %p")
        }
      }
      
      #reset rownames
      rownames(my_data) <- 1:nrow(my_data)
      
      
      
      while (is.na(my_data$DATE.TIME[nrow(my_data)])) {
        my_data <- my_data[-nrow(my_data),]
      }
      
      
      ## check how often is the cage sampled and if it is consistent through the file
      difference <- my_data$DATE.TIME[2] - my_data$DATE.TIME[1]
      for(r in 3:nrow(my_data)){
        difference <- c(difference, my_data$DATE.TIME[r] - my_data$DATE.TIME[r-1])
      }
      var(difference)
      unique(difference)
      # if the sampling is consistent there will be no variance (all values are the same), and unique() will return the only interval value in minutes
      
      # define light time and dark time
      # create check points for each specific exp
      days_in_series <- unique(as.Date(my_data$DATE.TIME))
      lights_off <- c()
      
      
      for (v in 1:nrow(my_data)){
        current_day <- as.Date(my_data$DATE.TIME[v])
        if(my_data$DATE.TIME[v] > strptime(paste(current_day, "07:00:00 EDT"), "%Y-%m-%d %H:%M:%S") & my_data$DATE.TIME[v] < strptime(paste(current_day, "19:00:00 EDT"), "%Y-%m-%d %H:%M:%S")){
          lights_off <- c(lights_off, FALSE)
        }else{
          lights_off <- c(lights_off, TRUE)
        }
      }
      
      my_data$LIGHTS.OFF <- lights_off
      
      
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
        if(my_data$DATE.TIME[k] > putative_dark_ranges[m,n]){
          if(n%%2 == 0){
            dark_ranges[m,n] <- my_data$INTERVAL[k-1]
            k <- k + 1
            n <- n - 1
            m <- m + 1
          }else{
            dark_ranges[m,n] <- my_data$INTERVAL[k]
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
      
      
      
      my_data$NEW.FEED1 <- compensate_negative_values(my_data$FEED1)
      my_data$NEW.FEED1[my_data$NEW.FEED1 >= threshold] <- threshold
      my_data$NEW.FEED.ACC <- cumsum(my_data$NEW.FEED1)
      
      
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = FEED1), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Food intake (g)") +
        theme_classic()
      tot_plots[[1 + (4*(j-1))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = FEED1.ACC), size = 0.1) +
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Cumulative food intake (g)") +
        theme_classic()
      tot_plots[[2 + (4*(j-1))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = NEW.FEED1), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Food intake fixed (g)") +
        theme_classic()
      tot_plots[[3 + (4*(j-1))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = NEW.FEED.ACC), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Cumulative food intake fixed (g)") +
        theme_classic()
      tot_plots[[4 + (4*(j-1))]] <- temp_plot
      
      # all_feeding_events <-c(all_feeding_events, my_data$NEW.FEED1)
      
      
      
      
      
      
      mouse_id <- read.table(mice[j],skip = 6, fill = T, sep = ",")[4,2]
      mouse_weight <- as.numeric(read.table(mice[j], skip = 6, fill = T, sep = ",")[5,2])
      mouse_cage <- as.numeric(read.table(mice[j], skip = 6, fill = T, sep = ",")[1,2])
      mouse_sex <- substring(mouse_id, 1, 1)
      if(mouse_id %in% raised_at_22){
        mouse_group <- "22C"
      }else{
        mouse_group <- "30C"
      }
      
      
      this_mouse <- list(mouse_id, mouse_weight, mouse_sex, mouse_group, mouse_cage, my_data)
      names(this_mouse) <- c("ID", "Weight", "Sex", "Group", "Cage", "Data")
      assign(mouse_id, this_mouse)
      if(mouse_sex == "F"){
        this_week[[j]] <- get(ls()[startsWith(ls(), "F")])
      }else{
        this_week[[j]] <- get(ls()[startsWith(ls(), "M")])
      }
      names(this_week[j]) <- mouse_id
      rm(list = c(mouse_id))
      
      
      
      
      
    }
    
    # now print the plots
    # ggexport(plotlist = tot_plots, filename = paste("Week", i, "plots food fixed final.pdf"),
    #          nrow = 2, ncol = 2)
    setwd("..")
  }
  
  assign(paste("Week", i, sep = "_"), this_week)
  rm(this_week)
  
}

CLAMS <- list(Week_1, Week_2, Week_3, Week_4, Week_5, Week_6, Week_7)

### let's set the names right

for(i in 1:7){
  setwd(folders[i])
  mouse_names <- c()
  if(dir.exists("./Females")){
    setwd("./Females")
    females_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    setwd("..")
    setwd("./Males")
    males_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    setwd("..")
    mice_number <- females_number + males_number
    #now let's loop first through the females folder
    setwd("./Females")
    females <- list.files(getwd(),pattern='2018')
    for(j in 1:length(females)){
      #store the number of females in a new variable for later
      bb <- j
      # all_feeding_events <-c(all_feeding_events, my_data$NEW.FEED1)
      mouse_id <- read.table(females[j], skip = 6, fill = T, sep = ",")[4,2]
      mouse_names <- c(mouse_names, mouse_id)
    }
    # now let's go back to the main folder of the week and loop through the males
    setwd("..")
    setwd("./Males")
    males <- list.files(getwd(),pattern='2018')
    for(j in 1:length(males)){
      mouse_id <- read.table(males[j], skip = 6, fill = T, sep = ",")[4,2]
      mouse_names <- c(mouse_names, mouse_id)
    }
    # now we get back to the week folder
    setwd("..")

    #and now back to the main folder CalR
    setwd("..")
  }else{
    # now the normal loop for folder where males and females are not split
    mice_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    mice <- list.files(getwd(),pattern='2018')
    for(j in 1:length(mice)){
      mouse_id <- read.table(mice[j],skip = 6, fill = T, sep = ",")[4,2]
      mouse_names <- c(mouse_names, mouse_id)
    }
    setwd("..")
  }
  names(CLAMS[[i]]) <- mouse_names
}


names(CLAMS) <- c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7")



save(CLAMS, file = "CLAMS.RData")

#remove all but CLAMS and the IDs of the mice per experiment
rm(list= ls()[!(ls() %in% c("CLAMS", "raised_at_30", "raised_at_22"))])

load("~/CLAMS.RData")


