# this is the script to evaluate the distribution of feeding bouts and to determine what is the reasonable threshold 
# for a feeding event


# we are goin to loop again through the whole dataset, like we did in the first script, 
# but this time we are only doing it to collect all the feeding events

source("~/compensate_negative_values.R")


library(readxl)
library(tidyverse)
library(stringr)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)

#this depends on your file organization
folders <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)[-1]   

all_feeding_events <-c()

for(i in 1:7){
  setwd(folders[i])
  
  #let's count how many samples we have for the current week
  if(dir.exists("./Females")){
    setwd("./Females")
    females_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    setwd("..")
    setwd("./Males")
    males_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    setwd("..")
    mice_number <- females_number + males_number
    #create the list to append all the plots of the week to (we are going to plot 8 parameters per mouse)
    tot_plots <- vector(mode = "list", length = mice_number*4)
    #now let's loop first through the females folder
    setwd("./Females")
    females <- list.files(getwd(),pattern='2018')
    for(j in 1:length(females)){
      my_data <- read.table(females[j], header = T, skip = 22, fill = T, sep = ",")
      my_data <- my_data[-1,]
      my_data <- my_data[-(dim(my_data)[1]),]
      my_data <- select(my_data, INTERVAL, CHAN, DATE.TIME, FEED1, FEED1.ACC, XTOT)
      
      #format the data
      for(o in colnames(my_data)){
        if(o %in% c("INTERVAL", "FEED1", "FEED1.ACC", "XTOT")){
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
      
      # my_data$DATE.TIME[my_data$LIGHTS.OFF]
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
        geom_line(aes(x = INTERVAL, y = XTOT), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("XTOT (counts)") +
        theme_classic()
      tot_plots[[3 + (4*(j-1))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = NEW.FEED1), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Food intake fixed (g)") +
        theme_classic()
      tot_plots[[4 + (4*(j-1))]] <- temp_plot
      

      #store the number of females in a new variable for later
      bb <- j
      all_feeding_events <-c(all_feeding_events, my_data$NEW.FEED1)
    }
    
    # now let's go back to the main folder of the week and loop through the males
    setwd("..")
    setwd("./Males")
    males <- list.files(getwd(),pattern='2018')
    for(j in 1:length(males)){
      my_data <- read.table(males[j], header = T, skip = 22, fill = T, sep = ",")
      my_data <- my_data[-1,]
      my_data <- my_data[-(dim(my_data)[1]),]
      my_data <- select(my_data, INTERVAL, CHAN, DATE.TIME, FEED1, FEED1.ACC, XTOT)
      
      #format the data
      for(o in colnames(my_data)){
        if(o %in% c("INTERVAL", "FEED1", "FEED1.ACC", "XTOT")){
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
      
      # my_data$DATE.TIME[my_data$LIGHTS.OFF]
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
        geom_line(aes(x = INTERVAL, y = XTOT), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("XTOT (counts)") +
        theme_classic()
      tot_plots[[3 + (4*(bb + (j-1)))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = NEW.FEED1), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Food intake fixed (g)") +
        theme_classic()
      tot_plots[[4 + (4*(bb + (j-1)))]] <- temp_plot
      
      all_feeding_events <-c(all_feeding_events, my_data$NEW.FEED1)
      # 
     
    }
    # now we get back to the week folder
    setwd("..")
    #print the plots
    # ggexport(plotlist = tot_plots, filename = paste("Week", i, "plots food fixed.pdf"),
             # nrow = 2, ncol = 2)
    #and now back to the main folder CalR
    setwd("..")
    
  }else{
    
    
    # now the normal loop for folder where males and females are not split
    mice_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    #create the list to append all the plots of the week to (we are going to plot 8 parameters per mouse)
    tot_plots <- vector(mode = "list", length = mice_number*4)
    
    mice <- list.files(getwd(),pattern='2018')
    for(j in 1:length(mice)){
      my_data <- read.table(mice[j], header = T, skip = 22, fill = T, sep = ",")
      my_data <- my_data[-1,]
      my_data <- my_data[-(dim(my_data)[1]),]
      my_data <- select(my_data, INTERVAL, CHAN, DATE.TIME, FEED1, FEED1.ACC, XTOT)
      
      #format the data
      for(o in colnames(my_data)){
        if(o %in% c("INTERVAL", "FEED1", "FEED1.ACC", "XTOT")){
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
      
      # my_data$DATE.TIME[my_data$LIGHTS.OFF]
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
        geom_line(aes(x = INTERVAL, y = XTOT), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("XTOT (counts)") +
        theme_classic()
      tot_plots[[3 + (4*(j-1))]] <- temp_plot
      
      temp_plot <- ggplot(my_data)+
        geom_line(aes(x = INTERVAL, y = NEW.FEED1), size = 0.1)+
        geom_rect(data = dark_ranges, mapping = aes(xmin = from - 1, xmax = to, ymin = ymin, ymax = ymax), alpha = 0.4) +
        xlab("Interval (time)") + 
        ylab("Food intake fixed (g)") +
        theme_classic()
      tot_plots[[4 + (4*(j-1))]] <- temp_plot
      
      all_feeding_events <-c(all_feeding_events, my_data$NEW.FEED1)
      
   
    }
    
    # now print the plots
    # ggexport(plotlist = tot_plots, filename = paste("Week", i, "plots food fixed.pdf"),
    #          nrow = 2, ncol = 2)
    setwd("..")
  }
}



###

dev.off()

#let's make some plots to visualize
hist(all_feeding_events,main="Histogram of observed data")
plot(density(all_feeding_events),main="Density estimate of data")

plot(ecdf(all_feeding_events),main="Empirical cumulative distribution function")

#there are some really high values that do not make much sense, 
#we can put a threshold at 15 to make sense out of the graphs


all_feeding_events2 <- all_feeding_events
all_feeding_events2[all_feeding_events2 > 15] <- 15




hist(all_feeding_events2,main="Histogram of all data points thresholded at 15", breaks = 1000)
plot(density(all_feeding_events2),main="Density estimate of all data points thresholded at 15")

plot(ecdf(all_feeding_events2),main="Empirical cumulative distribution function of all data points thresholded at 15")

quantile(all_feeding_events2, c(.95, .975, .99))

abline(v=quantile(all_feeding_events2, c(.95, .975, .99)), col=c("blue", "red", "black"), lty=c(1,1,1), lwd=c(1,1,1))



## let's remove the zeros

all_feeding_events3 <- all_feeding_events2[all_feeding_events2 != 0]
all_feeding_events3[all_feeding_events3 > 5] <- 5


hist(all_feeding_events3,main="Histogram of only feeding events", breaks = 1000)
plot(density(all_feeding_events3),main="Density estimate of only feeding events")

plot(ecdf(all_feeding_events3),main="Empirical cumulative distribution function of only feeding events")

quantile(all_feeding_events3, c(.95, .975, .98, .99))

abline(v=quantile(all_feeding_events3, c(.95, .975, .98, .99)), col=c("blue", "red", "green", "black"), lty=c(1,1,1,1), lwd=c(1,1,1,1))



## we decided to put a threshold at 97.5% 

threshold <- unname(quantile(all_feeding_events3, .975))


###### now let's redo everything putting a threshold and then redoing the cumulative plots
###### the following code is only executed to print the plots after fixing the feeding values.

########################################################################################
##################     SECOND ROUND OF LOOP, FINALIZING FEED1   ########################
##################     AND FIXING FEED1.ACC ACCORDINGLY         ########################
########################################################################################



folders <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)[-1]   # so it doesn't consider remove 3 days folder 

# all_feeding_events <-c()

for(i in 1:7){
  setwd(folders[i])
  
  #let's count how many samples we have for the current week
  if(dir.exists("./Females")){
    setwd("./Females")
    females_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    setwd("..")
    setwd("./Males")
    males_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    setwd("..")
    mice_number <- females_number + males_number
    #create the list to append all the plots of the week to (we are going to plot 8 parameters per mouse)
    tot_plots <- vector(mode = "list", length = mice_number*4)
    #now let's loop first through the females folder
    setwd("./Females")
    females <- list.files(getwd(),pattern='2018')
    for(j in 1:length(females)){
      my_data <- read.table(females[j], header = T, skip = 22, fill = T, sep = ",")
      my_data <- my_data[-1,]
      my_data <- my_data[-(dim(my_data)[1]),]
      my_data <- select(my_data, INTERVAL, CHAN, DATE.TIME, FEED1, FEED1.ACC, XTOT)
      
      #format the data
      for(o in colnames(my_data)){
        if(o %in% c("INTERVAL", "FEED1", "FEED1.ACC", "XTOT")){
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
      
      
      #store the number of females in a new variable for later
      bb <- j
      all_feeding_events <-c(all_feeding_events, my_data$NEW.FEED1)
    }
    
    # now let's go back to the main folder of the week and loop through the males
    setwd("..")
    setwd("./Males")
    males <- list.files(getwd(),pattern='2018')
    for(j in 1:length(males)){
      my_data <- read.table(males[j], header = T, skip = 22, fill = T, sep = ",")
      my_data <- my_data[-1,]
      my_data <- my_data[-(dim(my_data)[1]),]
      my_data <- select(my_data, INTERVAL, CHAN, DATE.TIME, FEED1, FEED1.ACC, XTOT)
      
      #format the data
      for(o in colnames(my_data)){
        if(o %in% c("INTERVAL", "FEED1", "FEED1.ACC", "XTOT")){
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
      
      # my_data$DATE.TIME[my_data$LIGHTS.OFF]
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
      
      all_feeding_events <-c(all_feeding_events, my_data$NEW.FEED1)
    }
    # now we get back to the week folder
    setwd("..")
    #print the plots
    ggexport(plotlist = tot_plots, filename = paste("Week", i, "plots food fixed final.pdf"),
             nrow = 2, ncol = 2)
    #and now back to the main folder CalR
    setwd("..")
    
  }else{
    
    
    # now the normal loop for folder where males and females are not split
    mice_number <- as.numeric(sapply(getwd(),function(dir){length(list.files(dir,pattern='2018'))}))
    #create the list to append all the plots of the week to (we are going to plot 8 parameters per mouse)
    tot_plots <- vector(mode = "list", length = mice_number*4)
    
    mice <- list.files(getwd(),pattern='2018')
    for(j in 1:length(mice)){
      my_data <- read.table(mice[j], header = T, skip = 22, fill = T, sep = ",")
      my_data <- my_data[-1,]
      my_data <- my_data[-(dim(my_data)[1]),]
      my_data <- select(my_data, INTERVAL, CHAN, DATE.TIME, FEED1, FEED1.ACC, XTOT)
      
      #format the data
      for(o in colnames(my_data)){
        if(o %in% c("INTERVAL", "FEED1", "FEED1.ACC", "XTOT")){
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
      
      all_feeding_events <-c(all_feeding_events, my_data$NEW.FEED1)
      
    }
    
    # now print the plots
    ggexport(plotlist = tot_plots, filename = paste("Week", i, "plots food fixed final.pdf"),
             nrow = 2, ncol = 2)
    setwd("..")
  }
}