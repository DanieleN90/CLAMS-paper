setwd("~/HFD challenge CLAMS/")

library(readxl)
library(tidyverse)
library(stringr)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)
library(data.table)
library(lubridate)
library(fpc)
library(readr)

library(lme4)
library(lmerTest)
library(multcomp)
library(MuMIn)


mice_names <- read_xlsx("GTT ITT.xlsx")[,1:2]
mice_names$Mouse

born_at_30 <- mice_names$Mouse[mice_names$...1 == "22C born"]
born_at_22 <- mice_names$Mouse[mice_names$...1 == "30C born"]

load("~/CLAMS_HFD_clean.RData")

#### we need to add the "food_intake" parameter

#==============================================================================================#
#let's start with calculating the average TEE for each week of each mouse
# now we need to find a way to have 1 value of heat per count/h
# the only problem is that the data is sampled every 30 min
# for heat I will average two consecutive value (day and night separately)
# for X.Amb I will sum two consecutive values (same thing as I did for the cold challenge CLAMS)

# the data on food intake during the first week is on chow, so we are going to skip it and go directly to week2


for (i in 2:7) {
  for(j in 1:length(CLAMS_HFD[[i]])){
    time_diff <- as.numeric(CLAMS_HFD[[i]][[j]][["Data"]][["Date.Time"]][2] - CLAMS_HFD[[i]][[j]][["Data"]][["Date.Time"]][1])
    
    CLAMS_HFD[[i]][[j]][["EE Analysis"]] <- NULL
    day_heat <- CLAMS_HFD[[i]][[j]][["Data"]][["Heat"]][CLAMS_HFD[[i]][[j]][["Data"]][["Light.Dark"]] == "Light"]
    day_heat <- .colMeans(day_heat, 2, length(day_heat) / 2) 
    night_heat <- CLAMS_HFD[[i]][[j]][["Data"]][["Heat"]][CLAMS_HFD[[i]][[j]][["Data"]][["Light.Dark"]] == "Dark"]
    night_heat <- .colMeans(night_heat, 2, length(night_heat) / 2)
    
    day_act <- CLAMS_HFD[[i]][[j]][["Data"]][["X.Ambulatory"]][CLAMS_HFD[[i]][[j]][["Data"]][["Light.Dark"]] == "Light"]
    day_act <- .colSums(day_act, 2, length(day_act) / 2) 
    night_act <- CLAMS_HFD[[i]][[j]][["Data"]][["X.Ambulatory"]][CLAMS_HFD[[i]][[j]][["Data"]][["Light.Dark"]] == "Dark"]
    night_act <- .colSums(night_act, 2, length(night_act) / 2)
    
    quadratic.model_day <-lm(day_heat ~ day_act + I(day_act^2))
    quadratic.model_night <-lm(night_heat ~ night_act + I(night_act^2))
    
    TEE_day <- mean(day_heat)
    TEE_night <- mean(night_heat)
    
    PAEE_day <- mean(day_heat) - as.numeric(quadratic.model_day$coefficients[1])
    PAEE_night <- mean(night_heat) - as.numeric(quadratic.model_night$coefficients[1])
    
    
    hours_in_week <- (time_diff * length(CLAMS_HFD[[i]][[j]][["Data"]][["Date.Time"]]))/60
    
    ## the hfd is D12492i from Research Diets Inc.
    ## it is actually the same HFD used in ABreu-Vieira 2015
    ## 5.24 kcal/g     so TEF is  0.419 kcal/g
    ## so for food intake we actually have only one number for each whole week, so we'll just split the food intake equally between day and night
    ## updated from 2/18/22: after looking at data from Marie's experiments with HFD (BIODAQ), we determined that 94.6% of the food is eaten during night time
    
    TEF_day <- ((CLAMS_HFD[[i]][[j]]$food_intake*0.419)*0.054)/(hours_in_week)
    TEF_night <- ((CLAMS_HFD[[i]][[j]]$food_intake*0.419)*0.946)/(hours_in_week)
    
    
    
    rest_of_TEE_day <- TEE_day - PAEE_day - TEF_day
    rest_of_TEE_night <- TEE_night - PAEE_night - TEF_night
    
    
    
    #regression plot for day
    my_data <- as.data.frame(cbind(day_heat, day_act))
    plot1 <-  ggplot(my_data, aes(x = day_act, y = day_heat)) +
      geom_point() +
      ylim(0,0.8) +
      xlim(0,8000) +
      ggtitle("regression plot to compute PAEE during the day", subtitle = paste(names(CLAMS_HFD[[i]])[j], " ", names(CLAMS_HFD)[i], " ", CLAMS_HFD[[i]][[j]][["Group"]])) +
      stat_smooth(method = "lm", col = "red",se=F) +
      stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
    
    #regression plot for the night
    my_data <- as.data.frame(cbind(night_heat, night_act))
    plot2 <-  ggplot(my_data, aes(x = night_act, y = night_heat)) +
      geom_point() +
      ylim(0,0.8) +
      xlim(0,8000) +
      ggtitle("regression plot to compute PAEE during the night", subtitle = paste(names(CLAMS_HFD[[i]])[j], " ", names(CLAMS_HFD)[i], " ", CLAMS_HFD[[i]][[j]][["Group"]])) +
      stat_smooth(method = "lm", col = "red",se=F) +
      stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
    
    
    
    
    EE_analysis <- list(quadratic.model_day, quadratic.model_night, TEE_day, TEE_night, PAEE_night, PAEE_day, TEF_night, TEF_day,rest_of_TEE_night, rest_of_TEE_day, plot1, plot2)
    names(EE_analysis) <- c("model day", "model night", "TEE day", "TEE night", "PAEE night", "PAEE day", "TEF night", "TEF day", "rest of TEE night", "rest of TEE day", "day regression plot", "night regression plot")
    
    EE_analysis <- list(EE_analysis)
    CLAMS_HFD[[i]][[j]] <- c(CLAMS_HFD[[i]][[j]], EE_analysis)
    names(CLAMS_HFD[[i]][[j]])[10] <- "EE Analysis"
  }
}


#### now let's add week 1
for(j in 1:length(CLAMS_HFD[[1]])){
  time_diff <- as.numeric(CLAMS_HFD[[1]][[j]][["Data"]][["Date.Time"]][2] - CLAMS_HFD[[1]][[j]][["Data"]][["Date.Time"]][1])
  
  CLAMS_HFD[[1]][[j]][["EE Analysis"]] <- NULL
  day_heat <- CLAMS_HFD[[1]][[j]][["Data"]][["Heat"]][CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]] == "Light"]
  day_heat <- .colMeans(day_heat, 2, length(day_heat) / 2) 
  night_heat <- CLAMS_HFD[[1]][[j]][["Data"]][["Heat"]][CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]] == "Dark"]
  night_heat <- .colMeans(night_heat, 2, length(night_heat) / 2)
  
  day_act <- CLAMS_HFD[[1]][[j]][["Data"]][["X.Ambulatory"]][CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]] == "Light"]
  day_act <- .colSums(day_act, 2, length(day_act) / 2) 
  night_act <- CLAMS_HFD[[1]][[j]][["Data"]][["X.Ambulatory"]][CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]] == "Dark"]
  night_act <- .colSums(night_act, 2, length(night_act) / 2)
  
  quadratic.model_day <-lm(day_heat ~ day_act + I(day_act^2))
  quadratic.model_night <-lm(night_heat ~ night_act + I(night_act^2))
  
  TEE_day <- mean(day_heat)
  TEE_night <- mean(night_heat)
  
  PAEE_day <- mean(day_heat) - as.numeric(quadratic.model_day$coefficients[1])
  PAEE_night <- mean(night_heat) - as.numeric(quadratic.model_night$coefficients[1])
  
  
  hours_in_week <- (time_diff * length(CLAMS_HFD[[1]][[j]][["Data"]][["Date.Time"]]))/60

  TEF_day <- ((CLAMS_HFD[[1]][[j]]$food_intake*0.342)*0.054)/(hours_in_week)
  TEF_night <- ((CLAMS_HFD[[1]][[j]]$food_intake*0.342)*0.946)/(hours_in_week)
  
  
  
  rest_of_TEE_day <- TEE_day - PAEE_day - TEF_day
  rest_of_TEE_night <- TEE_night - PAEE_night - TEF_night
  
  
  
  #regression plot for day
  my_data <- as.data.frame(cbind(day_heat, day_act))
  plot1 <-  ggplot(my_data, aes(x = day_act, y = day_heat)) +
    geom_point() +
    ylim(0,0.8) +
    xlim(0,8000) +
    ggtitle("regression plot to compute PAEE during the day", subtitle = paste(names(CLAMS_HFD[[1]])[j], " ", names(CLAMS_HFD)[1], " ", CLAMS_HFD[[1]][[j]][["Group"]])) +
    stat_smooth(method = "lm", col = "red",se=F) +
    stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
  
  #regression plot for the night
  my_data <- as.data.frame(cbind(night_heat, night_act))
  plot2 <-  ggplot(my_data, aes(x = night_act, y = night_heat)) +
    geom_point() +
    ylim(0,0.8) +
    xlim(0,8000) +
    ggtitle("regression plot to compute PAEE during the night", subtitle = paste(names(CLAMS_HFD[[1]])[j], " ", names(CLAMS_HFD)[1], " ", CLAMS_HFD[[1]][[j]][["Group"]])) +
    stat_smooth(method = "lm", col = "red",se=F) +
    stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
  
  
  
  
  EE_analysis <- list(quadratic.model_day, quadratic.model_night, TEE_day, TEE_night, PAEE_night, PAEE_day, plot1, plot2)
  names(EE_analysis) <- c("model day", "model night", "TEE day", "TEE night", "PAEE night", "PAEE day", "day regression plot", "night regression plot")
  
  EE_analysis <- list(EE_analysis)
  CLAMS_HFD[[1]][[j]] <- c(CLAMS_HFD[[1]][[j]], EE_analysis)
  names(CLAMS_HFD[[1]][[j]])[10] <- "EE Analysis"
}



save(CLAMS_HFD, file = "CLAMS_HFD_clean.RData")