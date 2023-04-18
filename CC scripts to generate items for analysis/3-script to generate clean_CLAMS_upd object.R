
library(readxl)
library(tidyverse)
library(stringr)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)
library(data.table)
library(lubridate)

raised_at_22 <- c("M744", "F745", "F746", "F747", "F749", "F750", "M751", "M754", "M755", "M756", "M757", "F758")
raised_at_30 <- c("M773", "M775", "F781", "M783", "M791", "M793", "M795", "F786", "F787", "F788", "F792", "F794")

#let's load up the object saved from script 1
load("~/CLAMS.RData")


## the goal of this scripts is to visualize problems in the recordings of the locomotor activity and heat
## the analysis that is done hre to calculate locomotor activity was inspired by Abreu-Vieira et al 2015
## however, the analysis here is only a temporary analysis to visualize results, the full components analysis is 
## reported in script 5

# The total energy expenditure (TEE) was fit by second order regression to the physical activity. 
# Taking the Y-axis intercept as the TEE in the absence of physical activity, 
# the energy expenditure of physical activity (PAEE) was defined as: PAEE = TEE - intercept.
# Each analysis was performed individually for each mouse, 
# at each ambient temperature and each of the light/dark phases (18 data sets per mouse)



#==============================================================================================#
#let's start with calculating the average TEE for each week of each mouse
# now we need to find a way to have 1 value of heat per count/h
# the problem is that there are four different sampling rates throughout the experiment: every 15, 26, 27, and 30 min
# I will average four consecutive datapoints when the sampling is every 15 min, and two datapoints when it's 26, 27, or 30

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    time_diff <- as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][2] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][1])
    if(time_diff == 15){
      CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
      day_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_heat <- .colMeans(day_heat, 4, length(day_heat) / 4) 
      night_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_heat <- .colMeans(night_heat, 4, length(night_heat) / 4)
      
      day_act <- CLAMS[[i]][[j]][["Data"]][["XTOT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_act <- .colSums(day_act, 4, length(day_act) / 4) 
      night_act <- CLAMS[[i]][[j]][["Data"]][["XTOT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_act <- .colSums(night_act, 4, length(night_act) / 4) 
    }else{
      CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
      day_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_heat <- .colMeans(day_heat, 2, length(day_heat) / 2) 
      night_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_heat <- .colMeans(night_heat, 2, length(night_heat) / 2)
      
      day_act <- CLAMS[[i]][[j]][["Data"]][["XTOT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_act <- .colSums(day_act, 2, length(day_act) / 2) 
      night_act <- CLAMS[[i]][[j]][["Data"]][["XTOT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_act <- .colSums(night_act, 2, length(night_act) / 2)
    }
    
    
    quadratic.model_day <-lm(day_heat ~ day_act + I(day_act^2))
    quadratic.model_night <-lm(night_heat ~ night_act + I(night_act^2))
    
    TEE_day <- mean(day_heat)
    TEE_night <- mean(night_heat)
    
    PAEE_day <- mean(day_heat) - as.numeric(quadratic.model_day$coefficients[1])
    PAEE_night <- mean(night_heat) - as.numeric(quadratic.model_night$coefficients[1])
    
    time_diff <- as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][2] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][1])
    hours_in_week <- (time_diff * length(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]]))/60
    
    
    TEF_day <- (sum(CLAMS[[i]][[j]][["Data"]][["NEW.FEED1"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]])*0.342)/hours_in_week
    TEF_night <- (sum(CLAMS[[i]][[j]][["Data"]][["NEW.FEED1"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]])*0.342)/hours_in_week

    rest_of_TEE_day <- TEE_day - PAEE_day - TEF_day
    rest_of_TEE_night <- TEE_night - PAEE_night - TEF_night
    
    
    
    #regression plot for day
    my_data <- as.data.frame(cbind(day_heat, day_act))
    plot1 <-  ggplot(my_data, aes(x = day_act, y = day_heat)) +
      geom_point() +
      ylim(0,0.8) +
      xlim(0,4000) +
      ggtitle("regression plot to compute PAEE during the day", subtitle = paste(names(CLAMS[[i]])[j], " ", names(CLAMS)[i], " ", CLAMS[[i]][[j]][["Group"]])) +
      stat_smooth(method = "lm", col = "red",se=F) +
      stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
    
    #regression plot for the night
    my_data <- as.data.frame(cbind(night_heat, night_act))
    plot2 <-  ggplot(my_data, aes(x = night_act, y = night_heat)) +
      geom_point() +
      ylim(0,0.8) +
      xlim(0,4000) +
      ggtitle("regression plot to compute PAEE during the night", subtitle = paste(names(CLAMS[[i]])[j], " ", names(CLAMS)[i], " ", CLAMS[[i]][[j]][["Group"]])) +
      stat_smooth(method = "lm", col = "red",se=F) +
      stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
    
    
    
    if(i == 1 & j>12){
      BMR_day <- TEE_day - PAEE_day - TEF_day
      BMR_night <- TEE_night - PAEE_night - TEF_night
      EE_analysis <- list(quadratic.model_day, quadratic.model_night, TEE_day, TEE_night, PAEE_night, PAEE_day, TEF_night, TEF_day, rest_of_TEE_night, rest_of_TEE_day, BMR_night, BMR_day, plot1, plot2)
      names(EE_analysis) <- c("model day", "model night", "TEE day", "TEE night", "PAEE night", "PAEE day", "TEF night", "TEF day", "rest of TEE night", "rest of TEE day", "BMR night", "BMR day", "day regression plot", "night regression plot")
      
      EE_analysis <- list(EE_analysis)
      CLAMS[[i]][[j]] <- c(CLAMS[[i]][[j]], EE_analysis)
      names(CLAMS[[i]][[j]])[7] <- "EE Analysis"
    }else{
      EE_analysis <- list(quadratic.model_day, quadratic.model_night, TEE_day, TEE_night, PAEE_night, PAEE_day, TEF_night, TEF_day,rest_of_TEE_night, rest_of_TEE_day, plot1, plot2)
      names(EE_analysis) <- c("model day", "model night", "TEE day", "TEE night", "PAEE night", "PAEE day", "TEF night", "TEF day", "rest of TEE night", "rest of TEE day", "day regression plot", "night regression plot")
      
      EE_analysis <- list(EE_analysis)
      CLAMS[[i]][[j]] <- c(CLAMS[[i]][[j]], EE_analysis)
      names(CLAMS[[i]][[j]])[7] <- "EE Analysis"
    }
  }
}


#### let's save the regression plots as a pdf

tot_plots <- vector(mode ="list", length = (12*7))
counter <- 1

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    if(CLAMS[[i]][[j]][["Group"]] == "22C"){
      tot_plots[[counter]] <- CLAMS[[i]][[j]][["EE Analysis"]][["night regression plot"]] 
      counter <- counter +1
    }
  }
  
}


tot_plots[sapply(tot_plots, is.null)] <- NULL

# now print the plots
ggexport(plotlist = tot_plots, filename = "__night time regression plots 22C.pdf",
         nrow = 3, ncol = 2, pointsize = 0.2)


# let's see the distribution of all these PAEE values

night_PAEE <- c()
day_PAEE <- c()
for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    night_PAEE <- c(night_PAEE, CLAMS[[i]][[j]][["EE Analysis"]][["PAEE night"]])
    day_PAEE <- c(day_PAEE, CLAMS[[i]][[j]][["EE Analysis"]][["PAEE day"]])
  }
}

hist(night_PAEE, breaks = 100)
plot(density(night_PAEE))
plot(night_PAEE)
abline(h=0, col = "red")

hist(day_PAEE, breaks = 100)
plot(density(day_PAEE))
plot(day_PAEE)
abline(h=0, col = "red")

## so it seems that there are at least 3 mice from the night_PAEE vector that have negative PAEE, which is fundamentally wrong
## let's see which one are these mice

which(night_PAEE <0)
# 43 144 162
43%/%24
#1
43%%24
#19
#so week2, mouse 19
CLAMS[[2]][[19]]
#week6 mouse 3
#week7 mouse 21




#ok so far so good, let's save the CLAMS file with the additional stuff
save(CLAMS, file = "CLAMS.RData")



#=======================================================================#
#let's try a cleaner analysis by removing outlier mice
# let's keep everything tidy and set up a new folder for this analysis


load("~CLAMS.RData")


#first, let's remove outlier mice
setdiff(names(CLAMS[[5]]), names(CLAMS[[6]]))
#"F750" "M754" "M757"
for (i in 1:5) {
  CLAMS[[i]][["F750"]] <- NULL
  CLAMS[[i]][["M754"]] <- NULL
  CLAMS[[i]][["M757"]] <- NULL
}

## now let's find which ones have 0 food intake due to an error of CLAMS
zero_food <- c()

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    if(sum(CLAMS[[i]][[j]][["Data"]][["NEW.FEED1"]]) == 0){
      zero_food <- c(zero_food, names(CLAMS[[i]])[[j]])
    }
  }
}
#F787 is missing data for food every week, whereas M795 has two weeks missing, but the on other weeks the values are extremely low (<5 for the week)
for (i in 1:7) {
  CLAMS[[i]][["F787"]] <- NULL
  CLAMS[[i]][["M795"]] <- NULL
}

#### now let's see which values are too high (>80)
too_much_food <- c()

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    if(sum(CLAMS[[i]][[j]][["Data"]][["NEW.FEED1"]]) >= 80){
      too_much_food <- c(too_much_food, names(CLAMS[[i]])[[j]])
    }
  }
}

unique(too_much_food)
#"M755" "F781" "F758" "F794" "M783"
#### let's see the course of cumulative food intake of each one of these animals to see if we can exclude just that one week
#### or if we need to completely eliminate these animals from the analysis

cumulative_food_intake_tot <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Animal ID", "food intake", "Week", "Group", "Sex"))))

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    cumulative_food_intake_tot[j + (24*(i-1)), 1] <- CLAMS[[i]][[j]][["ID"]]
    cumulative_food_intake_tot[j + (24*(i-1)), 2] <- tail(CLAMS[[i]][[j]][["Data"]][["NEW.FEED.ACC"]], n=1)
    cumulative_food_intake_tot[j + (24*(i-1)), 3] <- names(CLAMS)[i]
    cumulative_food_intake_tot[j + (24*(i-1)), 4] <- CLAMS[[i]][[j]][["Group"]]
    cumulative_food_intake_tot[j + (24*(i-1)), 5] <- CLAMS[[i]][[j]][["Sex"]]
  }
}

#NAs are due to outlier mice, let's remove them.
# we will also remove the 0s because they are errors
cumulative_food_intake_tot <- cumulative_food_intake_tot[complete.cases(cumulative_food_intake_tot), ]
cumulative_food_intake_tot <- cumulative_food_intake_tot[cumulative_food_intake_tot$food.intake != 0, ]

#let's keep only the mice we are interested in
cumulative_food_intake_tot <- cumulative_food_intake_tot[cumulative_food_intake_tot$Animal.ID %in% too_much_food, ]






temp_plot <- ggplot(cumulative_food_intake_tot, aes(x=Week, y=food.intake, color=Group)) +
  scale_x_discrete()+
  annotate("rect", xmin=1.1, xmax=4.1, ymin=0, ymax=80, alpha=0.2, fill="skyblue3") +
  annotate("rect", xmin=4.1, xmax=7, ymin=0, ymax=80, alpha=0.3, fill="skyblue4") +
  geom_line(aes(group = Animal.ID), alpha = 0.4)+
  geom_point(position = position_jitter(width = .1, height = 0), aes(size=Sex)) + 
  # geom_bernie() +
  stat_summary(aes(y = food.intake, group = Group), fun = mean, geom="line", size = 2) +
  xlab("Weeks") + 
  ylab("Cumulative food intake (g)") +
  ggtitle("Total cumulative food intake")+
  ylim(0,80)+
  annotate(geom="text", x=3, y=75, label="22C", color="black", size = 6) +
  annotate(geom="text", x=5.5, y=75, label="14C", color="black", size = 6) +
  scale_color_manual(values=c("steelblue", "darksalmon"))+
  stat_compare_means(aes(group = Group), label = "p.format", size = 3)+
  theme_bw()
temp_plot

CLAMS[[5]][["M755"]] <- NULL
CLAMS[[6]][["F781"]] <- NULL
CLAMS[[7]][["M783"]] <- NULL
CLAMS[[7]][["F758"]] <- NULL
CLAMS[[7]][["F794"]] <- NULL


### now there are still some mice that have a very low food intake on some weeks, let's see
too_little_food <- c()

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    if(sum(CLAMS[[i]][[j]][["Data"]][["NEW.FEED1"]]) < 12){
      too_little_food <- c(too_little_food, names(CLAMS[[i]])[[j]])
    }
  }
}

plot(ts(CLAMS[[4]][["F792"]]$Data$NEW.FEED.ACC))
plot(ts(CLAMS[[3]][["F792"]]$Data$NEW.FEED.ACC))
plot(ts(CLAMS[[5]][["F792"]]$Data$NEW.FEED.ACC))

## so this mouse seems to be normal on other weeks, I will just remove the data from week4
CLAMS[[4]][["F792"]] <- NULL



#### let's see now how does the heat look like

HEAT <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Animal ID", "HEAT", "Week", "Group", "Sex"))))

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    HEAT[j + (24*(i-1)), 1] <- CLAMS[[i]][[j]][["ID"]]
    HEAT[j + (24*(i-1)), 2] <- mean(CLAMS[[i]][[j]][["Data"]][["HEAT"]])
    HEAT[j + (24*(i-1)), 3] <- names(CLAMS)[i]
    HEAT[j + (24*(i-1)), 4] <- CLAMS[[i]][[j]][["Group"]]
    HEAT[j + (24*(i-1)), 5] <- CLAMS[[i]][[j]][["Sex"]]
  }
}


HEAT <- HEAT[complete.cases(HEAT), ]
HEAT <- HEAT[HEAT$HEAT != 0, ]


temp_plot <- ggplot(HEAT, aes(x=Week, y=HEAT, color=Group)) +
  scale_x_discrete()+
  annotate("rect", xmin=1.1, xmax=4.1, ymin=0, ymax=0.8, alpha=0.2, fill="skyblue3") +
  annotate("rect", xmin=4.1, xmax=7, ymin=0, ymax=0.8, alpha=0.3, fill="skyblue4") +
  geom_line(aes(group = Animal.ID), alpha = 0.4)+
  stat_summary(aes(y = HEAT, group = Group), fun = median, geom="line", size = 2) +
  geom_point(position = position_jitter(width = .1, height = 0), aes(size=Sex)) + 
  ylim(0,0.8) +
  xlab("Weeks") + 
  ylab("HEAT (kcal/hours)") +
  ggtitle("Mean Weekly HEAT")+
  scale_color_manual(values=c("steelblue", "darksalmon"))+
  annotate(geom="text", x=3, y=0.75, label="22C", color="black", size = 6) +
  annotate(geom="text", x=5.5, y=0.75, label="14C", color="black", size = 6) +
  stat_compare_means(aes(group = Group), label = "p.format", size = 3)+
  theme_bw()
temp_plot

#it seems like only one datapoint is clearly out of space
too_cold <- c()

for(j in 1:length(CLAMS[[6]])){
  if(mean(CLAMS[[6]][[j]][["Data"]][["HEAT"]]) < 0.4){
    too_cold <- c(too_cold, names(CLAMS[[6]])[[j]])
  }
}

plot(ts(CLAMS[[6]][["M744"]]$Data$HEAT))

CLAMS[[6]][["M744"]] <- NULL




#============================================================#
#ok now let's try again to redo the analysis
#I will save this version of CLAMS as clean_CLAMS
save(CLAMS, file = "clean_CLAMS_upd.RData")
load("/clean_CLAMS_upd.RData")


for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    time_diff <- as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][2] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][1])
    if(time_diff == 15){
      CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
      day_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_heat <- .colMeans(day_heat, 4, length(day_heat) / 4) 
      night_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_heat <- .colMeans(night_heat, 4, length(night_heat) / 4)
      
      day_act <- CLAMS[[i]][[j]][["Data"]][["XTOT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_act <- .colSums(day_act, 4, length(day_act) / 4) 
      night_act <- CLAMS[[i]][[j]][["Data"]][["XTOT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_act <- .colSums(night_act, 4, length(night_act) / 4) 
    }else{
      CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
      day_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_heat <- .colMeans(day_heat, 2, length(day_heat) / 2) 
      night_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_heat <- .colMeans(night_heat, 2, length(night_heat) / 2)
      
      day_act <- CLAMS[[i]][[j]][["Data"]][["XTOT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_act <- .colSums(day_act, 2, length(day_act) / 2) 
      night_act <- CLAMS[[i]][[j]][["Data"]][["XTOT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_act <- .colSums(night_act, 2, length(night_act) / 2)
    }
    
    
    quadratic.model_day <-rlm(day_heat ~ day_act + I(day_act^2))
    quadratic.model_night <-rlm(night_heat ~ night_act + I(night_act^2))
    
    # quadratic.model_night_test <-lm(night_heat ~ poly(night_act, 2, raw = T))
    
    TEE_day <- mean(day_heat)
    TEE_night <- mean(night_heat)
    
    PAEE_day <- mean(day_heat) - as.numeric(quadratic.model_day$coefficients[1])
    PAEE_night <- mean(night_heat) - as.numeric(quadratic.model_night$coefficients[1])
    
    time_diff <- as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][2] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][1])
    hours_in_week <- (time_diff * length(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]]))/60
    
    
    TEF_day <- (sum(CLAMS[[i]][[j]][["Data"]][["NEW.FEED1"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]])*0.342)/hours_in_week
    TEF_night <- (sum(CLAMS[[i]][[j]][["Data"]][["NEW.FEED1"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]])*0.342)/hours_in_week
    
    rest_of_TEE_day <- TEE_day - PAEE_day - TEF_day
    rest_of_TEE_night <- TEE_night - PAEE_night - TEF_night
    
    
    
    #regression plot for day
    my_data <- as.data.frame(cbind(day_heat, day_act))
    plot1 <-  ggplot(my_data, aes(x = day_act, y = day_heat)) +
      geom_point() +
      ylim(0,0.8) +
      xlim(0,10000) +
      ggtitle("regression plot to compute PAEE during the day", subtitle = paste(names(CLAMS[[i]])[j], " ", names(CLAMS)[i], " ", CLAMS[[i]][[j]][["Group"]])) +
      stat_smooth(method = "lm", col = "red",se=F) +
      stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
    
    #regression plot for the night
    my_data <- as.data.frame(cbind(night_heat, night_act))
    plot2 <-  ggplot(my_data, aes(x = night_act, y = night_heat)) +
      geom_point() +
      ylim(0,0.8) +
      xlim(0,10000) +
      ggtitle("regression plot to compute PAEE during the night", subtitle = paste(names(CLAMS[[i]])[j], " ", names(CLAMS)[i], " ", CLAMS[[i]][[j]][["Group"]])) +
      stat_smooth(method = "lm", col = "red",se=F) +
      stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
    
    
    
    if(i == 1 & j>12){
      BMR_day <- TEE_day - PAEE_day - TEF_day
      BMR_night <- TEE_night - PAEE_night - TEF_night
      EE_analysis <- list(quadratic.model_day, quadratic.model_night, TEE_day, TEE_night, PAEE_night, PAEE_day, TEF_night, TEF_day, rest_of_TEE_night, rest_of_TEE_day, BMR_night, BMR_day, plot1, plot2)
      names(EE_analysis) <- c("model day", "model night", "TEE day", "TEE night", "PAEE night", "PAEE day", "TEF night", "TEF day", "rest of TEE night", "rest of TEE day", "BMR night", "BMR day", "day regression plot", "night regression plot")
      
      EE_analysis <- list(EE_analysis)
      CLAMS[[i]][[j]] <- c(CLAMS[[i]][[j]], EE_analysis)
      names(CLAMS[[i]][[j]])[7] <- "EE Analysis"
    }else{
      EE_analysis <- list(quadratic.model_day, quadratic.model_night, TEE_day, TEE_night, PAEE_night, PAEE_day, TEF_night, TEF_day,rest_of_TEE_night, rest_of_TEE_day, plot1, plot2)
      names(EE_analysis) <- c("model day", "model night", "TEE day", "TEE night", "PAEE night", "PAEE day", "TEF night", "TEF day", "rest of TEE night", "rest of TEE day", "day regression plot", "night regression plot")
      
      EE_analysis <- list(EE_analysis)
      CLAMS[[i]][[j]] <- c(CLAMS[[i]][[j]], EE_analysis)
      names(CLAMS[[i]][[j]])[7] <- "EE Analysis"
    }
  }
}


#### let's save the regression plots as a pdf

tot_plots <- vector(mode ="list", length = (12*7))
counter <- 1

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    if(CLAMS[[i]][[j]][["Group"]] == "30C"){
      tot_plots[[counter]] <- CLAMS[[i]][[j]][["EE Analysis"]][["day regression plot"]] 
      counter <- counter +1
    }
  }
  
}


tot_plots[sapply(tot_plots, is.null)] <- NULL

# now print the plots
ggexport(plotlist = tot_plots, filename = "day time regression plots 30C.pdf",
         nrow = 3, ncol = 2, pointsize = 0.2)


# let's see the distribution of all these PAEE values

night_PAEE <- c()
day_PAEE <- c()
for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    night_PAEE <- c(night_PAEE, CLAMS[[i]][[j]][["EE Analysis"]][["PAEE night"]])
    day_PAEE <- c(day_PAEE, CLAMS[[i]][[j]][["EE Analysis"]][["PAEE day"]])
  }
}

hist(night_PAEE, breaks = 100)
plot(density(night_PAEE))
plot(night_PAEE)
abline(h=0, col = "red")

hist(day_PAEE, breaks = 100)
plot(density(day_PAEE))
plot(day_PAEE)
abline(h=0, col = "red")

## so it seems that there one mouse from the night_PAEE vector that have negative PAEE, which is fundamentally wrong
## let's see which one is

which(night_PAEE <0)
# 96

#so week6, mouse 3 (F747)



CLAMS[[6]][[3]] <- NULL





#ok so far so good, let's save the CLAMS file with the additional stuff
save(CLAMS, file = "clean_CLAMS_upd.RData")


