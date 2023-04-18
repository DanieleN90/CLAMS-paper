##### this analysis is inspired by Abreu-Vieira et al 2015 components analysis 



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

raised_at_22 <- c("M744", "F745", "F746", "F747", "F749", "F750", "M751", "M754", "M755", "M756", "M757", "F758")
raised_at_30 <- c("M773", "M775", "F781", "M783", "M791", "M793", "M795", "F786", "F787", "F788", "F792", "F794")
dead_mice <- c("F750", "M754", "M757")
other_mice_to_eliminate <- c("F786", "M793", "M795", "F747")


raised_at_22 <- setdiff(setdiff(raised_at_22, dead_mice), other_mice_to_eliminate)
raised_at_30 <- setdiff(setdiff(raised_at_30, dead_mice), other_mice_to_eliminate)
Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]
Nth.delete_vec<-function(vector, n)vector[-(seq(n,to=length(vector),by=n))]


load("~/clean_CLAMS_upd.RData")


# The total energy expenditure (TEE) was fit by second order regression to the physical activity. 
# Taking the Y-axis intercept as the TEE in the absence of physical activity, 
# the energy expenditure of physical activity (PAEE) was defined as: PAEE = TEE Y-intercept.
# Each analysis was performed individually for each mouse, 
# at each Ta and each of the light/dark phases (18 data sets per mouse)


#==============================================================================================#
#let's start with calculating the average TEE for each week of each mouse
# now we need to find a way to have 1 value of heat per count/h
# the problem is that there are four different sampling rates throughout the experiment: every 15, 26, 27, and 30 min
# I will average four consecutive datapoints when the sampling is every 15 min, and two datapoints when it's 26, 27, or 30
# update after 9/22/2021: for when the the interval is 26.25 min, I will take this strategy:
#  - for the heat: remove every 8th value
#  - for the locomotor activity: scale up by 26.25/30= 0.875, then remove every 8th value



for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    time_diff <- as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][2] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][1])
    if(time_diff == 15){
      CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
      day_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_heat <- .colMeans(day_heat, 4, length(day_heat) / 4) 
      night_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_heat <- .colMeans(night_heat, 4, length(night_heat) / 4)
      
      day_act <- CLAMS[[i]][[j]][["Data"]][["XAMB"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_act <- .colSums(day_act, 4, length(day_act) / 4) 
      night_act <- CLAMS[[i]][[j]][["Data"]][["XAMB"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_act <- .colSums(night_act, 4, length(night_act) / 4) 
    }else if(time_diff == 30){
      CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
      day_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_heat <- .colMeans(day_heat, 2, length(day_heat) / 2) 
      night_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_heat <- .colMeans(night_heat, 2, length(night_heat) / 2)
      
      day_act <- CLAMS[[i]][[j]][["Data"]][["XAMB"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_act <- .colSums(day_act, 2, length(day_act) / 2) 
      night_act <- CLAMS[[i]][[j]][["Data"]][["XAMB"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_act <- .colSums(night_act, 2, length(night_act) / 2)
    }else{
      CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
      day_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_heat <- Nth.delete_vec(day_heat, 8)
      day_heat <- .colMeans(day_heat, 2, length(day_heat) / 2) 
      night_heat <- CLAMS[[i]][[j]][["Data"]][["HEAT"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_heat <- Nth.delete_vec(night_heat, 8)
      night_heat <- .colMeans(night_heat, 2, length(night_heat) / 2)
      
      day_act <- CLAMS[[i]][[j]][["Data"]][["XAMB"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      day_act <- day_act/0.875
      day_act <- Nth.delete_vec(day_act, 8)
      day_act <- .colSums(day_act, 2, length(day_act) / 2) 
      night_act <- CLAMS[[i]][[j]][["Data"]][["XAMB"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]]
      night_act <- night_act/0.875
      night_act <- Nth.delete_vec(night_act, 8)
      night_act <- .colSums(night_act, 2, length(night_act) / 2)
    }
    
    
    quadratic.model_day <-lm(day_heat ~ day_act + I(day_act^2))
    quadratic.model_night <-lm(night_heat ~ night_act + I(night_act^2))
    
    TEE_day <- mean(day_heat)
    TEE_night <- mean(night_heat)
    
    PAEE_day <- mean(day_heat) - as.numeric(quadratic.model_day$coefficients[1])
    PAEE_night <- mean(night_heat) - as.numeric(quadratic.model_night$coefficients[1])
    
    
    time_diff <- as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][2] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][1])
    if(time_diff == 26 | time_diff == 27){
      time_diff == 26.25
    }
    hours_in_week <- (time_diff * length(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]]))/60
    
    TEF_day <- (sum(CLAMS[[i]][[j]][["Data"]][["NEW.FEED1"]][!CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]])*0.342)/(hours_in_week)
    TEF_night <- (sum(CLAMS[[i]][[j]][["Data"]][["NEW.FEED1"]][CLAMS[[i]][[j]][["Data"]][["LIGHTS.OFF"]]])*0.342)/(hours_in_week)
    
    rest_of_TEE_day <- TEE_day - PAEE_day - TEF_day
    rest_of_TEE_night <- TEE_night - PAEE_night - TEF_night
    
    
    
    #regression plot for day
    my_data <- as.data.frame(cbind(day_heat, day_act))
    plot1 <-  ggplot(my_data, aes(x = day_act, y = day_heat)) +
      geom_point() +
      ylim(0,0.8) +
      xlim(0,8000) +
      ggtitle("regression plot to compute PAEE during the day", subtitle = paste(names(CLAMS[[i]])[j], " ", names(CLAMS)[i], " ", CLAMS[[i]][[j]][["Group"]])) +
      stat_smooth(method = "lm", col = "red",se=F) +
      stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
    
    #regression plot for the night
    my_data <- as.data.frame(cbind(night_heat, night_act))
    plot2 <-  ggplot(my_data, aes(x = night_act, y = night_heat)) +
      geom_point() +
      ylim(0,0.8) +
      xlim(0,8000) +
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

tot_plots <- vector(mode ="list", length = (9*7))
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
         nrow = 4, ncol = 2, pointsize = 0.2)


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

## so it seems that there is at least 1 mouse from the night_PAEE vector that have negative PAEE, which is fundamentally wrong
## let's see which one are these mice

which(night_PAEE <0)
# 51
51%/%17
#3
51%%17
#0
#so week3, mouse 17 -> M791, week3, 30C
CLAMS[[3]][[17]]
# there are a bunch of high heat values with no locomotor activity that are messing up the whole data,
# so we need to remove at least those values
# I will use dbscan clustering on that dataset to identify the outliers

night_heat <- CLAMS[["Week3"]][["M791"]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]]
night_act <- CLAMS[["Week3"]][["M791"]][["EE Analysis"]][["model night"]][["model"]][["night_act"]]
my_data <- matrix(c(night_heat, night_act), nrow = length(night_act))

res2 <- dbscan(my_data, MinPts = 4, eps = 0.45)

my_data <- as.data.frame(cbind(my_data, res2$cluster))
my_data$V3 <- as.factor(my_data$V3)
plot2 <-  ggplot(my_data, aes(x = night_act, y = night_heat, color=res2$cluster)) +
  geom_point() +
  ylim(0,0.8) +
  xlim(0,4000) +
  ggtitle("regression plot to compute PAEE during the night", subtitle = paste(names(CLAMS[["Week3"]])[j], " ", names(CLAMS)[3], " ", CLAMS[["Week3"]][[j]][["Group"]])) +
  stat_smooth(method = "lm", col = "red",se=F) +
  stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))

plot2

#cluster 1 are the outliers, so I am keeping only those values in cluster 0
night_heat <- night_heat[res2$cluster == 0]
night_act <- night_act[res2$cluster == 0]

#recalculate the model only for the night
quadratic.model_night <-lm(night_heat ~ night_act + I(night_act^2))

TEE_night <- mean(night_heat)

PAEE_night <- mean(night_heat) - as.numeric(quadratic.model_night$coefficients[1])
rest_of_TEE_night <- TEE_night - PAEE_night - TEF_night
my_data <- as.data.frame(cbind(night_heat, night_act))
plot2 <-  ggplot(my_data, aes(x = night_act, y = night_heat)) +
  geom_point() +
  ylim(0,0.8) +
  xlim(0,4000) +
  ggtitle("regression plot to compute PAEE during the night", subtitle = paste(names(CLAMS[["Week3"]])[j], " ", names(CLAMS)[3], " ", CLAMS[["Week3"]][[j]][["Group"]])) +
  stat_smooth(method = "lm", col = "red",se=F) +
  stat_smooth(method = "lm", col = "blue",se=T,formula = y~poly(x,2))
plot2
#replace
CLAMS[["Week3"]][["M791"]][["EE Analysis"]][["model night"]] <- quadratic.model_night
CLAMS[["Week3"]][["M791"]][["EE Analysis"]][["TEE night"]] <- TEE_night
CLAMS[["Week3"]][["M791"]][["EE Analysis"]][["PAEE night"]] <- PAEE_night
CLAMS[["Week3"]][["M791"]][["EE Analysis"]][["rest of TEE night"]] <- rest_of_TEE_night
CLAMS[["Week3"]][["M791"]][["EE Analysis"]][["night regression plot"]] <- plot2

#now looks good, I will remake the plots for 30C night











#ok so far so good, let's save the CLAMS file with the additional stuff
save(CLAMS, file = "clean_CLAMS_upd.RData")

#=============================================================================#

#=============================================================================#
#let's make plots of the components  
# the dataset generate in the following lines will also be used for data analysis


EE_data <- data.frame(matrix(ncol=6, nrow = 0, dimnames = list(NULL, c("Value", "Component", "Group", "Week", "Mouse", "Sex"))))
counter <- 1
for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    for(k in 0:5){
      EE_data[counter + k, 1] <- CLAMS[[i]][[j]][["EE Analysis"]][[5 + k]]
      EE_data[counter + k, 2] <- names(CLAMS[[i]][[j]][["EE Analysis"]])[5 + k]
      EE_data[counter + k, 3] <- CLAMS[[i]][[j]][["Group"]]
      EE_data[counter + k, 4] <- names(CLAMS)[i]
      EE_data[counter + k, 5] <- names(CLAMS[[i]])[j]
      EE_data[counter + k, 6] <- CLAMS[[i]][[j]][["Sex"]]
    }
    counter <- counter + k + 1
  }
}
EE_data <- EE_data[EE_data$Mouse != "F787",]

EE_data$Group_Week <- paste(EE_data$Group, EE_data$Week, sep = "_")
EE_data$Component_Group_Week <- paste(EE_data$Component, EE_data$Group, EE_data$Week, sep = "_")
EE_data$Week_Sex <- paste(EE_data$Sex, EE_data$Week, sep = "_")

##let's divide the dataset into day and night
EE_data_day <- EE_data[EE_data$Component %like% "day",]
EE_data_night <- EE_data[EE_data$Component %like% "night",]




data_for_plot_night <- as.data.frame(group_by(EE_data_night, Component_Group_Week) %>% summarize(m = mean(Value), n=n(), se=sd(Value)/(sqrt(n))))
data_for_plot_day <- as.data.frame(group_by(EE_data_day, Component_Group_Week) %>% summarize(m = mean(Value), n=n(), se=sd(Value)/(sqrt(n))))

data_for_plot_night$Week <- rep(c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), 6)
data_for_plot_night$Component <- c(rep("PAEE", 14), rep("rest of TEE", 14), rep("TEF", 14))
data_for_plot_night$Group <- c(rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7))

data_for_plot_day$Week <- rep(c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), 6)
data_for_plot_day$Component <- c(rep("PAEE", 14), rep("rest of TEE", 14), rep("TEF", 14))
data_for_plot_day$Group <- c(rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7))

#need to make an adjusted dataframe for the error bars
data_for_plot_night_adj <- data_for_plot_night
for(i in 0:13){
  data_for_plot_night_adj[29 + i,2] <- data_for_plot_night[29 + i,2] + data_for_plot_night[15 + i,2]
}
for(i in 0:13){
  data_for_plot_night_adj[1 + i,2] <- data_for_plot_night[1 + i,2] + data_for_plot_night_adj[29 + i,2]
}

data_for_plot_day_adj <- data_for_plot_day
for(i in 0:13){
  data_for_plot_day_adj[29 + i,2] <- data_for_plot_day[29 + i,2] + data_for_plot_day[15 + i,2]
}
for(i in 0:13){
  data_for_plot_day_adj[1 + i,2] <- data_for_plot_day[1 + i,2] + data_for_plot_day_adj[29 + i,2]
}






ggplot(data_for_plot_night, aes(fill=factor(Component, levels = c("PAEE", "TEF", "rest of TEE")), y=m, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin= data_for_plot_night_adj$m - se, 
                    ymax= data_for_plot_night_adj$m + se ),
                width=.2,                    # Width of the error bars
                position="identity") +
  facet_grid(~ Week)+
  ylim(0,0.7) +
  ylab("Mean Weekly heat (KCal/h)") +
  scale_fill_manual(values = c("green3", "firebrick1", "steelblue4"))+ 
  ggtitle("Components of energy expenditure during night time") +
  guides(fill=guide_legend(title="Components"))


ggplot(data_for_plot_day, aes(fill=factor(Component, levels = c("PAEE", "TEF", "rest of TEE")), y=m, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin= data_for_plot_day_adj$m - se, 
                    ymax= data_for_plot_day_adj$m + se ),
                width=.2,                    # Width of the error bars
                position="identity") +
  facet_grid(~ Week)+
  ylim(0,0.7) +
  scale_fill_manual(values = c("green3", "firebrick1", "steelblue4"))+ 
  ggtitle("Components of energy expenditure during day time") +
  ylab("Mean Weekly heat (KCal/h)") +
  guides(fill=guide_legend(title="Components"))

