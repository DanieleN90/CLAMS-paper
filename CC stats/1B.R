#this script is to perform the analysis in Supplementary File 1B 

### to do so, we first need to consider how long is an interval, and then we can take the last or first n values

### this is not so straight forward because the length of an interval changes across the experiment because the
### number of mice change from week to week. However, we have this info under CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]]

setwd("~")

library(readxl)
library(tidyverse)
library(stringr)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)
library(lubridate)

raised_at_22 <- c("M744", "F745", "F746", "F747", "F749", "F750", "M751", "M754", "M755", "M756", "M757", "F758")
raised_at_30 <- c("M773", "M775", "F781", "M783", "M791", "M793", "M795", "F786", "F787", "F788", "F792", "F794")
dead_mice <- c("F750", "M754", "M757")
other_mice_to_eliminate <- c("F786", "M793", "M795", "F747")


raised_at_22 <- setdiff(setdiff(raised_at_22, dead_mice), other_mice_to_eliminate)
raised_at_30 <- setdiff(setdiff(raised_at_30, dead_mice), other_mice_to_eliminate)
Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]
Nth.delete_vec<-function(vector, n)vector[-(seq(n,to=length(vector),by=n))]
range01 <- function(x){(x-min(x))/(max(x)-min(x))}



load("~/CLAMS.RData")

### first we need to create a csv file for each hourly average of each week (week 1, 2, 4, and 5)
### this will help with downstream analysis 

#now let's print the hourly averages for plotting
# basically we will take the data from above and do an average every 2 consecutive intervals
# for the locomotor activity we will do the sum over two intervals, for the heat we will do the average
# for each value, we will assign the timestamp of the first of the two intervals involved

dff <- data.frame(matrix(ncol = 7, nrow = 0))
x <- c("Interval", "Mouse", "Heat", "Group", "Sex", "Time", "Lights_off")
colnames(dff) <- x



dff <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(dff) <- x
for(j in 1:length(CLAMS[[1]])){
  Interval <- CLAMS[[1]][[j]][["Data"]][["INTERVAL"]]
  enn <- length(CLAMS[[1]][[j]][["Data"]][["INTERVAL"]])
  Mouse <- rep(CLAMS[[1]][[j]][["ID"]], enn)
  Heat <- CLAMS[[1]][[j]][["Data"]][["HEAT"]]
  Time <- CLAMS[[1]][[j]][["Data"]][["DATE.TIME"]]
  Sex <- rep(CLAMS[[1]][[j]]$Sex, enn)
  Lights_off <- CLAMS[[1]][[j]][["Data"]][["LIGHTS.OFF"]]
  Group <- rep(CLAMS[[1]][[j]]$Group, enn)
  dffff <- data.frame(Interval,
                      Mouse,
                      Heat,
                      Group,
                      Sex,
                      Time,
                      Lights_off)
  if(dffff$Group[1] == "22C"){
    dffff <- dffff[seq(1, length(dffff[,1]), 2),]
    dffff <- dffff[14:334,]
    dffff$Interval <- seq(1, 321, 1)
  }else{
    dffff <- dffff[1:321,]
  }
  mean_heat <- colMeans(matrix(dffff$Heat, nrow=2))
  dffff <- dffff[seq(1, length(dffff[,1]), 2),]
  dffff$Heat <- mean_heat
  dffff$Interval <- seq(1, 161, 1)
  dff <- rbind(dff, dffff)
  colnames(dff) <- x
}
write.csv(dff,paste0("~", "Week 1 Heat Hourly Average.csv"), row.names = FALSE)

## now week 2

dff <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(dff) <- x
for(j in 1:length(CLAMS[[2]])){
  Interval <- CLAMS[[2]][[j]][["Data"]][["INTERVAL"]]
  enn <- length(CLAMS[[2]][[j]][["Data"]][["INTERVAL"]])
  Mouse <- rep(CLAMS[[2]][[j]][["ID"]], enn)
  Heat <- CLAMS[[2]][[j]][["Data"]][["HEAT"]]
  Time <- CLAMS[[2]][[j]][["Data"]][["DATE.TIME"]]
  Sex <- rep(CLAMS[[2]][[j]]$Sex, enn)
  Lights_off <- CLAMS[[2]][[j]][["Data"]][["LIGHTS.OFF"]]
  Group <- rep(CLAMS[[2]][[j]]$Group, enn)
  dffff <- data.frame(Interval,
                      Mouse,
                      Heat,
                      Group,
                      Sex,
                      Time,
                      Lights_off)
  if(dffff$Group[1] == "22C"){
    dffff <- dffff[seq(1, length(dffff[,1]), 2),]
    dffff <- dffff[1:332,] ## 332 is the shortest length among the mice
    dffff$Interval <- seq(1, 332, 1)
  }else{
    dffff <- dffff[1:332,]
  }
  mean_heat <- colMeans(matrix(dffff$Heat, nrow=2))
  dffff <- dffff[seq(1, length(dffff[,1]), 2),]
  dffff$Heat <- mean_heat
  dffff$Interval <- seq(1, 166, 1)
  dff <- rbind(dff, dffff)
  colnames(dff) <- x
}
write.csv(dff,paste0("~", "Week 2 Heat Hourly Average.csv"), row.names = FALSE)



## now week 4

dff <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(dff) <- x
for(j in 1:length(CLAMS[[4]])){
  Interval <- CLAMS[[4]][[j]][["Data"]][["INTERVAL"]]
  enn <- length(CLAMS[[4]][[j]][["Data"]][["INTERVAL"]])
  Mouse <- rep(CLAMS[[4]][[j]][["ID"]], enn)
  Heat <- CLAMS[[4]][[j]][["Data"]][["HEAT"]]
  Time <- CLAMS[[4]][[j]][["Data"]][["DATE.TIME"]]
  Sex <- rep(CLAMS[[4]][[j]]$Sex, enn)
  Lights_off <- CLAMS[[4]][[j]][["Data"]][["LIGHTS.OFF"]]
  Group <- rep(CLAMS[[4]][[j]]$Group, enn)
  dffff <- data.frame(Interval,
                      Mouse,
                      Heat,
                      Group,
                      Sex,
                      Time,
                      Lights_off)
  if(dffff$Group[1] == "22C"){
    dffff <- dffff[1:330,]
  }else{
    dffff <- Nth.delete(dffff, 8)
    dffff$Interval <- seq(1, length(dffff$Interval), 1)
    dffff <- dffff[1:330,]
  }
  mean_heat <- colMeans(matrix(dffff$Heat, nrow=2))
  dffff <- dffff[seq(1, length(dffff[,1]), 2),]
  dffff$Heat <- mean_heat
  dffff$Interval <- seq(1,330, 2)
  dff <- rbind(dff, dffff)
  colnames(dff) <- x
}
write.csv(dff,paste0("~", "Week 4 Heat Hourly Average.csv"), row.names = FALSE)



## now week 5

dff <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(dff) <- x
for(j in 1:length(CLAMS[[5]])){
  Interval <- CLAMS[[5]][[j]][["Data"]][["INTERVAL"]]
  enn <- length(CLAMS[[5]][[j]][["Data"]][["INTERVAL"]])
  Mouse <- rep(CLAMS[[5]][[j]][["ID"]], enn)
  Heat <- CLAMS[[5]][[j]][["Data"]][["HEAT"]]
  Time <- CLAMS[[5]][[j]][["Data"]][["DATE.TIME"]]
  Sex <- rep(CLAMS[[5]][[j]]$Sex, enn)
  Lights_off <- CLAMS[[5]][[j]][["Data"]][["LIGHTS.OFF"]]
  Group <- rep(CLAMS[[5]][[j]]$Group, enn)
  dffff <- data.frame(Interval,
                      Mouse,
                      Heat,
                      Group,
                      Sex,
                      Time,
                      Lights_off)
  if(dffff$Group[1] == "30C"){
    dffff <- Nth.delete(dffff, 8)
    dffff$Interval <- seq(1, length(dffff$Interval), 1)
  }
  mean_heat <- colMeans(matrix(dffff$Heat, nrow=2))
  dffff <- dffff[seq(1, length(dffff[,1]), 2),]
  dffff$Heat <- mean_heat
  dffff$Interval <- seq(1,165, 1)
  
  dff <- rbind(dff, dffff)
  colnames(dff) <- x
}
write.csv(dff,paste0("~", "Week 5 Heat Hourly Average.csv"), row.names = FALSE)

###########################################################################################################################
#################### now we can proceed with the rest of the analysis
load("~/clean_CLAMS_upd.RData")


### this is the to look at the 1h time intervals
week_1_heat <- read_csv("Week 1 Heat Hourly Average.csv")
week_2_heat <- read_csv("Week 2 Heat Hourly Average.csv")
week_4_heat <- read_csv("Week 4 Heat Hourly Average.csv")
week_5_heat <- read_csv("Week 5 Heat Hourly Average.csv")
week_1_heat$Week <- "Week1" 
week_2_heat$Week <- "Week2"
week_4_heat$Week <- "Week4" 
week_5_heat$Week <- "Week5"

# let's look at  end of week 1 vs start of week 2 (1B i and iii)


my_dataframe <- data.frame(matrix(ncol=5,nrow=0, 
                                  dimnames=list(NULL, c("Animal ID", "Heat", "Week", "Group", "Sex"))))
for(z in unique(week_2_heat$Mouse)){
  my_vector1 <- c(z,
                  mean(tail(week_1_heat[week_1_heat$Mouse == z,], n=1)$Heat),
                  "End of Week 1",
                  unique(week_1_heat$Group[week_1_heat$Mouse == z]),
                  unique(week_1_heat$Sex[week_1_heat$Mouse == z]))
  my_vector2 <- c(z,
                  mean(head(week_2_heat[week_2_heat$Mouse == z,], n=1)$Heat),
                  "Start of Week 2",
                  unique(week_2_heat$Group[week_2_heat$Mouse == z]),
                  unique(week_2_heat$Sex[week_2_heat$Mouse == z]))
  my_vector1 <- as.data.frame(t(my_vector1))
  my_vector2 <- as.data.frame(t(my_vector2))
  colnames(my_vector1) <- c("Animal ID", "Heat", "Week", "Group", "Sex")
  colnames(my_vector2) <- c("Animal ID", "Heat", "Week", "Group", "Sex")
  my_dataframe <- rbind(my_dataframe, my_vector1, my_vector2)
}

my_dataframe$Week<- factor(my_dataframe$Week, levels=c("End of Week 1", "Start of Week 2"))
my_dataframe$Subgroup <- interaction(my_dataframe$Group, my_dataframe$Sex, sep = "_")
my_dataframe$Heat <- as.numeric(my_dataframe$Heat)

#now the stats
one_hour <- my_dataframe
names(one_hour)[1] <- "Animal.ID"
one_hour$X <- NULL

one_hour$Animal.ID <- as.factor(one_hour$Animal.ID)
one_hour$Week <- as.factor(one_hour$Week)
one_hour$Heat <- as.numeric(one_hour$Heat)
one_hour$Group <- as.factor(one_hour$Group)
one_hour$Sex <- as.factor(one_hour$Sex)
one_hour$Subgroup <- as.factor(one_hour$Subgroup)

one_hour <- one_hour[one_hour$Week == "End of Week 1" | one_hour$Week == "Start of Week 2", ]

model_1h <- lmer(Heat ~ Week*Group*Sex + (1|Animal.ID), one_hour)
anova(model_1h)

one_hour$Time_Group <- interaction(one_hour$Group, one_hour$Week)

model <- lmer(Heat ~ Time_Group + (1|Animal.ID), one_hour)
posthoc <- glht(model, linfct=mcp(Time_Group = c("`30C.End of Week 1` - `22C.End of Week 1` == 0",
                                                 "`30C.Start of Week 2` - `22C.Start of Week 2` == 0")))

summary(posthoc)


# let's look at  end of week 4 vs start of week 5 (1B v)
my_dataframe <- data.frame(matrix(ncol=5,nrow=0, 
                                  dimnames=list(NULL, c("Animal ID", "Heat", "Week", "Group", "Sex"))))
for(z in unique(week_5_heat$Mouse)){
  my_vector1 <- c(z,
                  mean(tail(week_4_heat[week_4_heat$Mouse == z,], n=1)$Heat),
                  "End of Week 4",
                  unique(week_4_heat$Group[week_4_heat$Mouse == z]),
                  unique(week_4_heat$Sex[week_4_heat$Mouse == z]))
  my_vector2 <- c(z,
                  mean(head(week_5_heat[week_5_heat$Mouse == z,], n=1)$Heat),
                  "Start of Week 5",
                  unique(week_5_heat$Group[week_5_heat$Mouse == z]),
                  unique(week_5_heat$Sex[week_5_heat$Mouse == z]))
  my_vector1 <- as.data.frame(t(my_vector1))
  my_vector2 <- as.data.frame(t(my_vector2))
  colnames(my_vector1) <- c("Animal ID", "Heat", "Week", "Group", "Sex")
  colnames(my_vector2) <- c("Animal ID", "Heat", "Week", "Group", "Sex")
  my_dataframe <- rbind(my_dataframe, my_vector1, my_vector2)
}
my_dataframe$Week<- factor(my_dataframe$Week, levels=c("End of Week 4", "Start of Week 5"))
my_dataframe$Subgroup <- interaction(my_dataframe$Group, my_dataframe$Sex, sep = "_")
my_dataframe$Heat <- as.numeric(my_dataframe$Heat)

#now the stats
one_hour <- my_dataframe
names(one_hour)[1] <- "Animal.ID"
one_hour$X <- NULL

one_hour$Animal.ID <- as.factor(one_hour$Animal.ID)
one_hour$Week <- as.factor(one_hour$Week)
one_hour$Heat <- as.numeric(one_hour$Heat)
one_hour$Group <- as.factor(one_hour$Group)
one_hour$Sex <- as.factor(one_hour$Sex)
one_hour$Subgroup <- as.factor(one_hour$Subgroup)

one_hour <- one_hour[one_hour$Week == "End of Week 4" | one_hour$Week == "Start of Week 5", ]

model_1h <- lmer(Heat ~ Week*Group*Sex + (1|Animal.ID), one_hour)
anova(model_1h)

#################################################################
#and now let's look at the 4 hours analysis 

# let's look at  end of week 1 vs start of week 2 (1B ii and iv)


my_dataframe <- data.frame(matrix(ncol=5,nrow=0, 
                                  dimnames=list(NULL, c("Animal ID", "Heat", "Week", "Group", "Sex"))))
for(z in unique(week_2_heat$Mouse)){
  my_vector1 <- c(z,
                  mean(tail(week_1_heat[week_1_heat$Mouse == z,], n=4)$Heat),
                  "End of Week 1",
                  unique(week_1_heat$Group[week_1_heat$Mouse == z]),
                  unique(week_1_heat$Sex[week_1_heat$Mouse == z]))
  my_vector2 <- c(z,
                  mean(head(week_2_heat[week_2_heat$Mouse == z,], n=4)$Heat),
                  "Start of Week 2",
                  unique(week_2_heat$Group[week_2_heat$Mouse == z]),
                  unique(week_2_heat$Sex[week_2_heat$Mouse == z]))
  my_vector1 <- as.data.frame(t(my_vector1))
  my_vector2 <- as.data.frame(t(my_vector2))
  colnames(my_vector1) <- c("Animal ID", "Heat", "Week", "Group", "Sex")
  colnames(my_vector2) <- c("Animal ID", "Heat", "Week", "Group", "Sex")
  my_dataframe <- rbind(my_dataframe, my_vector1, my_vector2)
}

my_dataframe$Week<- factor(my_dataframe$Week, levels=c("End of Week 1", "Start of Week 2"))
my_dataframe$Subgroup <- interaction(my_dataframe$Group, my_dataframe$Sex, sep = "_")
my_dataframe$Heat <- as.numeric(my_dataframe$Heat)

#now the stats
one_hour <- my_dataframe
names(one_hour)[1] <- "Animal.ID"
one_hour$X <- NULL

one_hour$Animal.ID <- as.factor(one_hour$Animal.ID)
one_hour$Week <- as.factor(one_hour$Week)
one_hour$Heat <- as.numeric(one_hour$Heat)
one_hour$Group <- as.factor(one_hour$Group)
one_hour$Sex <- as.factor(one_hour$Sex)
one_hour$Subgroup <- as.factor(one_hour$Subgroup)

one_hour <- one_hour[one_hour$Week == "End of Week 1" | one_hour$Week == "Start of Week 2", ]

model_1h <- lmer(Heat ~ Week*Group*Sex + (1|Animal.ID), one_hour)
anova(model_1h)

one_hour$Time_Group <- interaction(one_hour$Group, one_hour$Week)

model <- lmer(Heat ~ Time_Group + (1|Animal.ID), one_hour)
posthoc <- glht(model, linfct=mcp(Time_Group = c("`30C.End of Week 1` - `22C.End of Week 1` == 0",
                                                 "`30C.Start of Week 2` - `22C.Start of Week 2` == 0")))

summary(posthoc)


# let's look at  end of week 4 vs start of week 5 (1B vi)
my_dataframe <- data.frame(matrix(ncol=5,nrow=0, 
                                  dimnames=list(NULL, c("Animal ID", "Heat", "Week", "Group", "Sex"))))
for(z in unique(week_5_heat$Mouse)){
  my_vector1 <- c(z,
                  mean(tail(week_4_heat[week_4_heat$Mouse == z,], n=4)$Heat),
                  "End of Week 4",
                  unique(week_4_heat$Group[week_4_heat$Mouse == z]),
                  unique(week_4_heat$Sex[week_4_heat$Mouse == z]))
  my_vector2 <- c(z,
                  mean(head(week_5_heat[week_5_heat$Mouse == z,], n=4)$Heat),
                  "Start of Week 5",
                  unique(week_5_heat$Group[week_5_heat$Mouse == z]),
                  unique(week_5_heat$Sex[week_5_heat$Mouse == z]))
  my_vector1 <- as.data.frame(t(my_vector1))
  my_vector2 <- as.data.frame(t(my_vector2))
  colnames(my_vector1) <- c("Animal ID", "Heat", "Week", "Group", "Sex")
  colnames(my_vector2) <- c("Animal ID", "Heat", "Week", "Group", "Sex")
  my_dataframe <- rbind(my_dataframe, my_vector1, my_vector2)
}
my_dataframe$Week<- factor(my_dataframe$Week, levels=c("End of Week 4", "Start of Week 5"))
my_dataframe$Subgroup <- interaction(my_dataframe$Group, my_dataframe$Sex, sep = "_")
my_dataframe$Heat <- as.numeric(my_dataframe$Heat)

#now the stats
one_hour <- my_dataframe
names(one_hour)[1] <- "Animal.ID"
one_hour$X <- NULL

one_hour$Animal.ID <- as.factor(one_hour$Animal.ID)
one_hour$Week <- as.factor(one_hour$Week)
one_hour$Heat <- as.numeric(one_hour$Heat)
one_hour$Group <- as.factor(one_hour$Group)
one_hour$Sex <- as.factor(one_hour$Sex)
one_hour$Subgroup <- as.factor(one_hour$Subgroup)

one_hour <- one_hour[one_hour$Week == "End of Week 4" | one_hour$Week == "Start of Week 5", ]

model_1h <- lmer(Heat ~ Week*Group*Sex + (1|Animal.ID), one_hour)
anova(model_1h)
