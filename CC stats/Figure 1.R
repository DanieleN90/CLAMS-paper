##############################
############ 1B  #############
##############################

setwd("~")

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
library(ggExtra)
library(openxlsx)
library(ggprism)
library(gridExtra)
library(aplot)


library(lme4)
library(lmerTest)
library(multcomp)
library(MuMIn)
library(emmeans)

################ remember we need to remove F787 because her food scale was broken
#this specific analysis is done on clea_CLAMS_upd object because it is the version that has the regression model
load("~/clean_CLAMS_upd.RData")

current_data <- data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("Mouse", 
                                                                        "Heat",
                                                                        "LA",
                                                                        "REE",
                                                                        "FI",
                                                                        "Bodyweight",
                                                                        "Group",
                                                                        "Week",
                                                                        "Temp",
                                                                        "Sex"))))

counter <- 1

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    # we will take the average of 4 intervals to determine how long is each interval
    average_interval_length <- mean(as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[2]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[1]]),
                                    as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[3]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[2]]),
                                    as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[4]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[3]]),
                                    as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[5]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[4]]))
    # now we can calculate how many hours are in each week
    hours_in_week <- length(CLAMS[[i]][[j]][["Data"]][["HEAT"]])*(average_interval_length/60)
    current_data[counter, 1] <- CLAMS[[i]][[j]][["ID"]]
    current_data[counter, 2] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]], CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]]))
    current_data[counter, 3] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]], CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]]))
    current_data[counter, 4] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]][CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]] <= 100], 
                                       CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]][CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]] <= 100]))
    #now we can calculate food intake in kcal/h by taking the total food intake in grams, multiplying by the caloric density (chow diet is 3.07kcal/g, and then dividing it by the total hours in the week)
    current_data[counter, 5] <- (tail(CLAMS[[i]][[j]][["Data"]][["NEW.FEED.ACC"]], n=1)*3.07)/hours_in_week
    current_data[counter, 6] <- CLAMS[[i]][[j]][["Weight"]]
    current_data[counter, 7] <- CLAMS[[i]][[j]][["Group"]]
    current_data[counter, 8] <- names(CLAMS)[i]
    current_data[counter, 9] <- 0
    current_data[counter, 10] <- CLAMS[[i]][[j]][["Sex"]]
    counter <- counter + 1
  }
}



current_data$Mouse <- as.factor(current_data$Mouse)
current_data$Heat <- as.numeric(current_data$Heat)
current_data$LA <- as.numeric(current_data$LA)
current_data$REE <- as.numeric(current_data$REE)
current_data$FI <- as.numeric(current_data$FI)
current_data$Bodyweight <- as.numeric(current_data$Bodyweight)
current_data$Group <- as.factor(current_data$Group)
current_data$Group <- factor(current_data$Group, levels = c("30C", "22C"))
current_data$Week <- as.factor(current_data$Week)
current_data$Temp <- as.factor(current_data$Temp)
current_data$Sex <- as.factor(current_data$Sex)


current_data$Subgroup <- interaction(current_data$Group, current_data$Sex, sep = "_")
current_data$Subgroup <- as.factor(current_data$Subgroup)


#### now we remove the weeks that we decided to exclude
current_data <- current_data[current_data$Week != "Week3" & current_data$Week != "Week6",]


anova(lmer(Bodyweight ~ Group*Sex*Week + (1|Mouse), current_data))




##############################
############ 1C  #############
##############################


MRI <- read_xlsx("~/metabolic parameters wk7.xlsx", col_names = T, skip = 1)
MRI <- MRI[1:21, 1:7]
names(MRI) <- c("Mouse", "bw week6", "fat mass week6", "lean mass week6", "bw week7", "fat mass week7", "lean mass week7")

MRI$Group <- as.character(MRI$Group)
for(i in 1:21){
  if(MRI[i,"Mouse"] %in% raised_at_22){
    MRI[i,"Group"] <- "RT"
  }else if(MRI[i,"Mouse"] %in% raised_at_30){
    MRI[i,"Group"] <- "TN"
  }else{
    MRI[i,"Group"] <- NA
  }
}

MRI <- MRI[complete.cases(MRI),]
MRI$Sex <- c("M", "F", "F", "F","M", "M", "M", "F", "M", "M", "F", "M", "F", "F", "M", "F","F")
MRI <- MRI[complete.cases(MRI),]

MRI$Subgroup <- interaction(MRI$Group, MRI$Sex)
MRI$Adiposity_w6 <- (MRI$`fat mass week6`/MRI$`bw week6`)*100
MRI$Adiposity_w7 <- (MRI$`fat mass week7`/MRI$`bw week7`)*100

### we are using the average of the values at week 6 and 7
MRI$mean_fatmass <- rowMeans(MRI[,c(3,6)])
MRI$mean_bodyweight <- rowMeans(MRI[,c(2,5)])
MRI$mean_leanmass <- rowMeans(MRI[,c(4,7)])
MRI$mean_adiposity <- (MRI$mean_fatmass/MRI$mean_bodyweight)*100


wilcox.test(MRI[MRI$Group == "RT" & MRI$Sex == "F",]$mean_leanmass,
            MRI[MRI$Group == "TN" & MRI$Sex == "F",]$mean_leanmass)


##############################
############ 1D  #############
##############################
wilcox.test(MRI[MRI$Group == "RT" & MRI$Sex == "M",]$mean_leanmass,
            MRI[MRI$Group == "TN" & MRI$Sex == "M",]$mean_leanmass)

##############################
############ 1E  #############
##############################
wilcox.test(MRI[MRI$Group == "RT" & MRI$Sex == "F",]$mean_adiposity,
            MRI[MRI$Group == "TN" & MRI$Sex == "F",]$mean_adiposity)

##############################
############ 1F  #############
##############################
wilcox.test(MRI[MRI$Group == "RT" & MRI$Sex == "M",]$mean_adiposity,
            MRI[MRI$Group == "TN" & MRI$Sex == "M",]$mean_adiposity)


##############################
############ 1I  #############
##############################


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


##############################
############ 1J  #############
##############################
load("~/clean_CLAMS_upd.RData")

current_data <- data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("Mouse", 
                                                                        "Heat",
                                                                        "LA",
                                                                        "REE",
                                                                        "FI",
                                                                        "Bodyweight",
                                                                        "Group",
                                                                        "Week",
                                                                        "Temp",
                                                                        "Sex"))))

counter <- 1

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    # we will take the average of 4 intervals to determine how long is each interval
    average_interval_length <- mean(as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[2]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[1]]),
                                    as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[3]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[2]]),
                                    as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[4]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[3]]),
                                    as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[5]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[4]]))
    # now we can calculate how many hours are in each week
    hours_in_week <- length(CLAMS[[i]][[j]][["Data"]][["HEAT"]])*(average_interval_length/60)
    current_data[counter, 1] <- CLAMS[[i]][[j]][["ID"]]
    current_data[counter, 2] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]], CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]]))
    current_data[counter, 3] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]], CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]]))
    current_data[counter, 4] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]][CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]] <= 100], 
                                       CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]][CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]] <= 100]))
    #now we can calculate food intake in kcal/h by taking the total food intake in grams, multiplying by the caloric density (chow diet is 3.07kcal/g, and then dividing it by the total hours in the week)
    current_data[counter, 5] <- (tail(CLAMS[[i]][[j]][["Data"]][["NEW.FEED.ACC"]], n=1)*3.07)/hours_in_week
    current_data[counter, 6] <- CLAMS[[i]][[j]][["Weight"]]
    current_data[counter, 7] <- CLAMS[[i]][[j]][["Group"]]
    current_data[counter, 8] <- names(CLAMS)[i]
    current_data[counter, 9] <- 0
    current_data[counter, 10] <- CLAMS[[i]][[j]][["Sex"]]
    counter <- counter + 1
  }
}



current_data$Mouse <- as.factor(current_data$Mouse)
current_data$Heat <- as.numeric(current_data$Heat)
current_data$LA <- as.numeric(current_data$LA)
current_data$REE <- as.numeric(current_data$REE)
current_data$FI <- as.numeric(current_data$FI)
current_data$Bodyweight <- as.numeric(current_data$Bodyweight)
current_data$Group <- as.factor(current_data$Group)
current_data$Group <- factor(current_data$Group, levels = c("30C", "22C"))
current_data$Week <- as.factor(current_data$Week)
current_data$Temp <- as.factor(current_data$Temp)
current_data$Sex <- as.factor(current_data$Sex)


current_data$Subgroup <- interaction(current_data$Group, current_data$Sex, sep = "_")
current_data$Subgroup <- as.factor(current_data$Subgroup)


#### now we remove the weeks that we decided to exclude
current_data <- current_data[current_data$Week != "Week3" & current_data$Week != "Week6",]
anova(lmer(Heat ~ Group*Sex*Week + (1|Mouse), current_data))
current_data$Interaction <- interaction(current_data$Group, current_data$Week, sep = "_")

model <- lmer(Heat ~ 0 + Interaction + (1|Mouse), current_data)
posthoc <- glht(model, linfct=mcp(Interaction = c("`30C_Week1` - `22C_Week1` == 0",
                                                  "`30C_Week2` - `22C_Week2` == 0",
                                                  "`30C_Week4` - `22C_Week4` == 0",
                                                  "`30C_Week5` - `22C_Week5` == 0",
                                                  "`30C_Week7` - `22C_Week7` == 0")))
summary(posthoc)


##############################
############ 1K  #############
##############################


load("~/clean_CLAMS_upd.RData")

EE_CalR<- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("Animal ID", 
                                                                 "Energy expenditure", 
                                                                 "Weight",
                                                                 "Week", 
                                                                 "Group", 
                                                                 "Sex"))))
for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    EE_CalR[j + (17*(i-1)), 1] <- CLAMS[[i]][[j]][["ID"]]
    EE_CalR[j + (17*(i-1)), 2] <- mean(CLAMS[[i]][[j]][["Data"]][["HEAT"]])
    EE_CalR[j + (17*(i-1)), 3] <- CLAMS[[i]][[j]][["Weight"]]
    EE_CalR[j + (17*(i-1)), 4] <- names(CLAMS)[i]
    EE_CalR[j + (17*(i-1)), 5] <- CLAMS[[i]][[j]][["Group"]]
    EE_CalR[j + (17*(i-1)), 6] <- CLAMS[[i]][[j]][["Sex"]]
  }
}

EE_CalR <- EE_CalR[complete.cases(EE_CalR), ]
# we need to get rid of M744 due to problems
EE_CalR <- EE_CalR[EE_CalR$Animal.ID != "M744",]

# let's remove week 3 and 6 since we are not inclulding them in the analysis

EE_CalR <- EE_CalR[EE_CalR$Week != "Week3",]
EE_CalR <- EE_CalR[EE_CalR$Week != "Week6",]

EE_CalR$Subgroup <- interaction(EE_CalR$Group, EE_CalR$Week)

model <- lmer(Energy.expenditure ~ Group + Weight + Sex + (1+Weight|Week) + (1|Animal.ID), data = EE_CalR)
model2 <- lmer(Energy.expenditure ~ Subgroup + (1+Weight|Week) + (1|Animal.ID), data = EE_CalR)
anova(model)

multcomp.sub <- glht(model2, linfct=mcp(Subgroup=c("`30C.Week1` - `22C.Week1` == 0",
                                                   "`30C.Week2` - `22C.Week2` == 0",
                                                   "`30C.Week4` - `22C.Week4` == 0",
                                                   "`30C.Week5` - `22C.Week5` == 0",
                                                   "`30C.Week7` - `22C.Week7` == 0")))
summary(multcomp.sub)


##############################
############ 1N  #############
##############################

load("~/clean_CLAMS_upd.RData")

current_data <- data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("Mouse", 
                                                                        "Heat",
                                                                        "LA",
                                                                        "REE",
                                                                        "FI",
                                                                        "Bodyweight",
                                                                        "Group",
                                                                        "Week",
                                                                        "Temp",
                                                                        "Sex"))))

counter <- 1

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    # we will take the average of 4 intervals to determine how long is each interval
    average_interval_length <- mean(as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[2]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[1]]),
                                    as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[3]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[2]]),
                                    as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[4]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[3]]),
                                    as.numeric(CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[5]] - CLAMS[[i]][[j]][["Data"]][["DATE.TIME"]][[4]]))
    # now we can calculate how many hours are in each week
    hours_in_week <- length(CLAMS[[i]][[j]][["Data"]][["HEAT"]])*(average_interval_length/60)
    current_data[counter, 1] <- CLAMS[[i]][[j]][["ID"]]
    current_data[counter, 2] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]], CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]]))
    current_data[counter, 3] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]], CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]]))
    current_data[counter, 4] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]][CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]] <= 100], 
                                       CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]][CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]] <= 100]))
    #now we can calculate food intake in kcal/h by taking the total food intake in grams, multiplying by the caloric density (chow diet is 3.07kcal/g, and then dividing it by the total hours in the week)
    current_data[counter, 5] <- (tail(CLAMS[[i]][[j]][["Data"]][["NEW.FEED.ACC"]], n=1)*3.07)/hours_in_week
    current_data[counter, 6] <- CLAMS[[i]][[j]][["Weight"]]
    current_data[counter, 7] <- CLAMS[[i]][[j]][["Group"]]
    current_data[counter, 8] <- names(CLAMS)[i]
    current_data[counter, 9] <- 0
    current_data[counter, 10] <- CLAMS[[i]][[j]][["Sex"]]
    counter <- counter + 1
  }
}



current_data$Mouse <- as.factor(current_data$Mouse)
current_data$Heat <- as.numeric(current_data$Heat)
current_data$LA <- as.numeric(current_data$LA)
current_data$REE <- as.numeric(current_data$REE)
current_data$FI <- as.numeric(current_data$FI)
current_data$Bodyweight <- as.numeric(current_data$Bodyweight)
current_data$Group <- as.factor(current_data$Group)
current_data$Group <- factor(current_data$Group, levels = c("30C", "22C"))
current_data$Week <- as.factor(current_data$Week)
current_data$Temp <- as.factor(current_data$Temp)
current_data$Sex <- as.factor(current_data$Sex)


current_data$Subgroup <- interaction(current_data$Group, current_data$Sex, sep = "_")
current_data$Subgroup <- as.factor(current_data$Subgroup)


#### now we remove the weeks that we decided to exclude
current_data <- current_data[current_data$Week != "Week3" & current_data$Week != "Week6",]
anova(lmer(LA ~ Group*Sex*Week + (1|Mouse), current_data))

current_data$Interaction <- interaction(current_data$Sex, current_data$Week, sep = "_")

model <- lmer(LA ~ 0 + Interaction + (1|Mouse), current_data)
posthoc <- glht(model, linfct=mcp(Interaction = c("`F_Week1` - `M_Week1` == 0",
                                                  "`F_Week2` - `M_Week2` == 0",
                                                  "`F_Week4` - `M_Week4` == 0",
                                                  "`F_Week5` - `M_Week5` == 0",
                                                  "`F_Week7` - `M_Week7` == 0")))
summary(posthoc)


##############################
############ 1O  #############
##############################

anova(lmer(FI ~ Group*Sex*Week + (1|Mouse), current_data))


##############################
############ 1P  #############
##############################

load("~/clean_CLAMS_upd+CITv2.RData")
EE_data <- data.frame(matrix(ncol=6, nrow = 0, dimnames = list(NULL, c("Value", "Component", "Group", "Week", "Mouse", "Sex"))))
counter <- 1
for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    if( i == 1 & j >8){
      for(k in 0:3){
        EE_data[counter + k, 1] <- CLAMS[[i]][[j]][["EE Analysis"]][[1 + k]]
        EE_data[counter + k, 2] <- names(CLAMS[[i]][[j]][["EE Analysis"]])[1 + k]
        EE_data[counter + k, 3] <- CLAMS[[i]][[j]][["Group"]]
        EE_data[counter + k, 4] <- names(CLAMS)[i]
        EE_data[counter + k, 5] <- names(CLAMS[[i]])[j]
        EE_data[counter + k, 6] <- CLAMS[[i]][[j]][["Sex"]]
      }
    }else{
      for(k in 0:4){
        EE_data[counter + k, 1] <- CLAMS[[i]][[j]][["EE Analysis"]][[1 + k]]
        EE_data[counter + k, 2] <- names(CLAMS[[i]][[j]][["EE Analysis"]])[1 + k]
        EE_data[counter + k, 3] <- CLAMS[[i]][[j]][["Group"]]
        EE_data[counter + k, 4] <- names(CLAMS)[i]
        EE_data[counter + k, 5] <- names(CLAMS[[i]])[j]
        EE_data[counter + k, 6] <- CLAMS[[i]][[j]][["Sex"]]
      }
    }
    counter <- counter + k + 1
  }
}
EE_data <- EE_data[EE_data$Mouse != "F787",]



EE_data$Group_Week <- paste(EE_data$Group, EE_data$Week, sep = "_")
EE_data$Component_Group_Week <- paste(EE_data$Component, EE_data$Group, EE_data$Week, sep = "_")
EE_data$Week_Sex <- paste(EE_data$Sex, EE_data$Week, sep = "_")


#### let's do the stats

EE_data$Value <- as.numeric(EE_data$Value)
EE_data$Component <- as.factor(EE_data$Component)
EE_data$Group <- as.factor(EE_data$Group)
EE_data$Week <- as.factor(EE_data$Week)
EE_data$Mouse <- as.factor(EE_data$Mouse)
EE_data$Group_Week <- as.factor(EE_data$Group_Week)
EE_data$Component_Group_Week <- as.factor(EE_data$Component_Group_Week)
EE_data$Week_Sex <- as.factor(EE_data$Week_Sex)

#PAEE 
PAEE <- EE_data[EE_data$Component == "PAEE",]
PAEE <- PAEE[PAEE$Week != "Week3",]
PAEE <- PAEE[PAEE$Week != "Week6",]


model_PAEE <- lmer(Value ~ Group*Week + (1|Mouse), PAEE)
anova(model_PAEE)


model <- lmer(Value ~ 0 + Week + (1|Mouse), PAEE)
posthoc <- glht(model, linfct=mcp(Week_Sex = c("`F_Week1` - `M_Week1` == 0",
                                               "`F_Week2` - `M_Week2` == 0",
                                               "`F_Week4` - `M_Week4` == 0",
                                               "`F_Week5` - `M_Week5` == 0",
                                               "`F_Week7` - `M_Week7` == 0")))

summary(posthoc)



#TEF 
TEF <- EE_data[EE_data$Component == "TEF",]
TEF <- TEF[TEF$Week != "Week3",]
TEF <- TEF[TEF$Week != "Week6",]



model_TEF <- lmer(Value ~ Group*Week + (1|Mouse), TEF)
anova(model_TEF)


model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), TEF)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week1` - `22C_Week1` == 0",
                                                 "`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week4` - `22C_Week4` == 0",
                                                 "`30C_Week5` - `22C_Week5` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))



summary(posthoc)


#CIT 
CIT <- EE_data[EE_data$Component == "CIT",]
CIT <- CIT[CIT$Week != "Week3",]
CIT <- CIT[CIT$Week != "Week6",]



model_CIT <- lmer(Value ~ Group*Week + (1|Mouse), CIT)
anova(model_CIT)


model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), CIT)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week4` - `22C_Week4` == 0",
                                                 "`30C_Week5` - `22C_Week5` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))



summary(posthoc)


##############################
############ 1Q  #############
##############################


EE_data_prop <- EE_data

### we need to fix the proportions here
### basically, because the of the mismatch problem of the BMR on week 1 of TN, if we just do the ratio of PAEE, TEF, and BMR over TEE we get results such that the sum is not 1
### the solution for this problem is to scale the proportion
### basically we will do the ratio between the TEE and the sum of the PAEE, TEF, and BMR, and then we will scale those three components by that value before doing the proportions
scale_factors <- data.frame(matrix(ncol=2, nrow = 0, dimnames = list(NULL, c("Mouse", "Value"))))
counter <- 1

for(j in 1:length(CLAMS[[1]])){
  if (CLAMS[[1]][[j]]$Group == "22C") next
  scale_factors[counter, "Mouse"] <- names(CLAMS[[1]])[j]
  scale_factors[counter, "Value"] <- CLAMS[[1]][[j]]$`EE Analysis`$TEE/(CLAMS[[1]][[j]]$`EE Analysis`$BMR + CLAMS[[1]][[j]]$`EE Analysis`$TEF + CLAMS[[1]][[j]]$`EE Analysis`$PAEE)
  counter <- counter + 1
}



for (i in 1:nrow(EE_data_prop)) {
  if(EE_data_prop[i, "Week"] == "Week1" & EE_data_prop[i, "Group"] == "30C" & EE_data_prop[i, "Component"] != "TEE"){
    EE_data_prop[i, 1] <- (EE_data_prop[i, 1]/CLAMS[[EE_data_prop[i,4]]][[as.character(EE_data_prop[i,5])]]$`EE Analysis`$`TEE`)*scale_factors[scale_factors$Mouse == as.character(EE_data_prop[i,5]), "Value"]
  }else{
    EE_data_prop[i, 1] <- EE_data_prop[i, 1]/CLAMS[[EE_data_prop[i,4]]][[as.character(EE_data_prop[i,5])]]$`EE Analysis`$`TEE`
  }
}
EE_data <- EE_data[EE_data$Mouse != "F787",]





data_for_plot_prop <- as.data.frame(group_by(EE_data_prop, Component_Group_Week) %>% dplyr::summarize(m = mean(Value), n=n(), se=sd(Value)/(sqrt(n))))


data_for_plot_prop$Week <- c(rep(c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), 3), c("Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), rep(c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), 6) )
data_for_plot_prop$Component <- c(rep("BMR", 14), rep("CIT", 13), rep("PAEE", 14), rep("TEE", 14), rep("TEF", 14))
data_for_plot_prop$Group <- c(rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 6), rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7))

#need to make an adjusted dataframe for the error bars
data_for_plot_prop_adj <- data_for_plot_prop
#adjust CIT which is the first one from the bottom after BMR
for(i in 0:6){
  data_for_plot_prop_adj[15 + i,2] <- data_for_plot_prop_adj[15 + i,2] + data_for_plot_prop_adj[1 + i,2]
}
for(i in 0:5){
  data_for_plot_prop_adj[22 + i,2] <- data_for_plot_prop_adj[22 + i,2] + data_for_plot_prop_adj[9 + i,2]
}
# now adjust TEF
for(i in 0:6){
  data_for_plot_prop_adj[56 + i,2] <- data_for_plot_prop_adj[56 + i,2] + data_for_plot_prop_adj[15 + i,2]
}
data_for_plot_prop_adj[63,2] <- data_for_plot_prop_adj[63,2] + data_for_plot_prop_adj[8,2]
for(i in 0:5){
  data_for_plot_prop_adj[64 + i,2] <- data_for_plot_prop_adj[64 + i,2] + data_for_plot_prop_adj[22 + i,2]
}
# and finally PAEE
for(i in 0:13){
  data_for_plot_prop_adj[28 + i,2] <- data_for_plot_prop_adj[28 + i,2] + data_for_plot_prop_adj[56 + i,2]
}





#### let's do the stats

EE_data_prop$Value <- as.numeric(EE_data_prop$Value)
EE_data_prop$Component <- as.factor(EE_data_prop$Component)
EE_data_prop$Group <- as.factor(EE_data_prop$Group)
EE_data_prop$Week <- as.factor(EE_data_prop$Week)
EE_data_prop$Mouse <- as.factor(EE_data_prop$Mouse)
EE_data_prop$Group_Week <- as.factor(EE_data_prop$Group_Week)
EE_data_prop$Component_Group_Week <- as.factor(EE_data_prop$Component_Group_Week)
EE_data_prop$Week_Sex <- as.factor(EE_data_prop$Week_Sex)


#PAEE 
PAEE <- EE_data_prop[EE_data_prop$Component == "PAEE",]
PAEE <- PAEE[PAEE$Week != "Week3",]
PAEE <- PAEE[PAEE$Week != "Week6",]


model_PAEE <- lmer(Value ~ Group*Week + (1|Mouse), PAEE)
anova(model_PAEE)


model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), PAEE)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week1` - `22C_Week1` == 0",
                                                 "`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week4` - `22C_Week4` == 0",
                                                 "`30C_Week5` - `22C_Week5` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))

summary(posthoc)



#TEF 
TEF <- EE_data_prop[EE_data_prop$Component == "TEF",]
TEF <- TEF[TEF$Week != "Week3",]
TEF <- TEF[TEF$Week != "Week6",]



model_TEF <- lmer(Value ~ Group*Week + (1|Mouse), TEF)
anova(model_TEF)


model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), TEF)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week1` - `22C_Week1` == 0",
                                                 "`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week4` - `22C_Week4` == 0",
                                                 "`30C_Week5` - `22C_Week5` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))


summary(posthoc)


#CIT 
CIT <- EE_data_prop[EE_data_prop$Component == "CIT",]
CIT <- CIT[CIT$Week != "Week3",]
CIT <- CIT[CIT$Week != "Week6",]



model_CIT <- lmer(Value ~ Group*Week + (1|Mouse), CIT)
anova(model_CIT)


model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), CIT)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week4` - `22C_Week4` == 0",
                                                 "`30C_Week5` - `22C_Week5` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))

summary(posthoc)
