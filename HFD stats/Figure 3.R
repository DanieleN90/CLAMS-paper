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
library(readr)

library(lme4)
library(lmerTest)
library(multcomp)
library(MuMIn)



load("~/CLAMS_HFD_clean.RData")

current_data <- data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, c("Mouse", 
                                                                       "Heat",
                                                                       "LA",
                                                                       "Bodyweight",
                                                                       "Lean_Mass",
                                                                       "Adiposity",
                                                                       "Group",
                                                                       "Week",
                                                                       "Sex"))))

counter <- 1
for (i in 1:7) {
  for(j in 1:length(CLAMS_HFD[[i]])){
    current_data[counter, 1] <- CLAMS_HFD[[i]][[j]][["ID"]]
    current_data[counter, 2] <- mean(CLAMS_HFD[[i]][[j]][["Data"]][["Heat"]])
    current_data[counter, 3] <- sum(CLAMS_HFD[[i]][[j]][["Data"]][["X.Ambulatory"]])
    current_data[counter, 4] <- CLAMS_HFD[[i]][[j]][["Weight"]]
    current_data[counter, 5] <- CLAMS_HFD[[i]][[j]][["Lean Mass"]]
    current_data[counter, 6] <- CLAMS_HFD[[i]][[j]][["Adiposity"]]
    current_data[counter, 7] <- CLAMS_HFD[[i]][[j]][["Group"]]
    current_data[counter, 8] <- names(CLAMS_HFD)[i]
    current_data[counter, 9] <- CLAMS_HFD[[i]][[j]][["Sex"]]
    counter <- counter + 1
  }
}



current_data$Mouse <- as.factor(current_data$Mouse)
current_data$Heat <- as.numeric(current_data$Heat)
current_data$LA <- as.numeric(current_data$LA)
current_data$Bodyweight <- as.numeric(current_data$Bodyweight)
current_data$Lean_Mass <- as.numeric(current_data$Lean_Mass)
current_data$Adiposity <- as.numeric(current_data$Adiposity)
current_data$Group <- as.factor(current_data$Group)
current_data$Group <- factor(current_data$Group, levels = c("30C", "22C"))
current_data$Week <- as.factor(current_data$Week)
current_data$Sex <- as.factor(current_data$Sex)


current_data$Subgroup1 <- interaction(current_data$Group, current_data$Sex, sep = "_")
current_data$Subgroup1 <- as.factor(current_data$Subgroup1)

current_data$Subgroup2 <- interaction(current_data$Sex, current_data$Week, sep = "_")
current_data$Subgroup2 <- as.factor(current_data$Subgroup2)

current_data$Subgroup3 <- interaction(current_data$Group, current_data$Subgroup2, sep = "_")
current_data$Subgroup3 <- as.factor(current_data$Subgroup3)



#### now we keep the weeks that we decided to include in the analysis
current_data <- current_data[current_data$Week == "Week1" | current_data$Week == "Week2" | current_data$Week == "Week7",]

##############################
############ 3C  #############
##############################

anova(lmer(Heat ~ Group*Sex + Week + (1|Mouse), current_data))

##############################
############ 3E  #############
##############################

load("~/CLAMS_HFD_clean.RData")


TEE_ANCOVA<- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("Animal ID", 
                                                                    "Energy expenditure", 
                                                                    "Weight",
                                                                    "Week", 
                                                                    "Group", 
                                                                    "Sex"))))
for (i in 1:7) {
  for(j in 1:length(CLAMS_HFD[[i]])){
    TEE_ANCOVA[j + (23*(i-1)), 1] <- CLAMS_HFD[[i]][[j]][["ID"]]
    TEE_ANCOVA[j + (23*(i-1)), 2] <- mean(CLAMS_HFD[[i]][[j]][["Data"]][["Heat"]])
    TEE_ANCOVA[j + (23*(i-1)), 3] <- CLAMS_HFD[[i]][[j]][["Weight"]]
    TEE_ANCOVA[j + (23*(i-1)), 4] <- names(CLAMS_HFD)[i]
    TEE_ANCOVA[j + (23*(i-1)), 5] <- CLAMS_HFD[[i]][[j]][["Group"]]
    TEE_ANCOVA[j + (23*(i-1)), 6] <- CLAMS_HFD[[i]][[j]][["Sex"]]
  }
}


model <- lmer(Energy.expenditure ~ Group*Weight + Sex + (1+Weight|Week) + (1|Animal.ID), data = TEE_ANCOVA)

anova(model)


##############################
############ 3F  #############
##############################


mice_names <- read_xlsx("GTT ITT.xlsx")[,1:2]
mice_names$Mouse

born_at_30 <- mice_names$Mouse[mice_names$...1 == "22C born"]
born_at_22 <- mice_names$Mouse[mice_names$...1 == "30C born"]

load("~/CLAMS_HFD_clean.RData")

# let's first define the days that I will need

chow_days_30 <- tail(unique(str_subset(unlist(strsplit(as.character(CLAMS_HFD[[1]][[1]]$Data$Date.Time), "\\s+")), "2019")), 3)
hfd_days_30 <- unique(str_subset(unlist(strsplit(as.character(CLAMS_HFD[[2]][[1]]$Data$Date.Time), "\\s+")), "2019"))
chow_days_22 <- tail(unique(str_subset(unlist(strsplit(as.character(CLAMS_HFD[[1]][[12]]$Data$Date.Time), "\\s+")), "2019")), 3)
hfd_days_22 <- unique(str_subset(unlist(strsplit(as.character(CLAMS_HFD[[2]][[12]]$Data$Date.Time), "\\s+")), "2019"))

days_30 <- union(chow_days_30, hfd_days_30)
days_22 <- union(chow_days_22, hfd_days_22)

all_mice <- c(born_at_22, born_at_30)

# I will take the 24h before HFD and the 24h after (consider that each time time point is 30 min, so 24h is 48 intervals)

Heat_24h <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("Animal ID", "Heat", "Time", "Diet", "Group", "Sex"))))

for (j in 1:23) {
  for (i in 1:2){
    if(i == 1 & all_mice[j] %in% born_at_22){
      Heat_24h[j + (23*(i-1)), 1] <- all_mice[j]
      Heat_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[1]][[all_mice[j]]][["Data"]][["Heat"]], n= 48), na.rm = T)
      Heat_24h[j + (23*(i-1)), 3] <- "Before"
      Heat_24h[j + (23*(i-1)), 4] <- "Chow"
      Heat_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      Heat_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]]
    }else if (i == 1 & all_mice[j] %in% born_at_30){
      Heat_24h[j + (23*(i-1)), 1] <- all_mice[j]
      Heat_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[1]][[all_mice[j]]][["Data"]][["Heat"]], n= 48), na.rm = T)
      Heat_24h[j + (23*(i-1)), 3] <- "Before"
      Heat_24h[j + (23*(i-1)), 4] <- "Chow"
      Heat_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      Heat_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]]
    }else if (i == 2 & all_mice[j] %in% born_at_22){
      Heat_24h[j + (23*(i-1)), 1] <- all_mice[j]
      Heat_24h[j + (23*(i-1)), 2] <- mean(head( CLAMS_HFD[[2]][[all_mice[j]]][["Data"]][["Heat"]], n= 48), na.rm = T)
      Heat_24h[j + (23*(i-1)), 3] <- "After"
      Heat_24h[j + (23*(i-1)), 4] <- "HFD"
      Heat_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      Heat_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]] 
    }else if (i == 2 & all_mice[j] %in% born_at_30){
      Heat_24h[j + (23*(i-1)), 1] <- all_mice[j]
      Heat_24h[j + (23*(i-1)), 2] <- mean(head( CLAMS_HFD[[2]][[all_mice[j]]][["Data"]][["Heat"]], n= 48), na.rm = T)
      Heat_24h[j + (23*(i-1)), 3] <- "After"
      Heat_24h[j + (23*(i-1)), 4] <- "HFD"
      Heat_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      Heat_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]] 
    }
  }
}

Heat_24h$Time <- factor(Heat_24h$Time, levels = c("Before", "After"))
Heat_24h$Subgroup <- interaction(Heat_24h$Group, Heat_24h$Sex, sep = "_")


anova(lmer(Heat ~ Group*Diet + (1|Animal.ID), Heat_24h))

# post-hoc for the analysis from above 
Heat_24h$Subgroup2 <- interaction(Heat_24h$Group, Heat_24h$Diet, sep = "_")
model <- lmer(Heat ~ 0 + Subgroup2 + (1|Animal.ID), Heat_24h)

posthoc <- glht(model, linfct=mcp(Subgroup2 = c("`22C_Chow` - `30C_Chow` == 0",
                                                "`22C_HFD` - `30C_HFD` == 0",
                                                "`22C_Chow` - `22C_HFD` == 0",
                                                "`30C_Chow` - `30C_HFD` == 0")))
summary(posthoc)

##and now let's test the heat gain between the two time points 
heat_gain <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Animal ID", "Heat", "Group", "Sex"))))

for (j in 1:23) {
  heat_gain[j, 1] <- Heat_24h$Animal.ID[j]
  heat_gain[j, 2] <- Heat_24h[Heat_24h$Animal.ID == Heat_24h$Animal.ID[j] & Heat_24h$Time == "After",]$Heat - Heat_24h[Heat_24h$Animal.ID == Heat_24h$Animal.ID[j] & Heat_24h$Time == "Before",]$Heat
  heat_gain[j, 3] <- Heat_24h$Group[j]
  heat_gain[j, 4] <- Heat_24h$Sex[j]
}

t.test(heat_gain[heat_gain$Group == "22C",]$Heat, heat_gain[heat_gain$Group == "30C",]$Heat)




##############################
############ 3G  #############
##############################
Heat_24h <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("Animal ID", "Heat", "Time", "Diet", "Group", "Sex"))))

for (j in 1:23) {
  for (i in 1:2){
    if(i == 1 & all_mice[j] %in% born_at_22){
      Heat_24h[j + (23*(i-1)), 1] <- all_mice[j]
      Heat_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[1]][[all_mice[j]]][["Data"]][["Heat"]], n= 48), na.rm = T)
      Heat_24h[j + (23*(i-1)), 3] <- "Before"
      Heat_24h[j + (23*(i-1)), 4] <- "Chow"
      Heat_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      Heat_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]]
    }else if (i == 1 & all_mice[j] %in% born_at_30){
      Heat_24h[j + (23*(i-1)), 1] <- all_mice[j]
      Heat_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[1]][[all_mice[j]]][["Data"]][["Heat"]], n= 48), na.rm = T)
      Heat_24h[j + (23*(i-1)), 3] <- "Before"
      Heat_24h[j + (23*(i-1)), 4] <- "Chow"
      Heat_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      Heat_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]]
    }else if (i == 2 & all_mice[j] %in% born_at_22){
      Heat_24h[j + (23*(i-1)), 1] <- all_mice[j]
      Heat_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[2]][[all_mice[j]]][["Data"]][["Heat"]], n= 48), na.rm = T)
      Heat_24h[j + (23*(i-1)), 3] <- "After"
      Heat_24h[j + (23*(i-1)), 4] <- "HFD"
      Heat_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[2]][[all_mice[j]]][["Group"]]
      Heat_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[2]][[all_mice[j]]][["Sex"]] 
    }else if (i == 2 & all_mice[j] %in% born_at_30){
      Heat_24h[j + (23*(i-1)), 1] <- all_mice[j]
      Heat_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[2]][[all_mice[j]]][["Data"]][["Heat"]], n= 48), na.rm = T)
      Heat_24h[j + (23*(i-1)), 3] <- "After"
      Heat_24h[j + (23*(i-1)), 4] <- "HFD"
      Heat_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[2]][[all_mice[j]]][["Group"]]
      Heat_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[2]][[all_mice[j]]][["Sex"]] 
    }
  }
}


Heat_24h$Time <- factor(Heat_24h$Time, levels = c("Before", "After"))
Heat_24h$Subgroup <- interaction(Heat_24h$Group, Heat_24h$Sex, sep = "_")

#### Supplementary File 2Cvii
anova(lmer(Heat ~ Group*Diet + (1|Animal.ID), Heat_24h))

# post-hoc analysis, Supplementary File 2Cix
Heat_24h$Subgroup2 <- interaction(Heat_24h$Group, Heat_24h$Diet, sep = "_")
model <- lmer(Heat ~ 0 + Subgroup2 + (1|Animal.ID), Heat_24h)

posthoc <- glht(model, linfct=mcp(Subgroup2 = c("`22C_Chow` - `30C_Chow` == 0",
                                                "`22C_HFD` - `30C_HFD` == 0",
                                                "`22C_Chow` - `22C_HFD` == 0",
                                                "`30C_Chow` - `30C_HFD` == 0")))
summary(posthoc)

#and now let's test the heat gain between the two time points (Supplemental File 2Cxi)
heat_gain <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Animal ID", "Heat", "Group", "Sex"))))

for (j in 1:23) {
  heat_gain[j, 1] <- Heat_24h$Animal.ID[j]
  heat_gain[j, 2] <- Heat_24h[Heat_24h$Animal.ID == Heat_24h$Animal.ID[j] & Heat_24h$Time == "After",]$Heat - Heat_24h[Heat_24h$Animal.ID == Heat_24h$Animal.ID[j] & Heat_24h$Time == "Before",]$Heat
  heat_gain[j, 3] <- Heat_24h$Group[j]
  heat_gain[j, 4] <- Heat_24h$Sex[j]
}

t.test(heat_gain[heat_gain$Group == "22C",]$Heat, heat_gain[heat_gain$Group == "30C",]$Heat)



##############################
############ 3I  #############
##############################
RER_24h <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("Animal ID", "RER", "Time", "Diet", "Group", "Sex"))))

for (j in 1:23) {
  for (i in 1:2){
    if(i == 1 & all_mice[j] %in% born_at_22){
      RER_24h[j + (23*(i-1)), 1] <- all_mice[j]
      RER_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[1]][[all_mice[j]]][["Data"]][["RER"]], n= 48), na.rm = T)
      RER_24h[j + (23*(i-1)), 3] <- "Before"
      RER_24h[j + (23*(i-1)), 4] <- "Chow"
      RER_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      RER_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]]
    }else if (i == 1 & all_mice[j] %in% born_at_30){
      RER_24h[j + (23*(i-1)), 1] <- all_mice[j]
      RER_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[1]][[all_mice[j]]][["Data"]][["RER"]], n= 48), na.rm = T)
      RER_24h[j + (23*(i-1)), 3] <- "Before"
      RER_24h[j + (23*(i-1)), 4] <- "Chow"
      RER_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      RER_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]]
    }else if (i == 2 & all_mice[j] %in% born_at_22){
      RER_24h[j + (23*(i-1)), 1] <- all_mice[j]
      RER_24h[j + (23*(i-1)), 2] <- mean(head( CLAMS_HFD[[2]][[all_mice[j]]][["Data"]][["RER"]], n= 48), na.rm = T)
      RER_24h[j + (23*(i-1)), 3] <- "After"
      RER_24h[j + (23*(i-1)), 4] <- "HFD"
      RER_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      RER_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]] 
    }else if (i == 2 & all_mice[j] %in% born_at_30){
      RER_24h[j + (23*(i-1)), 1] <- all_mice[j]
      RER_24h[j + (23*(i-1)), 2] <- mean(head( CLAMS_HFD[[2]][[all_mice[j]]][["Data"]][["RER"]], n= 48), na.rm = T)
      RER_24h[j + (23*(i-1)), 3] <- "After"
      RER_24h[j + (23*(i-1)), 4] <- "HFD"
      RER_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      RER_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]] 
    }
  }
}


RER_24h$Time <- factor(RER_24h$Time, levels = c("Before", "After"))
RER_24h$Subgroup <- interaction(RER_24h$Group, RER_24h$Sex, sep = "_")


#### Supplemental File 2Cii
anova(lmer(RER ~ Group*Diet + (1|Animal.ID), RER_24h))

# post-hoc for the analysis from above (Supplemental File 2Civ)
RER_24h$Subgroup2 <- interaction(RER_24h$Group, RER_24h$Diet, sep = "_")
model <- lmer(RER ~ 0 + Subgroup2 + (1|Animal.ID), RER_24h)

posthoc <- glht(model, linfct=mcp(Subgroup2 = c("`22C_Chow` - `30C_Chow` == 0",
                                                "`22C_HFD` - `30C_HFD` == 0",
                                                "`22C_Chow` - `22C_HFD` == 0",
                                                "`30C_Chow` - `30C_HFD` == 0")))
summary(posthoc)


#and now let's test the RER gain between the two time points (Supplemental File 2Cvi)
RER_gain <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Animal ID", "RER", "Group", "Sex"))))

for (j in 1:23) {
  RER_gain[j, 1] <- RER_24h$Animal.ID[j]
  RER_gain[j, 2] <- RER_24h[RER_24h$Animal.ID == RER_24h$Animal.ID[j] & RER_24h$Time == "After",]$RER - RER_24h[RER_24h$Animal.ID == RER_24h$Animal.ID[j] & RER_24h$Time == "Before",]$RER
  RER_gain[j, 3] <- RER_24h$Group[j]
  RER_gain[j, 4] <- RER_24h$Sex[j]
}

t.test(RER_gain[RER_gain$Group == "22C",]$RER, RER_gain[RER_gain$Group == "30C",]$RER)



##############################
############ 3J  #############
##############################
RER_24h <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("Animal ID", "RER", "Time", "Diet", "Group", "Sex"))))

for (j in 1:23) {
  for (i in 1:2){
    if(i == 1 & all_mice[j] %in% born_at_22){
      RER_24h[j + (23*(i-1)), 1] <- all_mice[j]
      RER_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[1]][[all_mice[j]]][["Data"]][["RER"]], n= 48), na.rm = T)
      RER_24h[j + (23*(i-1)), 3] <- "Before"
      RER_24h[j + (23*(i-1)), 4] <- "Chow"
      RER_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      RER_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]]
    }else if (i == 1 & all_mice[j] %in% born_at_30){
      RER_24h[j + (23*(i-1)), 1] <- all_mice[j]
      RER_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[1]][[all_mice[j]]][["Data"]][["RER"]], n= 48), na.rm = T)
      RER_24h[j + (23*(i-1)), 3] <- "Before"
      RER_24h[j + (23*(i-1)), 4] <- "Chow"
      RER_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[1]][[all_mice[j]]][["Group"]]
      RER_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[1]][[all_mice[j]]][["Sex"]]
    }else if (i == 2 & all_mice[j] %in% born_at_22){
      RER_24h[j + (23*(i-1)), 1] <- all_mice[j]
      RER_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[2]][[all_mice[j]]][["Data"]][["RER"]], n= 48), na.rm = T)
      RER_24h[j + (23*(i-1)), 3] <- "After"
      RER_24h[j + (23*(i-1)), 4] <- "HFD"
      RER_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[2]][[all_mice[j]]][["Group"]]
      RER_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[2]][[all_mice[j]]][["Sex"]] 
    }else if (i == 2 & all_mice[j] %in% born_at_30){
      RER_24h[j + (23*(i-1)), 1] <- all_mice[j]
      RER_24h[j + (23*(i-1)), 2] <- mean(tail( CLAMS_HFD[[2]][[all_mice[j]]][["Data"]][["RER"]], n= 48), na.rm = T)
      RER_24h[j + (23*(i-1)), 3] <- "After"
      RER_24h[j + (23*(i-1)), 4] <- "HFD"
      RER_24h[j + (23*(i-1)), 5] <- CLAMS_HFD[[2]][[all_mice[j]]][["Group"]]
      RER_24h[j + (23*(i-1)), 6] <- CLAMS_HFD[[2]][[all_mice[j]]][["Sex"]] 
    }
  }
}


RER_24h$Time <- factor(RER_24h$Time, levels = c("Before", "After"))
RER_24h$Subgroup <- interaction(RER_24h$Group, RER_24h$Sex, sep = "_")

#### Supplementary File 2Cviii
anova(lmer(RER ~ Group*Diet + (1|Animal.ID), RER_24h))

#post-hoc, Supplementary File 2Cx
RER_24h$Subgroup2 <- interaction(RER_24h$Group, RER_24h$Diet, sep = "_")
model <- lmer(RER ~ 0 + Subgroup2 + (1|Animal.ID), RER_24h)

posthoc <- glht(model, linfct=mcp(Subgroup2 = c("`22C_Chow` - `30C_Chow` == 0",
                                                "`22C_HFD` - `30C_HFD` == 0",
                                                "`22C_Chow` - `22C_HFD` == 0",
                                                "`30C_Chow` - `30C_HFD` == 0")))
summary(posthoc)

#and now let's test the RER gain between the two time points (Supplemental File 2Cxii)
RER_gain <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Animal ID", "RER", "Group", "Sex"))))

for (j in 1:23) {
  RER_gain[j, 1] <- RER_24h$Animal.ID[j]
  RER_gain[j, 2] <- RER_24h[RER_24h$Animal.ID == RER_24h$Animal.ID[j] & RER_24h$Time == "After",]$RER - RER_24h[RER_24h$Animal.ID == RER_24h$Animal.ID[j] & RER_24h$Time == "Before",]$RER
  RER_gain[j, 3] <- RER_24h$Group[j]
  RER_gain[j, 4] <- RER_24h$Sex[j]
}

t.test(RER_gain[RER_gain$Group == "22C",]$RER, RER_gain[RER_gain$Group == "30C",]$RER)



