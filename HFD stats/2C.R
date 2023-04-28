setwd("~/HFD challenge CLAMS")

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


#### Supplemental File 2Ci
anova(lmer(Heat ~ Group*Diet + (1|Animal.ID), Heat_24h))

# post-hoc for the analysis from above (Supplemental File 2Ciii)
Heat_24h$Subgroup2 <- interaction(Heat_24h$Group, Heat_24h$Diet, sep = "_")
model <- lmer(Heat ~ 0 + Subgroup2 + (1|Animal.ID), Heat_24h)

posthoc <- glht(model, linfct=mcp(Subgroup2 = c("`22C_Chow` - `30C_Chow` == 0",
                                                "`22C_HFD` - `30C_HFD` == 0",
                                                "`22C_Chow` - `22C_HFD` == 0",
                                                "`30C_Chow` - `30C_HFD` == 0")))
summary(posthoc)

##and now let's test the heat gain between the two time points (Supplemental File 2Cv)
heat_gain <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Animal ID", "Heat", "Group", "Sex"))))

for (j in 1:23) {
  heat_gain[j, 1] <- Heat_24h$Animal.ID[j]
  heat_gain[j, 2] <- Heat_24h[Heat_24h$Animal.ID == Heat_24h$Animal.ID[j] & Heat_24h$Time == "After",]$Heat - Heat_24h[Heat_24h$Animal.ID == Heat_24h$Animal.ID[j] & Heat_24h$Time == "Before",]$Heat
  heat_gain[j, 3] <- Heat_24h$Group[j]
  heat_gain[j, 4] <- Heat_24h$Sex[j]
}

t.test(heat_gain[heat_gain$Group == "22C",]$Heat, heat_gain[heat_gain$Group == "30C",]$Heat)



### now RER
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




##########################################################################
# and now let's do the same but comparing the 24h before HFD and last 24h of the first week under HFD


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




### now RER

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












