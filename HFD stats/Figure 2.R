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
############ 2B  #############
##############################

anova(lmer(Bodyweight ~ Group*Sex + Week + (1|Mouse), current_data))

### post-hoc analysis for bodyweight, first looking at the contrast between sexes across weeks
model <- lmer(Bodyweight ~ 0 + Subgroup2 + (1|Mouse), current_data)
posthoc <- glht(model, linfct=mcp(Subgroup2 = c("`F_Week1` - `M_Week1` == 0",
                                                "`F_Week2` - `M_Week2` == 0",
                                                "`F_Week7` - `M_Week7` == 0")))
summary(posthoc)

#and now the contrast between rearing temperatures within each sex, divided per week 
# 
model <- lmer(Bodyweight ~ 0 + Subgroup3 + (1|Mouse), current_data)

posthoc <- glht(model, linfct=mcp(Subgroup3 = c("`30C_M_Week1` - `22C_M_Week1` == 0",
                                                "`30C_F_Week1` - `22C_F_Week1` == 0",
                                                "`30C_M_Week2` - `22C_M_Week2` == 0",
                                                "`30C_F_Week2` - `22C_F_Week2` == 0",
                                                "`30C_M_Week7` - `22C_M_Week7` == 0",
                                                "`30C_F_Week7` - `22C_F_Week7` == 0")))
summary(posthoc)


##############################
############ 2C  #############
##############################

### now let's do lean mass
anova(lmer(Lean_Mass ~ Group*Sex + Week + (1|Mouse), current_data))

# post-hoc for lean mass, looking at contrasts between sexes per each week 
model <- lmer(Lean_Mass ~ 0 + Subgroup2 + (1|Mouse), current_data)
posthoc <- glht(model, linfct=mcp(Subgroup2 = c("`F_Week1` - `M_Week1` == 0",
                                                "`F_Week2` - `M_Week2` == 0",
                                                "`F_Week7` - `M_Week7` == 0")))
summary(posthoc)



##############################
############ 2D  #############
##############################

### now adiposity 
anova(lmer(Adiposity ~ Group*Sex + Week + (1|Mouse), current_data))

# post-hoc for lean mass, looking at contrasts between sexes per each week 
model <- lmer(Adiposity ~ 0 + Subgroup2 + (1|Mouse), current_data)
posthoc <- glht(model, linfct=mcp(Subgroup2 = c("`F_Week1` - `M_Week1` == 0",
                                                "`F_Week2` - `M_Week2` == 0",
                                                "`F_Week7` - `M_Week7` == 0")))
summary(posthoc)


##############################
############ 2E  #############
##############################

###### GTT
GTT <- read_xlsx("~/GTT and ITT rearranged.xlsx", col_names = T, skip = 0)

GTT$Subgroup <- interaction(GTT$Group, GTT$Sex, sep = "_")

GTT$Mouse <- as.character(GTT$Mouse)
GTT$Group <- as.factor(GTT$Group)
GTT$time <- as.factor(GTT$time)
GTT$value <- as.numeric(GTT$value)




#Welch t-tests
t.test(GTT[GTT$time == 0 & GTT$Group == "30C",]$value, GTT[GTT$time == 0 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 20 & GTT$Group == "30C",]$value, GTT[GTT$time == 20 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 40 & GTT$Group == "30C",]$value, GTT[GTT$time == 40 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 60 & GTT$Group == "30C",]$value, GTT[GTT$time == 60 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 90 & GTT$Group == "30C",]$value, GTT[GTT$time == 90 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 120 & GTT$Group == "30C",]$value, GTT[GTT$time == 120 & GTT$Group == "22C",]$value)



##############################
############ 2F  #############
##############################
#ANOVA
anova(lmer(value ~ Group +Sex + time + (1|Mouse), GTT))



##############################
############ 2G  #############
##############################

###### ITT
ITT <- read_xlsx("~/GTT and ITT rearranged.xlsx", col_names = T, skip = 0, sheet = "ITT")

ITT$Subgroup <- interaction(ITT$Group, ITT$Sex, sep = "_")

ITT$Mouse <- as.character(ITT$Mouse)
ITT$Group <- as.factor(ITT$Group)
ITT$time <- as.factor(ITT$time)
ITT$value <- as.numeric(ITT$value)


#ANOVA
anova(lmer(value ~ Group +Sex + time + (1|Mouse), ITT))

