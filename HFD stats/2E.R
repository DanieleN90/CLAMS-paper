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

setwd("~")

mice_names <- read_xlsx("GTT ITT.xlsx")[,1:2]
mice_names$Mouse

born_at_30 <- mice_names$Mouse[mice_names$...1 == "22C born"]
born_at_22 <- mice_names$Mouse[mice_names$...1 == "30C born"]

load("~/CLAMS_HFD_HFD-IT.RData")

### let's collect all the data for analysis in the EE_data dataframe
EE_data <- data.frame(matrix(ncol=6, nrow = 0, dimnames = list(NULL, c("Value", "Component", "Group", "Week", "Mouse", "Sex"))))
counter <- 1
for (i in 1:7) {
  for(j in 1:length(CLAMS_HFD[[i]])){
    if( i == 1){
      for(k in 0:3){
        EE_data[counter + k, 1] <- CLAMS_HFD[[i]][[j]][["EE Analysis"]][[1 + k]]
        EE_data[counter + k, 2] <- names(CLAMS_HFD[[i]][[j]][["EE Analysis"]])[1 + k]
        EE_data[counter + k, 3] <- CLAMS_HFD[[i]][[j]][["Group"]]
        EE_data[counter + k, 4] <- names(CLAMS_HFD)[i]
        EE_data[counter + k, 5] <- names(CLAMS_HFD[[i]])[j]
        EE_data[counter + k, 6] <- CLAMS_HFD[[i]][[j]][["Sex"]]
      }
    }else{
      for(k in 0:4){
        EE_data[counter + k, 1] <- CLAMS_HFD[[i]][[j]][["EE Analysis"]][[1 + k]]
        EE_data[counter + k, 2] <- names(CLAMS_HFD[[i]][[j]][["EE Analysis"]])[1 + k]
        EE_data[counter + k, 3] <- CLAMS_HFD[[i]][[j]][["Group"]]
        EE_data[counter + k, 4] <- names(CLAMS_HFD)[i]
        EE_data[counter + k, 5] <- names(CLAMS_HFD[[i]])[j]
        EE_data[counter + k, 6] <- CLAMS_HFD[[i]][[j]][["Sex"]]
      }
    }
    counter <- counter + k + 1
  }
}


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


# we want to keep only week 2,3 and 8, which in this case are week 1,2 and 7
EE_data <- EE_data[EE_data$Week == "Week1" | EE_data$Week == "Week2" | EE_data$Week == "Week7",]


#### PAEE Supplementary File 2Ei
PAEE <- EE_data[EE_data$Component == "PAEE",]

model_PAEE <- lmer(Value ~ Group*Week + (1|Mouse), PAEE)
anova(model_PAEE)

#post-hoc analysis, Supplementary File 2Eiv
model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), PAEE)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week1` - `22C_Week1` == 0",
                                                 "`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))
summary(posthoc)



#### TEF Supplementary File 2Eii
TEF <- EE_data[EE_data$Component == "TEF",]

model_TEF <- lmer(Value ~ Group*Week + (1|Mouse), TEF)
anova(model_TEF)


#### DIT Supplementary File 2Eiii, here called HFD-IT (old nomenclature)
HFD_IT <- EE_data[EE_data$Component == "HFD_IT" & EE_data$Week != "Week1",]

model_HFD_IT <- lmer(Value ~ Group*Week + (1|Mouse), HFD_IT)
anova(model_HFD_IT)