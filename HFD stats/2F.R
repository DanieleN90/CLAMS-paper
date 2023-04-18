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

#proportions
EE_data_prop <- EE_data

### we need to fix the proportions here
### basically, because the of the mismatch problem of the BMR on week 1 of TN, if we just do the ratio of PAEE, TEF, and BMR over TEE we get results such that the sum is not 1
### the solution for this problem is to scale the proportion
### basically we will do the ratio between the TEE and the sum of the PAEE, TEF, and BMR, and then we will scale those three components by that value before doing the proportions
scale_factors <- data.frame(matrix(ncol=2, nrow = 0, dimnames = list(NULL, c("Mouse", "Value"))))
counter <- 1

for(j in 1:length(CLAMS_HFD[[1]])){
  scale_factors[counter, "Mouse"] <- names(CLAMS_HFD[[1]])[j]
  scale_factors[counter, "Value"] <- CLAMS_HFD[[1]][[j]]$`EE Analysis`$TEE/(CLAMS_HFD[[1]][[j]]$`EE Analysis`$BMR + CLAMS_HFD[[1]][[j]]$`EE Analysis`$TEF + CLAMS_HFD[[1]][[j]]$`EE Analysis`$PAEE)
  counter <- counter + 1
}



for (i in 1:nrow(EE_data_prop)) {
  if(EE_data_prop[i, "Week"] == "Week1" & EE_data_prop[i, "Component"] != "TEE"){
    EE_data_prop[i, 1] <- (EE_data_prop[i, 1]/CLAMS_HFD[[EE_data_prop[i,4]]][[as.character(EE_data_prop[i,5])]]$`EE Analysis`$`TEE`)*scale_factors[scale_factors$Mouse == as.character(EE_data_prop[i,5]), "Value"]
  }else{
    EE_data_prop[i, 1] <- EE_data_prop[i, 1]/CLAMS_HFD[[EE_data_prop[i,4]]][[as.character(EE_data_prop[i,5])]]$`EE Analysis`$`TEE`
  }
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

# we want to keep only week 2,3 and 8, which in this case are week 1,2 and 7
EE_data_prop <- EE_data_prop[EE_data_prop$Week == "Week1" | EE_data_prop$Week == "Week2" | EE_data_prop$Week == "Week7",]


#### PAEE Supplementary File 2Fi
PAEE <- EE_data_prop[EE_data_prop$Component == "PAEE",]



model_PAEE <- lmer(Value ~ Group*Week + (1|Mouse), PAEE)
anova(model_PAEE)

#post-hoc analysis, Supplementary File 2Fiv
model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), PAEE)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week1` - `22C_Week1` == 0",
                                                 "`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))

summary(posthoc)



#### TEF Supplementary File 2Fii
TEF <- EE_data_prop[EE_data_prop$Component == "TEF",]

model_TEF <- lmer(Value ~ Group*Week + (1|Mouse), TEF)
anova(model_TEF)


#### DIT Supplementary File 2Fiii, here called HFD-IT (old nomenclature)
HFD_IT <- EE_data_prop[EE_data_prop$Component == "HFD_IT" & EE_data_prop$Week != "Week1",]

model_HFD_IT <- lmer(Value ~ Group*Week + (1|Mouse), HFD_IT)
anova(model_HFD_IT)

#post-hoc analysis, Supplementary File 2Fv
model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), HFD_IT)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))

summary(posthoc)





