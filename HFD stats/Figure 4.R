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
############ 4C  #############
##############################

anova(lmer(LA ~ Group*Sex + Week + (1|Mouse), current_data))



##############################
############ 4E  #############
##############################

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


#### PAEE
PAEE <- EE_data[EE_data$Component == "PAEE",]

model_PAEE <- lmer(Value ~ Group*Week + (1|Mouse), PAEE)
anova(model_PAEE)

#post-hoc analysis,
model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), PAEE)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week1` - `22C_Week1` == 0",
                                                 "`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))
summary(posthoc)



#### TEF 
TEF <- EE_data[EE_data$Component == "TEF",]

model_TEF <- lmer(Value ~ Group*Week + (1|Mouse), TEF)
anova(model_TEF)


#### DIT, here called HFD-IT (old nomenclature)
HFD_IT <- EE_data[EE_data$Component == "HFD_IT" & EE_data$Week != "Week1",]

model_HFD_IT <- lmer(Value ~ Group*Week + (1|Mouse), HFD_IT)
anova(model_HFD_IT)




##############################
############ 4F  #############
##############################

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


#### PAEE 
PAEE <- EE_data_prop[EE_data_prop$Component == "PAEE",]



model_PAEE <- lmer(Value ~ Group*Week + (1|Mouse), PAEE)
anova(model_PAEE)

#post-hoc analysis
model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), PAEE)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week1` - `22C_Week1` == 0",
                                                 "`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))

summary(posthoc)



#### TEF 
TEF <- EE_data_prop[EE_data_prop$Component == "TEF",]

model_TEF <- lmer(Value ~ Group*Week + (1|Mouse), TEF)
anova(model_TEF)


#### DIT, here called HFD-IT (old nomenclature)
HFD_IT <- EE_data_prop[EE_data_prop$Component == "HFD_IT" & EE_data_prop$Week != "Week1",]

model_HFD_IT <- lmer(Value ~ Group*Week + (1|Mouse), HFD_IT)
anova(model_HFD_IT)

#post-hoc analysis
model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), HFD_IT)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))

summary(posthoc)



