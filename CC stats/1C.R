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

################ remember we need to remove F787 because her food scale was broken
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


library(lme4)
library(lmerTest)
library(multcomp)
library(MuMIn)
library(emmeans)

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



#TEF day time
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


#CIT day time
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