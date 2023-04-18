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
    current_data[counter, 1] <- CLAMS[[i]][[j]][["ID"]]
    current_data[counter, 2] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]], CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]]))
    current_data[counter, 3] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]], CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]]))
    current_data[counter, 4] <- mean(c(CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]][CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]] <= 100], 
                                       CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]][CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]] <= 100]))
    current_data[counter, 5] <- tail(CLAMS[[i]][[j]][["Data"]][["NEW.FEED.ACC"]], n=1)
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


##########################################################################################
## now locomotor activity

anova(lmer(LA ~ Group*Sex*Week + (1|Mouse), current_data))

current_data$Interaction <- interaction(current_data$Sex, current_data$Week, sep = "_")

model <- lmer(LA ~ 0 + Interaction + (1|Mouse), current_data)
posthoc <- glht(model, linfct=mcp(Interaction = c("`F_Week1` - `M_Week1` == 0",
                                                  "`F_Week2` - `M_Week2` == 0",
                                                  "`F_Week4` - `M_Week4` == 0",
                                                  "`F_Week5` - `M_Week5` == 0",
                                                  "`F_Week7` - `M_Week7` == 0")))
summary(posthoc)

##########################################################################################
## now REE

anova(lmer(REE ~ Group*Sex*Week + (1|Mouse), current_data))

current_data$Interaction <- interaction(current_data$Group, current_data$Week, sep = "_")

model <- lmer(REE ~ 0 + Interaction + (1|Mouse), current_data)
posthoc <- glht(model, linfct=mcp(Interaction = c("`30C_Week1` - `22C_Week1` == 0",
                                                  "`30C_Week2` - `22C_Week2` == 0",
                                                  "`30C_Week4` - `22C_Week4` == 0",
                                                  "`30C_Week5` - `22C_Week5` == 0",
                                                  "`30C_Week7` - `22C_Week7` == 0")))
summary(posthoc)


##########################################################################################
## now bodyweight

anova(lmer(Bodyweight ~ Group*Sex*Week + (1|Mouse), current_data))

current_data$Interaction <- interaction(current_data$Sex, current_data$Week, sep = "_")

model <- lmer(Bodyweight ~ 0 + Interaction + (1|Mouse), current_data)
posthoc <- glht(model, linfct=mcp(Interaction = c("`F_Week1` - `M_Week1` == 0",
                                                  "`F_Week2` - `M_Week2` == 0",
                                                  "`F_Week4` - `M_Week4` == 0",
                                                  "`F_Week5` - `M_Week5` == 0",
                                                  "`F_Week7` - `M_Week7` == 0")))
summary(posthoc)

##########################################################################################
## now Food intake

anova(lmer(FI ~ Group*Sex*Week + (1|Mouse), current_data))

current_data$Interaction <- interaction(current_data$Sex, current_data$Week, sep = "_")

model <- lmer(FI ~ 0 + Interaction + (1|Mouse), current_data)
posthoc <- glht(model, linfct=mcp(Interaction = c("`F_Week1` - `M_Week1` == 0",
                                                  "`F_Week2` - `M_Week2` == 0",
                                                  "`F_Week4` - `M_Week4` == 0",
                                                  "`F_Week5` - `M_Week5` == 0",
                                                  "`F_Week7` - `M_Week7` == 0")))
summary(posthoc)
