#this analysis is to inspired by the ANCOVA from Alexander Banks, also present in CalR with the regression plots (Mina et al., 2018 Cell Metabolism)
setwd("~")

library(readxl)
library(tidyverse)
library(stringr)
library(ggpubr)
library(viridis)
library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(lubridate)
library(fpc)
library(ggprism)

library(lme4)
library(lmerTest)
library(multcomp)


raised_at_22 <- c("M744", "F745", "F746", "F747", "F749", "F750", "M751", "M754", "M755", "M756", "M757", "F758")
raised_at_30 <- c("M773", "M775", "F781", "M783", "M791", "M793", "M795", "F786", "F787", "F788", "F792", "F794")
dead_mice <- c("F750", "M754", "M757")
other_mice_to_eliminate <- c("F786", "M793", "M795", "F747")


raised_at_22 <- setdiff(setdiff(raised_at_22, dead_mice), other_mice_to_eliminate)
raised_at_30 <- setdiff(setdiff(raised_at_30, dead_mice), other_mice_to_eliminate)
Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]
Nth.delete_vec<-function(vector, n)vector[-(seq(n,to=length(vector),by=n))]


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