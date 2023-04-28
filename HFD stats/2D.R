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