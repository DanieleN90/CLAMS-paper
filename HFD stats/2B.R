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


mice_names <- read_xlsx("GTT ITT.xlsx")[,1:2]
mice_names$Mouse

born_at_30 <- mice_names$Mouse[mice_names$...1 == "22C born"]
born_at_22 <- mice_names$Mouse[mice_names$...1 == "30C born"]

load("~/CLAMS_HFD.RData")

###### GTT
GTT <- read_xlsx("~/GTT and ITT rearranged.xlsx", col_names = T, skip = 0)

GTT$Subgroup <- interaction(GTT$Group, GTT$Sex, sep = "_")

GTT$Mouse <- as.character(GTT$Mouse)
GTT$Group <- as.factor(GTT$Group)
GTT$time <- as.factor(GTT$time)
GTT$value <- as.numeric(GTT$value)


# just plotting the data for visualization 
temp_plot <- ggplot(GTT, aes(x=time, y=value, color=Group)) +
  geom_point(position = position_jitter(width = .1, height = 0), aes(shape=Sex, fill = Group), size = 6, colour = "black") + 
  xlab("Time (min)") + 
  ylab("Glycemia (mg/dL)") +
  ggtitle("GTT on HFD CLAMS")+
  stat_summary(aes(y = value, group = Group), fun = mean, geom="line", size = 2) +
  stat_summary(aes(y = value, group = Subgroup, linetype = Sex), fun = mean, geom="line", size = 1) +
  scale_color_manual(values=c("royalblue4", "darkred"),  guide = "none")+
  scale_fill_manual(values=c("steelblue", "darksalmon"), guide ="none")+
  scale_shape_manual(values=c(21,24))+
  ylim(0, 500)+
  guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list()))+
  theme_bw()
temp_plot

#ANOVA, Supplemental File 2Bi
anova(lmer(value ~ Group +Sex + time + (1|Mouse), GTT))

#Welch t-tests, Supplemental FIle 2Biii
t.test(GTT[GTT$time == 0 & GTT$Group == "30C",]$value, GTT[GTT$time == 0 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 20 & GTT$Group == "30C",]$value, GTT[GTT$time == 20 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 40 & GTT$Group == "30C",]$value, GTT[GTT$time == 40 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 60 & GTT$Group == "30C",]$value, GTT[GTT$time == 60 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 90 & GTT$Group == "30C",]$value, GTT[GTT$time == 90 & GTT$Group == "22C",]$value)
t.test(GTT[GTT$time == 120 & GTT$Group == "30C",]$value, GTT[GTT$time == 120 & GTT$Group == "22C",]$value)

##plot just time point 0 to show the difference
temp_plot <- ggplot(GTT[GTT$time == 0,], aes(x=Group, y=value, color=Group)) +
  geom_point(position = position_jitter(width = .1, height = 0), aes(shape=Sex, fill = Group), size = 6, colour = "black") + 
  xlab("Group") + 
  ylab("Glycemia (mg/dL)") +
  ggtitle("Resting Glycemia in HFD CLAMS", subtitle = "Glycemia at time 0")+
  stat_summary(aes(y = value, group = Group), fun = mean, size = 2) +
  scale_color_manual(values=c("royalblue4", "darkred"),  guide = "none")+
  scale_fill_manual(values=c("steelblue", "darksalmon"), guide ="none")+
  scale_shape_manual(values=c(21,24))+
  ylim(0, 100)+
  guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list()))+
  theme_bw()
temp_plot










###### ITT
ITT <- read_xlsx("~/GTT and ITT rearranged.xlsx", col_names = T, skip = 0, sheet = "ITT")

ITT$Subgroup <- interaction(ITT$Group, ITT$Sex, sep = "_")

ITT$Mouse <- as.character(ITT$Mouse)
ITT$Group <- as.factor(ITT$Group)
ITT$time <- as.factor(ITT$time)
ITT$value <- as.numeric(ITT$value)


# just plotting the data for visualization 
temp_plot <- ggplot(ITT, aes(x=time, y=value, color=Group)) +
  geom_point(position = position_jitter(width = .1, height = 0), aes(shape=Sex, fill = Group), size = 6, colour = "black") + 
  xlab("Time (min)") + 
  ylab("Glycemia (mg/dL)") +
  ggtitle("ITT on HFD CLAMS")+
  stat_summary(aes(y = value, group = Group), fun = mean, geom="line", size = 2) +
  stat_summary(aes(y = value, group = Subgroup, linetype = Sex), fun = mean, geom="line", size = 1) +
  scale_color_manual(values=c("royalblue4", "darkred"),  guide = "none")+
  scale_fill_manual(values=c("steelblue", "darksalmon"), guide ="none")+
  scale_shape_manual(values=c(21,24))+
  ylim(0, 200)+
  guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list()))+
  theme_bw()
temp_plot

#ANOVA, Supplemental File 2Bii
anova(lmer(value ~ Group +Sex + time + (1|Mouse), ITT))

#Welch t-tests, Supplemental FIle 2Biv
t.test(ITT[ITT$time == 0 & ITT$Group == "30C",]$value, ITT[ITT$time == 0 & ITT$Group == "22C",]$value)
t.test(ITT[ITT$time == 20 & ITT$Group == "30C",]$value, ITT[ITT$time == 20 & ITT$Group == "22C",]$value)
t.test(ITT[ITT$time == 40 & ITT$Group == "30C",]$value, ITT[ITT$time == 40 & ITT$Group == "22C",]$value)
t.test(ITT[ITT$time == 60 & ITT$Group == "30C",]$value, ITT[ITT$time == 60 & ITT$Group == "22C",]$value)
t.test(ITT[ITT$time == 90 & ITT$Group == "30C",]$value, ITT[ITT$time == 90 & ITT$Group == "22C",]$value)
t.test(ITT[ITT$time == 120 & ITT$Group == "30C",]$value, ITT[ITT$time == 120 & ITT$Group == "22C",]$value)

##plot just time point 0 to show the difference
temp_plot <- ggplot(ITT[ITT$time == 0,], aes(x=Group, y=value, color=Group)) +
  geom_point(position = position_jitter(width = .1, height = 0), aes(shape=Sex, fill = Group), size = 6, colour = "black") + 
  xlab("Group") + 
  ylab("Glycemia (mg/dL)") +
  ggtitle("Resting Glycemia in HFD CLAMS", subtitle = "Glycemia at time 0")+
  stat_summary(aes(y = value, group = Group), fun = mean, size = 2) +
  scale_color_manual(values=c("royalblue4", "darkred"),  guide = "none")+
  scale_fill_manual(values=c("steelblue", "darksalmon"), guide ="none")+
  scale_shape_manual(values=c(21,24))+
  ylim(0, 250)+
  guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list()))+
  theme_bw()
temp_plot

