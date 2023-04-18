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

## the EE_data item comes from the script 1C and we are using it as a starting point for the proportions
EE_data_prop <- EE_data

### we need to fix the proportions here
### basically, because the of the mismatch problem of the BMR on week 1 of TN, if we just do the ratio of PAEE, TEF, and BMR over TEE we get results such that the sum is not 1
### the solution for this problem is to scale the proportion
### basically we will do the ratio between the TEE and the sum of the PAEE, TEF, and BMR, and then we will scale those three components by that value before doing the proportions
scale_factors <- data.frame(matrix(ncol=2, nrow = 0, dimnames = list(NULL, c("Mouse", "Value"))))
counter <- 1

for(j in 1:length(CLAMS[[1]])){
  if (CLAMS[[1]][[j]]$Group == "22C") next
  scale_factors[counter, "Mouse"] <- names(CLAMS[[1]])[j]
  scale_factors[counter, "Value"] <- CLAMS[[1]][[j]]$`EE Analysis`$TEE/(CLAMS[[1]][[j]]$`EE Analysis`$BMR + CLAMS[[1]][[j]]$`EE Analysis`$TEF + CLAMS[[1]][[j]]$`EE Analysis`$PAEE)
  counter <- counter + 1
}



for (i in 1:nrow(EE_data_prop)) {
  if(EE_data_prop[i, "Week"] == "Week1" & EE_data_prop[i, "Group"] == "30C" & EE_data_prop[i, "Component"] != "TEE"){
    EE_data_prop[i, 1] <- (EE_data_prop[i, 1]/CLAMS[[EE_data_prop[i,4]]][[as.character(EE_data_prop[i,5])]]$`EE Analysis`$`TEE`)*scale_factors[scale_factors$Mouse == as.character(EE_data_prop[i,5]), "Value"]
  }else{
    EE_data_prop[i, 1] <- EE_data_prop[i, 1]/CLAMS[[EE_data_prop[i,4]]][[as.character(EE_data_prop[i,5])]]$`EE Analysis`$`TEE`
  }
}
EE_data <- EE_data[EE_data$Mouse != "F787",]





data_for_plot_prop <- as.data.frame(group_by(EE_data_prop, Component_Group_Week) %>% dplyr::summarize(m = mean(Value), n=n(), se=sd(Value)/(sqrt(n))))


data_for_plot_prop$Week <- c(rep(c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), 3), c("Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), rep(c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), 6) )
data_for_plot_prop$Component <- c(rep("BMR", 14), rep("CIT", 13), rep("PAEE", 14), rep("TEE", 14), rep("TEF", 14))
data_for_plot_prop$Group <- c(rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 6), rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7))

#need to make an adjusted dataframe for the error bars
data_for_plot_prop_adj <- data_for_plot_prop
#adjust CIT which is the first one from the bottom after BMR
for(i in 0:6){
  data_for_plot_prop_adj[15 + i,2] <- data_for_plot_prop_adj[15 + i,2] + data_for_plot_prop_adj[1 + i,2]
}
for(i in 0:5){
  data_for_plot_prop_adj[22 + i,2] <- data_for_plot_prop_adj[22 + i,2] + data_for_plot_prop_adj[9 + i,2]
}
# now adjust TEF
for(i in 0:6){
  data_for_plot_prop_adj[56 + i,2] <- data_for_plot_prop_adj[56 + i,2] + data_for_plot_prop_adj[15 + i,2]
}
data_for_plot_prop_adj[63,2] <- data_for_plot_prop_adj[63,2] + data_for_plot_prop_adj[8,2]
for(i in 0:5){
  data_for_plot_prop_adj[64 + i,2] <- data_for_plot_prop_adj[64 + i,2] + data_for_plot_prop_adj[22 + i,2]
}
# and finally PAEE
for(i in 0:13){
  data_for_plot_prop_adj[28 + i,2] <- data_for_plot_prop_adj[28 + i,2] + data_for_plot_prop_adj[56 + i,2]
}




ggplot(data_for_plot_prop[data_for_plot_prop$Component != "TEE",], aes(fill=factor(Component, levels = c("PAEE", "TEF", "CIT", "BMR")), y=m, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin= data_for_plot_prop_adj[data_for_plot_prop_adj$Component != "TEE",]$m - se,
                    ymax= data_for_plot_prop_adj[data_for_plot_prop_adj$Component != "TEE",]$m + se ),
                width=.2,                    # Width of the error bars
                position="identity") +
  facet_grid(~ Week)+
  ylim(0,1.05) +
  ylab("Mean Weekly heat (KCal/h)") +
  scale_fill_manual(values = c("green3", "firebrick2", "steelblue4", "grey33"))+ 
  ggtitle("Proportions of energy expenditure") +
  guides(fill=guide_legend(title="Components"))

write.csv(data_for_plot_prop, paste(getwd(), "CC components proportions.csv", sep = "/"))
write.csv(data_for_plot_prop_adj, paste(getwd(), "CC components proportions adjusted for stacked plot.csv", sep = "/"))



#### let's do the stats

EE_data_prop$Value <- as.numeric(EE_data_prop$Value)
EE_data_prop$Component <- as.factor(EE_data_prop$Component)
EE_data_prop$Group <- as.factor(EE_data_prop$Group)
EE_data_prop$Week <- as.factor(EE_data_prop$Week)
EE_data_prop$Mouse <- as.factor(EE_data_prop$Mouse)
EE_data_prop$Group_Week <- as.factor(EE_data_prop$Group_Week)
EE_data_prop$Component_Group_Week <- as.factor(EE_data_prop$Component_Group_Week)
EE_data_prop$Week_Sex <- as.factor(EE_data_prop$Week_Sex)


library(lme4)
library(lmerTest)
library(multcomp)
library(MuMIn)
library(emmeans)

#PAEE day time
PAEE <- EE_data_prop[EE_data_prop$Component == "PAEE",]
PAEE <- PAEE[PAEE$Week != "Week3",]
PAEE <- PAEE[PAEE$Week != "Week6",]


model_PAEE <- lmer(Value ~ Group*Week + (1|Mouse), PAEE)
anova(model_PAEE)


model <- lmer(Value ~ 0 + Group_Week + (1|Mouse), PAEE)
posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week1` - `22C_Week1` == 0",
                                                 "`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week4` - `22C_Week4` == 0",
                                                 "`30C_Week5` - `22C_Week5` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0")))

summary(posthoc)



#TEF day time
TEF <- EE_data_prop[EE_data_prop$Component == "TEF",]
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
CIT <- EE_data_prop[EE_data_prop$Component == "CIT",]
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
