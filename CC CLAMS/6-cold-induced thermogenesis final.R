### this is the update of cold induced thermogenesis,
### I need to recompute BMR just using day time BMR
### on top of that, I will also average day and night for the components, so that we will have only one graph for components and one for proportions
### this whole calculation might create some inconsistency with the fact that for the TN mice I am using only the day time BMR and I apply it to the average (day and night) TEE, and so this could
### end up in either small values of CIT in week 1 of TN mice or even negative values
### to avoid this issue, I will calculate the average BMR for TN and apply it as a flat number


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


mice_names <- read_xlsx("GTT ITT.xlsx")[,1:2]
mice_names$Mouse

born_at_30 <- mice_names$Mouse[mice_names$...1 == "22C born"]
born_at_22 <- mice_names$Mouse[mice_names$...1 == "30C born"]


Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]
Nth.delete_vec<-function(vector, n)vector[-(seq(n,to=length(vector),by=n))]


load("CLAMS_HFD.RData")


### since now we are putting together day and night, we need to recalculate the intercept and all that stuff for components
HFD_Week1 <- list()


for(j in 1:length(CLAMS_HFD[[1]])){
  
  # keep only TN mice
  if (CLAMS_HFD[[1]][[j]]$Group == "30C") next
  
  
  this_mouse <- list()
  mouse_id <- names(CLAMS_HFD[[1]])[j]
  
  Sex <- CLAMS_HFD[[1]][[j]][["Sex"]]
  
  
  time_diff <- as.numeric(CLAMS_HFD[[1]][[j]][["Data"]][["Date.Time"]][2] - CLAMS_HFD[[1]][[j]][["Data"]][["Date.Time"]][1])
  
  heat_day <- tail(CLAMS_HFD[[1]][[j]][["Data"]][["Heat"]], n=48)[tail(CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]], n=48) == "Light"]
  heat_day <- .colMeans(heat_day, 2, length(heat_day) / 2)
  heat_night <- tail(CLAMS_HFD[[1]][[j]][["Data"]][["Heat"]], n=48)[tail(CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]], n=48) == "Dark"]
  heat_night <- .colMeans(heat_night, 2, length(heat_night) / 2)
  
  act_day <- tail(CLAMS_HFD[[1]][[j]][["Data"]][["X.Ambulatory"]], n=48)[tail(CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]], n=48) == "Light"]
  act_day <- .colSums(act_day, 2, length(act_day) / 2)
  act_night <- tail(CLAMS_HFD[[1]][[j]][["Data"]][["X.Ambulatory"]], n=48)[tail(CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]], n=48) == "Dark"]
  act_night <- .colSums(act_night, 2, length(act_night) / 2)
  
  quadratic.model_night <-lm(heat_night ~ act_night + I(act_night^2))
  quadratic.model_day <-lm(heat_day ~ act_day + I(act_day^2))
  
  TEE_night <- mean(heat_night)
  TEE_day <- mean(heat_day)
  
  PAEE_night <- mean(heat_night) - as.numeric(quadratic.model_night$coefficients[1])
  PAEE_day <- mean(heat_day) - as.numeric(quadratic.model_day$coefficients[1])
  
  
  light_hours <- sum(tail(CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]], n=48) == "Light")/2
  dark_hours <- sum(tail(CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]], n=48) == "Dark")/2
  
  ## mice are on chow diet, which we estimated TEF of 0.342
  ##  
  
  food_night <- sum(tail(CLAMS_HFD[[1]][[j]][["Data"]][["Feed.Weight.1"]], n=48)[tail(CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]], n=48) == "Dark"])
  food_day <- sum(tail(CLAMS_HFD[[1]][[j]][["Data"]][["Feed.Weight.1"]], n=48)[tail(CLAMS_HFD[[1]][[j]][["Data"]][["Light.Dark"]], n=48) == "Light"])
  
  
  
  TEF_day <- (food_day*0.342)/(light_hours)
  TEF_night <- (food_night*0.342)/(dark_hours)
  
  rest_of_TEE_day <- TEE_day - PAEE_day - TEF_day
  rest_of_TEE_night <- TEE_night - PAEE_night - TEF_night
  # 
  # 
  
  this_mouse <- list(quadratic.model_day, TEE_day, PAEE_day, TEF_day, rest_of_TEE_day, Sex, quadratic.model_night, TEE_night, PAEE_night, TEF_night, rest_of_TEE_night)
  names(this_mouse) <- c("model_day", "TEE_day", "PAEE_day", "TEF_day", "rest of TEE_day", "Sex", "model_night", "TEE_night", "PAEE_night", "TEF_night", "rest of TEE_night")
  
  HFD_Week1[[j]] <- this_mouse
  names(HFD_Week1)[j] <- mouse_id
}

#clean everything, keep only HFD_Week1
rm(list=setdiff(ls(), "HFD_Week1"))
HFD_Week1 <- HFD_Week1[-which(sapply(HFD_Week1, is.null))]


#### let's plot it to inspect it
first_week <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c("Mouse", 
                                                                     "Data",
                                                                     "Sex",
                                                                     "Value"))))

j <- 1
# let's first loop on the HFD CLAMS
for(i in 1:96){
  first_week[i,"Mouse"] <- names(HFD_Week1)[j]
  first_week[i,"Sex"] <- HFD_Week1[[j]]$Sex
  if(i %in% seq(1,96,8)){
    first_week[i,"Value"] <- HFD_Week1[[j]][2]
    first_week[i,"Data"] <- names(HFD_Week1[[j]])[2]
  }else if (i %in% seq(2,96,8)){
    first_week[i,"Value"] <- HFD_Week1[[j]][3]
    first_week[i,"Data"] <- names(HFD_Week1[[j]])[3]
  }else if (i %in% seq(3,96,8)){
    first_week[i,"Value"] <- HFD_Week1[[j]][4]
    first_week[i,"Data"] <- names(HFD_Week1[[j]])[4]
  }else if (i %in% seq(4,96,8)){
    first_week[i,"Value"] <- HFD_Week1[[j]][5]
    first_week[i,"Data"] <- names(HFD_Week1[[j]])[5]
  }else if (i %in% seq(5,96,8)){
    first_week[i,"Value"] <- HFD_Week1[[j]][8]
    first_week[i,"Data"] <- names(HFD_Week1[[j]])[8]
  }else if (i %in% seq(6,96,8)){
    first_week[i,"Value"] <- HFD_Week1[[j]][9]
    first_week[i,"Data"] <- names(HFD_Week1[[j]])[9]
  }else if (i %in% seq(7,96,8)){
    first_week[i,"Value"] <- HFD_Week1[[j]][10]
    first_week[i,"Data"] <- names(HFD_Week1[[j]])[10]
  }else if (i %in% seq(8,96,8)){
    first_week[i,"Value"] <- HFD_Week1[[j]][11]
    first_week[i,"Data"] <- names(HFD_Week1[[j]])[11]
    j <-  1 + i%/%8
  }
}


first_week$Data <- factor(first_week$Data, levels=c("TEE_day", "TEE_night", "PAEE_day", "PAEE_night", "TEF_day", "TEF_night", "rest of TEE_day", "rest of TEE_night"))

plot1 <- ggplot(first_week, aes(y=Value, x= Data))+
  geom_point(position = position_jitter(width = .1, height = 0), aes(shape=Sex, fill = Data), size = 6, colour = "black") + 
  ylim(0, .3)+
  xlab("Component of TEE") +
  ylab("Heat (Kcal/H)") +
  ggtitle("Energy components in RT mice", subtitle = "last 24 of the first week of HFD CLAMS")+
  stat_summary(aes(y = Value, group = Data, fill = Data), fun = mean, geom="point", size = 3) +
  stat_summary(aes(y = Value, group = Data, fill = Data), fun.data = mean_se,  geom = "errorbar", width = 0.3) + 
  scale_shape_manual(values=c(21,24))+
  guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list()))+
  theme_bw()
plot1


#let's store the values
BEE_day <- mean(first_week$Value[first_week$Data == "rest of TEE_day"])
BEE_sd_day <- sd(first_week$Value[first_week$Data == "rest of TEE_day"])
BEE_n_day <- length(first_week$Value[first_week$Data == "rest of TEE_day"])
BEE_night <- mean(first_week$Value[first_week$Data == "rest of TEE_night"])
BEE_sd_night <- sd(first_week$Value[first_week$Data == "rest of TEE_night"])
BEE_n_night <- length(first_week$Value[first_week$Data == "rest of TEE_night"])
rm(list=setdiff(ls(), c("BEE_day", "BEE_n_day", "BEE_sd_day", "BEE_night", "BEE_n_night", "BEE_sd_night")))

RT_BEE <- list(BEE_day, BEE_n_day, BEE_sd_day, BEE_night, BEE_n_night, BEE_sd_night)
names(RT_BEE) <- c("BEE_day", "BEE_n_day", "BEE_sd_day", "BEE_night", "BEE_n_night", "BEE_sd_night")
rm(list=setdiff(ls(), c("RT_BEE")))


# ok so now we have this for RT mice, now we are going to take the average of the day BMR of week 1 in TN mice



######reusing code from components analysis.R


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


raised_at_22 <- c("M744", "F745", "F746", "F747", "F749", "F750", "M751", "M754", "M755", "M756", "M757", "F758")
raised_at_30 <- c("M773", "M775", "F781", "M783", "M791", "M793", "M795", "F786", "F787", "F788", "F792", "F794")
dead_mice <- c("F750", "M754", "M757")
other_mice_to_eliminate <- c("F786", "M793", "M795", "F747")


raised_at_22 <- setdiff(setdiff(raised_at_22, dead_mice), other_mice_to_eliminate)
raised_at_30 <- setdiff(setdiff(raised_at_30, dead_mice), other_mice_to_eliminate)
Nth.delete<-function(dataframe, n)dataframe[-(seq(n,to=nrow(dataframe),by=n)),]
Nth.delete_vec<-function(vector, n)vector[-(seq(n,to=length(vector),by=n))]


load("~/clean_CLAMS_upd.RData")

BMR_day <- c()
for(j in 1:length(CLAMS[[1]])){
  if (CLAMS[[1]][[j]]$Group == "22C") next
  BMR_day <- c(BMR_day, CLAMS[[1]][[j]][["EE Analysis"]][["rest of TEE day"]])
}
TN_BMR <- mean(BMR_day)



#we are going to loop through the CLAMS item, and for each mouse we are going to modify the EE analysis list to now account for BEE and cold-induced thermogenesis

for (i in 1:7) {
  for(j in 1:length(CLAMS[[i]])){
    if(CLAMS[[i]][[j]]$Group == "30C"){
      if(i == 1){
        BMR <- TN_BMR
        mouse_name <- names(CLAMS[[i]])[j]
        heat <- c(CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]], CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]])
        act <- c(CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]], CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]])
        TEE <- mean(heat)
        quadratic.model <-lm(heat ~ act + I(act^2))
        PAEE <- mean(heat) - as.numeric(quadratic.model$coefficients[1])
        TEF <- (tail(CLAMS[[i]][[j]][["Data"]][["NEW.FEED.ACC"]], n=1)*0.342)/(length(heat))
        this_mouse <- list(TEE, PAEE, TEF, BMR)
        names(this_mouse) <- c("TEE", "PAEE", "TEF", "BMR")
        CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
        CLAMS[[i]][[j]][["EE Analysis"]] <- this_mouse

      }else if (i >1){
        mouse_name <- names(CLAMS[[i]])[j]
        heat <- c(CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]], CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]])
        act <- c(CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]], CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]])
        TEE <- mean(heat)
        quadratic.model <-lm(heat ~ act + I(act^2))
        PAEE <- mean(heat) - as.numeric(quadratic.model$coefficients[1])
        TEF <- (tail(CLAMS[[i]][[j]][["Data"]][["NEW.FEED.ACC"]], n=1)*0.342)/(length(heat))
        BMR <- TN_BMR
        CIT <- TEE - PAEE - TEF - BMR
        this_mouse <- list(TEE, PAEE, TEF, BMR, CIT)
        names(this_mouse) <- c("TEE", "PAEE", "TEF", "BMR", "CIT")
        CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
        CLAMS[[i]][[j]][["EE Analysis"]] <- this_mouse
        
      }
    }else if (CLAMS[[i]][[j]]$Group == "22C"){
      mouse_name <- names(CLAMS[[i]])[j]
      heat <- c(CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_heat"]], CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_heat"]])
      act <- c(CLAMS[[i]][[j]][["EE Analysis"]][["model day"]][["model"]][["day_act"]], CLAMS[[i]][[j]][["EE Analysis"]][["model night"]][["model"]][["night_act"]])
      TEE <- mean(heat)
      quadratic.model <-lm(heat ~ act + I(act^2))
      PAEE <- mean(heat) - as.numeric(quadratic.model$coefficients[1])
      TEF <- (tail(CLAMS[[i]][[j]][["Data"]][["NEW.FEED.ACC"]], n=1)*0.342)/(length(heat))
      BMR <- RT_BEE[["BEE_day"]]
      CIT <- TEE - PAEE - TEF - BMR
      this_mouse <- list(TEE, PAEE, TEF, BMR, CIT)
      names(this_mouse) <- c("TEE", "PAEE", "TEF", "BMR", "CIT")
      CLAMS[[i]][[j]][["EE Analysis"]] <- NULL
      CLAMS[[i]][[j]][["EE Analysis"]] <- this_mouse
    }
  }
}


### let's save this new item perhaps under a new name to be safe
save(CLAMS, file = "clean_CLAMS_upd+CITv2.RData")


############### now we can proceed with the graphs and the stats
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





data_for_plot <- as.data.frame(group_by(EE_data, Component_Group_Week) %>% dplyr::summarize(m = mean(Value), n=n(), se=sd(Value)/(sqrt(n))))


data_for_plot$Week <- c(rep(c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), 3), c("Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), rep(c("Week1", "Week2", "Week3", "Week4", "Week5", "Week6", "Week7"), 6) )
data_for_plot$Component <- c(rep("BMR", 14), rep("CIT", 13), rep("PAEE", 14), rep("TEE", 14), rep("TEF", 14))
data_for_plot$Group <- c(rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 6), rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7), rep("22C", 7), rep("30C", 7))

#need to make an adjusted dataframe for the error bars
data_for_plot_adj <- data_for_plot
#adjust CIT which is the first one from the bottom after BMR
for(i in 0:6){
  data_for_plot_adj[15 + i,2] <- data_for_plot_adj[15 + i,2] + data_for_plot_adj[1 + i,2]
}
for(i in 0:5){
  data_for_plot_adj[22 + i,2] <- data_for_plot_adj[22 + i,2] + data_for_plot_adj[9 + i,2]
}
# now adjust TEF
for(i in 0:6){
  data_for_plot_adj[56 + i,2] <- data_for_plot_adj[56 + i,2] + data_for_plot_adj[15 + i,2]
}
data_for_plot_adj[63,2] <- data_for_plot_adj[63,2] + data_for_plot_adj[8,2]
for(i in 0:5){
  data_for_plot_adj[64 + i,2] <- data_for_plot_adj[64 + i,2] + data_for_plot_adj[22 + i,2]
}
# and finally PAEE
for(i in 0:13){
  data_for_plot_adj[28 + i,2] <- data_for_plot_adj[28 + i,2] + data_for_plot_adj[56 + i,2]
}

#let's rename BMR to BMR to keep it consistent 


ggplot(data_for_plot[data_for_plot$Component != "TEE",], aes(fill=factor(Component, levels = c("PAEE", "TEF", "CIT", "BMR")), y=m, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin= data_for_plot_adj[data_for_plot_adj$Component != "TEE",]$m - se,
                    ymax= data_for_plot_adj[data_for_plot_adj$Component != "TEE",]$m + se ),
                width=.2,                    # Width of the error bars
                position="identity") +
  facet_grid(~ Week)+
  ylim(0,0.7) +
  ylab("Mean Weekly heat (KCal/h)") +
  scale_fill_manual(values = c("green3", "firebrick2", "steelblue4", "grey33"))+ 
  ggtitle("Components of energy expenditure") +
  guides(fill=guide_legend(title="Components"))

write.csv(data_for_plot, paste(getwd(), "CC components.csv", sep = "/"))
write.csv(data_for_plot_adj, paste(getwd(), "CC components adjusted for stacked plot.csv", sep = "/"))


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

posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week4` - `22C_Week4` == 0",
                                                 "`30C_Week5` - `22C_Week5` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0",
                                                 "`30C_Week7` - `30C_Week5` == 0",
                                                 "`22C_Week7` - `22C_Week5` == 0")))


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

# posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week7` - `30C_Week5` == 0",
#                                                  "`22C_Week7` - `22C_Week5` == 0")))

posthoc <- glht(model, linfct=mcp(Group_Week = c("`30C_Week2` - `22C_Week2` == 0",
                                                 "`30C_Week4` - `22C_Week4` == 0",
                                                 "`30C_Week5` - `22C_Week5` == 0",
                                                 "`30C_Week7` - `22C_Week7` == 0",
                                                 "`30C_Week7` - `30C_Week5` == 0",
                                                 "`22C_Week7` - `22C_Week5` == 0")))

summary(posthoc)










#proportions
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















