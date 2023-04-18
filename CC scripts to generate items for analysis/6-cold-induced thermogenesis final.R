## to execute this code, it is important to have also worked on the high fat diet code
## and to have generated the "CLAMS_HFD_clean.RData" item

### this is the update of cold induced thermogenesis,
### I need to recompute Basal Metaboilic Rate (BMR) just using day time BMR
### on top of that, I will also average day and night for the components, so that we will have only one graph for components and one for proportions
### this whole calculation might create some inconsistency with the fact that for the Thermoneutral (TN, aka Tr30) mice I am using only the day time BMR and I apply it to the average (day and night) TEE, and so this could
### end up in either small values of Cold-Induced Thermogenesis (CIT) in week 1 of TN mice or even negative values
### to avoid this issue, I will calculate the average BMR for Tr30 and apply it as a flat number


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

load("CLAMS_HFD_clean.RData")

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


# ok so now we have this for Tr22 mice, now we are going to take the average of the day BMR of week 1 in Ts30 mice

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