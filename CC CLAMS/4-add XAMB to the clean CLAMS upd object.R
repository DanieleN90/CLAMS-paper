#### let's add the XAMB to the clean_CLAMS_updt

setwd("D:/OneDrive - cumc.columbia.edu/20181017-20181219 Cold challenge/updates from 9-15-2021")
setwd("..")
Weeks <- c("10.17.18", "10.24.18", "10.31.18", "11.7.18", "11.14.18", "11.21.18", "11.28.18", "12.5.18", "12.12.18")
# remember that only the first 7 weeks have the 22C group, and only the last 7 have the 30C group




library(readxl)
library(tidyverse)
library(stringr)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)

raised_at_22 <- c("M744", "F745", "F746", "F747", "F749", "F750", "M751", "M754", "M755", "M756", "M757", "F758")
raised_at_30 <- c("M773", "M775", "F781", "M783", "M791", "M793", "M795", "F786", "F787", "F788", "F792", "F794")
dead_mice <- c("F750", "M754", "M757")
other_mice_to_eliminate <- c("F786", "M793", "M795", "F747")


load("D:/OneDrive - cumc.columbia.edu/20181017-20181219 Cold challenge/updates from 9-15-2021/clean_CLAMS_upd.RData")




for (i in 1:7) {
  for(j in 1:17){
    if(CLAMS[[i]][[j]][["Group"]] == "22C"){
      setwd(paste(getwd(), Weeks[i], sep = "/"))
      week_data <- read.table(dir(pattern = "Oxymax"), header = T, fill = T, sep = ",")
      week_data <- select(week_data, Subject, X.Ambulatory)
      CLAMS[[i]][[j]][["Data"]]$XAMB <- as.numeric(week_data$X.Ambulatory[week_data$Subject == CLAMS[[i]][[j]][["ID"]]])
      setwd("..")
    }else{
      setwd(paste(getwd(), Weeks[i+2], sep = "/"))
      week_data <- read.table(dir(pattern = "Oxymax"), header = T, fill = T, sep = ",")
      week_data <- select(week_data, Subject, X.Ambulatory)
      CLAMS[[i]][[j]][["Data"]]$XAMB <- as.numeric(week_data$X.Ambulatory[week_data$Subject == CLAMS[[i]][[j]][["ID"]]])
      setwd("..")
    }
  }
}

correlations <- c()
for (i in 1:7) {
  for(j in 1:17){
    correlations <- c(correlations, cor(CLAMS[[i]][[j]][["Data"]]$XAMB, CLAMS[[i]][[j]][["Data"]]$XTOT))
  }
}

# the average correlation between XTOT and XAMB withing each subject is 0.98, lowest correlation is 0.9461
setwd("D:/OneDrive - cumc.columbia.edu/20181017-20181219 Cold challenge/updates from 9-15-2021")
save(CLAMS, file = "clean_CLAMS_upd.RData")
