#### let's add the XAMB to the CLAMS object

# setwd("~")

#these are the names of the folders with the weeks
Weeks <- c("10.17.18", "10.24.18", "10.31.18", "11.7.18", "11.14.18", "11.21.18", "11.28.18", "12.5.18", "12.12.18")
# remember that only the first 7 weeks have the 22C group, and only the last 7 have the 30C group

# in the previous scripts we used the files from the CalR folder, but now we are going to use the files from the folder which were named
# in chronological order
# the content of the files is identical, it just changes how the files are organized, and it makes it easier for me to extract data from them


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


load("~/clean_CLAMS_upd.RData")




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
# the average correlation between XTOT and XAMB withing each subject is 0.98, lowest correlation is 0.9461
setwd("D:/OneDrive - cumc.columbia.edu/20181017-20181219 Cold challenge/updates from 9-15-2021")
save(CLAMS, file = "clean_CLAMS_upd.RData")
