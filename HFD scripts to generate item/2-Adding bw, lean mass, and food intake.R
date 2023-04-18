setwd("~")

library(readxl)
library(tidyverse)
library(stringr)
library(ggpubr)
library(viridis)
library(ggplot2)
library(dplyr)

mice_names <- read_xlsx("GTT ITT.xlsx")[,1:2]
mice_names$Mouse

born_at_30 <- mice_names$Mouse[mice_names$...1 == "22C born"]
born_at_22 <- mice_names$Mouse[mice_names$...1 == "30C born"]

load("~/CLAMS_HFD.RData")

# let's load the weight of 30C mice
mice_at_30_bw <- read_excel("201904-06 CLAMS HFD challenge record.xlsx", 
                            sheet = "30C born, left side CLAMS - MRI", 
                            skip = 1, n_max = 12)
#fix the column names
names(mice_at_30_bw) <- c("Channel", "ID", "Week0", "Week1", "Week2", "Week3", "Week4", "Week5", "Week6","Week7", "Week8")
# mice_at_30_bw <- mice_at_30_bw %>% select(-contains(" "))
sapply(mice_at_30_bw, class)
mice_at_30_bw <- transform(mice_at_30_bw, Week3 = as.numeric(Week3))
mice_at_30_bw <- transform(mice_at_30_bw, Week4 = as.numeric(Week4))
mice_at_30_bw <- transform(mice_at_30_bw, Week5 = as.numeric(Week5))
mice_at_30_bw <- transform(mice_at_30_bw, Week6 = as.numeric(Week6))
mice_at_30_bw <- transform(mice_at_30_bw, Week7 = as.numeric(Week7))
mice_at_30_bw <- transform(mice_at_30_bw, Week8 = as.numeric(Week8))
mice_at_30_bw[is.na(mice_at_30_bw)] <- 1



#now the fat mass
mice_at_30_fm <- read_excel("201904-06 CLAMS HFD challenge record.xlsx", 
                            sheet = "30C born, left side CLAMS - MRI", 
                            skip = 33, n_max = 12)
names(mice_at_30_fm) <- c("Channel", "ID", "Week0", "Week1", "Week2", "Week3", "Week4", "Week5", "Week6","Week7", "Week8")
mice_at_30_fm <- transform(mice_at_30_fm, Week3 = as.numeric(Week3))
mice_at_30_fm <- transform(mice_at_30_fm, Week4 = as.numeric(Week4))
mice_at_30_fm <- transform(mice_at_30_fm, Week5 = as.numeric(Week5))
mice_at_30_fm <- transform(mice_at_30_fm, Week6 = as.numeric(Week6))
mice_at_30_fm <- transform(mice_at_30_fm, Week7 = as.numeric(Week7))
mice_at_30_fm <- transform(mice_at_30_fm, Week8 = as.numeric(Week8))
mice_at_30_fm[is.na(mice_at_30_fm)] <- 1


mice_at_30_adiposity <- mice_at_30_fm
for (i in 3:11) {
  mice_at_30_adiposity[,i] <- mice_at_30_fm[,i]/mice_at_30_bw[,i]
}

#now the lean mass
mice_at_30_lm <- read_excel("201904-06 CLAMS HFD challenge record.xlsx", 
                            sheet = "30C born, left side CLAMS - MRI", 
                            skip = 49, n_max = 12)
names(mice_at_30_lm) <- c("Channel", "ID", "Week0", "Week1", "Week2", "Week3", "Week4", "Week5", "Week6","Week7", "Week8")
mice_at_30_lm <- transform(mice_at_30_lm, Week3 = as.numeric(Week3))
mice_at_30_lm <- transform(mice_at_30_lm, Week4 = as.numeric(Week4))
mice_at_30_lm <- transform(mice_at_30_lm, Week5 = as.numeric(Week5))
mice_at_30_lm <- transform(mice_at_30_lm, Week6 = as.numeric(Week6))
mice_at_30_lm <- transform(mice_at_30_lm, Week7 = as.numeric(Week7))
mice_at_30_lm <- transform(mice_at_30_lm, Week8 = as.numeric(Week8))
mice_at_30_lm[is.na(mice_at_30_lm)] <- 1



#### now let's do the same for the 22C
mice_at_22_bw <- read_excel("201904-06 CLAMS HFD challenge record.xlsx", 
                            sheet = "22C born, right side CLAMS - MR", 
                            skip = 1, n_max = 12)
mice_at_22_bw <- mice_at_22_bw[,1:11]


#fix the column names
names(mice_at_22_bw) <- c("Channel", "ID", "Week0", "Week1", "Week2", "Week3", "Week4", "Week5", "Week6","Week7", "Week8")
sapply(mice_at_22_bw, class)

mice_at_22_fm <- read_excel("201904-06 CLAMS HFD challenge record.xlsx", 
                            sheet = "22C born, right side CLAMS - MR", 
                            skip = 33, n_max = 12)
mice_at_22_fm <- mice_at_22_fm[,1:11]


#fix the column names
names(mice_at_22_fm) <- c("Channel", "ID", "Week0", "Week1", "Week2", "Week3", "Week4", "Week5", "Week6","Week7", "Week8")
sapply(mice_at_22_fm, class)

mice_at_22_adiposity <- mice_at_22_fm
for (i in 3:11) {
  mice_at_22_adiposity[,i] <- mice_at_22_fm[,i]/mice_at_22_bw[,i]
}

mice_at_22_lm <- read_excel("201904-06 CLAMS HFD challenge record.xlsx", 
                            sheet = "22C born, right side CLAMS - MR", 
                            skip = 49, n_max = 12)
mice_at_22_lm <- mice_at_22_lm[,1:11]


#fix the column names
names(mice_at_22_lm) <- c("Channel", "ID", "Week0", "Week1", "Week2", "Week3", "Week4", "Week5", "Week6","Week7", "Week8")
sapply(mice_at_22_lm, class)

all_bw <- rbind(mice_at_22_bw, mice_at_30_bw)
all_bw <- select(all_bw, - Channel)  %>% remove_rownames %>% column_to_rownames(var="ID")
all_fm <- rbind(mice_at_22_fm, mice_at_30_fm)
all_fm <- select(all_fm, - Channel)  %>% remove_rownames %>% column_to_rownames(var="ID")
all_lm <- rbind(mice_at_22_lm, mice_at_30_lm)
all_lm <- select(all_lm, - Channel)  %>% remove_rownames %>% column_to_rownames(var="ID")




# finally let's add all this stuff to the CLAMS_HFD object
for(i in names(CLAMS_HFD)){
  for (j in names(CLAMS_HFD[[i]])) {
    CLAMS_HFD[[i]][[j]]$Weight <- all_bw[j, i]
    CLAMS_HFD[[i]][[j]]$`Fat Mass` <- all_fm[j, i]
    CLAMS_HFD[[i]][[j]]$`Lean Mass` <- all_lm[j, i]
    CLAMS_HFD[[i]][[j]]$Adiposity <- all_fm[j, i]/all_bw[j, i]
  }
}


### now go to the script to make group plots for each parameter and update the summary
#let's save
save(CLAMS_HFD, file = "CLAMS_HFD.RData")




## and now let's update the CLAMS_HFD_clean

for (i in 1:8) {
  CLAMS_HFD[[i]][["F980"]] <- NULL
}

### so now I removed the dead mouse


save(CLAMS_HFD, file = "CLAMS_HFD.RData")





### now let's do the chronoshift for the 22C
### this is because there was 1 week of difference between 22C and 30C, meaning that the recordings started at the same time, but 22C started HFD a week later than 30C
### so we need to adjust the week orders to aligh mice by experimental week
### thus, we are going to remove week 1 of 22C, shift all the other weeks one week back, and then remove week 8 from 30C
### we will have 7 weeks in total
### Additionally, the nomenclature in this item is shifted compared to the wording in the paper.
### what we refer to week 2,3, and 8 in the paper are here called 1,2,and 7. 
### This is just a nomenclature that I have used for data analysis, and it doesn't change 

for(i in 2:8){
  for(j in born_at_22){
    CLAMS_HFD[[i-1]][[j]] <- CLAMS_HFD[[i]][[j]]
  }
}
## and now let's remove the last week for everybody
CLAMS_HFD[[8]] <- NULL






####  now let's add the food intake for the HFD

food_intake <- read_excel("BW composition and food intake HFD.xlsx", 
                          sheet = "food intake", range = "B1:H27") 
food_intake <- rbind.data.frame(food_intake[1:11,], food_intake[15:26,])


for(i in 2:7){
  for(j in 1:length(CLAMS_HFD[[i]])){
    CLAMS_HFD[[i]][[j]]$food_intake <- as.numeric(food_intake[food_intake$Mouse == CLAMS_HFD[[i]][[j]][["ID"]], i])
  }
}

#and now the food intake for week 1, which is chow and we can use the data from CLAMS

for(j in 1:length(CLAMS_HFD[[1]])){
  CLAMS_HFD[[1]][[j]]$food_intake <- sum(CLAMS_HFD[[1]][[j]][["Data"]][["Feed.Weight.1"]])
}



save(CLAMS_HFD, file = "CLAMS_HFD_clean.RData")

rm(list= ls()[!(ls() %in% c("CLAMS_HFD", "born_at_22", "born_at_30"))])
