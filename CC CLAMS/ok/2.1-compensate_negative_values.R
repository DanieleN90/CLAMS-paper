compensate_negative_values <- function(series){
  for(i in 1:length(series)){
    if(series[i] < 0){
      if(!is.na(series[i+1]) & series[i] + series[i+1] == 0){
        series[i] <- 0
        series[i+1] <- 0
      }else if(!is.na(series[i+1]) & series[i] + series[i+1] > 0){
        series[i+1] <- series[i] + series[i+1]
        series[i] <- 0
      }else if(!is.na(series[i+2]) & series[i] + series[i+2] > 0){
        series[i+2] <- series[i] + series[i+2]
        series[i] <- 0
      }else{
        series[i] <- 0
      }
    }
  }
return(series)
}
