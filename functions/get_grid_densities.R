# Determine density of Searching behaviour on a grid
# based on model predictions
# Function by David Pritchard

get_grid_densities <- function(grid, expdata, landmarks){

#  expdata <- exp2Rotated
#  landmarks <- "N"  
#  grid1 <- makeGrid(6,r)
  
#  grid1 <- cbind(grid1,rep(0,times = nrow(grid)),
#                rep(0,times = nrow(grid)),
#                rep(0,times = nrow(grid)))
  # Indsum is sum of probabilities for index
  # Indmean is mean probability, and 
  # Indprop is proportion relative to maximum sum.
#  names(grid1)[3:5] <- c('Indsum','Indmean','Indprop') 
  
  birdid <- expdata %>% filter(LM==landmarks) %>% distinct(ID)
  birdid <- birdid$ID
  
  for(bird in birdid){
    
    # set up empty vectors
    count <- rep(0,times = nrow(grid))
    sum <- rep(0,times = nrow(grid))
    mean <- rep(0,times = nrow(grid))
    
    for(i in 1:length(expdata$X[which(expdata$ID==bird)])){
      # count through the points on the path
      x <- expdata$X[which(expdata$ID==bird)][i]
      y <- expdata$Y[which(expdata$ID==bird)][i]
      # find the closest grid point
      ind <- which(abs(grid$X-x)==min(abs(grid$X-x))&abs(grid$Y-y)==min(abs(grid$Y-y)))
      # count the number of points falling in each grid point per bird
      count[ind] <- count[ind]+1 
      # sum the investigation probabilities
      sum[ind] <- sum[ind]+expdata$Search[which(expdata$ID==bird)][i] 
    }
    mean <- sum/count 
    prop <- sum/max(sum) #scale probabilities to a max of 1
    mean[which(is.na(mean))] <- 0
    for(j in 1:length(mean)){
      grid$Indsum[j] <- grid$Indsum[j]+sum[j]
      grid$Indmean[j] <- grid$Indmean[j]+mean[j]
      grid$Indprop[j] <- grid$Indprop[j]+prop[j] # add birds together
      
    }
  }  
  return(grid)
}
