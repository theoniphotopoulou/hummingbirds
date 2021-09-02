# This function acts to rotate and resize data so they can be 
# plotted on top of one another
# David Pritchard


#-----------------------------------------------------------------------
rotate_resize <- function(expdata, exp_birdid){
  
  for(birdid in exp_birdid){ # loop through birds

    # First step is to standardise the data to a common frame of reference. 
    # We do this by making the flower coordinates the same for all birds.
    # Convert all coordinates so Flower is at 0,0
    relativex <- expdata$X[which(expdata$ID==birdid)]-expdata$Flowerx[which(expdata$ID==birdid)][1]
    relativey <- expdata$Z[which(expdata$ID==birdid)]-expdata$Flowerz[which(expdata$ID==birdid)][1]
    relativeLLMx <- expdata$LeftLMx[which(expdata$ID==birdid)][1]-expdata$Flowerx[which(expdata$ID==birdid)][1]
    relativeLLMy <- expdata$LeftLMz[which(expdata$ID==birdid)][1]-expdata$Flowerz[which(expdata$ID==birdid)][1]
    relativeRLMx <- expdata$RightLMx[which(expdata$ID==birdid)][1]-expdata$Flowerx[which(expdata$ID==birdid)][1]
    relativeRLMy <- expdata$RightLMz[which(expdata$ID==birdid)][1]-expdata$Flowerz[which(expdata$ID==birdid)][1]
    
    # For all birds the landmarks are at different angles to the camera. To standardise 
    # for plotting we need to rotate the data to a common orientation. We do this by 
    # setting a fixed vector to be parallel to the x axis. Here we close the flower-right 
    # landmark vector.
    # Set up empty vectors to be filled with rotated values
    rotatedx <- c()
    rotatedy <- c()
    
    # Get coordinates for the Right landmark
    RLMx <- relativeRLMx
    RLMy <- relativeRLMy
    
    # Calculate the angle of the vector between the right Landmark the flower (at 0,0)
    if((RLMx == 0) & (RLMy > 0)){ angle <- rad(90)
    } else {
      if((RLMx == 0)&(RLMy < 0)){ angle = rad(270)
      } else {
        if(RLMx > 0){ angle=(atan(RLMy/RLMx))
        } else {
          if(RLMx < 0){ angle <- rad(180 + deg((atan(RLMy/RLMx))))
          }
        }  
      }
    }  
    
    # Using the angle derived above, rotate the points so that angle of 
    # the Right Landmarks-Flower vector is 0
    for( i in 1:(length(relativex))){
      xi <- relativex[i]
      yi <- relativey[i]
      if(is.na(relativex[i])){
        x <- NA 
        y <- NA
      } else {
        x <- (xi*cos(-angle)) - (yi*sin(-angle))
        y <- (xi*sin(-angle)) + (yi*cos(-angle))
      }
      rotatedx <- append(rotatedx, x)
      rotatedy <- append(rotatedy, y)
    }
    
    # Rotate the Left Landmark
    xi <- relativeLLMx
    yi <- relativeLLMy
    if(is.na(relativeLLMx)){
      x=NA
      y=NA
    } else {
      x <- (xi*cos(-angle)) - (yi*sin(-angle))
      y <- (xi*sin(-angle)) + (yi*cos(-angle))
    }
    rotatedLLMx <- x
    rotatedLLMy <- y
    
    # Rotate the Right Landmark
    lmxi <- RLMx
    lmyi <- RLMy
    lmx <- (lmxi*cos(-angle)) - (lmyi*sin(-angle))
    lmy <- (lmxi*sin(-angle)) + (lmyi*cos(-angle))
    rotatedRLMx <- lmx
    rotatedRLMy <- lmy
    
    # Having centred on the flower and rotated to match orientation, we then resize 
    # the data to correct for slight differences between sites. We know that the 
    # right landmark was 0.3m from the flower. As it has been rotated, the X 
    # coordinate should be 0.3. By dividing the observed X coordinate from 0.3 
    # we get a scaling factor that can be used to resize all of the points for a bird
    
    # Get scaling factor
    scale <- (rotatedRLMx)/0.30
    
    # Resize data by dividing by scaling factor
    resizedx <- rotatedx/scale
    resizedy <- rotatedy/scale
    
    # Make up empty vectors the same length as the number of points in the flightpath
    resizedRLMx <- rep(NA, times = length(resizedx))
    resizedRLMy <- rep(NA, times = length(resizedx))
    resizedLLMx <- rep(NA, times = length(resizedx))
    resizedLLMy <- rep(NA, times = length(resizedx))
    
    # Resize landmark coordinates using scaling factor
    resizedRLMx[1] <- rotatedRLMx/scale
    resizedRLMy[1] <- rotatedRLMy/scale
    resizedLLMx[1] <- rotatedLLMx/scale
    resizedLLMy[1] <- rotatedLLMy/scale  
    
    # Set up dataframe of standardised data
    
    # make bird ID vector
    ID <- rep(birdid, times = length(resizedx))
    # combine vectors
    current <- cbind(ID,resizedx,resizedy,resizedLLMx,resizedLLMy,resizedRLMx,resizedRLMy) 
    # if first loop then make a dataframe
    if(birdid==1){
      rotatedBirds <- data.frame(current)
      # if not first loop then add to dataframe
    } else{
      rotatedBirds <- rbind(rotatedBirds,current)
    }
  }
  return(rotatedBirds)
}