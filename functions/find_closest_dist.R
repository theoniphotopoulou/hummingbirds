# The find_closest_dist() function is for linking estimated transition probabilities
# at fixed distances from a flower to observed, irregular distances  
#
# Theoni Photopoulou
# 20200220

# I have estimated transition probabilities between states at regular intervals
# in the range of the CurrFlowerDist variable. This function is for finding 
# which distance to extract transition probabilities at, since the observed 
# distances are irregular.
# obsdist is a scalar and distvector a vector
find_closest_dist <- function(obsdist,distvector, ...){
  tprobs <- distvector[which.min(abs(obsdist-distvector))]
  return(tprobs)
}