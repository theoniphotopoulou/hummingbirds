# This function acts to produce a grid of x,y points of a given diameter 
# in inter-point distance
# David Pritchard


#-----------------------------------------------------------------------
makeGrid = function(diam,dist){
  sc = seq(from = -(diam/2), to = (diam/2), by  = dist)
  hor = rep(sc, times = length(sc)) # get X coordinates
  ver = c()
  for(i in sc){
    bb = rep(i, times = length(sc))
    ver = append(ver,bb) # For every value of X there is every value of Y
  }
  grid = data.frame(cbind(hor,ver))
  names(grid) = c('X','Y')
  return(grid) # returns 2 colums: X & Y
}