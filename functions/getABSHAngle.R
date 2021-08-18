# This function extracts the horizontal change in angle (Yaw) of absolute directions 
# between stops cross 3 points: A B and C
# David Pritchard


#-----------------------------------------------------------------------
getABSHAngle = function(A,B,C){
  V = c(B[1]-A[1], B[2]-A[2], B[3]-A[3]) # difference between points A & B
  Hangle = atan2(V[2], V[1]) # difference in horizontal direction between vectors 
  while(Hangle<=(-pi)){Hangle = Hangle+2*pi}
  while(Hangle>(pi)){Hangle = Hangle-2*pi}
  return(Hangle)
}