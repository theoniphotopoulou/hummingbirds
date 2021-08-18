# This function extracts the vertical change in angle (Pitch) of absolute directions 
# between stops cross 3 points: A B and C
# David Pritchard


#-----------------------------------------------------------------------
getABSVAngle = function(A,B){
  V = c(B[1]-A[1], B[2]-A[2], B[3]-A[3]) # difference between points A & B
  Vangle = atan2(V[3], sqrt((V[1]^2) + (V[2]^2))) # difference in vertical directions between vectors
  while(Vangle<=(-pi)){Vangle = Vangle+2*pi}
  while(Vangle>(pi)){Vangle = Vangle-2*pi}
  VertAngle = append(VertAngle,Vangle)
  return(Vangle)
}