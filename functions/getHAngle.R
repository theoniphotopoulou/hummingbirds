# This function extracts horizontal change in angle (Yaw) across three
# 3D points: A B and C
# David Pritchard


#-----------------------------------------------------------------------
getHAngle <- function(A,B,C){
  V = c(B[1]-A[1],B[2]-A[2],B[3]-A[3]) # difference between points A & B
  W = c(C[1]-B[1],C[2]-B[2],C[3]-B[3]) # difference between points B & C
  Hangle = atan2(W[2],W[1])-atan2(V[2],V[1]) # difference in horizontal direction between vectors 
  while(Hangle<=(-pi)){Hangle = Hangle+2*pi}
  while(Hangle>(pi)){Hangle = Hangle-2*pi}
  return(Hangle)
}