# This function extracts vertical change in angle (Pitch) across 
# three 3D points: A B and C
# David Pritchard


#-----------------------------------------------------------------------
getVAngle <- function(A,B,C){
  V = c(B[1]-A[1],B[2]-A[2],B[3]-A[3]) # difference between points A & B
  W = c(C[1]-B[1],C[2]-B[2],C[3]-B[3]) # difference between points B & C
  Vangle = atan2(W[3],sqrt((W[1]^2)+(W[2]^2)))-atan2(V[3],sqrt((V[1]^2)+(V[2]^2))) # difference in vertical directions between vectors
  while(Vangle<=(-pi)){Vangle = Vangle+2*pi}
  while(Vangle>(pi)){Vangle = Vangle-2*pi}
  VertAngle = append(VertAngle,Vangle)
  return(Vangle)
}