#-----------------------------------------------------------------------
# Read in processed data
bird <- read.csv("data/processed-data.csv")


#-----------------------------------------------------------------------
# Set up empty vectors to fill in
id <- c() # Bird ID
expt <- c() # Experiment
LM <- c() # Landmarks present of not

firstStop <- c() # Distance of first stop from the flower
closestStop <- c() # Distance of closest stop to the flower
overall <- c() # Closest flown overall


#-----------------------------------------------------------------------
# Loop over experiments (1 to 3)
for(e in c(1:3)){
  
  # loop over birds
  for(i in unique(bird$ID[which(bird$Exp==e & bird$Section == 'Test')])){ 
    
    id <- append(id,i) # get ID
    expt <- append(expt,e) # get Experiment
    
    # Get Landmark code (in processed data sheet)
    LM <- append(LM,head(bird$LM[which(bird$ID==i & bird$Exp==e & bird$Section=='Test')], n=1)) 
    # Identify distance of bird when first moving less than 5mm
    a <- head(bird$CurrFlowerDist[which(bird$ID==i & bird$Exp==e & bird$Section=='Test' & bird$step<5)], n=1) 
    # What was the closest distance when moving less than 5mm
    b <- min(bird$CurrFlowerDist[which(bird$ID==i & bird$Exp==e & bird$Section=='Test' & bird$step<5)])
    # What was the closest distance overall?
    c <- min(bird$CurrFlowerDist[which(bird$ID==i & bird$Exp==e & bird$Section=='Test')]) 
    
    # Add to vectors
    firstStop = append(firstStop,a)
    closestStop = append(closestStop,b)
    overall = append(overall,c)
  } # closes i loop
} # closes e loop


#-----------------------------------------------------------------------
# Bind vectors to data.frame
routeSummary <- data.frame(cbind(id, expt, LM, firstStop, closestStop, overall))


#-----------------------------------------------------------------------
# Make variables numeric
routeSummary$firstStop <- as.numeric(as.character(routeSummary$firstStop))
routeSummary$closestStop <- as.numeric(as.character(routeSummary$closestStop))
routeSummary$overall <- as.numeric(as.character(routeSummary$overall))



#-----------------------------------------------------------------------
# Save summary data
write.csv(routeSummary,file = "data/summary_data.csv", row.names = FALSE)






