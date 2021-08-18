# Load packages
require(circular)
require(CircStats)

#-----------------------------------------------------------------------
# Read in raw data from file
sheet <- read.csv("data/allData_Aug2021.csv",header=TRUE) 
head(sheet)


#-----------------------------------------------------------------------
# Set up variables
All <- c(1:14) # Bird IDs
LM <- c(2,4,6,7,8,9,13) # IDs of birds with landmarks present
noLM <- c(1,3,5,10,11,12,14) # IDs of birds with landmarks removed


#-----------------------------------------------------------------------
# Make a data.frame of test data only
testt <- data.frame(cbind(sheet$Expt,
                          sheet$ID,
                          rep('NA',times = length(sheet$Expt)),
                          rep('NA',times = length(sheet$Expt)),
                          rep('NA',times = length(sheet$Expt)),
                          rep('Test',times = length(sheet$Expt)),
                          sheet$Testx,
                          sheet$Testy,
                          sheet$Testz))

names(testt) = c('Exp','ID', 'Site','LM','NearFar', 'Section','X','Y','Z')
allData <- testt
allData$X <- as.numeric(as.character(allData$X))
allData$Y <- as.numeric(as.character(allData$Y))
allData$Z <- as.numeric(as.character(allData$Z))
allData$LM <- factor(allData$LM, levels = c('Y','N',NA))


#-----------------------------------------------------------------------
# Load function to extract horizontal change in angle (Yaw) across 
# three 3D points: A B and C
source("functions/getHAngle.R")


#-----------------------------------------------------------------------
# Load function to extract vertical change in angle (Pitch) across 
# three 3D points: A B and C
source("functions/getVAngle.R")


#-----------------------------------------------------------------------
# Consolidate across tests into single data.frame

# loop over experiments (1 to 3) 
for(e in 1:3){
  
  currsheet <- allData[which(allData$Exp==e),]
  
  sheet.col.names <- c('Exp','ID', 'Section', 'LM','NearFar', 'Site', 
                       'X','Y','Z', 'step', 'yaw', 'pitch', 'CurrFlowerDist')
  output <- matrix(rep(NA,times=(nrow(currsheet)*length(sheet.col.names))),
                   nrow=nrow(currsheet), ncol=length(sheet.col.names))
  # Make an empty data frame
  output <- as.data.frame(output)
  names(output) <- sheet.col.names
  
  # Fill in variables from older sheets
  output$Exp <- currsheet$Exp           # Experiment
  output$ID <- currsheet$ID             # Bird ID
  output$Site <- currsheet$Site         # Location
  output$Section <- currsheet$Section   # Test
  output$LM <- currsheet$LM             # Landmark present or not
  output$NearFar <- currsheet$NearFar   # We did not find an effect of moving near and 
                                        # far so did not end up using this
  # Coordinates
  output$X <- currsheet$X
  output$Y <- currsheet$Y
  output$Z <- currsheet$Z  
  
  # loop over individual birds (1 to 14)
  for(i in 1:14){ 
    for(k in c('Test')){
      
      # Bird coordinates
      Tx <- na.omit(currsheet$X[which(currsheet$ID==i & currsheet$Section==k)])
      Ty <- na.omit(currsheet$Y[which(currsheet$ID==i & currsheet$Section==k)])
      Tz <- na.omit(currsheet$Z[which(currsheet$ID==i & currsheet$Section==k)])
      
      # Flower coordinates
      Flower <- c(na.omit(currsheet$Flowerx[which(currsheet$ID==i)]),
                  na.omit(currsheet$Flowery[which(currsheet$ID==i)]),
                  na.omit(currsheet$Flowerz[which(currsheet$ID==i)]))
      
      # Set up empty vectors to fill in
      stepl <- c()        # Step length
      Fcurr <- c()
      Fprev <- c()
      VertAngle <- c(NA)  # Pitch angle
      HorAngle <- c(NA)   # Yaw angle
      
      
      for(j in 2:length(Tx)-1){
        
        X = c(Tx[j],Ty[j],Tz[j])         # previous point, starting from 1
        Z = c(Tx[j+2],Ty[j+2],Tz[j+2])   # next point
        Y = c(Tx[j+1],Ty[j+1],Tz[j+1])   # current point starting from 2
        if(is.na(Z[1])== TRUE|is.na(X[1]==TRUE)){ 
          stepl = append(stepl,NA)
          HorAngle = append(HorAngle,NA)
          VertAngle = append(VertAngle,NA)
          
        } # if furthest forward not NA...
        else{
          
          # get distance between points
          ac = sqrt((Y[1]-X[1])^2+(Y[2]-X[2])^2+(Y[3]-X[3])^2) 
          
          stepl = append(stepl,ac)
          
          # Yaw & Pitch
          Hangle = getHAngle(X,Y,Z)
          HorAngle = append(HorAngle,Hangle) 
          
          Vangle = getVAngle(X,Y,Z)
          VertAngle = append(VertAngle,Vangle)
          
        }
      }
      # fill in empty cells in output sheet
      output$step[which(output$ID==i & output$Section==k & is.na(output$X)==FALSE)][1:length(stepl)] = stepl
      output$yaw[which(output$ID==i & output$Section==k & is.na(output$X)==FALSE)][1:length(HorAngle)] = HorAngle
      output$pitch[which(output$ID==i & output$Section==k & is.na(output$X)==FALSE)][1:length(VertAngle)] = VertAngle
      
    }
  }
  if(e==1){Output1 = output}
  if(e==2){Output2 = output}
  if(e==3){Output3 = output}
  
  rm(output) # remove used sheet
}

# convert Output to a single data.frame "bird"
bird <- rbind(Output1, Output2, Output3) 
# remove Output sheets
rm(Output1, Output2, Output3) 

# remove unrealistic values (based on max dive speed of 22 m/s)
bird <- bird[which(bird$step < 880),] 


#-----------------------------------------------------------------------
# Calculate current flower distance at each time point

# loop over experiments (1 to 3)
for(e in c(1:3)){
  for(i in unique(bird$ID)){
    
    flowerx = sheet$Flowerx[which(sheet$Expt==e&sheet$ID==i)][1]
    flowery = sheet$Flowery[which(sheet$Expt==e&sheet$ID==i)][1]
    flowerz = sheet$Flowerz[which(sheet$Expt==e&sheet$ID==i)][1]
    
    bird$CurrFlowerDist[which(bird$Exp==e & bird$ID==i)] <- sqrt((bird$X[which(bird$Exp==e&bird$ID==i)] - flowerx)^2 + 
                                                                 (bird$Y[which(bird$Exp==e&bird$ID==i)] - flowery)^2 + 
                                                                 (bird$Z[which(bird$Exp==e&bird$ID==i)] - flowerz)^2)
  }
}


#-----------------------------------------------------------------------
# Add stops

stops <- rep(NA,times = nrow(bird))
stops[which(bird$step<0.005)] <- 'Y'
stops[which(bird$step>=0.005)] <- 'N'
stops <- factor(stops)
bird <- cbind(bird,stops)


#-----------------------------------------------------------------------
# Add if landmarks present

bird$LM[which(bird$Section=='Test' & bird$ID %in% noLM)] <- "N"

bird$LM[which(bird$Section=='Test' & bird$ID %in% LM)] <- "Y"

bird$LM[which(bird$Section=='Test' & bird$Exp==1)] <- "Y"

# remove other used sheets
rm(allData,currsheet,testt,sheet) 


#-----------------------------------------------------------------------
# Save processed data
write.csv(bird, file = "data/processed-data.csv",row.names = FALSE)

