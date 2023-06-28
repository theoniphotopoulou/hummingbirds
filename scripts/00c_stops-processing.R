######## PROCESSING OF STOPS ######## 
# This is for comparison between experiments
# David Pritchard, August 2021

# Load packages

library(gridExtra) # arranging plots
library(circular) # circular data functions
library(CircStats) # circular statistics functions
library(lme4) # linear mixed models
library(tidyverse) # tidyverse functions
library(MuMIn) # model averaging
require(RColorBrewer) # colour palette
library(reshape2) # rearranging datasets


## READ IN FUNCTIONS 

# source function to extract horizontal change in angle (Yaw) cross 3 points: A B and C
source("functions/getHAngle.R")

# source function to extract vertical change in angle (Pitch) across three 3D points: A B and C
source("functions/getVAngle.R")

# source function to extract absolute horizontal change in angle (Yaw) across three 3D points: A B and C
#  *** WHAT DOES THIS FUNCTION DO? ***
source("functions/getABSHAngle.R")

# source function to extract absolute vertical change in angle (Pitch) across three 3D points: A B and C
#  *** WHAT DOES THIS FUNCTION DO? ***
source("functions/getABSVAngle.R")

# source summarising function
source("functions/summarySE.R")


## READ IN PROCESSED TRAJECTORY DATA
bird <- read.csv("data/processed-data.csv")


##############################
## 1. Identify stops and make stops data frame
##############################

# Set up empty vectors to fill in
id <- c() # Bird ID
expt <- c() # Experiment
LM <- c() # Landmarks present of not

# Coordinates
X <- c() # X coordinates of stops
Y <- c() # Y coordinates of stops
Z <- c() # Z coordinates of stops
Dist <- c() # Distance to the flower
Dur <- c() # Duration of stop

# Cycle through experiments and identify discrete stops through run-length encoding

for(e in c(1:3)){ # e loops over experiments 1 to 3
  for(i in unique(bird$ID[which(bird$Exp==e & bird$Section=='Test')])){ # cycle through birds
    
    # get TRUE/FALSE if data in test where moved less than 1cm
    stops <- bird$ID==i & bird$Exp==e & bird$Section=='Test' & bird$step<10 
    # run-length encoding to identify continuous runs of TRUE (step<1cm) or FALSE (step>1cm)
    aa <- rle(stops) 
    # start indexes of runs
    ind <- c(1, lag(cumsum(aa$lengths))[-1] + 1)
    # start indexes of TRUE values (i.e. stops)
    start <- ind[which(aa$values=='TRUE')] 
    # use length of TRUE runs to extract duration
    Dur <- append(Dur,aa$lengths[which(aa$values=='TRUE')]/25) 
    
    x <- bird$X[start] # X coordinates of stops
    y <- bird$Y[start] # Y coordinates of stops
    z <- bird$Z[start] # Z coordinates of stops
    
    # distance of stops from flower location
    dist <- bird$CurrFlowerDist[start] 
    # get ID
    id <- append(id,rep(i,times = length(start))) 
    # get Experiment
    expt <- append(expt,rep(e,times = length(start)))
    # get Landmark code (in processed data sheet)
    LM <- append(LM,
                 rep(head(bird$LM[which(bird$ID==i & bird$Exp==e & bird$Section=='Test')],n=1),
                     times = length(start))) 
    
    # add to vectors
    X <- append(X,x)
    Y <- append(Y,y)
    Z <- append(Z,z)
    Dist <- append(Dist,dist)
  }
}

# Bind vectors to data.frame
stops <- cbind.data.frame(id, expt, LM, X, Y, Z, Dist, Dur)

# Ensure numeric data are numeric
stops$X <- as.numeric(stops$X)
stops$Y <- as.numeric(stops$Y)
stops$Z <- as.numeric(stops$Z)
stops$Dist <- as.numeric(stops$Dist)
stops$Dur <- as.numeric(stops$Dur)



##############################
## 2. Get distances and directions between stops
##############################

# Set up empty vectors to fill in
stepl <- c() #step length
ID <- c()
EXP <- c()
VertAngle <- c() #Pitch
HorAngle <- c() #Yaw
VertTurn <- c() #Pitch
HorTurn <- c() #Yaw

for(exp in unique(stops$expt)){ # exp loops over experiments in the stops data frame
  for(id in unique(stops$id[which(stops$expt==exp)])){ # id loops over individual birds within eacg experiment
    
    # distances and directions are between pairs of stops so start with NA
    stepl <- append(stepl,NA)
    HorAngle <- append(HorAngle,NA)
    VertAngle <- append(VertAngle,NA)
    HorTurn <- append(HorTurn,c(NA))
    VertTurn <- append(VertTurn,c(NA))
    # append id and experiment
    ID <- append(ID,id)
    EXP <- append(EXP,exp)
    
    # get xyz coordinates for stops by focal bird in focal experiment
    Tx <- as.numeric(stops$X[which(stops$expt==exp&stops$id==id)])
    Ty <- as.numeric(stops$Y[which(stops$expt==exp&stops$id==id)])
    Tz <- as.numeric(stops$Z[which(stops$expt==exp&stops$id==id)])
    
    if(length(Tx)==1){next} # if there is only 1 stops in a bird/experiment combo then skip to next iteration
    # otherwise loop through stops
    for(j in 1:(length(Tx)-1)){
      
      X <- c(Tx[j], Ty[j], Tz[j]) # previous stop, starting from 1
      Y <- c(Tx[j+1], Ty[j+1], Tz[j+1]) # current stop starting from 2
      
      ac <- sqrt((Y[1]-X[1])^2 + (Y[2]-X[2])^2 + (Y[3]-X[3])^2) # get distance between stops
      
      stepl <- append(stepl,ac) # append between-stop-distance to vector
      ID <- append(ID,id) # append ID
      EXP <- append(EXP,exp) # append experiment
      
      # Yaw & Pitch of absolute directions between stops
      Hangle <- getABSHAngle(X,Y)
      HorAngle <- append(HorAngle,Hangle) 
      
      Vangle <- getABSVAngle(X,Y)
      VertAngle <- append(VertAngle,Vangle)
      
    }
    
    # Yaw and Pitch of turning direction based on angle between 3 stops
    HorTurn <- append(HorTurn,c(NA))
    VertTurn <- append(VertTurn,c(NA))
    
    if(length(Tx)==2){next} # if 2 or less stops in bird/experiment combo then skip to next iteration
    # otherwise loop through stops
    for(k in 2:(length(Tx)-1)){
      
      X <- c(Tx[k-1], Ty[k-1], Tz[k-1]) # previous stop, starting from 1
      Y <- c(Tx[k], Ty[k], Tz[k]) # current stop starting from 2
      Z <- c(Tx[k+1], Ty[k+1], Tz[k+1]) # next stop starting from 3
      
      # append turning angles
      Hangle2 <- getHAngle(X,Y,Z)
      HorTurn <- append(HorTurn,Hangle2) 
      
      Vangle2 <- getVAngle(X,Y,Z)
      VertTurn <- append(VertTurn,Vangle2)
      
    }
  }
}

# Make data frame of between-stop distances and directions
stopVector <- cbind.data.frame(stops, stepl, 
                               HorAngle, VertAngle, 
                               HorTurn, VertTurn)



##############################
## 3. Compare first & closest stops and closest flown overall between experiments
##############################

# Read in data on first & closest stops and closest flown overall
allData <- read.csv('data/summary_data.csv')

# Make sure landmark is a factor
allData$LM <- factor(allData$LM,label = c('Absent','Present')) 

# Create data frames for comparing experiments 1 & 2 and 2 & 3
summary12 <- allData[which(allData$expt %in% c(1,2)),]
summary23 <- allData[which(allData$expt %in% c(2,3)),]
summary12$LM[which(summary12$id %in% summary12$id[which(summary12$LM=='Absent')])]='Absent' 
# convert experiment 1 birds in landmark removed group to "Absent" 
# (even though present in expt 1)

# Calculate the difference between closest stop and overall closest flown
summary12 <- summary12 %>% 
  mutate(diffStop = closestStop-overall)

summary23 <- summary23 %>% 
  mutate(diffStop = closestStop-overall) 

# Make sure experiment and id are factors
summary12$expt <- factor(summary12$expt)
summary23$expt <- factor(summary23$expt)

summary12$id <- factor(summary12$id)
summary23$id <- factor(summary23$id)

# Scale the firstStop, closestStop, overall and diffStop variables 
# in the summary dataframes. This is done to move the values away 
# from zero so that the Gamma distribution can be used for modelling, 
# below. 

# We use a constant of 3 to offset all variables in data sets, 
# from all experiments
scaling.offset <- 3

summary12 <- summary12 %>% 
  mutate(firstScaled = scale(firstStop,
                             center = TRUE,
                             scale = TRUE) + scaling.offset,
         closeScaled = scale(closestStop,
                             center = TRUE,
                             scale = TRUE) + scaling.offset,
         overallScaled = scale(overall,
                               center = TRUE,
                               scale = TRUE) + scaling.offset, 
         diffScaled = scale(diffStop,
                            center = TRUE,
                            scale = TRUE) + scaling.offset) 

summary23 <- summary23 %>% 
  mutate(firstScaled = scale(firstStop,
                             center = TRUE,
                             scale = TRUE) + scaling.offset,
         closeScaled = scale(closestStop, 
                             center = TRUE,
                             scale = TRUE) + scaling.offset,
         overallScaled = scale(overall,
                               center = TRUE,
                               scale = TRUE) + scaling.offset,
         diffScaled = scale(diffStop,
                            center = TRUE,
                            scale = TRUE) + scaling.offset)



# Make data frames for comparing all flight distances in 1 & 2 and 2 & 3
bird12 <- bird[which(bird$Exp %in% c(1,2)),] # paths in experiments 1 & 2
bird23 <-  bird[which(bird$Exp %in% c(2,3)),] # paths in experiments 2 & 3
bird12$LM[which(bird12$ID %in% bird12$ID[which(bird12$LM=='N')])]='N' 
# convert experiment 1 birds in landmark removed 
# group to "N" (even though present in expt 1)

# Scale the CurrFlowerDist variable in the summary dataframes
# This is done to move the values away from zero so that the Gamma
# distribution can be used for modelling, below. 

bird12 <- bird12 %>% 
  mutate(distScaled = scale(CurrFlowerDist, 
                            center = TRUE, 
                            scale = TRUE) + scaling.offset) 

bird23 <- bird23 %>% 
  mutate(distScaled = scale(CurrFlowerDist, 
                            center = TRUE, 
                            scale = TRUE) + scaling.offset) 

# Make sure experiment and id are factors
bird12$Exp <- factor(bird12$Exp)
bird23$Exp <- factor(bird23$Exp)

bird12$ID <- factor(bird12$ID)
bird23$ID <- factor(bird23$ID)

# Make data frames for comparing stop distances in 1 & 2 and 2 & 3

stops12 <- stops[which(stops$expt %in% c(1,2)),] # stops in experiments 1 & 2
stops23 <- stops[which(stops$expt %in% c(2,3)),] # stops in experiments 2 & 3
stops12$LM[which(stops12$id %in% stops12$id[which(stops12$LM=='N')])]='N' 
# convert experiment 1 birds in landmark removed group to "N" 
# (even though present in expt 1)

# Scale the Dist variable in the summary dataframes
# This is done to move the values away from zero so that the Gamma
# distribution can be used for modelling, below. 

stops12 <- stops12 %>% 
  mutate(distScaled = scale(Dist, 
                            center = TRUE, 
                            scale = TRUE) + scaling.offset)
stops23 <- stops23 %>% 
  mutate(distScaled = scale(Dist, 
                            center = TRUE, 
                            scale = TRUE) + scaling.offset) 


# Make sure experiment and id are factors
stops12$expt <- factor(stops12$expt)
stops23$expt <- factor(stops23$expt)

stops12$id <- factor(stops12$id)
stops23$id <- factor(stops23$id)

# stops12$id <- factor(stops12$LM)
# stops23$id <- factor(stops23$LM)

# Make data frames for comparing between-stop distances in 1 & 2 and 2 & 3

stopsV12 <- stopVector[which(stopVector$expt %in% c(1,2)),] # between-stops distance in experiments 1 & 2
stopsV23 <- stopVector[which(stopVector$expt %in% c(2,3)),] # between-stops distance in experiments 2 & 3
stopsV12$LM[which(stopsV12$id %in% stopsV12$id[which(stopsV12$LM=='N')])]='N' 
# convert experiment 1 birds in landmark removed group to "N" 
# (even though present in expt 1)

# Make sure experiment and id are factors and remove NAs and 0 distances
stopsV12$stepl[which(stopsV12$stepl==0)] = NA
stopsV23$stepl[which(stopsV23$stepl==0)] = NA

stopsV12NA <- stopsV12[is.na(stopsV12$stepl)==F,]
stopsV23NA <- stopsV23[is.na(stopsV23$stepl)==F,]

stopsV12NA$expt <- factor(stopsV12NA$expt)
stopsV23NA$expt <- factor(stopsV23NA$expt)

stopsV12NA$id <- factor(stopsV12NA$id)
stopsV23NA$id <- factor(stopsV23NA$id)

# Scale the stepl variable in the summary dataframes
# This is done to move the values away from zero so that the Gamma
# distribution can be used for modelling, below. 

stopsV12NA <- stopsV12NA %>% 
  mutate(stepScaled = scale(stepl,
                            center = TRUE,
                            scale = TRUE) + scaling.offset) 

stopsV23NA <- stopsV23NA %>% 
  mutate(stepScaled = scale(stepl,
                            center = TRUE,
                            scale = TRUE) + scaling.offset) 

# Make data frames for comparing number of stops
stopNum <- aggregate(X ~ expt + id + LM, data = stops, FUN = 'length')
stopNum$LM[which(stopNum$id %in% stopNum$id[which(stopNum$LM=='N')])]='N'

# Make data frames of paths and stops for each experiment and landmarks group
exp1data <- bird[which(bird$Exp==1),]
exp2data <- bird[which(bird$Exp==2),]
exp3data <- bird[which(bird$Exp==3),]

  # Path for landmarks present group
ALM1 <- exp1data$CurrFlowerDist[which(exp1data$Exp==1 & exp1data$ID %in% 
                                        unique(exp2data$ID[which(exp2data$Exp==2 & exp2data$LM=='Y')]))]
ALM2 <- exp2data$CurrFlowerDist[which(exp2data$LM=='Y' & exp2data$Exp==2)]
ALM3 <- exp3data$CurrFlowerDist[which(exp3data$LM=='Y' & exp3data$Exp==3)]

  # Stops for landmarks present group
LM1 <- stops$Dist[which(stops$LM=='Y' & stops$expt==1 & stops$id %in% 
                          unique(stops$id[which(stops$expt==2 & stops$LM=='Y')]))]
LM2 <- stops$Dist[which(stops$LM=='Y' & stops$expt==2)]
LM3 <- stops$Dist[which(stops$LM=='Y' & stops$expt==3)]

  # Path for landmarks removed group
AnLM1 <- exp1data$CurrFlowerDist[which(exp1data$Exp==1 & exp1data$ID %in% 
                                         unique(exp2data$ID[which(exp2data$Exp==2 & exp2data$LM=='N')]))]
AnLM2 <- exp2data$CurrFlowerDist[which(exp2data$LM=='N' & exp2data$Exp==2)]
AnLM3 <- exp3data$CurrFlowerDist[which(exp3data$LM=='N' & exp3data$Exp==3)]

  # Stops for landmarks removed group
noLM1 <- stops$Dist[which(stops$LM=='Y' & stops$expt==1 & stops$id %in% 
                            unique(stops$id[which(stops$expt==2 & stops$LM=='N')]))]
noLM2 <- stops$Dist[which(stops$LM=='N' & stops$expt==2)]
noLM3 <- stops$Dist[which(stops$LM=='N' & stops$expt==3)]

# Make data frames of between-landmarks distances for each experiment and landmarks group

  # Landmarks removed
nodLM1 <- stopVector$stepl[which(stopVector$LM=='Y' & stopVector$expt==1 & stopVector$id %in% 
                                   unique(stopVector$id[which(stopVector$expt==2 & stopVector$LM=='N')]))]
nodLM2 <- stopVector$stepl[which(stopVector$LM=='N' & stopVector$expt==2)]
nodLM3 <- stopVector$stepl[which(stopVector$LM=='N' & stopVector$expt==3)]

  # Landmarks present
dLM1 <- stopVector$stepl[which(stopVector$LM=='Y' & stopVector$expt==1 & stopVector$id %in% 
                                 unique(stopVector$id[which(stopVector$expt==2 & stopVector$LM=='Y')]))]
dLM2 <- stopVector$stepl[which(stopVector$LM=='Y' & stopVector$expt==2)]
dLM3 <- stopVector$stepl[which(stopVector$LM=='Y' & stopVector$expt==3)]


# SAVE PROCESSED STOPS DATA
save(summary12, summary23,
     bird12, bird23,
     stops12, stops23,
     stopsV12NA, stopsV23NA,
     ALM1, ALM2, ALM3,
     AnLM1, AnLM2, AnLM3,
     LM1, LM2, LM3,
     noLM1, noLM2, noLM3,
     dLM1, dLM2, dLM3,
     nodLM1, nodLM2, nodLM3,
     stopNum,
     file = "data/processed-stops.RData")









