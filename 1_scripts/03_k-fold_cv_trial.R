library(rjags)
library(nnet) 
library(MASS)
library(pscl)
library(dplyr)
library(lubridate)
library(stringi)
#packages to run a zero-truncated poisson in R
library(VGAM)
library(boot)
select <- dplyr::select

setwd("C:/Users/Emily/Desktop/R_WORKING/Data/Songscope/HM_ARU_data_06June2018")
BSSC <- read.csv("BSSC_sampleV3.csv")

#####################################################################
### Clean the song count - breeding status data ###
BSSC <- BSSC %>% 
  mutate(BreedingStatus = factor(BreedingStatus, levels = c("Single", 
                                                            "Paired",
                                                            "Feeding Young"))) %>%
  mutate(Date = as.POSIXct(Date)) %>% 
  mutate(ODate = yday(Date)) %>% 
  rename(SongCount = Song_Count) %>% 
  mutate(Cluster = stri_sub(OSFL_ID, from = 1, length = 2)) %>% 
  mutate(Prov = ifelse(Cluster == "AB", "AB", "NT"))

lookup <- data.frame(Sunrise_sample = unique(BSSC$Sunrise_sample),
                     sunrise.time = c(-15, -10, -5, 0, 5, 10, 15, 20, 25))

BSSC<- full_join(BSSC, lookup)

Daily_BSSC_derived <- BSSC %>% 
  group_by(ARU_ID, ODate, BreedingStatus) %>% 
  summarize(DailyMeanSC = mean(SongCount),
            DailyMaxSC = max(SongCount),
            DailyProp0s = sum(SongCount == 0)/9,
            DailyMedian = median(SongCount))

Daily_BSSC_derived <- Daily_BSSC_derived %>% 
  mutate(DailyMeanSC_rounded = round(DailyMeanSC))%>%
  select(DailyMeanSC_rounded, ARU_ID, BreedingStatus, ODate) %>% 
  rename(SongCount = DailyMeanSC_rounded)

#Intermediate step to make a clean data table from the dplyr grouped/summarized data
To_join_ARUs <- Daily_BSSC_derived$ARU_ID
TO_join_ODate <- Daily_BSSC_derived$ODate
To_join_BreedingStatus <- Daily_BSSC_derived$BreedingStatus
To_join_SongCount <- Daily_BSSC_derived$SongCount

Daily_BSSC_derived2 <- data.frame(ARU_ID = To_join_ARUs,
                                  ODate = TO_join_ODate,
                                  BreedingStatus = To_join_BreedingStatus,
                                  SongCount = To_join_SongCount)

Daily_BSSC_derived2_no0s <- Daily_BSSC_derived2 %>% 
  filter(SongCount > 0)

#format ARU data of interest to be used in Bayesian model

cleanit <- function(ds){
  
  cleandata <- ds %>% 
    select(ARU_ID, SongCount, BreedingStatus, ODate) %>% 
    rename(birdID = ARU_ID, 
           songCount = SongCount, 
           breedingStatus = BreedingStatus,
           JDate = ODate)
  
  scaledDate = scale(cleandata$JDate, center = TRUE, scale = TRUE)
  cleandata$JDate <- as.vector(scaledDate)
  
  return(cleandata)
  
}

Daily_BSSC_derived2.clean <- cleanit(Daily_BSSC_derived2)
Daily_BSSC_derived2_no0s.clean <- cleanit(Daily_BSSC_derived2_no0s)

#########First, model zero-tuncated dataset#########
### set MCMC specifications ###

n.ch = 5  # number of chains (5)
n.adpt = 1000  # number of adaptation iterations (1000)
burn.val = 1000 # length of burn in (1000)
n.itr = 1000 # chain length (1000)
dta.cln = Daily_BSSC_derived2_no0s.clean

#Check that the order of BS factor values is FY, P, S
print(levels(dta.cln$breedingStatus))

#If not in correct order, fix order
dta.cln <- dta.cln %>% 
  mutate(breedingStatus = factor(breedingStatus, levels = c("Feeding Young","Paired", "Single")))

# #Test out models to use the 0-truncated poisson (package VGAM, need R>= 3.4.0)
# 
# m <- vglm(formula = songCount ~ 0+breedingStatus, family = pospoisson(), data = dta.cln)
# summ <- summary(m)
# AICvlm(m)
# model2_pois <- glm(songCount ~ 0+breedingStatus, data=dta.cln, family=poisson)
# summary(model2_pois)
# #We can see AIC is lower with the zero-truncated Poisson (AIC = 892), vs the original glm (AIC = 905)

###########################################################################
### Hierarchical model function for zero-truncated data ###

BayesianModel.zT <- function(n.ch, dta, testDat, ntestDat){ 
  trainDat <- dta[-((nrow(dta)-ntestDat+1):nrow(dta)),]
  
  ### get calibration priors from GLM models ###
  { 
    # 1. song count as a function of bs (using zero-tuncated Poisson distribution)
    m <- vglm(formula = songCount ~ 0+breedingStatus, family = pospoisson(), data = trainDat)
    summ <- summary(m)
    
    # mean prior values
    beta0_FY <- summ@coef3[1,1]
    beta0_P <- summ@coef3[2,1]
    beta0_S <- summ@coef3[3,1]

    # SE for prior values
    beta0_FY_SE <- summ@coef3[1,2]
    beta0_P_SE <- summ@coef3[2,2]
    beta0_S_SE <- summ@coef3[3,2]
    
    # 2. bs as a function of date
    trainDat$breedingStatus <- relevel(trainDat$breedingStatus, ref = "Paired")
    model_4 <- multinom(breedingStatus ~ JDate, data=trainDat) # multinom Model
    summ2 <- summary(model_4)
    
    # mean prior values
    eta_FY <- summ2$coefficients[1,1]
    eta_S <- summ2$coefficients[2,1]
    eta_FY_date <- summ2$coefficients[1,2]
    eta_S_date <- summ2$coefficients[2,2]
    # SE prior values
    eta_FY_SE <- summ2$standard.errors[1,1]
    eta_S_SE <- summ2$standard.errors[2,1]
    eta_FY_date_SE <- summ2$standard.errors[1,2]
    eta_S_date_SE <- summ2$standard.errors[2,2]
  }
  
  jagsmodel <- jags.model("model4_0truncpoisson.bug",
                          data = list('songcount' = dta$songCount,
                                      'breedingstatus' = dta$breedingStatus,
                                      'N' = nrow(dta),
                                      'JDate' = dta$JDate,
                                      'beta0_P'= beta0_P, 'beta0_P_SE' = beta0_P_SE,
                                      'beta0_FY' = beta0_FY, 'beta0_FY_SE' = beta0_FY_SE,
                                      'beta0_S' = beta0_S, 'beta0_S_SE' = beta0_S_SE,
                                      'eta_FY' = eta_FY, 'eta_FY_SE' = eta_FY_SE,
                                      'eta_S' = eta_S, 'eta_S_SE' = eta_S_SE,
                                      'eta_FY_date' = eta_FY_date, 'eta_FY_date_SE' = eta_FY_date_SE,
                                      'eta_S_date' = eta_S_date, 'eta_S_date_SE' = eta_S_date_SE),
                          n.chains = n.ch, # number of parallel chains
                          n.adapt = n.adpt) # number of iterations
  update(jagsmodel, 1000) # burn in (gets rid of initial transient behaviour)
  
  out <- coda.samples(jagsmodel,
                      variable.names=c(paste('breedingstatus[',nrow(dta)-ntestDat+1,':',nrow(dta),']',sep=""),
                                       "beta0","eta"),  # parameter vector (length 5,for each breeding status)
                      n.iter=10000)  # could increase chain length from default (10000) if I want
  
  Props <- matrix(nrow = ntestDat, ncol = 3) # to save proportions from traces
  Predicted <- matrix(nrow = ntestDat) # to save predicted from highest proportion
  
  out.vals = numeric()
  for(j in 1:n.ch){  # put values from each chain all together into one vector
    out.vals = rbind(out.vals,out[[j]][,1:ntestDat + 3])
  }
  
  for(k in 1:ntestDat){  # for each node we are monitoring (for the bird we set aside as a test bird)
    for(l in 1:3){  # for each breeding status
      if(ntestDat == 1){ # if ntestDat = 1, dimensions of vectors just require changed syntax
        Props[k,l] = sum(out.vals == l)/length(out.vals)
      } else {
        Props[k,l] = sum(out.vals[,k] == l)/nrow(out.vals)
      }
    }
    Predicted[k] = which(Props[k,] == max(Props[k,]))
  }
  return(list(Predicted, Props))
}

###########################################################################
### k-fold validation (for each bird) ###

# create prediction storage matrix
df <- data.frame(matrix(ncol = 8, nrow = 0)) # data frame to store all predictions
colnames(df) <- c("birdID","songCount", "breedingStatus","JDate", "Predicted","Prob.FY","Prob.P","Prob.S")

for (testBird in unique(dta.cln$birdID)){ 
  print(testBird) # for troubleshooting purposes
  dta <- dta.cln # temporary dta file to play with
  testDat <- dta[which(dta$birdID == testBird),]
  ntestDat <- nrow(testDat)
  dta <- rbind(dta[-which(dta$birdID == testBird),], testDat)  # put test bird at bottom of dataframe
  dta$breedingStatus[(nrow(dta)-ntestDat+1):nrow(dta)] <- NA   # set test bird breeding statuses to NA 
  
  Mod.output <- BayesianModel.zT(n.ch, dta, testDat, ntestDat)
  Predicted <- Mod.output[[1]]
  Props <- Mod.output[[2]]
  
  #put into data frame of all predictions ######
  Predicted[which(Predicted == 1)] <- "Feeding Young"
  Predicted[which(Predicted == 2)] <- "Paired"
  Predicted[which(Predicted == 3)] <- "Single"
  Predicted <- factor(Predicted, levels = c("Feeding Young", "Paired", "Single")) # reorder factor levels
  
  Predicted <- as.factor(Predicted)
  df.temp <- cbind(testDat,Predicted,Props)
  colnames(df.temp) <- c("birdID","songCount", "breedingStatus", "JDate","Predicted","Prob.FY","Prob.P","Prob.S")
  df <- rbind(df,df.temp)
}

##############################################################################################################
### output total proportion of breeding statuses correctly predicted ###

print(sum(df$breedingStatus == df$Predicted)/nrow(df))
table(df$breedingStatus,df$Predicted)
library("caret")
confusionMatrix(df$Predicted, df$breedingStatus)

write.csv(df, "BS_preds_Daily_meanSC_no0s.csv")

###############Second, model dataset with song counts including zeros (regular Poisson)###########################
### set MCMC specifications ###

n.ch = 5  # number of chains (5)
n.adpt = 1000  # number of adaptation iterations (1000)
burn.val = 1000 # length of burn in (1000)
n.itr = 1000 # chain length (1000)
dta.cln = Daily_BSSC_derived2.clean

#Check that the order of BS factor values is FY, P, S
print(levels(dta.cln$breedingStatus))

#If not in correct order, fix
dta.cln <- dta.cln %>% 
  mutate(breedingStatus = factor(breedingStatus, levels = c("Feeding Young","Paired", "Single")))

###########################################################################
### Hierarchical model function for zero-truncated data ###

BayesianModel.notime <- function(n.ch, dta, testDat, ntestDat){ 
  trainDat <- dta[-((nrow(dta)-ntestDat+1):nrow(dta)),]
  
  ### get calibration priors from GLM models ###
  { 
    # 1. song count as a function of bs and time
    model2_pois <- glm(songCount ~ 0+breedingStatus, data=trainDat, family=poisson)
    summ <- summary(model2_pois)
    
    # mean prior values
    beta0_FY <- summ$coefficients[1,1]
    beta0_P <- summ$coefficients[2,1]
    beta0_S <- summ$coefficients[3,1]
    # SE for prior values
    beta0_FY_SE <- summ$coefficients[1,2]
    beta0_P_SE <- summ$coefficients[2,2]
    beta0_S_SE <- summ$coefficients[3,2]
    
    # 2. bs as a function of date
    trainDat$breedingStatus <- relevel(trainDat$breedingStatus, ref = "Paired")
    model_4 <- multinom(breedingStatus ~ JDate, data=trainDat) # multinom Model
    summ2 <- summary(model_4)
    
    # mean prior values
    eta_FY <- summ2$coefficients[1,1]
    eta_S <- summ2$coefficients[2,1]
    eta_FY_date <- summ2$coefficients[1,2]
    eta_S_date <- summ2$coefficients[2,2]
    # SE prior values
    eta_FY_SE <- summ2$standard.errors[1,1]
    eta_S_SE <- summ2$standard.errors[2,1]
    eta_FY_date_SE <- summ2$standard.errors[1,2]
    eta_S_date_SE <- summ2$standard.errors[2,2]
  }
  
  jagsmodel <- jags.model("model4_poisson_notime.bug",
                          data = list('songcount' = dta$songCount,
                                      'breedingstatus' = dta$breedingStatus,
                                      'N' = nrow(dta),
                                      'JDate' = dta$JDate,
                                      'beta0_P'= beta0_P, 'beta0_P_SE' = beta0_P_SE,
                                      'beta0_FY' = beta0_FY, 'beta0_FY_SE' = beta0_FY_SE,
                                      'beta0_S' = beta0_S, 'beta0_S_SE' = beta0_S_SE,
                                      'eta_FY' = eta_FY, 'eta_FY_SE' = eta_FY_SE,
                                      'eta_S' = eta_S, 'eta_S_SE' = eta_S_SE,
                                      'eta_FY_date' = eta_FY_date, 'eta_FY_date_SE' = eta_FY_date_SE,
                                      'eta_S_date' = eta_S_date, 'eta_S_date_SE' = eta_S_date_SE),
                          n.chains = n.ch, # number of parallel chains
                          n.adapt = n.adpt) # number of iterations
  update(jagsmodel, 1000) # burn in (gets rid of initial transient behaviour)
  
  out <- coda.samples(jagsmodel,
                      variable.names=c(paste('breedingstatus[',nrow(dta)-ntestDat+1,':',nrow(dta),']',sep=""),
                                       "beta0","eta"),  # parameter vector (length 5,for each breeding status)
                      n.iter=10000)  # could increase chain length from default (10000) if I want
  
  Props <- matrix(nrow = ntestDat, ncol = 3) # to save proportions from traces
  Predicted <- matrix(nrow = ntestDat) # to save predicted from highest proportion
  
  out.vals = numeric()
  for(j in 1:n.ch){  # put values from each chain all together into one vector
    out.vals = rbind(out.vals,out[[j]][,1:ntestDat + 3])
  }
  
  for(k in 1:ntestDat){  # for each node we are monitoring (for the bird we set aside as a test bird)
    for(l in 1:3){  # for each breeding status
      if(ntestDat == 1){ # if ntestDat = 1, dimensions of vectors just require changed syntax
        Props[k,l] = sum(out.vals == l)/length(out.vals)
      } else {
        Props[k,l] = sum(out.vals[,k] == l)/nrow(out.vals)
      }
    }
    Predicted[k] = which(Props[k,] == max(Props[k,]))
  }
  return(list(Predicted, Props))
}

###########################################################################
### k-fold validation (for each bird) ###

# create prediction storage matrix
df <- data.frame(matrix(ncol = 8, nrow = 0)) # data frame to store all predictions
colnames(df) <- c("birdID","songCount", "breedingStatus","JDate", "Predicted","Prob.FY","Prob.P","Prob.S")

for (testBird in unique(dta.cln$birdID)){ 
  print(testBird) # for troubleshooting purposes
  dta <- dta.cln # temporary dta file to play with
  testDat <- dta[which(dta$birdID == testBird),]
  ntestDat <- nrow(testDat)
  dta <- rbind(dta[-which(dta$birdID == testBird),], testDat)  # put test bird at bottom of dataframe
  dta$breedingStatus[(nrow(dta)-ntestDat+1):nrow(dta)] <- NA   # set test bird breeding statuses to NA 
  
  Mod.output <- BayesianModel.notime(n.ch, dta, testDat, ntestDat)
  Predicted <- Mod.output[[1]]
  Props <- Mod.output[[2]]
  
  #put into data frame of all predictions ######
  Predicted[which(Predicted == 1)] <- "Feeding Young"
  Predicted[which(Predicted == 2)] <- "Paired"
  Predicted[which(Predicted == 3)] <- "Single"
  Predicted <- factor(Predicted, levels = c("Feeding Young", "Paired", "Single")) # reorder factor levels
  
  Predicted <- as.factor(Predicted)
  df.temp <- cbind(testDat,Predicted,Props)
  colnames(df.temp) <- c("birdID","songCount", "breedingStatus", "JDate","Predicted","Prob.FY","Prob.P","Prob.S")
  df <- rbind(df,df.temp)
}

##############################################################################################################
### output total proportion of breeding statuses correctly predicted ###

print(sum(df$breedingStatus == df$Predicted)/nrow(df))
table(df$breedingStatus,df$Predicted)
library("caret")
confusionMatrix(df$Predicted, df$breedingStatus)

write.csv(df, "BS_preds_Daily_meanSC.csv")

