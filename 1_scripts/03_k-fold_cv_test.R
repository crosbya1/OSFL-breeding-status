# R code to run k-fold cross validation for modeling 

# Load required package
library(dplyr)
library(tidyr)
library(data.table)
library(jagsUI)


load("0_data/processed/osfl_data_package.Rdata")


### set MCMC specifications ###

nc = 3  # number of chains (5)
nb = 500  # number of adaptation iterations (1000)
nt = 1
ni = 1000 # chain length (1000)

#Check that the order of BS factor values is FY, P, S
print(levels(as.factor(daily_sc$standardized_bs)))


k <- 5 # the number of folds


bird_day <- nrow(daily_sc)
nk <- 1/k
leavout <- sample(rep(1:(bird_day/(bird_day*nk)), ceiling(bird_day*nk))[1:bird_day])

for(i in 1:k){
  train <- which(leavout != k)
  test <- which(leavout == k)
  
  sc.train <- daily_sc[train, ]
  sc.test <- daily_sc[test, ]
  sc.test$standardized_bs <- NA
  sc.dat <- rbind(sc.train, sc.test)
  
  y.train <- y[train, ]
  y.dat <- rbind(y.train, y[test, ])
  
  j.train <- Jday[train, ]
  j.dat <- rbind(j.train, Jday[test, ])
  
  t.train <- TimeRel2Sun[train, ]
  t.dat <- rbind(t.train, TimeRel2Sun[test, ])
  
  
  
  
  # Run the basic version of the model
  data <- list(N = nrow(y.dat), breedingstatus = as.numeric(as.factor(sc.dat$standardized_bs)), 
               nsamp = sc.dat$n_samples_per_day, y = y.dat[, -1], JDate = j.dat$jday, JDate2 = j.dat$jday2, 
               TimeRel2Sun = t.dat[, -1]/60/60, ss = as.numeric(as.factor(sc.dat$SS)), n.ss = nlevels(as.factor(sc.dat$SS)))
  
  params <- c("eta", "beta0", "beta1", "breedingstatus")
  
  system.time({
    out.poisson.basic <- jags.basic(data = data, parameters.to.save = params, 
                              model.file = "1_scripts/model-scripts/poisson_basic.txt", n.chains = 3, n.iter = 1000, 
                              n.burnin = 500, n.thin = 1, parallel = T)
  })
  
  Props <- matrix(nrow = ntestDat, ncol = 3) # to save proportions from traces
  Predicted <- matrix(nrow = ntestDat) # to save predicted from highest proportion
  Likelihood <- matrix(nrow = ntestDat) # to save predicted from highest proportion
  
  out.vals = numeric()
  for(j in 1:n.ch){  # put values from each chain all together into one vector
    out.vals = rbind(out.vals,out[[j]][,1:ntestDat + 4])
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
    b <- rep(0, 3)
    b[as.numeric(testDat$breedingStatus)[k]] <- 1
    Likelihood[k] = dmultinom(b, prob = Props[k, ])
    modlike <- dmultinom(table(testDat$breedingStatus), prob = table(as.factor(out.vals)))
    lprob <- -2*log(modlike)
  }
  
  
}



