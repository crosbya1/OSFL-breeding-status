# R code to run k-fold cross validation for modeling 

# Load required packages
library(tidyverse)
library(data.table)
library(jagsUI)
library(coda)

load("0_data/processed/osfl_data_package.Rdata")

### set MCMC specifications ###

nc = 3  # number of chains
nb = 500  # number of iterations at the beginning of the chain to discard (i.e., the burn-in)
nt = 1  # thinning rate
ni = 1000 # number of iterations per chain (including burn-in)

#Check that the order of BS factor values is FY, P, S
print(levels(as.factor(daily_sc$standardized_bs)))

k <- 5 # the number of folds

#Create k-fold cross validation groupings (i.e., leave out) by randomly assigning bird-days to 1 of 5 folds
bird_day <- nrow(daily_sc)
nk <- 1/k
leavout <- sample(rep(1:(bird_day/(bird_day*nk)), ceiling(bird_day*nk))[1:bird_day])

for(i in 1:k){
  train <- which(leavout != i)
  test <- which(leavout == i)
  
  sc.train <- daily_sc[train, ]
  sc.test <- daily_sc[test, ]
  sc.test$data_bs <- sc.test$standardized_bs
  sc.test$standardized_bs <- NA
  sc.dat <- rbind(sc.train, sc.test)
  
  y.train <- y[train, ]
  y.dat <- rbind(y.train, y[test, ])
  
  j.train <- Jday[train, ]
  j.dat <- rbind(j.train, Jday[test, ])
  
  t.train <- TimeRel2Sun[train, ]
  t.dat <- rbind(t.train, TimeRel2Sun[test, ])
  
  ntest <- length(test)
  
  # Run the basic version of the model
  data <- list(N = nrow(y.dat), breedingstatus = as.numeric(as.factor(sc.dat$standardized_bs)), 
               nsamp = sc.dat$n_samples_per_day, y = y.dat[, -1], JDate = j.dat$jday, JDate2 = j.dat$jday2, 
               TimeRel2Sun = t.dat[, -1]/60/60, ss = as.numeric(as.factor(sc.dat$SS)), n.ss = nlevels(as.factor(sc.dat$SS)))
  
  params <- c("eta", "beta0", "beta1", paste('breedingstatus[',nrow(sc.dat)-ntest+1,':',nrow(sc.dat),']',sep=""))  # model parameters to be monitored
  
  params <- c("eta", "beta0", "beta1", "p")

  params <- c("eta", "beta0", "beta1", paste("p[", nrow(sc.dat)-ntest+1, ":", nrow(sc.dat), ",1:3]", sep = ""))
  
  
  
  system.time({
    out <- jags(data = data, parameters.to.save = params, 
                              model.file = "1_scripts/model-scripts/poisson_basic.txt", n.chains = nc, n.iter = ni, 
                              n.burnin = nb, n.thin = nt, parallel = T)
  })
  
  Props <- matrix(nrow = ntest, ncol = 3) # to save proportions from traces
  Predicted <- matrix(nrow = ntest) # to save predicted from highest proportion
  Likelihood <- matrix(nrow = ntest) # to save likelihood of the model given the data
  
  out.vals = numeric()
  for(j in 1:nc){  # put values from each chain all together into one vector
    out.vals = rbind(out.vals,out.poisson.basic[[j]][,1:ntest + 4])
  }
  
  for(t in 1:ntest){  # for each test bird-day
    for(l in 1:3){  # for each breeding status (FY, P, S)
      if(ntest == 1){ # if ntestDat = 1, dimensions of vectors just require changed syntax
        Props[t,l] = sum(out.vals == l)/length(out.vals)
      } else {
        Props[t,l] = sum(out.vals[,t] == l)/nrow(out.vals)
      }
    }
    Predicted[t] = which(Props[t,] == max(Props[t,]))
    
    b <- rep(0, 3)
    b[as.numeric(as.factor(sc.test$data_bs))[t]] <- 1
    Likelihood[t] = dmultinom(b, prob = Props[t, ])

    df.temp <- cbind(sc.test$SS[t], sc.test$JDate[t], sc.test$data_bs[t], Predicted[t], Props[t,1], Props[t,2], Props[t,3], Likelihood[t]) 
    df.temp <- data.frame(df.temp)
    colnames(df.temp) <- c("SS","JDate","breedingStatus","Predicted","Prob.FY","Prob.P","Prob.S","Likelihood")
    
    if (t==1) {
      df = df.temp
    } else {
      df = rbind(df,df.temp)
    }
    
  }

  if (i==1) {
    df_predictions = df
  } else {
    df_predictions = rbind(df_predictions,df)
  }

  
modlike <- dmultinom(table(sc.test$data_bs), prob = table(as.factor(out.vals))) # likelihood of all the data, given the model
lprob <- -2*log(modlike) # -2 log likelihood of all the data (i.e., deviance)

df_lik.temp <- cbind(i, round(modlike,digits = 5), round(lprob, digits = 5))
df_lik.temp <- data.frame(df_lik.temp)
colnames(df_lik.temp) <- c("k-fold","Modlike", "lprob")

if (i==1) {
  df_lik = df_lik.temp
} else {
  df_lik = rbind(df_lik, df_lik.temp)
}

print(paste0("Fold ",i," of ",k, " Completed"))
  
}



