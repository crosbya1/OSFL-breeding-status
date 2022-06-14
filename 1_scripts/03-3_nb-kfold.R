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


# ======================================================
# Find the number of iterations needed for convergence on a single model run 
i = 1
train <- which(leavout != i)
test <- which(leavout == i)

sc.train <- daily_sc[train, ]
sc.train$data_bs <- sc.train$standardized_bs
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



params <- c("eta", "beta0", "beta1", paste("p[", nrow(sc.dat)-ntest+1, ":", nrow(sc.dat), ",1:3]", sep = ""), 
            paste('breedingstatus[',nrow(sc.dat)-ntest+1,':',nrow(sc.dat),']',sep=""))


ni = 50000
nb = 30000
nt = 5

system.time({
  out <- jags(data = data, parameters.to.save = params, 
              model.file = "1_scripts/model-scripts/negbinom_basic.txt", n.chains = nc, n.iter = ni, 
              n.burnin = nb, n.thin = nt, parallel = T)
})
#Full convergence in ~15 minutes with 20k iterations

negbinom_basic_cv <- lapply(1:k, function(i){
  
  print(paste0("************************* Fold ", i, " *********************"))
  
  train <- which(leavout != i)
  test <- which(leavout == i)
  
  sc.train <- daily_sc[train, ]
  sc.train$data_bs <- sc.train$standardized_bs
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
  ntrain <- length(train)
  
  # Run the basic version of the model
  data <- list(N = nrow(y.dat), breedingstatus = as.numeric(as.factor(sc.dat$standardized_bs)), 
               nsamp = sc.dat$n_samples_per_day, y = y.dat[, -1], JDate = j.dat$jday, JDate2 = j.dat$jday2, 
               TimeRel2Sun = t.dat[, -1]/60/60, ss = as.numeric(as.factor(sc.dat$SS)), n.ss = nlevels(as.factor(sc.dat$SS)))
  
  
  
  params <- c("eta", "beta0", "beta1", paste("p[", nrow(sc.dat)-ntest+1, ":", nrow(sc.dat), ",1:3]", sep = ""), 
              paste('breedingstatus[',nrow(sc.dat)-ntest+1,':',nrow(sc.dat),']',sep=""))
  
  
  
  system.time({
    out <- jags(data = data, parameters.to.save = params, 
                model.file = "1_scripts/model-scripts/negbinom_basic.txt", n.chains = nc, n.iter = ni, 
                n.burnin = nb, n.thin = nt, parallel = T)
  })
  
  # Retrieve the values we need 
  loglik.birdday <- lapply(1:nrow(out$sims.list$p), function(iter){
    # Calculate the probabilities of being in each breeding status
    p.test <- (nrow(sc.dat)-ntest+1):nrow(sc.dat)
    pb <- out$sims.list$p[iter, p.test, ]
    
    # Calculate expected song count at each survey for each breeding status, based on the model
    lam.test <- lapply(p.test, function(x){
      ns <- 1:data$nsamp[x]
      lam <- do.call(rbind, lapply(ns, function(i){
        unlist(sapply(1:3, function(j) exp(out$sims.list$beta0[iter, j] + out$sims.list$beta1[iter]*setDF(data$TimeRel2Sun)[x, i]))) 
      }))
      return(lam)
    })
    
    lik.bd <- lapply(1:length(lam.test), function(t){
      # likelihood of survey-specific predictions, given the data and model predictions
      lik <- do.call(rbind, lapply(1:nrow(lam.test[[t]]), function(k){
        # Probability of the song count data, given expected song counts from the model, for each breeding status
        pl <- dpois(setDF(data$y)[p.test[t], k], lam.test[[t]][k, ])
        # Probability of each breeding status, given model predictions (Equation 1 in the manuscript)
        pbs <- pl*pb[t]/sum(pl*pb[t])
        return(pbs)
      }))
      b <- rep(0, 3)
      b[as.numeric(as.factor(sc.test$data_bs))[t]] <- 1
      
      # Sum the log likelihoods of data|model for each survey in a day 
      lik.day <- sum(sapply(1:nrow(lik), function(l) dmultinom(b, prob = lik[l, ], log = TRUE)))
      # For each test bird-day, return the joint breeding status/singing probabilities (lik), the actual breeding status data (b), and the bird-day log likelihood (lik.day)
      return(list(lik = lik, b = b, lik.day = lik.day))
    })
    return(lik.bd)
  })
  
  # Get the full model likelihood as the sum of bird-day likelihoods (lik.day) averaged over all iterations
  modlik.birdday <- mean(sapply(1:length(loglik.birdday), function(i){
    -2*sum(sapply(1:length(loglik.birdday[[i]]), function(j) loglik.birdday[[i]][[j]]$lik.day))
  }))
  
  # Get the data for each test bird-year
  dta.bs <- do.call(rbind, lapply(1:length(loglik.birdday[[i]]), function(j) loglik.birdday[[i]][[j]]$b))
  
  # Calculate the expected proportion of bird-days in each category, according to the model 
  bs.pred <- lapply(1:length(loglik.birdday), function(i){
    d.p <- do.call(rbind, lapply(1:length(loglik.birdday[[i]]), function(j){
      t <- loglik.birdday[[i]][[j]]$lik
      t1 <- do.call(cbind, lapply(1:nrow(t), function(k) rmultinom(1, 1, t[k, ])))
      return(t(rmultinom(1, 1, rowMeans(t1))))
    }))
  })
  
  loglik.prop <- mean(-2*sapply(1:length(bs.pred), function(i) dmultinom(colSums(dta.bs), prob = colMeans(bs.pred[[i]]), log = TRUE)))
  
  Props <- matrix(nrow = ntest, ncol = 3) # to save proportions from traces
  Predicted <- matrix(nrow = ntest) # to save predicted from highest proportion
  
  for(t in 1:ntest){  # for each test bird-day
    for(l in 1:3){  # for each breeding status (FY, P, S)
      if(ntest == 1){ # if ntestDat = 1, dimensions of vectors just require changed syntax
        Props[t,l] = sum(out.vals == l)/length(out.vals)
      }else{
        Props[t,l] = sum(out$sims.list$breedingstatus[, (t + ntrain)] == l)/nrow(out$sims.list$breedingstatus)
      }
    }
    Predicted[t] = which(Props[t,] == max(Props[t,]))
  } 
  
  #put into data frame of all predictions ######
  Predicted[which(Predicted == 1)] <- "Feeding Young"
  Predicted[which(Predicted == 2)] <- "Paired"
  Predicted[which(Predicted == 3)] <- "Single"
  Predicted <- factor(Predicted, levels = c("Feeding Young", "Paired", "Single")) # reorder factor levels
  
  # Predicted <- as.factor(Predicted)
  df.temp <- data.frame(sc.test[, c("SS", "total_songs_per_day", "data_bs", "JDate")], Predicted, Props)
  colnames(df.temp) <- c("birdID", "totalSongsPerDay", "breedingStatus", "JDate","Predicted","Prob.FY","Prob.P","Prob.S")
  
  
  return(list(out = out, df = df.temp, loglik.birdday = loglik.birdday, modlik.birdday = modlik.birdday, dta.bs = dta.bs, bs.pred = bs.pred, loglik.prop = loglik.prop))
  
})

modsum <- do.call(rbind, lapply(1:length(poisson_basic_cv), function(x) data.frame(modlik = poisson_basic_cv[[x]]$modlik.birdday, proplik = poisson_basic_cv[[x]]$loglik.prop)))
score <- colMeans(modsum)
