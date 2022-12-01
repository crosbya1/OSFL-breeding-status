# R code to run best model selected through cross validation 

# Load required packages
library(tidyverse)
library(data.table)
library(jagsUI)
library(coda)

load("0_data/processed/osfl_data_package.Rdata")

### set MCMC specifications ###

#Check that the order of BS factor values is FY, P, S
print(levels(as.factor(daily_sc$standardized_bs)))


# ======================================================
# Find the number of iterations needed for convergence on a single model run 


sc.dat <- daily_sc
sc.dat$data_bs <- sc.dat$standardized_bs

y.dat <- y

j.dat <- Jday

t.dat <- TimeRel2Sun


# Run the basic version of the model
data <- list(N = nrow(y.dat), breedingstatus = as.numeric(as.factor(sc.dat$standardized_bs)), 
             nsamp = sc.dat$n_samples_per_day, y = y.dat[, -1], JDate = j.dat$jday, JDate2 = j.dat$jday2, 
             TimeRel2Sun = t.dat[, -1]/60/60, ss = as.numeric(as.factor(sc.dat$SS)), n.ss = nlevels(as.factor(sc.dat$SS)))



params <- c("eta", "beta0", "beta1", "r", "p")

nc = 3
ni = 50000
nb = 30000
nt = 5

system.time({
  out <- jags(data = data, parameters.to.save = params, 
              model.file = "1_scripts/model-scripts/negbinom_basic.txt", n.chains = nc, n.iter = ni, 
              n.burnin = nb, n.thin = nt, parallel = T)
})
#Full convergence with 50k iterations


# The observed breeding status in matrix format
b <- data.frame(matrix(0, nrow(sc.dat), 3))

for(x in 1:nrow(b)){
  b[x, as.numeric(as.factor(sc.dat$standardized_bs))[x]] <- 1
} 

b$day <- sc.dat$JDate

status.day <- data.frame(b %>% group_by(day) %>% summarize(FY = sum(X1), P = sum(X2), S = sum(X3)))
status.prop <- ((status.day %>% select(-day))/rowSums(status.day %>% select(-day))) %>% filter(table(sc.dat$JDate) >= 5)

ggplot(data = status.day) + geom_line(aes(x = day, y = FY), col = "blue") + geom_line(aes(x = day, y = P), col = "red") + geom_line(aes(x = day, y = S), col = "green")

pr <- out$sims.list$r


loglik.birdday <- lapply(1:nrow(out$sims.list$p), function(iter){
  # Calculate the probabilities of being in each breeding status
  #p.test <- (nrow(sc.dat)-ntest+1):nrow(sc.dat)
  pb <- out$sims.list$p[iter, , ]
  
  # Calculate expected song count at each survey for each breeding status, based on the model
  lam.test <- lapply(1:nrow(pb), function(x){
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
      pl <- dnbinom(setDF(data$y)[t, k], size = pr[iter], mu = lam.test[[t]][k, ])
      # Probability of each breeding status, given model predictions (Equation 1 in the manuscript)
      pbs <- pl*pb[t, ]/sum(pl*pb[t, ])
      return(pbs)
    }))
    b <- rep(0, 3)
    b[as.numeric(as.factor(sc.dat$data_bs))[t]] <- 1
    
    # Sum the log likelihoods of data|model for each survey in a day 
    lik.day <- sum(sapply(1:nrow(lik), function(l) dmultinom(b, prob = lik[l, ], log = TRUE)))
    # For each test bird-day, return the joint breeding status/singing probabilities (lik), the actual breeding status data (b), and the bird-day log likelihood (lik.day)
    return(list(lik = lik, b = b, lik.day = lik.day))
  })
  return(lik.bd)
})

save(out, loglik.birdday, file = "2_pipeline/store/nb-final-output.Rdata")

