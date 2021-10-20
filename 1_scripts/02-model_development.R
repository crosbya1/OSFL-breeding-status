# R code for developing different models of OSFL breeding status as a function of Julian day and 
# singing rate as a funciton of breeding status

rm(list=ls(all=T))
getwd()

# Load required package
library(dplyr)
library(tidyr)
library(data.table)
library(jagsUI)


load("0_data/processed/osfl_data_package.Rdata")

# The basic Poison model (no random effects or quadratics)

sink("poisson_basic.txt")
cat("
     
# note: the order of breeding states is: FY, P, S

model{
  ### Likelihood
  for (i in 1:N){

    breedingstatus[i] ~ dcat(p[i,1:3]) 
    for (j in 1:3){
      p[i,j] <- delta[i,j]/sum(delta[i,])
      log(delta[i,j]) = alpha[i,j]
    }	
    alpha[i,1] <- eta[1] + eta[4]*JDate[i]
    alpha[i,2] <- eta[2]
    alpha[i,3] <- eta[3] + eta[5]*JDate[i]
	
	  for(j in 1:nsamp[i]){
	    y[i, j] ~ dpois(lambda[i, j]) 
	    log(lambda[i, j]) <- beta0[breedingstatus[i]] + beta1*TimeRel2Sun[i, j]
	    }
    }
  
  # Priors
  # etas are the log odds of each state (with paired as the reference here)
  eta[1] ~ dnorm(0, 0.01)
  eta[2] ~ dnorm(0, 1000)
  eta[3] ~ dnorm(0, 0.01)
  eta[4] ~ dnorm(0, 0.01)
  eta[5] ~ dnorm(0, 0.01)
  
  # exp(beta0) is the expected song count (for each state)
  beta0[1] ~ dnorm(0, 0.01)
  beta0[2] ~ dnorm(0, 0.01)
  beta0[3] ~ dnorm(0, 0.01) 

  beta1 ~ dnorm(0, 0.01)
}

", fill = TRUE)
sink()



# Run the basic version of the model
data <- list(N = nrow(y), breedingstatus = as.numeric(as.factor(daily_sc$standardized_bs)), 
             nsamp = daily_sc$n_samples_per_day, y = y[, -1], JDate = Jday$jday, JDate2 = Jday$jday2, 
             TimeRel2Sun = TimeRel2Sun[, -1]/60/60, ss = as.numeric(as.factor(daily_sc$SS)), n.ss = nlevels(as.factor(daily_sc$SS)))

params <- c("eta", "beta0", "beta1")

system.time({
  out.poisson.basic <- jags(data = data, parameters.to.save = params, 
                            model.file = "1_scripts/model-scripts/poisson_basic.txt", n.chains = 3, n.iter = 1000, 
                            n.burnin = 500, n.thin = 1, parallel = T)
})



# Add a random effect of bird-year
sink("poisson_re.txt")
cat("
     
# note: the order of breeding states is: FY, P, S

model{
  ### Likelihood
  for (i in 1:N){

    breedingstatus[i] ~ dcat(p[i,1:3]) 
    for (j in 1:3){
      p[i,j] <- delta[i,j]/sum(delta[i,])
      log(delta[i,j]) = alpha[i,j]
    }	
    alpha[i,1] <- eta[1] + eta[4]*JDate[i]
    alpha[i,2] <- eta[2]
    alpha[i,3] <- eta[3] + eta[5]*JDate[i]
	
	  for(j in 1:nsamp[i]){
	    y[i, j] ~ dpois(lambda[i, j]) 
	    log(lambda[i, j]) <- beta0[breedingstatus[i]] + beta1*TimeRel2Sun[i, j] + eps.ss[ss[i]]
	    }
    }
  
  # Priors
  # etas are the log odds of each state (with paired as the reference here)
  eta[1] ~ dnorm(0, 0.01)
  eta[2] ~ dnorm(0, 1000)
  eta[3] ~ dnorm(0, 0.01)
  eta[4] ~ dnorm(0, 0.01)
  eta[5] ~ dnorm(0, 0.01)
  
  # exp(beta0) is the expected song count (for each state)
  beta0[1] ~ dnorm(0, 0.01)
  beta0[2] ~ dnorm(0, 0.01)
  beta0[3] ~ dnorm(0, 0.01) 

  beta1 ~ dnorm(0, 0.01)
  
  sd.ss ~ dunif(0, 5)
  tau.ss <- pow(sd.ss, -2)
  for(i in 1:n.ss){
    eps.ss[i] ~ dnorm(0, tau.ss)
	}
}


", fill = TRUE)
sink()

params.re <- c("eta", "beta0", "beta1", "sd.ss")


system.time({
  out.poisson.re <- jags(data = data, parameters.to.save = params.re, model.file = "1_scripts/model-scripts/poisson_re.txt", 
                 n.chains = 3, n.iter = 1000, n.burnin = 500, n.thin = 1, parallel = T)
})

params.zip <- c("eta", "beta0", "beta1", "psi")

zinit <- as.matrix(y %>% dplyr::select(-pkey))
zinit[zinit > 1] <- 1
inits <- function(){list(z = zinit)}

system.time({
  out.zip <- jags(data = data, inits = inits, parameters.to.save = params.zip, model.file = "1_scripts/model-scripts/zip.txt", 
                  n.chains = 3, n.iter = 1000, n.burnin = 500, n.thin = 1, parallel = T)
})

system.time({
  out.zip.re <- jags(data = data, inits = inits, parameters.to.save = params.zip, model.file = "1_scripts/model-scripts/zip_re.txt", 
                  n.chains = 3, n.iter = 1000, n.burnin = 500, n.thin = 1, parallel = T)
})

# It seems like the zero inflation and the random effect are both trying to measure the same thing, 
# so if we put them together in the same model it will confuse the effects and not be able to estimate the betas 

params.nb <- c("eta", "beta0", "beta1", "r")

system.time({
  out.nb.basic <- jags(data = data, inits = inits, parameters.to.save = params.nb, model.file = "1_scripts/model-scripts/negbinom_basic.txt", 
                           n.chains = 3, n.iter = 1000, n.burnin = 500, n.thin = 1, parallel = T)
})


system.time({
  out.poisson.quad <- jags(data = data, inits = inits, parameters.to.save = params.zip, model.file = "1_scripts/model-scripts/poisson_quad.txt", 
                     n.chains = 3, n.iter = 1000, n.burnin = 500, n.thin = 1, parallel = T)
})

system.time({
  out.zip.quad <- jags(data = data, inits = inits, parameters.to.save = params.zip, model.file = "1_scripts/model-scripts/zip_quad.txt", 
                           n.chains = 3, n.iter = 1000, n.burnin = 500, n.thin = 1, parallel = T)
})


system.time({
  out.negbinom.quad <- jags(data = data, inits = inits, parameters.to.save = params.zip, model.file = "1_scripts/model-scripts/negbinom_quad.txt", 
                       n.chains = 3, n.iter = 1000, n.burnin = 500, n.thin = 1, parallel = T)
})

system.time({
  out.nb.quad <- jags(data = data, inits = inits, parameters.to.save = params.nb, model.file = "1_scripts/model-scripts/negbinom_quad.txt", 
                       n.chains = 3, n.iter = 1000, n.burnin = 500, n.thin = 1, parallel = T)
})
