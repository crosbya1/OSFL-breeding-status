
library(jagsUI)
library(tidyverse)
library(ggplot2)


load("2_pipeline/store/nb-final-output.Rdata")
load("0_data/processed/osfl_data_package.Rdata")

# Get the observed breeding status at each bird day
dta.bs <- matrix(0, nrow(sc.dat), 3)
for(i in 1:nrow(dta.bs)){
  dta.bs[i, as.numeric(as.factor(sc.dat$data_bs))[i]] <- 1
}

# Break bird-days into 13 5-day bins
dnum <- data.frame(dnum = as.numeric(as.factor(daily_sc$JDate)))

bins <- dnum %>% mutate(cuts = cut(dnum, c(seq(0, 60, 5), Inf))) %>% 
  group_by(cuts)

bins$bin <- as.numeric(as.factor(bins$cuts))

bs.pred <- lapply(1:nlevels(as.factor(bins$bin)), function(i){
  d.p <- do.call(rbind, lapply(1:length(loglik.birdday), function(j){
    d <- loglik.birdday[[j]][which(bins$bin == i)]
    p <- do.call(rbind, lapply(1:length(d), function(k) t(rmultinom(1, 1,  colMeans(d[[k]]$lik)))))
    return(colMeans(p))
  }))
})

save(bs.pred, dta.bs, dnum, file = "2_pipeline/store/nb-bins-prediction.Rdata")

# Format for plot 

source("1_scripts/functions/plot-functions.R")

bs.plot <- do.call(rbind, lapply(1:length(bs.pred), function(x){
  data.frame(bs.pred[[x]], x)
}))
colnames(bs.plot) <- c("fy", "p", "s", "bin")

head(bs.plot)

bs.plot.long <- gather(bs.plot, status, prob, fy:s, factor_key = TRUE)

bs.plot.sum <- summaryCRI(bs.plot.long, measurevar = "prob", cri = 0.95, groupvars = c("bin", "status"))

data.bs.wide <- data.frame(dta.bs, bins$bin)
colnames(data.bs.wide) <- c("fy", "p", "s", "bin")
data.bs.long <- gather(data.bs.wide, status, obs, fy:s, factor_key = TRUE)
data.bs.sum <- data.bs.long %>% group_by(status, bin) %>% dplyr::summarise(prop = mean(obs))


bs_prob_plot <- ggplot(bs.plot.sum, aes(x = bin, y = prob, colour = status)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_line() +
  geom_point() + geom_point(data = data.bs.sum, aes(x = bin, y = prop, colour = status), size = 5) + 
  theme_bw() +labs(x = "Time period", y = "Proportion of \nsamples") + 
  scale_colour_discrete(name = "Breeding \nstatus", labels = c("Feeding young", "Paired", "Single")) +
  scale_x_continuous(breaks = seq(1,13,1)) +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 15))






# -----------------------------------------------------------------------------------------------
 
bs.pred <- lapply(1:length(loglik.birdday), function(i){
  # Calculate the expected proportion of bird-days in each category, according to the model
  d.p <- do.call(rbind, lapply(1:length(loglik.birdday[[i]]), function(j){
    t <- loglik.birdday[[i]][[j]]$lik
    return(t(rmultinom(1, 1, colMeans(t))))
  }))
  prop.pred <- do.call(rbind, lapply(1:nlevels(as.factor(bins$bin)), function(j){
    b <- d.p[which(bins$bin == j), ]
    b.p <- colMeans(b)
    return(b.p)
  }))
  return(prop.pred)
})

dnum <- data.frame(dnum = as.numeric(as.factor(daily_sc$JDate)))

bins <- dnum %>% mutate(cuts = cut(dnum, c(seq(0, 60, 5), Inf))) %>% 
  group_by(cuts)

bins$bin <- as.numeric(as.factor(bins$cuts))

for(i in 1:nlevels(as.factor(bins$bin))){
  b <- d.p[which(bins$bin == i), ]
  b.p <- colMeans(b)
}

b.dat <- dta.bs[which(bins$bin == j), ]

b <- matrix(0, nrow(sc.dat), 3)
for(i in 1:nrow(b)){
  b[i, as.numeric(as.factor(sc.dat$data_bs))[i]] <- 1
}

b[as.numeric(as.factor(sc.dat$data_bs))[t]] <- 1




