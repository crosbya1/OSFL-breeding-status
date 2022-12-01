
library(jagsUI)

# Get outputs from cross-validation model runs 



# List of files in the cross-validation objects (in each list element)
#   - out: the outputs from the model run
#   - df: the data frame of Predicted and Props for creating the confusion matrix 
#   - loglik.birdday: list of the likelihoods for each test bird-day at each iteration (saved draws from the posterior)
#     - lik: the likelihood of survey-specific predictions, given the data and model predictions
#     - b: the observed breeding status
#     - lik.day: sum of log likelihoods of data|model for each survey in a day
#   - modlik.birdday: the full model likelihood as the sum of bird-day likelihoods (lik.day) averaged over all iterations 
#   - dta.bs: the observed breeding status for each test bird-day 
#   - bs.pred: the expected proportion of bird-days in each category, according to the model
#   - loglik.prop: of the observed proportion of bird-days in each breeding status, given the expected proportion of bird-days in each breeding status from the model


# For likelihood-based model comparison, take modlike.birdday and loglik.prop averaged across each fold 
mod.comp <- function(model){
  ln <- length(model)
  ml <- do.call(rbind, lapply(1:ln, function(x){
    mod.lik <- model[[x]]$modlik.birdday
    prop.lik <- model[[x]]$loglik.prop
    return(c(mod.lik, prop.lik))
  }))
  r <- sapply(1:ncol(ml), function(x) mean(ml[, x][is.finite(ml[, x])]))
  names(r) <- c("mod.lik", "prop.lik")
  return(r)
  gc()
}

d <- "2_pipeline/store"
ls <- list.files(d, pattern = ".Rdata")

mod.comp.table <- do.call(rbind, lapply(1:length(ls), function(x){
  load(file.path(d, ls[x]))
  t <- get(strsplit(ls[x], "\\.")[[1]][1])
  m <- mod.comp(t)
  return(m)
  gc()
}))
rownames(mod.comp.table) <- sapply(1:length(ls), function(x){
  n <- strsplit(strsplit(ls[x], "\\.")[[1]][1], "\\_")[[1]][1:2]
  n1 <- paste0(n[1], "_", n[2])
})

print(mod.comp.table, dig = 3)

saveRDS(mod.comp.table, file = "3_outputs/tables/mod_comp_table.rds")





