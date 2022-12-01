

# The estimated negative binomial dispersion parameter (estimated from the model at iteration i)
r <- 0.225 

# The expected song count for each breeding status as a function of time since sunrise (estimated from the model at iteration i)
lam <- c(2.836685,  8.303503,  5.398654)

# The observed song count
y <- 12


# Probability of the observed song count data, given expected song count from the model, for each breeding status
pl <- dpois(y, lam)  # If the poisson model
pl <- pl/sum(pl)

pn <- dnbinom(y, size = r, mu = lam)  # If the negative binomial model
pn <- pn/sum(pn)

# The estimated probability of being in each breeding status as a function of Julian day (estimated from the model at iteration i)
pb <- c(0.04073711, 0.75050228, 0.20876061)


# The probability of being in each breeding status as a function of Julian day and song rate (Equation 1 in the manuscript)
pbl <- pl*pb/sum(pl*pb)  # For the poisson model
pbn <- pn*pb/sum(pn*pb)  # For the negative binomial model

# The observed breeding status in vector format
b <- c(0, 1, 0)

# The log likelihood of the estimated probabilities, given the data
lik.l <- dmultinom(b, prob = pbl, log = TRUE)  # For the negative binomial model
lik.n <- dmultinom(b, prob = pbn, log = TRUE)  # For the negative binomial model


