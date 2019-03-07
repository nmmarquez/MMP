rm(list=ls())

library(tidyverse)
library(rstan)
library(shinystan)

# load in the cleaned dataset from the MMP
persFullDF <- read_rds("./cleaned_data/persFullDF.rds") %>%
    as_tibble()

# subset for testing purposes
sampleN <- Inf

ageFuncs <- list(
    function(x) x == 16,
    function(x) x == 17,
    function(x) x == 18,
    function(x) x == 19,
    function(x) x == 20,
    function(x) x == 21,
    function(x) x == 22,
    function(x) x == 23,
    function(x) x >= 24
)

# Subset data
subPersFullDF <- filter(persFullDF, id <= sampleN)

# Create Covariates
X <- cbind(
    rep(1, nrow(subPersFullDF)),
    sapply(ageFuncs, function(f) as.integer(f(subPersFullDF$age))))

Y <- cbind(
    as.numeric(!subPersFullDF$migration),
    as.numeric(!subPersFullDF$lpr),
    as.numeric(subPersFullDF$lpr)
)

Y[is.na(Y)] <- 0

cidx <- as.integer(as.factor(subPersFullDF$commun))

migData <- list(
    Y = Y,
    N = nrow(Y),
    X = X,
    K = ncol(X),
    cidx = cidx,
    C = max(cidx)
)

fit <- stan(
    file = './analysis/multinomial.stan', 
    data = migData,
    chains=4,
    cores=4,
    iter=20000,
    pars=c("beta", "zeta", "sigma", "rho"))

launch_shinystan(fit)
