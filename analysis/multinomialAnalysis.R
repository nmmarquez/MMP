rm(list=ls())

library(tidyverse)
library(rstan)
library(shinystan)
library(TMB)
library(tmbstan)

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
    # Intercept
    rep(1, nrow(subPersFullDF)),
    ## IRCA
    as.numeric(subPersFullDF$year >= 1987),
    # Dummy for ages
    sapply(ageFuncs, function(f) as.integer(f(subPersFullDF$age))),
    # Year variable starting at zero
    subPersFullDF$year - min(subPersFullDF$year),
    # year Trend squared
    (subPersFullDF$year - min(subPersFullDF$year))^2,
    ## Personal Effects
    persFullDF$edyrs,
    persFullDF$edyrs^2,
    ## Comunity effects
    log(persFullDF$COMPOP),
    persFullDF$LFPM,
    persFullDF$MINX2,
    persFullDF$pHeadUS,
    ## Policy in US
    # unemployment
    persFullDF$UR,
    # deportations
    log(persFullDF$deportations),
    # cohort numbers
    persFullDF$cohort,
    persFullDF$cohort^2
    )

Y <- cbind(
    as.numeric(!subPersFullDF$migration),
    as.numeric(!subPersFullDF$lpr),
    as.numeric(subPersFullDF$lpr)
)

Y[is.na(Y)] <- 0

Yvec <- apply(Y, 1, function(x) which(x == 1))

cidx <- as.integer(as.factor(subPersFullDF$commun))

migData <- list(
    Y = Y,
    N = nrow(Y),
    X = X,
    K = ncol(X),
    cidx = cidx,
    C = max(cidx)
)


estimateMNLR <- function(Y, X, cidx, silent=T, mcmc=F, re=FALSE, ...){
    model_name <- "./analysis/mnlr"
    hideme <- capture.output(compile(paste0(model_name, ".cpp")))
    
    # load the model
    dyn.load(dynlib(model_name))
    
    random <- "zetas"
    Map <- list()
    if(!re){
        random <- NULL
        Map[["zetas"]] <- factor(rep(NA, max(cidx)*2))
        dim(Map[["zetas"]]) <- c(max(cidx),2)
        Map[["sigmas"]] <- factor(rep(NA, 2))
        Map[["logrho"]] <- factor(rep(NA, 1))
    }
    
    # set the starting parameter points and data
    Params <- list(
        betas=matrix(0, nrow=ncol(X), ncol=2), 
        sigmas=c(0,0),
        logrho=0,
        zetas=matrix(0, nrow=max(cidx), ncol=2)
        )
    Data <- list(
        group=Y-1, 
        covs=X,
        cidx=cidx-1,
        C=max(cidx)
        )

    # build and optimize the objective function
    Obj <- MakeADFun(
        data=Data, parameters=Params, DLL="mnlr", silent=silent,
        random=random, map=Map)
    if(mcmc){
        return(tmbstan(Obj, ...))
    }
    start.time <- Sys.time()
    Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr, ...)
    runtime <- Sys.time() - start.time
    sdrep <- TMB::sdreport(Obj, getJointPrecision=TRUE)

    list(
        obj = Obj,
        opt = Opt,
        runtime = runtime,
        sd = sdrep)
}

predictMNLR <- function(X, model, draws=NULL){
    if(is.null(model$par)){
        allPars <- model$opt$par
        model$par <- allPars[names(allPars) == "betas"]
    }
    G <- 3
    if(is.null(X)){
        # lets not redo this for every draw
        X <- model.matrix(as.formula(model$call$formula[-2]), data)
    }
    if(!is.null(draws)){
        # recursive process for simulations
        betas_ <- mvtnorm::rmvnorm(draws, model$par, model$sd$cov.fixed)
        y_prob <- array(0, dim=c(nrow(X), G, draws))
        for(d in 1:draws){
            m <- model
            m$par <- betas_[d,] # treat the sims like true values
            y_prob[,,d] <- predictMNLR(X, m) # rerun function with sim
        }
    }
    else{
        # transform vector of beta paramters into matrix for ease of use
        betas_ <- matrix(model$par, ncol=G-1, nrow=ncol(X))
        # get the ratios for the groups vs ref
        ratio_ref <- exp(X %*% betas_)
        # the response var is transformed to probability space
        p1 <- apply(ratio_ref, 1, function(x) 1 / (1 + sum(x)))
        y_prob <- unname(cbind(p1, ratio_ref * p1))
    }
    return(y_prob)
}


testAge <- estimateMNLR(Yvec, X, cidx, F, 
                        control=list(eval.max=1e9, iter.max=1e9))
testAge$opt$convergence
testAge$opt$message
testPred <- predictMNLR(X, testAge, 100)

