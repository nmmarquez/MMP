rm(list=ls())

library(tidyverse)
library(rstan)
library(shinystan)
library(TMB)
library(tmbstan)
library(brms)
library(INLA)

load("./cleaned_data/inlaRuns.Rdata")

samplePars <- function(model, n){
    drawList <- inla.posterior.sample(n, model)
    npars <- nrow(drawList[[1]]$latent)
    betaDraws <- sapply(drawList, function(x){
        x$latent[(npars+1-length(model$names.fixed)):npars,]
    })
    commIDX <- startsWith(row.names(drawList[[1]]$latent), "commun")
    communDraws <- sapply(drawList, function(x){
        x$latent[commIDX,]
    })
    if(class(communDraws) == "list"){
        communDraws <- NULL
    }
    list(beta=betaDraws, commun=communDraws)
}

if(!file.exists("./cleaned_data/fullPars.Rds")){
    fullPars <- samplePars(fullList$interHier, 1000)
    saveRDS(fullPars, file="./cleaned_data/fullPars.Rds")
}

fullPars <- read_rds("./cleaned_data/fullPars.Rds")

predictNew <- function(model, DF, n=1000, prob=T, pars=NULL){
    if(is.null(pars)){
        pars <- samplePars(model, n)
    }
    ff <- update(model$.args$formula[-2], ~ . - f(commun, model="iid"))
    X <- model.matrix(ff, DF)
    linPred <- X %*% pars$beta
    if("commun" %in% names(DF) & !is.null(pars$commun)){
        commIDX <- paste0("commun:", DF$commun)
        linPred <- linPred + pars$commun[commIDX,]
    }
    if(prob){
        linPred <- arm::invlogit(linPred)
    }
    linPred
}

agegroups <- c(unique(persFullDF$ageC), rep(unique(persFullDF$ageC)[10], 11))

yearDF <- persFullDF %>%
    select(yearScale, lnDeport, IRCA, UR, year) %>%
    unique %>%
    group_by(yearScale) %>%
    summarize_all(mean) %>%
    mutate(yearScaleSq=yearScale^2)

meanDF <- persFullDF %>% 
    distinct(id, .keep_all=T) %>%
    select(-ageC, -statebrn) %>%
    summarize_all(mean) %>%
    select(-commun) %>%
    mutate(edyrsSq=edyrs^2, cohortSq=cohort^2)

estDF <- rbind(
    meanDF %>%
        mutate(year=1978, yearScale=1, IRCA=F) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=F) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1978, yearScale=1, IRCA=T) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=T) %>%
        cbind(ageC=agegroups, age=15:35)) %>%
    cbind(predictNew(fullList$interHier, ., 1000, pars=fullPars))

estYearDF <- rbind(
    meanDF %>%
        mutate(year=1978, yearScale=1, IRCA=F) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=F) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1978, yearScale=1, IRCA=T) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=T) %>%
        cbind(ageC=agegroups, age=15:35)) %>%
    select(-yearScale, -lnDeport, -UR, -yearScaleSq) %>%
    left_join(select(yearDF, -IRCA), by=c("year")) %>%
    cbind(predictNew(fullList$interHier, ., 1000, pars=fullPars))

estDF %>%
    gather("Draw", "p", `1`:`1000`) %>%
    select(year, IRCA, Draw, p, age) %>%
    arrange(year, IRCA, Draw, age) %>%
    group_by(year, IRCA, Draw) %>%
    mutate(Survival=cumprod(1-p)) %>%
    group_by(year, IRCA, age) %>%
    summarise(
        mu=median(Survival), 
        lo=quantile(Survival, probs=.025),
        hi=quantile(Survival, probs=.975)) %>%
    ungroup %>%
    mutate(Year=as.character(year)) %>%
    ggplot(aes(x=age, y=mu, ymin=lo, ymax=hi, group=Year)) +
    geom_line(aes(color=Year)) +
    geom_ribbon(aes(fill=Year), alpha=.3) +
    facet_wrap(~IRCA)+ 
    theme_classic() +
    labs(x="Survival Till First Migration", y="Porbability")

estYearDF %>%
    gather("Draw", "p", `1`:`1000`) %>%
    select(year, IRCA, Draw, p, age) %>%
    arrange(year, IRCA, Draw, age) %>%
    group_by(year, IRCA, Draw) %>%
    mutate(Survival=cumprod(1-p)) %>%
    group_by(year, IRCA, age) %>%
    summarise(
        mu=median(Survival), 
        lo=quantile(Survival, probs=.025),
        hi=quantile(Survival, probs=.975)) %>%
    ungroup %>%
    mutate(Year=as.character(year)) %>%
    ggplot(aes(x=age, y=mu, ymin=lo, ymax=hi, group=Year)) +
    geom_line(aes(color=Year)) +
    geom_ribbon(aes(fill=Year), alpha=.3) +
    facet_wrap(~IRCA)+ 
    theme_classic() +
    labs(x="Survival Till First Migration", y="Porbability")

estYearDF %>%
    gather("Draw", "p", `1`:`1000`) %>%
    select(year, IRCA, Draw, p, age) %>%
    arrange(year, IRCA, Draw, age) %>%
    group_by(year, IRCA, Draw) %>%
    mutate(Survival=cumprod(1-p)) %>%
    group_by(year, age, Draw) %>%
    summarise(diffSurv=diff(Survival)) %>%
    group_by(year, age) %>%
    summarise(
        mu=median(diffSurv), 
        lo=quantile(diffSurv, probs=.025),
        hi=quantile(diffSurv, probs=.975)) %>%
    ungroup %>%
    mutate(Year=as.character(year)) %>%
    ggplot(aes(x=age, y=mu, ymin=lo, ymax=hi, group=Year)) +
    geom_line(aes(color=Year)) +
    geom_ribbon(aes(fill=Year), alpha=.3) +
    theme_classic() +
    labs(x="Difference In Survival For IRCA", y="Porbability")

