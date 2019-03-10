rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(rstan)
library(TMB)
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

cfPlots <- list()

cfPlots$inconsis <- estDF %>%
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
    labs(x="Age", y="Survival Porbability till First Migration") +
    ggtitle("Counter Factual IRCA Effects for 2 Years (Consistent Controls)")

cfPlots$year <- estYearDF %>%
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
    labs(x="Age", y="Survival Porbability till First Migration") +
    ggtitle("Counter Factual IRCA Effects for 2 Years (Time Appropriate Controls)")

cfPlots$yearDiff <- estYearDF %>%
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
    labs(y="Increased Probability", x="Age") +
    ggtitle("Increased Survival For IRCA") +
    geom_hline(yintercept=0, linetype=2)

eduYearDF <- rbind(
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=F, edyrs=6, edyrsSq=6) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=F, edyrs=16, edyrsSq=16^2) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=T, edyrs=6, edyrsSq=6) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=T, edyrs=16, edyrsSq=16^2) %>%
        cbind(ageC=agegroups, age=15:35)) %>%
    select(-yearScale, -lnDeport, -UR, -yearScaleSq) %>%
    left_join(select(yearDF, -IRCA), by=c("year")) %>%
    cbind(predictNew(fullList$interHier, ., 1000, pars=fullPars))

cfPlots$edu <- eduYearDF %>%
    gather("Draw", "p", `1`:`1000`) %>%
    select(edyrs, IRCA, Draw, p, age) %>%
    arrange(edyrs, IRCA, Draw, age) %>%
    group_by(edyrs, IRCA, Draw) %>%
    mutate(Survival=cumprod(1-p)) %>%
    group_by(edyrs, IRCA, age) %>%
    summarise(
        mu=median(Survival), 
        lo=quantile(Survival, probs=.025),
        hi=quantile(Survival, probs=.975)) %>%
    ungroup %>%
    mutate(Years=as.character(edyrs)) %>%
    ggplot(aes(x=age, y=mu, ymin=lo, ymax=hi, group=Years)) +
    geom_line(aes(color=Years)) +
    geom_ribbon(aes(fill=Years), alpha=.3) +
    facet_wrap(~IRCA)+ 
    theme_classic() +
    labs(x="Age", y="Survival Porbability till First Migration") +
    ggtitle("Counter Factual IRCA Effects by Education Level")

cfPlots$eduDiff <- eduYearDF %>%
    gather("Draw", "p", `1`:`1000`) %>%
    select(edyrs, IRCA, Draw, p, age) %>%
    arrange(edyrs, IRCA, Draw, age) %>%
    group_by(edyrs, IRCA, Draw) %>%
    mutate(Survival=cumprod(1-p)) %>%
    group_by(edyrs, age, Draw) %>%
    summarise(diffSurv=diff(Survival)) %>%
    group_by(edyrs, age) %>%
    summarise(
        mu=median(diffSurv), 
        lo=quantile(diffSurv, probs=.025),
        hi=quantile(diffSurv, probs=.975)) %>%
    ungroup %>%
    mutate(Year=as.character(edyrs)) %>%
    ggplot(aes(x=age, y=mu, ymin=lo, ymax=hi, group=Year)) +
    geom_line(aes(color=Year)) +
    geom_ribbon(aes(fill=Year), alpha=.3) +
    theme_classic() +
    labs(y="Increased Probability", x="Age") +
    ggtitle("Increased Survival For IRCA by Education Level") +
    geom_hline(yintercept=0, linetype=2)

survDF <- read.csv("./data/pers161.csv") %>%
    group_by(commun) %>%
    summarize(surveyyr=min(surveyyr))
    

cfPlots$re <- fullList$interHier$summary.random$commun %>%
    rename(mu=`0.5quant`, lo=`0.025quant`, hi=`0.975quant`) %>%
    arrange(mu) %>%
    mutate(ID2=1:n(), commun=ID) %>%
    left_join(survDF, by="commun") %>%
    ggplot(aes(x=ID2, y=mu, ymin=lo, ymax=hi, color=surveyyr)) +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    theme_classic() +
    geom_hline(yintercept=0, linetype=2) + 
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
    labs(x="", y="Random Effect Estimate", color="Survey\nYear") +
    scale_color_distiller(palette = "Spectral") +
    ggtitle("Random Effect Estimates")
    
communTest <- fullList$interHier$summary.random$commun %>%
    arrange(mean) %>%
    filter(mean == max(mean) | ID == 135) %>%
    pull(ID)


communDF <- persFullDF %>%
    filter(commun %in% communTest) %>%
    select(commun, edyrs, lnPop, LFPM, MINX2, pHeadUS) %>%
    group_by(commun) %>%
    summarize_all(mean) %>%
    mutate(year=1997, cohort=10, cohortSq=100, edyrsSq=edyrs^2) %>%
    left_join(yearDF) %>%
    rbind(mutate(., IRCA=0)) %>%
    mutate(key=1) %>%
    left_join(tibble(ageC=agegroups, age=15:35, key=1)) %>%
    cbind(predictNew(fullList$interHier, ., 1000, pars=fullPars))
    
cfPlots$commun <- communDF %>%
    gather("Draw", "p", `1`:`1000`) %>%
    select(commun, IRCA, Draw, p, age) %>%
    arrange(commun, IRCA, Draw, age) %>%
    group_by(commun, IRCA, Draw) %>%
    mutate(Survival=cumprod(1-p)) %>%
    group_by(commun, IRCA, age) %>%
    summarise(
        mu=median(Survival), 
        lo=quantile(Survival, probs=.025),
        hi=quantile(Survival, probs=.975)) %>%
    ungroup %>%
    mutate(Community=as.character(commun)) %>%
    ggplot(aes(x=age, y=mu, ymin=lo, ymax=hi, group=Community)) +
    geom_line(aes(color=Community)) +
    geom_ribbon(aes(fill=Community), alpha=.3) +
    facet_wrap(~IRCA)+ 
    theme_classic() +
    labs(x="Age", y="Survival Porbability till First Migration") +
    ggtitle("Counter Factual IRCA Effects by Community")

cfPlots$communDiff <- communDF %>%
    gather("Draw", "p", `1`:`1000`) %>%
    select(commun, IRCA, Draw, p, age) %>%
    arrange(commun, IRCA, Draw, age) %>%
    group_by(commun, IRCA, Draw) %>%
    mutate(Survival=cumprod(1-p)) %>%
    group_by(commun, age, Draw) %>%
    summarise(diffSurv=diff(Survival)) %>%
    group_by(commun, age) %>%
    summarise(
        mu=median(diffSurv), 
        lo=quantile(diffSurv, probs=.025),
        hi=quantile(diffSurv, probs=.975)) %>%
    ungroup %>%
    mutate(Community=as.character(commun)) %>%
    ggplot(aes(x=age, y=mu, ymin=lo, ymax=hi, group=Community)) +
    geom_line(aes(color=Community)) +
    geom_ribbon(aes(fill=Community), alpha=.3) +
    theme_classic() +
    labs(y="Increased Probability", x="Age") +
    ggtitle("Increased Survival For IRCA by Community") +
    geom_hline(yintercept=0, linetype=2)

saveRDS(cfPlots, file="./plots/cfPlots.Rds")