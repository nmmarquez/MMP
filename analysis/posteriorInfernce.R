.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())

library(dplyr)
library(tidyr)
library(tibble)
library(TMB)
library(INLA)
library(glmmTMB)
library(ggplot2)
library(plotly)
library(stringr)

load("./cleaned_data/tmbRuns.Rdata")

simulateCommRandomEffects <- function(model, n=1000){
    nCommun <- nrow(ranef(model)$cond$commun)
    nRE <- ncol(ranef(model)$cond$commun)
    effects <- c(as.matrix(ranef(model)$cond$commun))
    Sigma <- matrix(0, nrow=nCommun*nRE, ncol=nCommun*nRE)
    diag(Sigma) <- model$sdr$diag.cov.random[1:(nCommun*nRE)]
    if(nRE>1){
        rho <- attr(summary(model)$varcor$cond$commun, "correlation")[1,2]
        for(i in 1:nCommun){
            Sigma[i, i+nCommun] <- sqrt(Sigma[i, i]) * 
                sqrt(Sigma[i+nCommun, i+nCommun]) * rho
            Sigma[i+nCommun, i] <- Sigma[i,i+nCommun] 
        }
    }
    draws <- mvtnorm::rmvnorm(n, effects, Sigma)
    
    reDraws <- lapply(1:nRE, function(i){
        z <- t(draws[,(nCommun*(i-1)+1):(nCommun*i)])
        row.names(z) <- paste0("commun:", unique(persFullDF$commun))
        z
    })
    names(reDraws) <- names(ranef(model)$cond$commun)
    reDraws
}

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

# if(!file.exists("./cleaned_data/fullPars.Rds")){
#     fullPars <- samplePars(reIntList$perscommun, 1000)
#     saveRDS(fullPars, file="./cleaned_data/fullPars.Rds")
# }

# fullPars <- read_rds("./cleaned_data/fullPars.Rds")

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

predictNew2 <- function(model, DF, n=1000, prob=T){
    ff <- update(
        model$call$formula[-2], 
        ~ . -(1 | hhid) - (1 | commun) - (1 + IRCA | commun))
    X <- model.matrix(ff, DF)
    idx <- names(model$fit$par) == "beta"
    b_ <- model$fit$par[idx]
    S_ <- model$sdr$cov.fixed[idx,idx]
    linPred <- X %*% t(mvtnorm::rmvnorm(n, b_, S_))
    if("commun" %in% names(DF)){
        reDraws <- simulateCommRandomEffects(model, n)
        commIDX <- paste0("commun:", DF$commun)
        for(j in names(reDraws)){
            linPred <- linPred + reDraws[[j]][commIDX,] * 
                sapply(1:n, function(i) X[,j])
        }
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
    select(-ageC, -statebrn, -hhid, -commun2, -eduC) %>%
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
    cbind(predictNew2(reSlopeList$full, ., 1000))

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
    cbind(predictNew2(reSlopeList$full, ., 1000))

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
        mutate(year=1991, yearScale=14, IRCA=F, edyrs=3, edyrsSq=3^2) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=F, edyrs=16, edyrsSq=16^2) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=T, edyrs=3, edyrsSq=3^2) %>%
        cbind(ageC=agegroups, age=15:35),
    meanDF %>%
        mutate(year=1991, yearScale=14, IRCA=T, edyrs=16, edyrsSq=16^2) %>%
        cbind(ageC=agegroups, age=15:35)) %>%
    select(-yearScale, -lnDeport, -UR, -yearScaleSq) %>%
    left_join(select(yearDF, -IRCA), by=c("year")) %>%
    cbind(predictNew2(reSlopeList$full, ., 1000))

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

# eduYearDF %>%
#     gather("Draw", "p", `1`:`1000`) %>%
#     select(edyrs, IRCA, Draw, p, age) %>%
#     arrange(edyrs, IRCA, Draw, age) %>%
#     group_by(edyrs, IRCA, Draw) %>%
#     mutate(Survival=cumprod(1-p)) %>%
#     group_by(edyrs, age, Draw) %>%
#     summarise(diffSurv=diff(Survival)) %>%
#     arrange(edyrs, age, Draw, edyrs) %>%
#     group_by(age, Draw) %>%
#     summarise(diffSurv=diff(diffSurv)) %>%
#     group_by(age) %>%
#     summarise(
#         mu=median(diffSurv), 
#         lo=quantile(diffSurv, probs=.025),
#         hi=quantile(diffSurv, probs=.975)) %>%
#     ungroup %>%
#     ggplot(aes(x=age, y=mu, ymin=lo, ymax=hi)) +
#     geom_line() +
#     geom_ribbon(alpha=.3) +
#     theme_classic() +
#     labs(y="Increased Probability", x="Age") +
#     ggtitle("Increased Survival For IRCA by Education Level") +
#     geom_hline(yintercept=0, linetype=2)

survDF <- read.csv("./data/pers161.csv") %>%
    group_by(commun) %>%
    summarize(surveyyr=min(surveyyr))

reEffects <- simulateCommRandomEffects(reSlopeList$full, 1000)
b0 <- reSlopeList$full$fit$par[1]
bIRCA <- reSlopeList$full$fit$par[2]

reDF <- bind_rows(
    tibble(
        mu = apply(reEffects$`(Intercept)`, 1, median) + b0,
        lo = apply(reEffects$`(Intercept)`, 1, quantile, probs=.025) + b0,
        hi = apply(reEffects$`(Intercept)`, 1, quantile, probs=.975) + b0,
        commun = row.names(reEffects$IRCA),
        ov=b0,
        effect="Intercept"),
    tibble(
        mu = apply(reEffects$IRCA, 1, median) + bIRCA,
        lo = apply(reEffects$IRCA, 1, quantile, probs=.025) + bIRCA,
        hi = apply(reEffects$IRCA, 1, quantile, probs=.975) + bIRCA,
        commun = row.names(reEffects$IRCA),
        ov=bIRCA,
        effect="IRCA")) %>%
    mutate(commun=as.numeric(str_replace(commun, "commun:", "")))
    
cfPlots$re <- reDF %>%
    group_by(effect) %>%
    arrange(effect, mu) %>%
    mutate(ID2=1:n()) %>%
    ungroup() %>%
    left_join(survDF, by="commun") %>%
    ggplot(aes(x=ID2, y=mu, ymin=lo, ymax=hi, color=surveyyr)) +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    theme_classic() +
    geom_hline(aes(yintercept=ov), linetype=2) + 
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
    labs(x="", y="Random Effect Estimate", color="Survey\nYear") +
    scale_color_distiller(palette = "Spectral") +
    ggtitle("Random Effect Estimates") +
    facet_wrap(~effect)
    
cfPlots$reBivar <- reDF %>%
    select(-ov) %>%
    nest(mu, lo, hi, .key = "temp") %>%
    spread(effect, temp) %>%
    unnest(Intercept, IRCA, .sep = '_') %>%
    mutate(IRCAsig=IRCA_lo > 0 | IRCA_hi < 0) %>%
    ggplot(aes(
        x=Intercept_mu, xmin=Intercept_lo, xmax=Intercept_hi,
        y=IRCA_mu, ymin=IRCA_lo, ymax=IRCA_hi, color=IRCAsig, text=commun)) +
    geom_point() +
    geom_errorbar(alpha=.3, size=.2) +
    geom_errorbarh(alpha=.3, size=.2) +
    theme_classic() +
    geom_vline(xintercept=b0, linetype=2) +
    geom_hline(yintercept=bIRCA, linetype=2) +
    labs(x="Intercept", y="IRCA", title="Random Effects Estimates",
         color="Signficant\nIRCA Effect")

communTest <- c(22, 60)

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
    cbind(predictNew2(reSlopeList$full, ., 1000))

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
    mutate(IRCA=as.logical(IRCA)) %>%
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

cfPlots$covPlot <- summary(reSlopeList$full)$coefficients$cond %>%
    as.data.frame() %>%
    mutate(Covariate=row.names(.), mu=Estimate) %>%
    mutate(lo=mu-1.96*`Std. Error`, hi=mu+1.96*`Std. Error`) %>%
    select(-Estimate:-`Pr(>|z|)`) %>%
    mutate_if(is.numeric, exp) %>%
    arrange(mu) %>%
    mutate(Covariate=factor(Covariate, levels=Covariate)) %>%
    mutate(tx=paste0(round(mu,2), "(", round(lo,2), ",", round(hi,2), ")")) %>%
    ggplot(aes(x=Covariate, y=mu, ymin=lo, ymax=hi, label=tx)) +
    geom_point() +
    geom_errorbar() +
    theme_classic() +
    coord_flip() +
    geom_hline(yintercept=1, linetype=2) +
    geom_text(aes(y=ifelse(mu<5, hi+.9, lo-.9))) +
    labs(y="Estimate")

saveRDS(cfPlots, file="./plots/cfPlots.Rds")
ggsave(covPlot, filename="./plots/covariateEst.png")
ggsave(cfPlots$edu + theme(text = element_text(size=20)), 
       filename="./plots/edu.png")
ggsave(cfPlots$eduDiff + theme(text = element_text(size=20)), 
       filename="./plots/eduDiff.png")
ggsave(cfPlots$commun + theme(text = element_text(size=20)), 
       filename="./plots/commun.png")
ggsave(cfPlots$communDiff + theme(text = element_text(size=20)), 
       filename="./plots/communDiff.png")
ggsave(cfPlots$re + theme(text = element_text(size=20)), 
       filename="./plots/wrrrrry.png")

