.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())

library(dplyr)
library(tidyr)
library(tibble)
library(TMB)
library(INLA)
library(readr)
library(BMA)
library(glmmTMB)

# load in the cleaned dataset from the MMP
set.seed(123)
persFullDF <- read_rds("./cleaned_data/persFullDF.rds") %>%
    as_tibble() %>%
    mutate(migration=as.numeric(migration)) %>%
    mutate(ageC=case_when(
        age == 15 ~ "age 15",
        age == 16 ~ "age 16",
        age == 17 ~ "age 17",
        age == 18 ~ "age 18",
        age == 19 ~ "age 19",
        age == 20 ~ "age 20",
        age == 21 ~ "age 21",
        age == 22 ~ "age 22",
        age == 23 ~ "age 23",
        TRUE ~ "age 24+"
    )) %>%
    mutate(eduC=case_when(
        edyrs <= 2 ~ "0-2 year Edu",
        edyrs <= 5 ~ "3-5 years Edu",
        edyrs <= 8 ~ "6-8 years Edu",
        edyrs <= 11 ~ "9-11 years Edu",
        TRUE ~ "12+ years Edu"
    )) %>%
    select(-age) %>%
    mutate(IRCA=as.numeric(year >= 1987)) %>%
    mutate(yearScale=year-min(year), yearScaleSq=yearScale^2) %>%
    mutate(edyrsSq=edyrs^2, lnPop=log(COMPOP)) %>%
    mutate(lnDeport=log(deportations), cohortSq=cohort^2) %>%
    # randomly sample one individual for each household
    # group_by(commun, hhnum) %>%
    # filter(id == sample(unique(id), 1)) %>%
    ungroup() %>%
    mutate(commun2=commun, hhid=paste0(commun, "_", hhnum))

# This model takes about 35 min for a much smaller portion of the covariates
# user   system  elapsed 
# 71.271   18.840 2157.156 
# system.time (b1 <- brm (
#     Y ~ 1 + IRCA + ageC,
#     data=persFullDF, family="bernoulli", chains=3, iter=3000, warmup=600,
#     cores=4, prior=c(set_prior ("normal (0, 8)"))))
# saveRDS(b1, file="cleaned_data/stanrun.Rds")

ffFullModelsEduGroup <-list(
    base = migration ~ 1 + IRCA + ageC,
    person = migration ~ 1 + IRCA * eduC + ageC,
    commun = migration ~ 1 + IRCA + ageC + lnPop + LFPM + MINX2 + pHeadUS,
    perscommun = migration ~ 1 + IRCA * eduC + ageC + lnPop + LFPM +
        MINX2 + pHeadUS,
    full = migration ~ 1 + IRCA * eduC + ageC +lnPop + LFPM +
        MINX2 + pHeadUS + UR + lnDeport + cohort + yearScale
)

ffFullModels <-list(
    base = migration ~ 1 + IRCA + ageC,
    person = migration ~ 1 + IRCA * edyrs + edyrsSq + ageC,
    commun = migration ~ 1 + IRCA + ageC + lnPop + LFPM + MINX2 + pHeadUS,
    perscommun = migration ~ 1 + IRCA * edyrs + edyrsSq + ageC + lnPop + LFPM +
        MINX2 + pHeadUS,
    full = migration ~ 1 + IRCA * edyrs + edyrsSq + ageC +lnPop + LFPM +
        MINX2 + pHeadUS + UR + lnDeport + cohort + yearScale
)

ffMigModels <- lapply(ffFullModels, function(x)
    update(x, lpr ~ .))

# lets check out the BMA
bmaFit <- bic.glm(
    x = model.matrix(ffFullModels$full[-2], persFullDF)[,-1],
    y = persFullDF$migration,
    glm.family="binomial")

imageplot.bma(bmaFit)

# lets remove the year and cohort effect as they dont seem to be providing much
#ffFullModels$full <- update(ffFullModels$full, . ~ . - yearScale - cohort)
#ffFullModelsEduGroup$full <- update(
#    ffFullModelsEduGroup$full, . ~ . - yearScale - cohort)

glmList <- lapply(ffFullModels, function(f){
    glm(f, family="binomial", data=persFullDF)
})

glmListGroup <- lapply(ffFullModelsEduGroup, function(f){
    glm(f, family="binomial", data=persFullDF)
})

sapply(glmList, BIC)
sapply(glmListGroup, BIC)

ffRandModels <- lapply(ffFullModels[4:5], function(f_){
    update(f_, . ~ . + (1|commun) + (1|hhid))
})

ffRandSlope <- lapply(ffFullModels[4:5], function(f_){
    update(f_, . ~ . + (1+IRCA|commun) + (1|hhid))
})

reIntList <- lapply(ffRandModels, function(f_){
    glmmTMB(f_, family="binomial", data=persFullDF, verbose=T)
})

reSlopeList <- lapply(ffRandSlope, function(f_){
    glmmTMB(f_, family="binomial", data=persFullDF, verbose=T)
})

# run the same models but conditional on having migrated what is type of 
# migration

# migDF <- persFullDF %>%
#     filter(!is.na(lpr)) %>%
#     mutate(lpr=as.numeric(lpr))
# 
# migList <- lapply(ffMigModels, function(f){
#     inla(f, family="binomial", data=migDF, control.compute=list(config=T))
# })
# 
# migList$persInterHier <- inla(
#     lpr ~ 1 + IRCA + ageC + edyrs + edyrsSq + IRCA*edyrs + IRCA*edyrsSq + 
#         f(commun, model="iid"), 
#     family="binomial", data=migDF, control.compute=list(config=T))

inlaBIC <- function(inlaModel){
    k <- length(inlaModel$names.fixed) + inlaModel$nhyper -1
    n <- length(inlaModel$misc$linkfunctions$link)
    lnL <- inlaModel$mlik[1,1]
    unname(log(n) * k  - 2 * lnL)
}


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
        t(draws[,(nCommun*(i-1)+1):(nCommun*i)])
    })
    names(reDraws) <- names(ranef(model)$cond$commun)
    reDraws
}

modelBICList <- list(
    reInt = lapply(reIntList, BIC),
    reSlope = lapply(reSlopeList, BIC),
    fixed = sapply(glmList, BIC),
    fixedGroupEdu = sapply(glmListGroup, BIC))

save(persFullDF, reIntList, modelBICList, reSlopeList,
     file = "./cleaned_data/tmbRuns.Rdata")
