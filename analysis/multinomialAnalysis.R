rm(list=ls())

library(tidyverse)
library(rstan)
library(shinystan)
library(TMB)
library(tmbstan)
library(brms)
library(INLA)

# load in the cleaned dataset from the MMP
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
    select(-age) %>%
    mutate(IRCA=as.numeric(year >= 1987)) %>%
    mutate(yearScale=year-min(year), yearScaleSq=yearScale^2) %>%
    mutate(edyrsSq=edyrs^2, lnPop=log(COMPOP)) %>%
    mutate(lnDeport=log(deportations), cohortSq=cohort^2)

# This model takes about 35 min for a much smaller portion of the covariates
# user   system  elapsed 
# 71.271   18.840 2157.156 
# system.time (b1 <- brm (
#     Y ~ 1 + IRCA + ageC,
#     data=persFullDF, family="bernoulli", chains=3, iter=3000, warmup=600,
#     cores=4, prior=c(set_prior ("normal (0, 8)"))))
# saveRDS(b1, file="cleaned_data/stanrun.Rds")

ffFullModels <-list(
    base = migration ~ 1 + IRCA + ageC,
    person = migration ~ 1 + IRCA + ageC + edyrs + edyrsSq,
    commun = migration ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop +
        LFPM + MINX2 + pHeadUS,
    full = migration ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM +
        MINX2 + pHeadUS + UR + lnDeport + cohort + cohortSq + yearScale + 
        yearScaleSq,
    full = migration ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM +
        MINX2 + pHeadUS + UR + lnDeport + cohort + cohortSq + yearScale + 
        yearScaleSq,
    inter = migration ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM +
        MINX2 + pHeadUS + UR + lnDeport + cohort + cohortSq + yearScale + 
        yearScaleSq + IRCA*edyrs + IRCA*edyrsSq,
    personHier = migration ~ 1 + IRCA + ageC + edyrs + edyrsSq +
        f(commun, model="iid"),
    hier = migration ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM +
        MINX2 + pHeadUS + UR + lnDeport + cohort + cohortSq + yearScale + 
        yearScaleSq + f(commun, model="iid"),
    interHier = migration ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM +
        MINX2 + pHeadUS + UR + lnDeport + cohort + cohortSq + yearScale + 
        yearScaleSq + IRCA*edyrs + IRCA*edyrsSq +
        f(commun, model="iid")
)

ffMigModels <- lapply(ffFullModels, function(x)
    update(x, lpr ~ .))

fullList <- lapply(ffFullModels, function(f){
    inla(f, family="binomial", data=persFullDF, control.compute=list(config=T))
})

# run the same models but conditional on having migrated what is type of 
# migration

migDF <- persFullDF %>%
    filter(!is.na(lpr)) %>%
    mutate(lpr=as.numeric(lpr))

migList <- lapply(ffMigModels, function(f){
    inla(f, family="binomial", data=migDF, control.compute=list(config=T))
})

migList$persInterHier <- inla(
    lpr ~ 1 + IRCA + ageC + edyrs + edyrsSq + IRCA*edyrs + IRCA*edyrsSq + 
        f(commun, model="iid"), 
    family="binomial", data=migDF, control.compute=list(config=T))

inlaBIC <- function(inlaModel){
    k <- length(inlaModel$names.fixed) + inlaModel$nhyper -1
    n <- length(inlaModel$misc$linkfunctions$link)
    lnL <- inlaModel$mlik[1,1]
    unname(log(n) * k  - 2 * lnL)
}

migListBIC <- sapply(migList, inlaBIC)
fullListBIC <- sapply(fullList, inlaBIC)

save(migList, fullList, persFullDF, migDF, migListBIC, fullListBIC,
     file = "./cleaned_data/inlaRuns.Rdata")
