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

fullList <- list()

fullList$baseInla <- inla(
    migration ~ 1 + IRCA + ageC,
    family = "binomial",
    data = persFullDF,  control.compute = list(dic = TRUE, waic = TRUE))

fullList$persInla <- inla(
    migration ~ 1 + IRCA + ageC + edyrs + edyrsSq,
    family = "binomial",
    data = persFullDF,  control.compute = list(dic = TRUE, waic = TRUE))

fullList$pcomInla <- inla(
    migration ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM + MINX2 + 
        pHeadUS,
    family = "binomial",
    data = persFullDF,  control.compute = list(dic = TRUE, waic = TRUE)
)

fullList$fixedInla <- inla(
    migration ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM + MINX2 + 
        pHeadUS + UR + lnDeport + cohort + cohortSq + yearScale + yearScaleSq,
    family = "binomial",
    data = persFullDF,  control.compute = list(dic = TRUE, waic = TRUE)
)

fullList$hierINLA <- inla(
    migration ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM + MINX2 + 
        pHeadUS + UR + lnDeport + cohort + cohortSq + yearScale + yearScaleSq +
        f(commun, model="iid"),
    family = "binomial",
    data = persFullDF,  control.compute = list(dic = TRUE, waic = TRUE)
)

# run the same models but conditional on having migrated what is type of 
# migration

migDF <- persFullDF %>%
    filter(!is.na(lpr)) %>%
    mutate(lpr=as.numeric(lpr))

migList <- list()

migList$baseInla <- inla(
    lpr ~ 1 + IRCA + ageC,
    family = "binomial",
    data = migDF,  control.compute = list(dic = TRUE, waic = TRUE))

migList$persInla <- inla(
    lpr ~ 1 + IRCA + ageC + edyrs + edyrsSq,
    family = "binomial",
    data = migDF,  control.compute = list(dic = TRUE, waic = TRUE))

migList$pcomInla <- inla(
    lpr ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM + MINX2 + 
        pHeadUS,
    family = "binomial",
    data = migDF,  control.compute = list(dic = TRUE, waic = TRUE)
)

migList$fixedInla <- inla(
    lpr ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM + MINX2 + 
        pHeadUS + UR + lnDeport + cohort + cohortSq + yearScale + yearScaleSq,
    family = "binomial",
    data = migDF,  control.compute = list(dic = TRUE, waic = TRUE)
)

migList$hierINLA <- inla(
    lpr ~ 1 + IRCA + ageC + edyrs + edyrsSq + lnPop + LFPM + MINX2 + 
        pHeadUS + UR + lnDeport + cohort + cohortSq + yearScale + yearScaleSq +
        f(commun, model="iid"),
    family = "binomial",
    data = migDF,  control.compute = list(dic = TRUE, waic = TRUE)
)

save(migList, fullList, persFullDF, migDF, 
     file = "./cleaned_data/inlaRuns.Rdata")