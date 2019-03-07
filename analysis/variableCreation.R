# variables we need to create
# - education at time of survey (time invariant) x
# - migrant cohort (birth year - 1950) (time invariant) x
# - community share of household heads who were in the US (time invariant) x
# - Dummies for state of birth x
# - Community Population (time varying)  x
# - Male Labor Force Participation (time varying) x 
# - % municipality labor force earning more than 2x minimum wage(time varying) x 
# - 1 period lagged annual number of mexican deported from US
# - 1 period lagged unemployment rate of hispanics in us (over age 16 males)
# - time of first migration X
# - type of first migration X

## Additional household characteristics
# - family of orgin migration experience
# - size of household
# - lfp of household

rm(list=ls())
library(XML)
library(RCurl)
library(rlist)
library(readxl)
library(parallel)
library(tidyverse)

# grab the deporation data from the webs
theurl <- "https://www.dhs.gov/immigration-statistics/yearbook/2016/table39" %>%
    getURL(httpheader = c("User-Agent"="Mozilla/5.0 (Windows NT 6.1; WOW64)"))
deportDF <- readHTMLTable(theurl, header=TRUE) %>%
    .[[1]] %>%
    mutate_all(as.character) %>%
    # we add 1 year beacuse this is for lagged years
    mutate(year = as.numeric(str_sub(Year, 1, 4))+1) %>%
    mutate(deportations=as.numeric(str_remove(Removals1, ","))) %>%
    select(year, deportations) %>%
    as_tibble

# pool labor data for hispanic males here (unadjusted males over 16)
# https://data.bls.gov/data/
unemployDF <- "./data/SeriesReport-20190305222617_fdd8ae.xlsx" %>%
    # remove a bunch of header data
    read_xlsx(skip=12) %>%
    # Create the average unemployment across the year
    gather("Month", "UR", -Year) %>%
    group_by(Year) %>%
    summarize(UR=mean(UR)) %>%
    # lag the years for our data set
    mutate(year=Year+1) %>%
    select(year, UR)

vars <- c("COMMUN", "COMPOP", "MINX2", "LFPM")

# Calculate head of households currently migrated
migDF <- read.csv("./data/mig161.csv") %>%
    as_tibble() %>%
    select(commun, uscurtrp) %>%
    # all of the data values here are either 1 or 2
    mutate(curMigrate=uscurtrp==1) %>%
    group_by(commun) %>%
    summarize(pHeadUS=mean(curMigrate))

# get community level variables
commDF <- read_csv("./data/commun161.csv") %>%
    select(matches( paste0("^(", paste(vars, collapse="|"), ")"))) %>%
    gather("Measure", "Value", -COMMUN) %>%
    mutate(Year=as.numeric(str_sub(Measure, -2, -1))) %>%
    mutate(Year=ifelse(Year < 40, Year + 2000, Year + 1900)) %>%
    mutate(Measure = str_sub(Measure, 1, -3)) %>%
    spread("Measure", "Value") %>%
    rename(year=Year, commun=COMMUN) %>%
    filter(year > 1970) %>%
    mutate(MINX2=ifelse(MINX2==9999 | MINX2==8888, NA, MINX2)) %>%
    mutate(LFPM=ifelse(LFPM==8888, NA, LFPM)) %>%
    # it looks like in the original analyis they removed any commun that
    # was missing data for any metric after 1980
    group_by(commun) %>%
    mutate(fullData=all(!is.na(LFPM) & !is.na(MINX2))) %>%
    ungroup %>%
    filter(fullData) %>%
    select(-fullData) %>%
    rename(yeard=year) %>%
    # hack for cleaner merge later
    mutate(yeard=yeard-.1)

# We need to expand out the datset to include all years and the authors do this
# by using the closest year that has data we should maybe also try linear
# interpolation later as well
commFullDF <- expand.grid(commun=unique(commDF$commun), year=1977:2015) %>%
    as_tibble %>%
    left_join(commDF) %>%
    mutate(diffabs=abs(yeard-year)) %>%
    group_by(commun, year) %>%
    mutate(minDiff=min(diffabs) == diffabs) %>%
    filter(minDiff) %>%
    ungroup %>%
    select(-diffabs, -minDiff, -yeard)

# lets get our state ids in order
stateID <- c(
    "AGU", "BCN", "BCS", "CAM", "COA", "COL", "CHP", "CHH", "CMX", "DUR", "GUA", 
    "GRU", "HID", "JAL", "MEX", "MIC", "MOR", "NAY", "NLE", "OAX", "PUE", "QUE", 
    "ROO", "SLP", "SIN", "SON", "TAB", "TAM", "TLA", "VER", "YUC", "ZAC")

persDF <- read.csv("./data/pers161.csv") %>%
    # Only males for this analysis
    filter(sex == 1) %>%
    # If an individual is deceased education data is not collected
    filter(edyrs != 8888) %>%
    # remove not reported education years
    filter(edyrs != 9999) %>%
    # Remove Individuals that migration data is not reported for
    filter(usyr1 != 9999) %>%
    # Remove Individuals that migration documentation is not reported for
    filter(usdoc1 != 9999) %>%
    # Remove Individuals whom state data is not reported for
    filter(statebrn != 9999) %>%
    # Only include individuals born in Mexico
    filter(statebrn <= 32) %>%
    # Change states to characters
    mutate(statebrn = stateID[statebrn]) %>%
    # change 8888 to never migrated NA (right censored)
    mutate(usyr1=ifelse(usyr1 == 8888, NA, usyr1)) %>%
    # only observe individuals who turned 14 between 1976-1985  
    filter(yrborn >= 1962 & yrborn <= 1971) %>%
    # recrete age based on surveyyr
    mutate(age=surveyyr - yrborn) %>%
    # calculate migration age if available
    mutate(migAge = usyr1-yrborn) %>%
    # remove individuals who first migrated before the age of 15
    filter(migAge >= 15 | is.na(migAge)) %>%
    mutate(lpr = case_when(
        usdoc1 == 8 ~ FALSE,
        usdoc1 == 1 ~ TRUE,
        TRUE ~ NA
    )) %>%
    # remove migrations that arent of concern
    filter(!(usdoc1 %in% c(2:7, 9:12))) %>%
    # Filter where we have an observed migration but uknown type
    filter(!(is.na(lpr) & !is.na(usyr1))) %>%
    # filter for anything below age 
    mutate(cohort=yrborn-1950) %>%
    as_tibble() %>%
    # remove places where we are missing community level data
    filter(commun %in% unique(commDF$commun)) %>%
    # IF I DO THE NEXT LINE I MATCH UP EXACTLY WITH THE OBS NUMBER
    # IN THE ORIGINAL PAPER (14,580)
    # filter(commun <= 154) %>%
    # i tried checking year dead here but like i mentioned before all deaths
    # have no education data colected so they are completely removed from
    # our datset unfortunately
    select(commun, hhnum, surveyyr, edyrs, cohort, yrborn,
           statebrn, usyr1, lpr, migAge) %>%
    mutate(id=1:n())

# We want to create a person year data set from the original persons data set
# that covers that covers the age of each individual who migrated from 15 to 35
# up until either death or migration

persFullDF <- bind_rows(mclapply(1:nrow(persDF), function(i){
    subDF <- filter(persDF, id == i) %>%
        select(-surveyyr)
    migAge <- ifelse(is.na(subDF$migAge), 36, subDF$migAge)
    tibble(id=i, age=15:migAge) %>%
        left_join(subDF, by="id") %>%
        mutate(year=yrborn+age) %>%
        select(-migAge, -yrborn, -usyr1) %>%
        mutate(lpr=ifelse(max(year) == year, lpr, NA)) %>%
        mutate(migration=!is.na(lpr)) %>%
        filter(age <= 35) %>%
        as.data.frame}, mc.cores = 5)) %>%
    left_join(commFullDF) %>%
    left_join(migDF) %>%
    left_join(unemployDF) %>%
    left_join(deportDF)

saveRDS(persFullDF, file = "./cleaned_data/persFullDF.rds")