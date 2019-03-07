rm(list=ls())

# Approx 98 percent of Mexican FB users are 15 and over 
# https://www.slideshare.net/wearesocialsg/digital-in-2017-central-america #66
fbprop15above <- .98

# Social media penetration
# https://www.slideshare.net/wearesocialsg/digital-in-2017-central-america #63
socialMediaPenetration <- .59

# total number of monthly FB users in 2017 76 million
# https://www.slideshare.net/wearesocialsg/digital-in-2017-central-america #64
totalFB <- 76000000

# FB Penetration
# https://www.emarketer.com/Article/Facebook-Dominates-Social-Media-Market-Mexico/1013828
fbPenetration <- .97

# estimated population percentage age 15 and above as of 2017
# https://www.cia.gov/library/publications/the-world-factbook/fields/2010.html
popprop15above <- 1 - .2693

# estimated population size Mexico 2017
# https://esa.un.org/unpd/wpp/DVD/Files/1_Indicators%20(Standard)/EXCEL_FILES/1_Population/WPP2017_POP_F01_1_TOTAL_POPULATION_BOTH_SEXES.xlsx
popTotal <- 129163000

# estimated number of users age 15 and above
fbnum15above <- fbprop15above * popTotal * socialMediaPenetration * fbPenetration

# estimated population age 15 and above
popnum15above <- popTotal * popprop15above

# estimate penetration of FB users 15 and above
round(fbnum15above / popnum15above * 100, 2)

round(fbPenetration * socialMediaPenetration * 100, 2)

round(totalFB/popTotal * 100, 2)
