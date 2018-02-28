library(readxl)
library(tidyverse)

envdat <- read_excel('C:/Users/Marcus.SCCWRP2K/Desktop/MeanSurface100m2016WCOA.xlsx') %>% 
  select(`CTD Station`, `'pCO2'`, `'Aragonite'`, `'Oxygen'`, `'Temperature'`, `'Salinity'`, `'Alkalinity'`) %>% 
  rename(
    CTD = `CTD Station`, 
    Ara = `'Aragonite'`, 
    O2 = `'Oxygen'`,
    pCO2 = `pCO2`,
    Temp = `'Temperature'`, 
    Sal = `'Salinity'`, 
    TA = `'Alkalinity'`
  )
biodat <- read_excel('C:/Users/Marcus.SCCWRP2K/Desktop/PteropodIntegrated_oxibiomarkers.xlsx', 
                     sheet = 'Pteropod30mbis') %>% 
  select(CTD, GSHonGSSG, GST, GR, LPX, CAT, ORAC, SOC)