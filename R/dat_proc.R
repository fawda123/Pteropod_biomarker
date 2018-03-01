library(readxl)
library(tidyverse)

envdat <- read_excel('raw/MeanSurface100m2016WCOA.xlsx') %>% 
  select(`CTD Station`, `'pCO2'`, `'pH'`, `'CO3'`, `'Fluorescence'`, `'Aragonite'`, `'Oxygen'`, `'Temperature'`, `'Salinity'`, `'Alkalinity'`) %>% 
  rename(
    CTD = `CTD Station`, 
    Ara = `'Aragonite'`, 
    O2 = `'Oxygen'`,
    pH = `'pH'`,
    CO3 = `'CO3'`,
    Fluor = `'Fluorescence'`,
    pCO2 = `'pCO2'`,
    Temp = `'Temperature'`, 
    Sal = `'Salinity'`, 
    TA = `'Alkalinity'`
  )
biodat <- read_excel('raw/PteropodIntegrated_oxibiomarkers.xlsx', 
                     sheet = 'Pteropod30mbis') %>% 
  select(CTD, GSHonGSSG, GST, GR, LPX, CAT, ORAC, SOD, Latitude) %>% 
  mutate(
    LPX = ifelse(LPX > 4, NA, LPX),
    LPX = ifelse(is.na(LPX), mean(LPX, na.rm = T), LPX)
    ) %>% 
  mutate(
    ORACvLPX = ORAC / LPX
  )

save(envdat, file = 'data/envdat.RData', compress = 'xz')
save(biodat, file = 'data/biodat.RData', compress = 'xz')
