library(readxl)
library(tidyverse)

# environmental data
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
  ) %>% 
  mutate_if(is.character, as.numeric)

# biomarker data
biodat <- read_excel('raw/PteropodIntegrated_oxibiomarkers.xlsx', 
                     sheet = 'Pteropod30mbis') %>% 
  select(CTD, GSHonGSSG, GST, GR, LPX, CAT, ORAC, SOD, Latitude) %>% 
  mutate(
    ORACvLPX = ORAC / LPX
  ) %>% 
  gather('var', 'val', -CTD) %>% 
  group_by(CTD, var) %>% 
  summarize(val = mean(val, na.rm = T)) %>% 
  group_by(var) %>% 
  mutate(val = ifelse(is.na(val), mean(val, na.rm = T), val)) %>% 
  spread(var, val) %>% 
  ungroup

# abundance
abudat <- read_excel('raw/Copy of abundances to be usef_Marcus2.xlsx') %>% 
  rename(
    CTD = `CTD station`,
    abund = `real abundance`
  ) %>% 
  mutate(
    abu = log10(1 + abund)
  ) %>% 
  group_by(CTD) %>% 
  summarize(
    abu = mean(abu, na.rm = T)
  ) %>% 
  ungroup

# shell dissolution
disdat <- read_excel('raw/parameters for WCOA 2016 combined.xlsx', sheet = 'dissolution') %>% 
  fill(station) %>% 
  rename( 
    CTD = station, 
    dis = `dissolution extent`, 
    ty2 = `type II`,
    ty3 = `type III`,
    scr = scarring
    ) %>% 
  mutate(dis = asin(dis / 100)) %>% 
  group_by(CTD) %>% 
  summarise(
    dis = mean(dis, na.rm = T), 
    ty2 = mean(ty2, na.rm = T),
    ty3 = mean(ty3, na.rm = T),
    scr = mean(scr, na.rm = T)
  ) %>% 
  ungroup

# shell length
diadat <- read_excel('raw/parameters for WCOA 2016 combined.xlsx', sheet = 'diameter') %>% 
  rename(
    CTD = Individual, 
    len = `Diameter (mm)`
  ) %>% 
  group_by(CTD) %>% 
  summarise(
    len = mean(len, na.rm = T)
  ) %>% 
  ungroup

# full_join of biodat with abundance, dissolution, and length data
ptedat <- biodat %>% 
  full_join(abudat, by = 'CTD') %>% 
  full_join(disdat, by = 'CTD') %>% 
  full_join(diadat, by = 'CTD')

# save envdat, ptedat
save(envdat, file = 'data/envdat.RData', compress = 'xz')
save(ptedat, file = 'data/ptedat.RData', compress = 'xz')
