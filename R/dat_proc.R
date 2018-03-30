library(readxl)
library(tidyverse)

source('R/funcs.R')

# environmental data
envdat <- read_excel('raw/MeanSurface100m2016WCOA.xlsx') %>% 
  select(`CTD Station`, `Latitude`, `'pCO2'`, `'pH'`, `'CO3'`, `'Fluorescence'`, `'Aragonite'`, `'Oxygen'`, `'Temperature'`, `'Salinity'`, `'Alkalinity'`) %>% 
  rename(
    CTD = `CTD Station`,
    Lat = `Latitude`,
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
  mutate(
    dis = asin(dis / 100),
    ty2 = asin(ty2),
    ty3 = asin(ty3),
    scr = asin(scr)
    ) %>% 
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

# exposure data
expdat <- read_excel('raw/experimental treatment for Marcus.xlsx') %>% 
  rename(
    temp_trt = T, 
    pco2_trt = PCO2, 
    dead = `%Dead`, 
    alive = `%Alive`
  ) %>% 
  mutate(
    temp_trt = ifelse(temp_trt < 10, 'lo', 'hi'),
    pco2_trt = ifelse(pco2_trt < 800, 'lo', 'hi')
  ) %>% 
  select(CTD, temp_trt, pco2_trt, alive) %>% 
  group_by(CTD, temp_trt, pco2_trt) %>% 
  summarize(alive = mean(alive)) %>% 
  ungroup

save(expdat, file = 'data/expdat.RData', compress = 'xz')

######
# all stressor covar models

##
# cellular

# data
data(envdat)
data(ptedat)

envchr <- c('pCO2', 'pH', 'CO3', 'Ara', 'O2', 'Temp', 'Fluor')
ptechr <- c('CAT', 'GR', 'GSHonGSSG', 'GST', 'LPX', 'ORAC', 'SOD', 'ORACvLPX')

dat_frm <- ptedat %>% 
  select(one_of('CTD', ptechr)) %>% 
  na.omit %>% 
  inner_join(envdat, by = 'CTD') %>% 
  select(one_of(c('CTD', envchr, ptechr))) %>% 
  data.frame %>% 
  remove_rownames %>% 
  column_to_rownames('CTD')

# models
env_cmb <- envchr %>% 
  combn(2) %>% 
  t %>% 
  data.frame(stringsAsFactors = F) %>%
  crossing(ptechr, .) %>% 
  t %>% 
  data.frame(stringsAsFactors = F) %>% 
  as.list

# models
mod_all <- env_cmb %>%
  map(function(x, dat = dat_frm){
    
    # formula
    frm <- paste(x[1], '~', x[2], '*', x[3]) %>% 
      as.formula
    
    # model only
    lm(frm, data = dat)
    
  }) %>% 
  enframe('Model', 'Modobj') %>% 
  mutate(
    Model = gsub('^X', 'mod', Model), 
    data = map(Modobj, function(x){
      
      # response
      yvar <- all.vars(formula(x))[1]
      
      # rsq
      rsq <- summary(x) %>% 
        .$r.squared
      
      # coeff summary
      out <- x %>% 
        summary %>% 
        .$coefficients %>% 
        .[-1, c(1, 4)] %>% 
        data.frame %>% 
        rownames_to_column('env_lab') %>% 
        mutate(
          pte_lab = yvar,
          Rsq = round(rsq, 2),
          Est = round(Estimate, 2), 
          Pvl = p_ast(Pr...t..)
        ) %>% 
        select(pte_lab, env_lab, Rsq, Est, Pvl)
      
      return(out)
      
    })
  ) 

biomod <- mod_all
save(biomod, file = 'data/biomod.RData', compress = 'xz')

##
# physiological/abundance

# data
data(envdat)
data(ptedat)

envchr <- c('pCO2', 'pH', 'CO3', 'Ara', 'O2', 'Temp', 'Fluor')
ptechr <- c('abu', 'dis', 'len', 'ty2', 'ty3', 'scr')

dat_frm <- ptedat %>% 
  select(one_of('CTD', ptechr)) %>% 
  left_join(envdat, by = 'CTD') %>% 
  data.frame %>% 
  remove_rownames %>% 
  column_to_rownames('CTD')

# models

# variables to model, chr strings only
env_cmb <- envchr %>% 
  combn(2) %>% 
  t %>% 
  data.frame(stringsAsFactors = F) %>%
  crossing(ptechr, .) %>% 
  t %>% 
  data.frame(stringsAsFactors = F) %>% 
  as.list

# # back-transform abu for models
# dat_frm <- dat_frm %>% 
#   mutate(
#     abu = (10^abu) - 1,
#     abu = as.integer(round(abu, 0))
#   )
#   
# models
mod_all <- env_cmb %>%
  
  map(function(x, dat = dat_frm){

    # formula
    frm <- paste(x[1], '~', x[2], '*', x[3]) %>% 
      as.formula
    
    # if(x[1] == 'abu'){
    #   
    #   # zinf model
    #   out <- zeroinfl(frm, data = dat_frm)
    # 
    # } else {
    
    # model only
    out <- lm(frm, data = dat)
    
    # }
    
    return(out)
    
  }) %>% 
  enframe('Model', 'Modobj') %>% 
  mutate(
    Model = gsub('^X', 'mod', Model), 
    data = map(Modobj, function(x){
      
      # response
      yvar <- all.vars(formula(x))[1]
      
      # rsq
      rsq <- summary(x) %>% 
        .$r.squared
      
      # coeff summary
      out <- x %>% 
        summary %>% 
        .$coefficients %>% 
        .[-1, c(1, 4)] %>% 
        data.frame %>% 
        rownames_to_column('env_lab') %>% 
        mutate(
          pte_lab = yvar,
          Rsq = round(rsq, 2),
          Est = round(Estimate, 2), 
          Pvl = p_ast(Pr...t..)
        ) %>% 
        select(pte_lab, env_lab, Rsq, Est, Pvl)
      
      return(out)
      
    })
  ) 

phymod <- mod_all
save(phymod, file = 'data/phymod.RData', compress = 'xz')
