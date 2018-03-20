---
output:
  html_document:
    keep_md: yes
    self_contained: no
    code_folding: hide
---

# Pteropod stressor interaction: figures
  

```r
library(tidyverse)
library(vegan)
library(ggord)
library(scales)
library(effects)
library(knitr)
library(gridExtra)

source("R/funcs.R")

opts_chunk$set(fig.align = 'center', message = F, echo = T, cache = F, dev = 'png', dev.args = list(family = 'serif'), dpi = 400, fig.pos = '!h', warning = F, background = 'white', out.width = '100%',
               fig.process = function(x) {
                 x2 = sub('-\\d+([.][a-z]+)$', '\\1', x)
                 if (file.rename(x, x2)) x2 else x
               })

data(envdat)
data(ptedat)

envchr <- c('Lat', 'pCO2', 'pH', 'CO3', 'Ara', 'O2', 'Temp', 'Fluor')
biochr <- c('CAT', 'GR', 'GSHonGSSG', 'GST', 'LPX', 'ORAC', 'SOD', 'ORACvLPX')
phychr <- c('abu', 'dis', 'len', 'ty2', 'ty3', 'scr')
```


```r
# biomarker rda model
dat_bio <- ptedat %>% 
  select(one_of('CTD', biochr)) %>% 
  na.omit %>% 
  inner_join(envdat, by = 'CTD') %>% 
  select(one_of(c('CTD', envchr, biochr))) %>% 
  data.frame %>% 
  remove_rownames %>% 
  column_to_rownames('CTD')

envphy <- select(dat_bio, one_of(envchr)) %>% 
  select(-Lat) %>% 
  decostand(method = 'range')
ptephy <- select(dat_bio, one_of(biochr)) %>% 
  decostand(method = 'range')

biomod <- rda(ptephy, envphy)

# physiology rda model
dat_phy <- ptedat %>% 
  select(one_of('CTD', phychr)) %>% 
  na.omit %>% 
  inner_join(envdat, by = 'CTD') %>% 
  select(one_of(c('CTD', envchr, phychr))) %>% 
  data.frame %>% 
  remove_rownames %>% 
  column_to_rownames('CTD')

envphy <- select(dat_phy, one_of(envchr)) %>% 
  select(-Lat) %>% 
  decostand(method = 'range')
ptephy <- select(dat_phy, one_of(phychr)) %>% 
  decostand(method = 'range')

phymod <- rda(ptephy, envphy)
```

RDA plots:

```r
# globals
coord_fix <- F
repel <- F
addsize <- 3
alpha <- 0.8
arrow <- 0.3
cols <- 'lightblue'

p1a <- ggord(biomod, ptslab = T, parse = T, repel = repel, coord_fix = coord_fix, addsize = addsize, size = dat_bio$Lat, sizelab = 'Latitude', alpha = alpha, arrow = arrow, ellipse = F) + 
  ggtitle('(a) Cellular stressors')
p2a <- ggord(phymod, ptslab = T, parse = T, repel = repel, coord_fix = coord_fix, addsize = addsize, size = dat_phy$Lat, sizelab = 'Latitude', alpha = alpha, arrow = arrow, ellipse = F) + 
  ggtitle('(b) Physiological stressors')

grid.arrange(p1a, p2a, ncol = 2)
```

<img src="figures_files/figure-html/rda_plt.png" width="100%" style="display: block; margin: auto;" />

