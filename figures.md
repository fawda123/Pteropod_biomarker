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
library(grid)

source("R/funcs.R")

opts_chunk$set(fig.align = 'center', message = F, echo = T, cache = F, dev = 'png', dev.args = list(family = 'serif'), dpi = 400, fig.pos = '!h', warning = F, background = 'white', out.width = '100%',
               fig.process = function(x) {
                 x2 = sub('-\\d+([.][a-z]+)$', '\\1', x)
                 if (file.rename(x, x2)) x2 else x
               })

data(envdat)
data(ptedat)
data(phymod)
data(biomod)

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

biorda<- rda(ptephy, envphy)

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

phyrda <- rda(ptephy, envphy)
```


```r
# globals
coord_fix <- F
repel <- F
addsize <- 3
alpha <- 0.8
arrow <- 0.3
cols <- 'lightblue'

p1a <- ggord(biorda, ptslab = T, parse = T, repel = repel, coord_fix = coord_fix, addsize = addsize, size = dat_bio$Lat, sizelab = 'Latitude', alpha = alpha, arrow = arrow, ellipse = F) + 
  ggtitle('(a) Cellular endpoints')
p2a <- ggord(phyrda, ptslab = T, parse = T, repel = repel, coord_fix = coord_fix, addsize = addsize, size = dat_phy$Lat, sizelab = 'Latitude', alpha = alpha, arrow = arrow, ellipse = F) + 
  ggtitle('(b) Physiological endpoints')

grid.arrange(p1a, p2a, ncol = 2)
```

<img src="figures_files/figure-html/rdaplo.png" width="100%" style="display: block; margin: auto;" />

Fig. 1 Results of redundancy analyses for environmental variables with (a) cellular and (b) physiological endpoints of pteropod response to OA stressors.  Points are site locations in multivariate space with the size proportional to latitude. Separate RDAs were created for cellular and physiological endpoints because not all data were available across all stations.   


```r
# physiology rda model
dat_cor <- ptedat %>% 
  select(one_of('CTD', phychr, biochr)) %>% 
  # na.omit %>% 
  left_join(envdat, by = 'CTD') %>% 
  data.frame %>% 
  remove_rownames %>% 
  column_to_rownames('CTD')

# all correlations
crs <- crossing(var1 = names(dat_cor), var2 = names(dat_cor)) %>% 
  filter(var1 != var2) %>% 
  rownames_to_column() %>% 
  group_by(rowname) %>% 
  nest %>% 
  mutate(
    crs = map(data, function(x){
      
      # variables
      vr1 <- dat_cor[[x$var1]]
      vr2 <- dat_cor[[x$var2]]
      
      # pearson
      pr_ts <- cor.test(vr1, vr2, method = 'pearson')
      pr_cr <- round(pr_ts$estimate, 2)
      pr_pv <- p_ast(pr_ts$p.value)
      pr <- paste(pr_cr, pr_pv)
    
      out <- data.frame(pr = pr, stringsAsFactors = F)
      return(out)
      
    })
  ) %>% 
  unnest %>% 
  select(-rowname)
```



```r
levs <- c(sort(envchr), sort(biochr), sort(phychr))
labs <- c('Omega[ar]', 'CO[3]^-2', 'chla', 'Lat', 'O[2]', 'pCO[2]', 'pH', 'Temp', 'CAT', 'GR', 'GSHonGSSG', 'GST', ' LPX', 'ORAC', 'ORACvLPX', 'SOD', 'abundance', 'dissolution', 'growth', 'scarring', 'typeII', 'typeIII')
prplo <- crs %>% 
  separate(pr, c('cor', 'sig'), sep = ' ') %>% 
  filter(var1 %in% levs & var2 %in% levs) %>%  
  mutate(
    cor = as.numeric(cor), 
    var1 = factor(var1, levels = rev(levs), labels = rev(labs)), 
    var2 = factor(var2, levels = rev(levs), labels = rev(labs)), 
    sig = gsub('ns', '', sig)
  )

pbase <- theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), 
  axis.text.y = element_text(size = 8),
  legend.position = c(0.5, 1.12),
  legend.direction = 'horizontal',
  plot.margin = unit(c(4,4,0,0), "lines"),
  strip.background = element_blank(), 
  strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5), 
  panel.background = element_rect(fill = 'black')
  ) 

outlab <- data.frame(
  y = c(3.5, 10.5, 18.5), 
  lab = c('Population/\nPhysiology', 'Cellular', 'Environment')
)

p <- ggplot(prplo) + 
  geom_tile(aes(y = var1, x = var2, fill = cor), colour = 'black') + 
  geom_text(aes(y = var1, x = var2, label = sig)) +
  annotation_custom(grob = textGrob(label = outlab$lab[1], hjust = 0, gp = gpar(cex = 0.7)),
                    ymin = outlab$y[1], ymax = outlab$y[1], xmin = 23, xmax = 23) +
  annotation_custom(grob = textGrob(label = outlab$lab[2], hjust = 0, gp = gpar(cex = 0.7)),
                    ymin = outlab$y[2], ymax = outlab$y[2], xmin = 23, xmax = 23) +  
  annotation_custom(grob = textGrob(label = outlab$lab[3], hjust = 0, gp = gpar(cex = 0.7)),
                    ymin = outlab$y[3], ymax = outlab$y[3], xmin = 23, xmax = 23) +
  annotation_custom(grob = textGrob(label = outlab$lab[1], hjust = 0.5, gp = gpar(cex = 0.7)),
                    xmin = outlab$y[1], xmax = outlab$y[1], ymin = 23.5, ymax = 23.5) +
  annotation_custom(grob = textGrob(label = outlab$lab[2], hjust = 0.5, gp = gpar(cex = 0.7)),
                    xmin = outlab$y[2], xmax = outlab$y[2], ymin = 23.5, ymax = 23.5) +  
  annotation_custom(grob = textGrob(label = outlab$lab[3], hjust = 0.5, gp = gpar(cex = 0.7)),
                    xmin = outlab$y[3], xmax = outlab$y[3], ymin = 23.5, ymax = 23.5) +
  pbase +
  scale_y_discrete('', expand = c(0, 0), labels = parse(text = rev(labs))) + 
  scale_x_discrete('', expand = c(0, 0), labels = parse(text = rev(labs))) +
  scale_fill_gradient2('Correlation', low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0) +
  guides(fill = guide_colourbar(barheight = 0.5, barwidth = 5, label.theme = element_text(size = 6, angle = 0))) +
  geom_hline(yintercept = 6.5, size = 1.5) +
  geom_hline(yintercept = 14.5, size = 1.5) +
  geom_vline(xintercept = 6.5, size = 1.5) +
  geom_vline(xintercept = 14.5, size = 1.5) 

# Code to override clipping
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
```

<img src="figures_files/figure-html/corplo.png" width="100%" style="display: block; margin: auto;" />
Fig. 2 Correlation matrix of environmental variables, cellular response endpoints, and population (abundance) and physiological response endpoints for pteropods.  Darker red values are strong positive correlations and darker purple values are strong negative correlations.  Significance values at alpha = 0.05 are shown by stars (p < 0.05 \*, p < 0.005 \*).


```r
# cellular 92 (lpx, ara, ph - syn), 101 (lpx, ara, t - neg add), 164 (sod, ara, t, neg but no thrsh) 
biosel <- biomod %>% 
  filter(Model %in% paste0('mod', c(92, 101, 164))) %>% 
  dplyr::select(-data)

cols <- RColorBrewer::brewer.pal(9, 'RdBu')

##
# plots

# mod 92 cell (lpx, ara, ph - syn)
pl1 <- filter(biosel, grepl('92$', Model)) %>% 
  .$Modobj %>% 
  .[[1]] %>% 
  get_pldat(., 'pH')

p1 <- ggplot() +
  geom_line(data = pl1[[1]], aes(x = Ara, y = LPX, group = pH, colour = pH), size = 1) + 
  geom_text(data = pl1[[2]], aes(x= x, y = y, label = lab), hjust = 0) +
  theme_bw() +
  scale_x_continuous(expression(Omega[Ar])) +
  scale_colour_gradientn(colours = cols) + 
  scale_y_continuous('LPX') +
  ggtitle("(a) Synergistic effect of pH on aragonite\nsaturation")

# mod 101 cell (lpx, ara, t - neg add), 164 (sod, ara, t, neg but no thrsh) 
pl2 <- filter(biosel, grepl('101$', Model)) %>% 
  .$Modobj %>% 
  .[[1]] %>% 
  get_pldat(., 'Temp')

p2 <- ggplot() +
  geom_line(data = pl2[[1]], aes(x = Ara, y = LPX, group = Temp, colour = Temp), size = 1) + 
  geom_text(data = pl2[[2]], aes(x= x, y = y, label = lab), hjust = 0) +
  theme_bw() +
  scale_x_continuous(expression(Omega[Ar])) +
  scale_colour_gradientn(colours = cols) + 
  scale_y_continuous('LPX') +
  ggtitle("(b) Negative additive effect below threshold of\n temperature on aragonite saturation")

# mod 164 cell (sod, ara, t, neg but no thrsh) 
pl3 <- filter(biosel, grepl('164$', Model)) %>% 
  .$Modobj %>% 
  .[[1]] %>% 
  get_pldat(., 'Temp')

p3 <- ggplot() +
  geom_line(data = pl3[[1]], aes(x = Ara, y = SOD, group = Temp, colour = Temp), size = 1) + 
  geom_text(data = pl3[[2]], aes(x= x, y = y, label = lab), hjust = 0) +
  theme_bw() +
  scale_x_continuous(expression(Omega[Ar])) +
  scale_colour_gradientn(colours = cols) + 
  scale_y_continuous('SOD') +
  ggtitle("(c) Negative additive effect of temperature on \naragonite saturation")

grid.arrange(p1, p2, p3, ncol = 2)
```

<img src="figures_files/figure-html/effbio.png" width="100%" style="display: block; margin: auto;" />
Fig. 3 Examples of model interactions of co-occuring environmental variables on cellular response measures.  Each subplot shows a different relationship as either additive or synergistic effects between the variables. All y-axes are transformed to conform to model output. Covarying environmental variables were held constant at the minimum, 25th, median, 75th, and maximum values in the observed data.


```r
# physio 4 (abu, o2, pco2 - pos additive), 38 (dis, ara, t - syn)
physel <- phymod %>% 
  filter(Model %in% paste0('mod', c(4, 38))) %>% 
  dplyr::select(-data)

cols <- RColorBrewer::brewer.pal(9, 'RdBu')

# mod 4 physio (abu, o2, pco2 - pos additive)
pl4 <- filter(physel, grepl('4$', Model)) %>% 
  .$Modobj %>% 
  .[[1]] %>% 
  get_pldat(., 'pCO2', fct = 0.9)

p4 <- ggplot() +
  geom_line(data = pl4[[1]], aes(x = O2, y = abu, group = pCO2, colour = pCO2), size = 1) + 
  geom_text(data = pl4[[2]], aes(x= x, y = y, label = lab), hjust = 0) +
  theme_bw() +
  scale_colour_gradientn(colours = cols) + 
  scale_y_continuous('Abundance') +
  ggtitle("(a) Synergistic effect of O2 on\npCO2 stress")

# mod 38 physio (dis, ara, t - syn)
pl5 <- filter(physel, grepl('38$', Model)) %>% 
  .$Modobj %>% 
  .[[1]] %>% 
  get_pldat(., 'Temp', 'right', fct = 1.05)

p5 <- ggplot() +
  geom_line(data = pl5[[1]], aes(x = Ara, y = dis, group = Temp, colour = Temp), size = 1) + 
  geom_text(data = pl5[[2]], aes(x= x, y = y, label = lab), hjust = 1) +
  theme_bw() +
  scale_colour_gradientn(colours = cols) + 
  scale_y_continuous('Shell dissolution') +
  scale_x_continuous(expression(Omega[Ar])) +
  ggtitle("(b) Synergistic effect of temperature on \naragonite saturation")

grid.arrange(p4, p5, ncol = 2)
```

<img src="figures_files/figure-html/effphy.png" width="100%" style="display: block; margin: auto;" />
Fig. 4 Examples of model interactions of co-occuring environmental variables on abundance and shell dissolution. Each subplot shows a different relationship as either additive or synergistic effects between the variables. All y-axes are transformed to conform to model output. Covarying environmental variables were held constant at the minimum, 25th, median, 75th, and maximum values in the observed data.


```r
mods <- readxl::read_excel('raw/info for table.xlsx') %>% 
  .[-1, -1] %>% 
  gather('var', 'Model', everything()) %>% 
  mutate(
    org = ifelse(var %in% c('abundance', 'dissolution', 'length'), 'phy', 'bio'),
    Model = paste0('mod', Model)
    ) %>% 
  split(.$org)

biotab <- mods$bio %>% 
  left_join(biomod, by = 'Model') %>% 
  select(Model, data) %>%
  mutate(data = map(data, function(x){
    x %>% mutate(vr = seq(1:nrow(.)))
  })) %>% 
  unnest %>% 
  mutate(
    Model = gsub('mod', '', Model),
    Model = as.numeric(Model)
    ) %>% 
  group_by(Model) %>% 
  mutate(
    Pvl = gsub('ns', '', Pvl), 
    Rsq = ifelse(duplicated(Rsq), '', Rsq)
    ) %>% 
  ungroup %>% 
  unite('Est', Est, Pvl, sep = '') %>% 
  gather('estnm', 'estvl', Rsq, Est) %>% 
  unite('est', pte_lab, estnm, sep = ', ', remove = F) %>% 
  select(-estnm) %>% 
  split(.$pte_lab) %>%
  lapply(., function(x){
    
    out <- spread(x, est, estvl) %>% 
      select(-pte_lab) %>%
      arrange(Model, vr) %>% 
      select(-vr, -Model)
    
    if(!any(grepl('^LPX', names(out))))
      out <- out %>% 
        select(-env_lab)
    
    return(out)  
    
  }) %>% 
  do.call('cbind', .) %>% 
  mutate(
    Model = rep(letters[1:(nrow(.)/3)], each = 3),
    Model = paste0('(', Model, ')'),
    Model = ifelse(duplicated(Model), '', Model)
  ) %>% 
  select(Model, everything())

names(biotab) <- gsub('^.*\\.', '', names(biotab))
vrs <- gsub('(^.*),\\s.*$', '\\1', names(biotab))
vrs[duplicated(vrs)] <- ''
vrs <- paste0(vrs, '<br>')
vrs[vrs %in% 'env_lab<br>'] <- 'Parameter<br>'
vrs2 <- gsub('^.*,\\s|Model|env_lab', '', names(biotab))
vrs <- paste(vrs, vrs2, sep = '')
  
names(biotab) <- vrs

knitr::kable(biotab, caption = "Table 1: Model results for pteropod cellular response to pairs of co-occurring environmental variables. The estimated joint effects of variables in each model are shown.  The overall R-squred value for each model is also shown.  Significance of each effect is noted as * p < 0.05, ** p < 0.005.")
```



Table: Table 1: Model results for pteropod cellular response to pairs of co-occurring environmental variables. The estimated joint effects of variables in each model are shown.  The overall R-squred value for each model is also shown.  Significance of each effect is noted as * p < 0.05, ** p < 0.005.

Model<br>   Parameter<br>   LPX<br>Est   <br>Rsq   ORAC<br>Est   <br>Rsq   ORACvLPX<br>Est   <br>Rsq   SOD<br>Est   <br>Rsq 
----------  --------------  -----------  --------  ------------  --------  ----------------  --------  -----------  --------
(a)         pCO2            0.05*        0.65      0.31          0.57      -0.86             0.51      0.19         0.52    
            Ara             17.71                  246.29                  -271.39                     72.58                
            pCO2:Ara        -0.03                  0.06                    0.43                        -0.03                
(b)         pCO2            0            0.87      -0.61         0.15      0.51              0.91      -0.23        0.57    
            O2              -0.04                  -1.49                   2.66*                       -0.96                
            pCO2:O2         0                      0                       0                           0                    
(c)         pCO2            0.09         0.62      -0.5          0.47      -0.98             0.55      -0.23        0.66    
            Temp            4.77                   5.68                    -58.87                      -5.39                
            pCO2:Temp       -0.01                  0.06                    0.08                        0.03                 
(d)         pH              -101.03*     0.72      -810.38       0.58      1876.68           0.6       -475.74      0.59    
            Ara             -406.15*               555                     6496.4                      -619.83              
            pH:Ara          51.79*                 -25.27                  -835.38                     90.07                
(e)         pH              11.41        0.94      1794.69       0.29      -1034.57          0.91      584.91       0.68    
            O2              -1.18*                 33.71                   7.48                        7.4                  
            pH:O2           0.13*                  -4.46                   -0.53                       -1.06                
(f)         pH              -163.87      0.68      824.67        0.45      2025.87           0.62      240.63       0.65    
            Temp            -123.79                842.19                  1314.98                     358.65               
            pH:Temp         15.63                  -101.23                 -168.02                     -43.4                
(g)         Ara             -14.85       0.88      383.89        0.51      134.38            0.8       95.75        0.71    
            O2              -0.19**                -0.1                    3.32*                       -0.41                
            Ara:O2          0.09*                  -0.71                   -1.22                       -0.14                
(h)         Ara             -62.93       0.72      206.6         0.4       681.44            0.68      59.34        0.61    
            Temp            -8.68                  79.92                   36.91                       46.48                
            Ara:Temp        5.94                   -26.52                  -49.42                      -15.21               
(i)         Ara             -4.72        0.58      66.14         0.28      137.19            0.69      -19.6        0.27    
            Fluor           -5.35                  -133.22                 345.5                       170.06               
            Ara:Fluor       -2.21                  149.85                  -60.2                       -132.48              
(j)         O2              -0.33*       0.81      1.64          0.43      3.98              0.76      1.13         0.73    
            Temp            -6.68                  77.43                   53.74                       46.45                
            O2:Temp         0.03                   -0.2                    -0.32                       -0.16                
(k)         O2              -0.04*       0.79      0.13          0.08      1.06**            0.91      -0.27        0.52    
            Fluor           -14.83                 -256.69                 656.19                      -85.16               
            O2:Fluor        0.03                   1.25                    -1.96                       0.26                 
(l)         Temp            -0.67        0.22      21.36         0.61      21.32             0.26      1.35         0.04    
            Fluor           47.63                  -2603.6                 -756.16                     610.65               
            Temp:Fluor      -5.75                  294.24                  101.91                      -63.93               


```r
phytab <- mods$phy %>% 
  left_join(phymod, by = 'Model') %>% 
  select(Model, data) %>%
  mutate(data = map(data, function(x){
    x %>% mutate(vr = seq(1:nrow(.)))
  })) %>% 
  unnest %>% 
  mutate(
    Model = gsub('mod', '', Model),
    Model = as.numeric(Model)
  ) %>% 
  group_by(Model) %>% 
  mutate(
    Pvl = gsub('ns', '', Pvl), 
    Rsq = ifelse(duplicated(Rsq), '', Rsq)
  ) %>% 
  ungroup %>% 
  unite('Est', Est, Pvl, sep = '') %>% 
  gather('estnm', 'estvl', Rsq, Est) %>% 
  unite('est', pte_lab, estnm, sep = ', ', remove = F) %>% 
  select(-estnm) %>% 
  split(.$pte_lab) %>%
  lapply(., function(x){
    
    out <- spread(x, est, estvl) %>% 
      select(-pte_lab) %>%
      arrange(Model, vr) %>% 
      select(-vr, -Model)
    
    if(!any(grepl('^abu', names(out))))
      out <- out %>% 
        select(-env_lab)
    
    return(out)  
    
  }) %>% 
  do.call('cbind', .) %>% 
  mutate(
    Model = rep(letters[1:(nrow(.)/3)], each = 3),
    Model = paste0('(', Model, ')'),
    Model = ifelse(duplicated(Model), '', Model)
  ) %>% 
  select(Model, everything())

names(phytab) <- gsub('^.*\\.', '', names(phytab))
vrs <- gsub('(^.*),\\s.*$', '\\1', names(phytab))
vrs[duplicated(vrs)] <- ''
vrs <- paste0(vrs, '<br>')
vrs[vrs %in% 'env_lab<br>'] <- 'Parameter<br>'
vrs2 <- gsub('^.*,\\s|Model|env_lab', '', names(phytab))
vrs <- paste(vrs, vrs2, sep = '')
vrs <- gsub('abu', 'Abundance', vrs)
vrs <- gsub('len', 'Length', vrs)
vrs <- gsub('dis', 'Shell dissolution', vrs)

names(phytab) <- vrs

knitr::kable(phytab, caption = "Table 2: Model results for pteropod physiological and population response to pairs of co-occurring environmental variables. The estimated joint effects of variables in each model are shown.  The overall R-squred value for each model is also shown.  Significance of each effect is noted as * p < 0.05, ** p < 0.005.")
```



Table: Table 2: Model results for pteropod physiological and population response to pairs of co-occurring environmental variables. The estimated joint effects of variables in each model are shown.  The overall R-squred value for each model is also shown.  Significance of each effect is noted as * p < 0.05, ** p < 0.005.

Model<br>   Parameter<br>   Abundance<br>Est   <br>Rsq   Shell dissolution<br>Est   <br>Rsq   Length<br>Est   <br>Rsq 
----------  --------------  -----------------  --------  -------------------------  --------  --------------  --------
(a)         pCO2            -0.01*             0.17      0*                         0.83      0               0.08    
            Ara             -3.1*                        0.4*                                 -1.53                   
            pCO2:Ara        0                            0*                                   0                       
(b)         pCO2            0                  0.32      0                          0.91      0.01            0.19    
            O2              0.02                         0                                    0.02                    
            pCO2:O2         0*                           0**                                  0                       
(c)         pCO2            0                  0.44      0*                         0.82      0               0.11    
            Temp            -0.23                        0.12*                                0                       
            pCO2:Temp       0                            0                                    0                       
(d)         pH              14.66*             0.29      -2.18*                     0.84      6.61            0.08    
            Ara             39.44                        -9.36*                               4.03                    
            pH:Ara          -5.25                        1.19*                                -0.74                   
(e)         pH              4.44               0.14      -0.44                      0.93      -16.71          0.39    
            O2              0.36                         -0.06**                              0.1                     
            pH:O2           -0.04                        0.01**                               -0.01                   
(f)         pH              3.1                0.48      -3.7*                      0.86      -1.93           0.1     
            Temp            -0.34                        -2.57                                -3.8                    
            pH:Temp         0                            0.33                                 0.45                    
(g)         Ara             1.63               0.24      -0.41                      0.9       -2.24           0.23    
            O2              0.03*                        -0.01**                              0.02                    
            Ara:O2          -0.01                        0*                                   0                       
(h)         Ara             1.47               0.54      -1.27*                     0.91      0.52            0.1     
            Temp            -0.43                        -0.11                                -0.41                   
            Ara:Temp        0                            0.1*                                 0.05                    
(i)         Ara             0.49               0.03      0.06                       0.69      0.26            0.06    
            Fluor           8.24                         4.93                                 1.73                    
            Ara:Fluor       -5.38                        -3.59                                -0.3                    
(j)         O2              0                  0.46      -0.01*                     0.91      -0.01           0.16    
            Temp            -0.41                        -0.15                                -0.65                   
            O2:Temp         0                            0*                                   0                       
(k)         O2              0.01               0.12      0                          0.8       0.01            0.15    
            Fluor           14.87                        3.62                                 23.3                    
            O2:Fluor        -0.07                        -0.02                                -0.1                    
(l)         Temp            -0.24              0.18      0.02                       0.3       0               0.04    
            Fluor           -7.4                         8.15                                 1.15                    
            Temp:Fluor      0.78                         -0.88                                -0.02                   

