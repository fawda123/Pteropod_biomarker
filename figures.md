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

Heatmaps of correlations:


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
prplo <- crs %>% 
  separate(pr, c('cor', 'sig'), sep = ' ') %>% 
  filter(var1 %in% levs & var2 %in% levs) %>% 
  mutate(
    cor = as.numeric(cor), 
    var1 = factor(var1, levels = rev(levs)), 
    var2 = factor(var2, levels = rev(levs)), 
    sig = gsub('ns', '', sig)
  )

pbase <- theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), 
  axis.text.y = element_text(size = 8),
  legend.position = 'top', 
  plot.margin = unit(c(0,4,0,0), "lines"),
  strip.background = element_blank(), 
  strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5), 
  panel.background = element_rect(fill = 'black')
  ) 

outlab <- data.frame(
  y = c(3.5, 10.5, 18.5), 
  lab = c('Physiology', 'Cellular', 'Environment')
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
  pbase +
  scale_y_discrete('', expand = c(0, 0)) + 
  scale_x_discrete('', expand = c(0, 0)) +
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

