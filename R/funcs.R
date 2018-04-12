# function for formatting p-values in tables
p_ast <- function(x){
  
  sig_cats <- c('**', '*', 'ns')
  sig_vals <- c(-Inf, 0.005, 0.05, Inf)
  
  out <- cut(x, breaks = sig_vals, labels = sig_cats, right = FALSE)
  out <- as.character(out)
  
  return(out)
  
}

# vif function
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  library(fmsb)
  
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

# Data to plot for effects
#
# modin lm model
# cvar chr string of variable to hold constant
# pos is where the labels are, left or right of effects line
# fct is scaling factor for labels from end of lines
get_pldat <- function(modin, cvar, pos = c('left', 'right'), fct = NULL){
  
  pos <- match.arg(pos)
  
  # crossing of model data by range
  x <- modin$model %>% 
    .[, -1] %>% 
    data.frame %>% 
    as.list %>% 
    map(range) %>%
    map(function(x) seq(x[1], x[2], length = 100))
  
  # quantiles for cvar
  x[[cvar]] <- modin$model[[cvar]]%>% quantile(., c(0, 1))
  
  # make data frame
  nms <- names(x) 
  x <- crossing(x[[1]], x[[2]])
  names(x) <- nms
  x <- x[, c(names(x)[!names(x) %in% cvar], cvar)]

  # get predictions, combine with exp vars
  prd_vl <- predict(modin, newdata = x, se = T) %>% 
    data.frame(., x) %>% 
    dplyr::select(-df, -residual.scale) %>% 
    mutate(
      hi = fit + se.fit, 
      lo = fit - se.fit
      )
  names(prd_vl)[1] <- all.vars(formula(modin))[1]
  
  # min x axis values for quantile labels
  yvar <- names(prd_vl)[1]
  xvar <- all.vars(formula(modin))
  xvar <- xvar[!xvar %in% c(yvar, cvar)]

  locs <- prd_vl %>% 
    group_by(.dots = list(cvar))   
  if(pos == 'right'){
    if(is.null(fct)) fct <- 1.05
    locs <- filter(locs, row_number() == n())
  } else {
    if(is.null(fct)) fct <- 0.95
    locs <- filter(locs, row_number() == 1)
  }
  
  yval <- locs[[yvar]]
  xval <- locs[[xvar]] %>% unique %>% `*`(fct)
  xlab <- data.frame(
    lab = c('Max', 'Min'), 
    x = xval, y = yval,
    stringsAsFactors = F)
  dr <- locs[[cvar]] %>% range %>% diff %>% sign
  if(dr == 1) xlab$lab <- rev(xlab$lab)
  
  # output
  out <- list(prd_vl = prd_vl, xlab = xlab)
  return(out)
  
}