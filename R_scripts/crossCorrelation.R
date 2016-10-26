#################################
# Cross-correlation coefficient #
#################################

crossCorrelate <- function(t, s, lag=0){
  t <- as.numeric(t)
  s <- as.numeric(s)
  
  # use the cross-correlation function, return max.lag
  # the R function returns the normalised cross-correlation
  # coefficient, unlike the Numpy version <- no need to 
  # pre-normalise the time series values
  if(lag == 0){
    x.cor <- ccf(t, s, lag.max=0, plot=F)$acf[1]
  }
  else{
    x.cor <- ccf(t, s, lag.max=lag, plot=F)$acf[(lag*2) + 1]
  }
  return(x.cor)
}

