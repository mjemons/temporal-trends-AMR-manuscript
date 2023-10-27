set.seed(1234)

AICc = function(m){
  k = length(coefficients(m))
  n = length(residuals(m))
  return(2*k - 2*logLik(m)+(2*k^2 + 2*k)/(n-k-1))
  
}


fitting0 = function(tmp1){
  # move timeseries downwards (fitting more robust)
  tmp1$Year = tmp1$Year - 2000
  # fitting
  m = tryCatch({m = lm(p~1,data = tmp1, weights = tmp1$N)},
               warning = function(m){m},
               error = function(m){NA},
               finally = function(m){m})
  # detect fails
  if(length(m)==1){to.return = c(NA,NA,NA)}
  #write output
  else{
    mean.abs.error = mean(abs(coef(m)-tmp1$p))
    # c(intercept, AIC, AICc, mean.abs.error)
    to.return = c(coef(m)[1],AICc(m),mean.abs.error)}
    to.return = as.data.frame(t(to.return))
    colnames(to.return) = c('m0.plateau',"m0.AICc","m0.mean.abs.error")
    return(to.return)
  }

fitting1 = function(tmp1){
  # move timeseries downwards (fitting more robust)
  tmp1$Year = tmp1$Year - 2000
  meanp <- mean(tmp1$p, na.rm = T)
  # fitting
  m = tryCatch({m = nlsLM(p~1/(1+exp(b0-b1*Year)),data = tmp1,start = list(b0 = log((1-meanp)/meanp), b1 = 0),
                          weights = tmp1$N,
                          control = nls.lm.control(maxiter = 1000))},
               warning = function(m){m},
               error = function(m){NA},
               finally = function(m){m})
  # detect fails
  if(length(m)==1){to.return = c(NA,NA,NA,NA,NA,NA,NA,NA)}
  #write output
  else{
    w = tryCatch({test = c(summary(m)$parameters[2,1]-qt(0.975, df = summary(m)$df[2])*summary(m)$parameters[2,2], summary(m)$parameters[2,1] + qt(0.975, df = summary(m)$df[2])*summary(m)$parameters[2,2])},
                 warning = function(test){test},
                 error = function(test){c(NA,NA)},
                 finally = function(test){test})
    r.squared = 1 - sum(summary(m)$residuals^2)/sum((tmp1$p-sum(tmp1$p*m$weights/sum(m$weights)))^2*m$weights)
    mean.abs.error = mean(abs(1/(1+exp(coef(m)[1]-coef(m)[2]*tmp1$Year))-tmp1$p))
    # correct b0 because of shift of timeseries:
    b0 = coef(m)[1] + coef(m)[2]*2000
    # c(intercept,slope,slope.high,slope.low,pval,AIC,AICc,abs.error,r.squared)
    to.return = c(b0,coef(m)[2],w,summary(m)$parameters[2,4],AICc(m),mean.abs.error,r.squared)}
    to.return = as.data.frame(t(to.return))
    colnames(to.return) = c("m1.intercept",'m1.slope',"m1.slope.lower","m1.slope.upper","m1.slope.pval","m1.AICc","m1.mean.abs.error","m1.r2")
    return(to.return)
}

fitting2 =function(tmp1,startval,startval2){
  # move timeseries downwards (fitting more robust)
  tmp1$Year = tmp1$Year - 2000
  A = NA
  xx = 1
  initial.val.multipliers = c(1,1.25,1.5,2,3,5,10,15)
  # fitting
  while(length(A) == 1 & xx < 9){
    
    A = tryCatch({m = nlsLM(p~K/(1+exp(b0-b1*Year)),data = tmp1,start = list(b0 = 0,b1 = startval*initial.val.multipliers[xx], K= startval2),
                            weights = tmp1$N,
                            control = nls.lm.control(maxiter = 1000),upper = c(Inf,Inf,1),lower = c(-Inf,0,0))},
                 warning = function(m){m},
                 error = function(m){c(NA)},
                 finally = function(m){m})
    if(length(A) == 1){xx = xx+1}
  }
  if(length(A)==1){to.return = c(NA,NA,NA,NA,NA,NA,NA,NA,NA)}
  else{
    w = tryCatch({test = confint(A)[2,]},
                 warning = function(test){test},
                 error = function(test){c(NA,NA)},
                 finally = function(test){test})
    # if both error an warning, does not work
    if(typeof(w) == 'list'){w = c(NA,NA)}
    mean.abs.error = mean(abs(coefficients(A)[3]/(1+exp(coefficients(A)[1]-coefficients(A)[2]*tmp1$Year))-tmp1$p))
    r.squared = 1 - sum(summary(A)$residuals^2)/sum((tmp1$p-sum(tmp1$p*A$weights/sum(A$weights)))^2*A$weights)
    # correct b0 because of shift of timeseries:
    b0 = coefficients(A)[1] + coefficients(A)[2]*2000
    # c(intercept,plateau,slope,slope.high,slope.low,AIC,AICc,abs.error,r.squared)
    to.return = c(b0,coefficients(A)[3],coefficients(A)[2],w,summary(A)$parameters[2,4],AICc(m),mean.abs.error,r.squared)}
    to.return = as.data.frame(t(to.return))
    colnames(to.return) = c("m2.intercept","m2.plateau",'m2.slope',"m2.slope.lower","m2.slope.upper","m2.slope.pval","m2.AICc","m2.mean.abs.error","m2.r2")
    return(to.return)
}