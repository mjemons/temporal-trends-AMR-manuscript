source("code/Utils.R")
source("code/logistic_fitting_functions.R")

# This script generates estimates of uncertainty for the plateau and slope parameters of stabilising trajectories.
# It then compares the magnitude of the uncertainty to the magnitude of the estimates.


#### Defining functions ####

fitting2_uncertainty =function(tmp1,startval,startval2){
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
  if(length(A)==1){to.return = c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)}
  else{
    w = tryCatch({test = confint(A)[2,]},
                 warning = function(test){test},
                 error = function(test){c(NA,NA)},
                 finally = function(test){test})
    # if both error an warning, does not work
    if(typeof(w) == 'list'){w = c(NA,NA)}
    w2 = tryCatch({test = confint(A)[3,]},
                 warning = function(test){test},
                 error = function(test){c(NA,NA)},
                 finally = function(test){test})
    # if both error an warning, does not work
    if(typeof(w2) == 'list'){w = c(NA,NA)}
    mean.abs.error = mean(abs(coefficients(A)[3]/(1+exp(coefficients(A)[1]-coefficients(A)[2]*tmp1$Year))-tmp1$p))
    r.squared = 1 - sum(summary(A)$residuals^2)/sum((tmp1$p-sum(tmp1$p*A$weights/sum(A$weights)))^2*A$weights)
    # correct b0 because of shift of timeseries:
    b0 = coefficients(A)[1] + coefficients(A)[2]*2000
    # c(intercept,plateau,plateau.high,plateau.low,slope,slope.high,slope.low,AIC,AICc,abs.error,r.squared)
    to.return = c(b0,coefficients(A)[3],w2,coefficients(A)[2],w,summary(A)$parameters[2,4],AICc(m),mean.abs.error,r.squared)}
  to.return = as.data.frame(t(to.return))
  colnames(to.return) = c("m2.intercept","m2.plateau","m2.plateau.lower","m2.plateau.upper",'m2.slope',"m2.slope.lower","m2.slope.upper","m2.slope.pval","m2.AICc","m2.mean.abs.error","m2.r2")
  return(to.return)
}

fitting0_uncertainty = function(tmp1){
  # move timeseries downwards (fitting more robust)
  tmp1$Year = tmp1$Year - 2000
  # fitting
  m = tryCatch({m = lm(p~1,data = tmp1, weights = tmp1$N)},
               warning = function(m){m},
               error = function(m){NA},
               finally = function(m){m})
  # detect fails
  if(length(m)==1){to.return = c(NA,NA,NA,NA,NA)}
  #write output
  else{
    w2 = tryCatch({test = confint(m)},
                  warning = function(test){test},
                  error = function(test){c(NA,NA)},
                  finally = function(test){test})
    # if both error an warning, does not work
    if(typeof(w2) == 'list'){w = c(NA,NA)}
    mean.abs.error = mean(abs(coef(m)-tmp1$p))
    # c(intercept, lower, upper, AIC, AICc, mean.abs.error)
    to.return = c(coef(m)[1],w2,AICc(m),mean.abs.error)}
  to.return = as.data.frame(t(to.return))
  colnames(to.return) = c('m0.plateau',"m0.plateau.lower","m0.plateau.upper","m0.AICc","m0.mean.abs.error")
  return(to.return)
}


#### Get stabilising trajectories ####

# load the results of temporal trend fits and clean

summary_fits <- read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")
# rename for plotting + get rid of poor fits
summary_fits_toplot <- summary_fits[summary_fits$category!='poor fit',]
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns decreasing',]$category = "decreasing (ns)"
summary_fits_toplot[summary_fits_toplot$category == 'logistic s decreasing',]$category = 'decreasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns increasing',]$category = 'increasing (ns)' 
summary_fits_toplot[summary_fits_toplot$category == 'logistic s increasing',]$category = 'increasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'flat ',]$category = 'stable'
summary_fits_toplot[summary_fits_toplot$category == 'plateauing',]$category = 'stabilising'
summary_fits_toplot$category = factor(summary_fits_toplot$category, levels = c('decreasing (s)',"decreasing (ns)",'stable','stabilising','increasing (ns)','increasing (s)'))

# get stalising categories:

stabilising = summary_fits_toplot[which(summary_fits_toplot$category == 'stabilising'),]$combR

# get trajectories
amr_summary <- read.csv(file = "data/summary_AMR_filtered.csv")
amr_summary <- amr_summary %>% arrange(Year)

amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Country, amr_summary$Antibiotic, sep = "|")

### Get stabilising fits with uncertainty estimates

# create data frame for analysis
fit.flat = data.frame()
fit.log = data.frame()
fit.plateau = data.frame()
median.r = data.frame()

time_points = 5
observations = 30
minimum_N_R = 10

# fit flat, logistic and logistic with plateau models
for(i in 1:length(stabilising)){
  if(i/100 == floor(i/100)) print(i)
  mycomb <- stabilising[i]
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  if((nrow(tmp1) > time_points) & (all(tmp1$N >= observations)) & (sum(tmp1$N_I_R) >= minimum_N_R)){ # FB added condition on minimum number R isolates
    print(tmp1)
    f1 = cbind(mycomb,fitting0(tmp1))
    fit.flat = rbind(fit.flat,f1)
    f2 = cbind(mycomb,fitting1(tmp1))
    fit.log = rbind(fit.log,f2)
    f3 = cbind(mycomb,fitting2_uncertainty(tmp1,f2$m1.slope,mean(tmp1$p,na.rm = T)))
    fit.plateau = rbind(fit.plateau,f3) 
    f4 = ddply(tmp1,.(combR),summarise,median.amr = median(p,na.rm = T), min.year = min(Year, na.rm = T), max.year = max(Year, na.rm = T))
    f4=cbind(f4,Country = unique(tmp1$Country), Pathogen = unique(tmp1$Pathogen), Antibiotic = unique(tmp1$Antibiotic)) #ME added country, pathogen and antibiotic for forest plots
    median.r = rbind(median.r,f4)
  }
}

# compare uncertainty stabilising trajectories

fit.plateau$plateau.interval= fit.plateau$m2.plateau.upper - fit.plateau$m2.plateau.lower
fit.plateau$slope.interval= fit.plateau$m2.slope.upper - fit.plateau$m2.slope.lower

estimate = fit.plateau[which(!is.na(fit.plateau$plateau.interval)),]$m2.plateau
interval = fit.plateau[which(!is.na(fit.plateau$plateau.interval)),]$plateau.interval

plateau = data.frame(cbind(estimate,interval))

estimate = fit.plateau[!is.na(fit.plateau$slope.interval),]$m2.slope
interval = fit.plateau[!is.na(fit.plateau$slope.interval),]$slope.interval

slope = data.frame(cbind(estimate,interval))

p1 = ggplot(data = plateau) +
  geom_histogram(aes(x = estimate, fill = "estimate"), alpha = 0.5) +
  geom_histogram(aes(x = interval, fill = "interval"), alpha = 0.5) +
  labs(title = "Stabilising: plateau parameter", x = "Value", y = "Count")


p2 = ggplot(data = slope) +
  geom_histogram(aes(x = estimate, fill = "estimate"), alpha = 0.5) +
  geom_histogram(aes(x = interval, fill = "interval"), alpha = 0.5) +
  labs(title = "Stabilising: slope parameter", x = "Value", y = "Count")


# get uncertainty associated with slope for increasing trajectories

increasing = summary_fits_toplot[which(summary_fits_toplot$category == 'increasing (s)'),]$combR

# create data frame for analysis
fit.flat = data.frame()
fit.log = data.frame()
fit.plateau = data.frame()
median.r = data.frame()

# fit flat, logistic and logistic with plateau models
for(i in 1:length(increasing)){
  if(i/100 == floor(i/100)) print(i)
  mycomb <- increasing[i]
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  if((nrow(tmp1) > time_points) & (all(tmp1$N >= observations)) & (sum(tmp1$N_I_R) >= minimum_N_R)){ # FB added condition on minimum number R isolates
    print(tmp1)
    f1 = cbind(mycomb,fitting0(tmp1))
    fit.flat = rbind(fit.flat,f1)
    f2 = cbind(mycomb,fitting1(tmp1))
    fit.log = rbind(fit.log,f2)
    f3 = cbind(mycomb,fitting2_uncertainty(tmp1,f2$m1.slope,mean(tmp1$p,na.rm = T)))
    fit.plateau = rbind(fit.plateau,f3) 
    f4 = ddply(tmp1,.(combR),summarise,median.amr = median(p,na.rm = T), min.year = min(Year, na.rm = T), max.year = max(Year, na.rm = T))
    f4=cbind(f4,Country = unique(tmp1$Country), Pathogen = unique(tmp1$Pathogen), Antibiotic = unique(tmp1$Antibiotic)) #ME added country, pathogen and antibiotic for forest plots
    median.r = rbind(median.r,f4)
  }
}

fit.log$slope.interval= fit.log$m1.slope.upper - fit.log$m1.slope.lower

estimate = fit.log[!is.na(fit.log$slope.interval),]$m1.slope
interval = fit.log[!is.na(fit.log$slope.interval),]$slope.interval

increasing.slope = data.frame(cbind(estimate,interval))

p3 = ggplot(data = increasing.slope) +
  geom_histogram(aes(x = estimate, fill = "estimate"), alpha = 0.5) +
  geom_histogram(aes(x = interval, fill = "interval"), alpha = 0.5) +
  labs(title = "Increasing: slope parameter", x = "Value", y = "Count")

# get uncertainty associated with slope for stable trajectories

stable = summary_fits_toplot[which(summary_fits_toplot$category == 'stable'),]$combR

# create data frame for analysis
fit.flat = data.frame()
fit.log = data.frame()
fit.plateau = data.frame()
median.r = data.frame()

# fit flat, logistic and logistic with plateau models
for(i in 1:length(stable)){
  if(i/100 == floor(i/100)) print(i)
  mycomb <- stable[i]
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  if((nrow(tmp1) > time_points) & (all(tmp1$N >= observations)) & (sum(tmp1$N_I_R) >= minimum_N_R)){ # FB added condition on minimum number R isolates
    print(tmp1)
    f1 = cbind(mycomb,fitting0_uncertainty(tmp1))
    fit.flat = rbind(fit.flat,f1)
    f2 = cbind(mycomb,fitting1(tmp1))
    fit.log = rbind(fit.log,f2)
    f3 = cbind(mycomb,fitting2_uncertainty(tmp1,f2$m1.slope,mean(tmp1$p,na.rm = T)))
    fit.plateau = rbind(fit.plateau,f3) 
    f4 = ddply(tmp1,.(combR),summarise,median.amr = median(p,na.rm = T), min.year = min(Year, na.rm = T), max.year = max(Year, na.rm = T))
    f4=cbind(f4,Country = unique(tmp1$Country), Pathogen = unique(tmp1$Pathogen), Antibiotic = unique(tmp1$Antibiotic)) #ME added country, pathogen and antibiotic for forest plots
    median.r = rbind(median.r,f4)
  }
}

fit.flat$plateau.interval= fit.flat$m0.plateau.upper - fit.flat$m0.plateau.lower

estimate = fit.flat[!is.na(fit.flat$plateau.interval),]$m0.plateau
interval = fit.flat[!is.na(fit.flat$plateau.interval),]$plateau.interval

stable.plateau = data.frame(cbind(estimate,interval))

p4 = ggplot(data = stable.plateau) +
  geom_histogram(aes(x = estimate, fill = "estimate"), alpha = 0.5) +
  geom_histogram(aes(x = interval, fill = "interval"), alpha = 0.5) +
  labs(title = "Stable: plateau parameter", x = "Value", y = "Count")


### Plot everything
 
pdf("output/supporting/Uncertainty_estimates.pdf",width = 12,height = 10)
grid.arrange(p1, p2, p4, p3, ncol=2)
dev.off()
