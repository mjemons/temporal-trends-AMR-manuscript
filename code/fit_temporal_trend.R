# This script takes to temporal trend data, fits 3 possible models and picks best based on AICc

logistic_fits_wrapper = function(patient_type = 'INPAT', time_points = 5, observations = 30, minimum_N_R = 10){

#patient_type = 'INPAT'; time_points = 5; observations = 30; minimum_N_R = 10
amr_summary <- read.csv(file = "data/summary_AMR.csv")
amr_summary <- amr_summary %>% arrange(Year)
# print info on AMR files
range(amr_summary$Year)
table(amr_summary$Country)
for(pp in unique(amr_summary$Pathogen)) cat(pp, sort(unique(amr_summary$Antibiotic[amr_summary$Pathogen==pp])), "\n")

# create combinations of pathogens, countries, antibiotic:
amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Country, amr_summary$Antibiotic, sep = "|")
all_combR <- unique(amr_summary$combR)
all_combR <- sort(all_combR)

# create data frame for analysis
fit.flat = data.frame()
fit.log = data.frame()
fit.plateau = data.frame()
median.r = data.frame()

# fit flat, logistic and logistic with plateau models
for(i in 1:length(all_combR)){
  if(i/100 == floor(i/100)) print(i)
  mycomb <- all_combR[i]
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType==patient_type),]
  if((nrow(tmp1) > time_points) & (all(tmp1$N >= observations)) & (sum(tmp1$N_R + tmp1$N_I) >= minimum_N_R)){ # FB added condition on minimum number R isolates
    f1 = cbind(mycomb,fitting0(tmp1))
    fit.flat = rbind(fit.flat,f1)
    f2 = cbind(mycomb,fitting1(tmp1))
    fit.log = rbind(fit.log,f2)
    f3 = cbind(mycomb,fitting2(tmp1,f2$m1.slope,mean(tmp1$p,na.rm = T)))
    fit.plateau = rbind(fit.plateau,f3) 
    f4 = ddply(tmp1,.(combR),summarise,median.amr = median(p,na.rm = T), min.year = min(Year, na.rm = T), max.year = max(Year, na.rm = T))
    f4=cbind(f4,Country = unique(tmp1$Country), Pathogen = unique(tmp1$Pathogen), Antibiotic = unique(tmp1$Antibiotic)) #ME added country, pathogen and antibiotic for forest plots
    median.r = rbind(median.r,f4)
  }
}

# change column name that codes bug-drug-country data
colnames(fit.plateau)[1] = 'combR'
colnames(fit.flat)[1] = 'combR'
colnames(fit.log)[1] = 'combR'

colnames(median.r) = c('combR', "median.amr", "min.year", "max.year", "Country", "Pathogen", "Antibiotic")

# merge all the fits

all.fits = merge(fit.flat,fit.log,by = 'combR')
all.fits = merge(all.fits,fit.plateau,by='combR')
all.fits = merge(all.fits,median.r,by='combR')
print(median.r)
write.csv(all.fits, paste("data/results/all_fit_details", patient_type,".csv",sep=''), row.names = FALSE)

# ok, then find best fit by comparing AICc

minAICc = c()
set.seed(1234)

for(i in 1:dim(fit.flat)[1]){  
  AICcs = c(fit.flat$m0.AICc[i],fit.log$m1.AICc[i],fit.plateau$m2.AICc[i])
  AICcs[is.na(AICcs)] == Inf
  minAICc = c(minAICc,which.min(AICcs))
}

# add information about which fit is best to data frames
fit.flat$best.AICc = minAICc == 1
fit.log$best.AICc = minAICc == 2
fit.plateau$best.AICc = minAICc == 3

# For AICc: same as above but with AICc
meta.data = amr_summary[amr_summary$patientType == patient_type,c('combR', 'Country', 'Pathogen', 'Antibiotic', 'patientType', 'Antibiotic_long', 'Antibiotic_class', 'Pathogen_long')]
meta.data = unique(meta.data)
meta.data = meta.data[meta.data$combR %in% fit.flat$combR,]

# sort into same order as other data frames
meta.data = meta.data[match(fit.flat$combR,meta.data$combR),]

# get best fit
meta.data$best.model = ''
meta.data[minAICc ==1,]$best.model = 'flat'
meta.data[minAICc ==2,]$best.model = 'logistic'
meta.data[minAICc ==3,]$best.model = 'plateau'

# add info r2, error and median
meta.data$r2.best.fit = NA
meta.data[minAICc ==1,]$r2.best.fit = NA
meta.data[minAICc ==2,]$r2.best.fit = fit.log[minAICc ==2,]$m1.r2
meta.data[minAICc ==3,]$r2.best.fit = fit.plateau[minAICc ==3,]$m2.r2

meta.data$error.best.fit = NA
meta.data[minAICc ==1,]$error.best.fit = fit.flat[minAICc ==1,]$m0.mean.abs.error
meta.data[minAICc ==2,]$error.best.fit = fit.log[minAICc ==2,]$m1.mean.abs.error
meta.data[minAICc ==3,]$error.best.fit = fit.plateau[minAICc ==3,]$m2.mean.abs.error

meta.data$median.amr = median.r$median.amr

# get trend for plateau and logistic

meta.data$trend = ''

meta.data[minAICc ==2 & fit.log$m1.slope>0 & fit.log$m1.slope.pval < 0.05,]$trend = 's increasing'
meta.data[minAICc ==2 & fit.log$m1.slope>0 & fit.log$m1.slope.pval > 0.05,]$trend = 'ns increasing'
meta.data[minAICc ==2 & fit.log$m1.slope<0 & fit.log$m1.slope.pval < 0.05,]$trend = 's decreasing'
meta.data[minAICc ==2 & fit.log$m1.slope<0 & fit.log$m1.slope.pval > 0.05,]$trend = 'ns decreasing'

meta.data[minAICc ==3 & fit.plateau$m2.slope>0 & fit.plateau$m2.slope.pval < 0.05,]$trend = 's increasing'
meta.data[minAICc ==3 & fit.plateau$m2.slope>0 & fit.plateau$m2.slope.pval > 0.05,]$trend = 'ns increasing'

meta.data$category = paste(meta.data$best.model,meta.data$trend)
meta.data$full.category = paste(meta.data$best.model,meta.data$trend)
meta.data$category = meta.data$full.category
meta.data[meta.data$category == "plateau ns increasing",]$category = 'plateauing'
meta.data[meta.data$category == "plateau s increasing",]$category = 'plateauing'
meta.data[meta.data$error.best.fit > 0.05,]$category = 'poor fit'


AICc.best.fit = meta.data

write.csv(AICc.best.fit,paste("data/results/temporal_trend_AICc",patient_type,".csv",sep=''), row.names = FALSE)

}
