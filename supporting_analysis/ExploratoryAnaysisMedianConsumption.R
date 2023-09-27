### does median antibiotic consumption predict category of trend?

rm(list = ls())

# get category data:
summary_fits = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")

# rename for plotting + get rid of poor fits
summary_fits_toplot <- summary_fits[summary_fits$category!='poor fit',]
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns decreasing',]$category = "decreasing (ns)"
summary_fits_toplot[summary_fits_toplot$category == 'logistic s decreasing',]$category = 'decreasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns increasing',]$category = 'increasing (ns)' 
summary_fits_toplot[summary_fits_toplot$category == 'logistic s increasing',]$category = 'increasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'flat ',]$category = 'stable'
summary_fits_toplot[summary_fits_toplot$category == 'plateauing',]$category = 'stabilising'

# get amc data
amc_summary = read.csv(file = "data/summary_AMC_byclass_improved.csv")

# look at correlation between year (not doing linear regression because of large variation in magnitude of DDDs by class)
amc_summary$Year = as.integer(amc_summary$Year)
amc_summary = amc_summary[amc_summary$Antimicrobial.Type == 'AllNew',] 


median.amc = ddply(amc_summary,.(Country,Antibiotic_class,Sector), function(x)median(x$DDD))

colnames(median.amc)[4] = 'median.DDD'

# trend by category is indeed what we expect, but not enough power to really say anything conclusive

mean_CI = function(x)
  {c(mean(x,na.rm = T), mean(x,na.rm = T)-qnorm(0.975) * sd(x,na.rm = T)/sqrt(sum(!is.na(x))),mean(x,na.rm = T)+qnorm(0.975) * sd(x,na.rm = T)/sqrt(sum(!is.na(x))))
}

median.amc = merge(median.amc,summary_fits_toplot,by = c("Country","Antibiotic_class"))
ddply(median.amc,.(Sector,category), function(x) mean_CI(x$median.DDD))
ddply(median.amc,.(Sector,category,Pathogen), function(x) mean(x$median.DDD))

m = aov(median.DDD~category+Pathogen,data = median.amc[median.amc$Sector == 'Hospital Sector',])
summary(m)
