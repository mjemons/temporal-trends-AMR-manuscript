# Looking at the speed of increase in stabalising vs increasing models
source("code/Utils.R")

rm(list = ls())

# get data
all.fits = read.csv(file = "data/results/all_fit_detailsINPAT.csv")
summary_fits = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")

summary_fits_toplot <- summary_fits[summary_fits$category!='poor fit',]
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns decreasing',]$category = "decreasing (ns)"
summary_fits_toplot[summary_fits_toplot$category == 'logistic s decreasing',]$category = 'decreasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns increasing',]$category = 'increasing (ns)' 
summary_fits_toplot[summary_fits_toplot$category == 'logistic s increasing',]$category = 'increasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'flat ',]$category = 'stable'
summary_fits_toplot[summary_fits_toplot$category == 'plateauing',]$category = 'stabilising'

summary_fits_toplot$category = factor(summary_fits_toplot$category, levels = c("decreasing (ns)",'decreasing (s)','stabilising','stable','increasing (s)','increasing (ns)'))

all.fits = merge(all.fits,summary_fits_toplot[,c('combR','category','Pathogen','Antibiotic_class','Country')],by = c("combR","Pathogen","Country"))

# to make models comparable, calculate maximal predicted rate of change in observed range

# speed of increase in year x: (b1 E^(b0 - b1 x))/(1 + E^(b0 - b1 x))^2

max.speed.m1 = function(x){
  # max speed in standard logistic model
  b0 = x$m1.intercept
  b1 = x$m1.slope
  year.range = c(x$min.year:x$max.year)
  increases = b1 * exp(b0 - b1 * year.range)/(1+exp(b0 - b1 * year.range))^2
  max(increases)
}

min.speed.m1 = function(x){
  # min speed in standard logistic model
  b0 = x$m1.intercept
  b1 = x$m1.slope
  year.range = c(x$min.year:x$max.year)
  increases = b1 * exp(b0 - b1 * year.range)/(1+exp(b0 - b1 * year.range))^2
  min(increases)
}

max.speed.m2 = function(x){
  # max speed in plateauing logistic model
  b0 = x$m2.intercept
  b1 = x$m2.slope
  k = x$m2.plateau
  year.range = c(x$min.year:x$max.year)
  increases = k*b1 * exp(b0 - b1 * year.range)/(1+exp(b0 - b1 * year.range))^2
  max(increases)
}

### calculate relevant measure

# increasing
s.fits = all.fits[all.fits$category == 'increasing (s)',]
speed.increasing = ddply(s.fits,.(combR), max.speed.m1)
s.fits = merge(s.fits,speed.increasing)

# decreasing
d.fits = all.fits[all.fits$category == 'decreasing (s)',]
speed.decreasing = ddply(d.fits,.(combR), min.speed.m1)
d.fits = merge(d.fits,speed.decreasing)

# stabilising
p.fits = all.fits[all.fits$category == 'stabilising',]
speed.plateau = ddply(p.fits,.(combR), max.speed.m2)
p.fits = merge(p.fits,speed.plateau)

# compare increasing vs stablising
fits = rbind(s.fits,p.fits)

print(ddply(fits,.(category),summarise,round(mean(V1),3)))

m = glm(V1 ~ category,data = fits)
temp = summary(m)
print(paste('p val effect stabalising vs increasing (uncorrected)',round(temp$coefficients[2,4],7)))

m = glm(V1 ~ category + Pathogen + Antibiotic+ Country,data = fits)
temp = summary(m)
print(paste('p val effect stabalising vs increasing (corrected)',round(temp$coefficients[2,4],7)))


ggplot(fits, aes(x=V1, fill=category)) + geom_histogram(alpha=0.7, position="identity", binwidth = 0.01) + 
  theme_classic()+
  scale_fill_manual(values=c("wheat4", "slategray1"))+
  labs(x = 'maximum rate of change')

ggsave("output/speed_of_change.pdf", width=7.5, height=5)

## testing speed vs consumption

# get consumption data and compute overall 2000 consumption

amc_summary <- read.csv("data/summary_AMC_byclass_improved.csv")
subset = amc_summary[amc_summary$Antimicrobial.Type == 'AllNew',]
total.consumption = ddply(subset,.(Country,Year,Sector),summarise,DDD = sum(DDD))

# format data for analysis
mean.speed = ddply(p.fits,.(Country),summarise,mean.speed = mean(V1))
to.test = merge(total.consumption[total.consumption$Year == 2000,],mean.speed, by = 'Country')

correlation = cor.test(to.test[to.test$Sector=='Community',]$DDD,to.test[to.test$Sector=='Community',]$mean.speed,method = 'spearman')

print(paste('Correlation between consumption 2000 and speed of initial increase rho:',round(correlation$estimate,2)))
print(paste('Correlation between consumption 2000 and speed of initial increase p:',round(correlation$p.value,4)))
