# Looking at the speed of increase in stabalising vs increasing models

# slope
rm(list = ls())

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

all.fits = merge(all.fits,summary_fits_toplot[,c('combR','category')],by = "combR")

# speed of increase in year x: (b1 E^(b0 - b1 x))/(1 + E^(b0 - b1 x))^2

max.speed.m1 = function(x){
  b0 = x$m1.intercept
  b1 = x$m1.slope
  year.range = c(x$min.year:x$max.year)
  increases = b1 * exp(b0 - b1 * year.range)/(1+exp(b0 - b1 * year.range))^2
  max(increases)
}

min.speed.m1 = function(x){
  b0 = x$m1.intercept
  b1 = x$m1.slope
  year.range = c(x$min.year:x$max.year)
  increases = b1 * exp(b0 - b1 * year.range)/(1+exp(b0 - b1 * year.range))^2
  min(increases)
}

max.speed.m2 = function(x){
  b0 = x$m2.intercept
  b1 = x$m2.slope
  k = x$m2.plateau
  year.range = c(x$min.year:x$max.year)
  increases = k*b1 * exp(b0 - b1 * year.range)/(1+exp(b0 - b1 * year.range))^2
  max(increases)
}

s.fits = all.fits[all.fits$category == 'increasing (s)',]
speed.increasing = ddply(s.fits,.(combR), max.speed.m1)
s.fits = merge(s.fits,speed.increasing)

d.fits = all.fits[all.fits$category == 'decreasing (s)',]
speed.decreasing = ddply(d.fits,.(combR), min.speed.m1)
d.fits = merge(d.fits,speed.decreasing)

p.fits = all.fits[all.fits$category == 'stabilising',]
speed.plateau = ddply(p.fits,.(combR), max.speed.m2)
p.fits = merge(p.fits,speed.plateau)



fits = rbind(s.fits,p.fits)

m = glm(V1 ~ category + Pathogen + Antibiotic+ Country,data = fits)
summary(m)

m = glm(V1 ~ category,data = fits)
summary(m)

fits.temp = fits
fits.temp = fits.temp[fits.temp$V1 < 0.4,]

ggplot(fits.temp, aes(x=V1, fill=category)) + geom_histogram(alpha=0.2, position="identity")


# Similar analysis looking at maximal predicted rate of change (i.e. not necessarily observed)- but this does not reflect actual observed range.

s.fits = all.fits[all.fits$category == 'increasing (s)',]
s.fits$max.change = s.fits$m1.slope/4

p.fits = all.fits[all.fits$category == 'stabilising',]
p.fits$max.change = p.fits$m2.slope/4 * p.fits$m2.plateau

fits = rbind(s.fits,p.fits)

m = glm(max.change ~ category + Pathogen + Antibiotic_class + Country,data = fits)
summary(m)

m = glm(max.change ~ category,data = fits)
summary(m)

