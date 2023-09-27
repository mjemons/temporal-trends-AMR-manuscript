# Does temporal trend in consumption reflect temporal trend in resistance?

rm(list = ls())
library(ggplot2)
library(plyr)

# get categoriesed temporal trends
summary_fits = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")


# rename for plotting + get rid of poor fits
summary_fits_toplot <- summary_fits[summary_fits$category!='poor fit',]
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns decreasing',]$category = "decreasing (ns)"
summary_fits_toplot[summary_fits_toplot$category == 'logistic s decreasing',]$category = 'decreasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns increasing',]$category = 'increasing (ns)' 
summary_fits_toplot[summary_fits_toplot$category == 'logistic s increasing',]$category = 'increasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'flat ',]$category = 'stable'
summary_fits_toplot[summary_fits_toplot$category == 'plateauing',]$category = 'stabilising'

summary_fits_toplot$category = factor(summary_fits_toplot$category, levels = c('decreasing (s)',"decreasing (ns)",'stable','stabilising','increasing (ns)','increasing (s)', 'poor fit'))

#### Trends in consumption data ### 

# get consumption data
amc_summary = read.csv(file = "data/summary_AMC_byclass_improved.csv")

# look at correlation between year (not doing linear regression because of large variation in magnitude of DDDs by class)
amc_summary$Year = as.integer(amc_summary$Year)
trends_in_amc= ddply(amc_summary,.(Country,Antibiotic_class,Sector,Antimicrobial.Type), function(x) cor(x$DDD,x$Year,method = 'spearman'))
colnames(trends_in_amc)[5] = 'rho'

# median

amc_median <- ddply(amc_summary,.(Country,Sector,Antimicrobial.Type,Antibiotic_class), summarise, Median_AMC = median(DDD,na.rm = T), mad_AMC = mad(DDD, na.rm = T))
amc_median = amc_median[trends_in_amc$Antimicrobial.Type == 'AllNew',] 
amc_median = merge(amc_median,summary_fits_toplot,by = c("Country","Antibiotic_class"))

# trend by category is indeed what we expect, but not enough power to really say anything conclusive

trends_in_amc = trends_in_amc[trends_in_amc$Antimicrobial.Type == 'AllNew',] 
trends_in_amc = merge(trends_in_amc,summary_fits_toplot,by = c("Country","Antibiotic_class"))
ddply(trends_in_amc,.(Sector,category), function(x) c(mean(x$rho,na.rm = T), mean(x$rho,na.rm = T)-qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho))),mean(x$rho,na.rm = T)+qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho)))))

# 1. make plot nicer

data_summary <- function(x) {
  m <- median(x)
  meds = numeric(5000)
  for(i in 1:5000){
    x_t = sample(x, length(x), replace = T)
    meds[i] = median(x_t)}
  x_t = sort(meds)
  ymin = quantile(x_t, 0.975)
  names(ymin) = NULL
  ymax = quantile(x_t,0.025)
  names(ymax) = NULL
  return(c(y=m,ymin=ymin,ymax=ymax))
}

h = ggplot(trends_in_amc[trends_in_amc$Sector == 'Hospital Sector',], aes(x=category, y=rho, fill = category,alpha = 0.5)) + 
  geom_violin() +
  stat_summary(fun.data=data_summary)+
  scale_fill_manual(values = rev(c("slategray", "slategray1", "wheat", "wheat4", "pink", "pink4")))+
  theme_classic()+
  theme(legend.position = 'none', text = element_text(size = 15))+
  ylab('Trend in antibiotic consumption (correlation coefficient)')+
  xlab('Trend in resistance')


c = ggplot(trends_in_amc[trends_in_amc$Sector == 'Community',], aes(x=category, y=rho, fill = category,alpha = 0.5)) + 
  geom_violin() +
  stat_summary(fun.data=data_summary)+
  scale_fill_manual(values = rev(c("slategray", "slategray1", "wheat", "wheat4", "pink", "pink4")))+
  theme_classic()+
  theme(legend.position = 'none', text = element_text(size = 15))+
  ylab('Trend in antibiotic consumption (correlation coefficient)')+
  xlab('Trend in resistance')


# Save plots

ggsave(filename = "output/trend_in_antibiotic_consumption_by_category_hospital.pdf", h,width = 200, height = 200, units = "mm")
ggsave(filename = "output/trend_in_antibiotic_consumption_by_category_community.pdf", c,width = 200, height = 200, units = "mm")

# 2. compare increasing vs decreasing 

trends_in_amc$rough.trend = (trends_in_amc$category == "increasing (s)" | trends_in_amc$category == "increasing (ns)")*1
trends_in_amc[which(trends_in_amc$category == "decreasing (s)" | trends_in_amc$category == "decreasing (ns)"),]$rough.trend = -1

# remove stable and stabalising
trends_to_test = trends_in_amc[trends_in_amc$rough.trend != 0,]

c_stat = wilcox.test(trends_to_test[trends_to_test$Sector == 'Community'&trends_to_test$rough.trend==-1,]$rho,trends_to_test[trends_to_test$Sector == 'Community'&trends_to_test$rough.trend==1,]$rho)
h_stat = wilcox.test(trends_to_test[trends_to_test$Sector == 'Hospital Sector'&trends_to_test$rough.trend==-1,]$rho,trends_to_test[trends_to_test$Sector == 'Hospital Sector'&trends_to_test$rough.trend==1,]$rho)

print(paste('p val difference in consumption trend between increasing and decreasing resistance - community',round(c_stat$p.value,3)))
print(paste('p val difference in consumption trend between increasing and decreasing resistance - hospital',round(h_stat$p.value,3)))

# Looking at median consumption: no clear insights

ggplot(amc_median[amc_median$Sector == 'Hospital Sector',], aes(x=category, y=Median_AMC, fill = category,alpha = 0.5)) + 
  geom_violin() +
  stat_summary(fun.data=data_summary)+
  scale_fill_manual(values = rev(c("slategray", "slategray1", "wheat", "wheat4", "pink", "pink4")))+
  theme_classic()+
  theme(legend.position = 'none')+
  ylab('Antibiotic Consumption')+
  xlab('Trend in resistance')

ggplot(amc_median[amc_median$Sector == 'Community',], aes(x=category, y=Median_AMC, fill = category,alpha = 0.5)) + 
  geom_violin() +
  stat_summary(fun.data=data_summary)+
  scale_fill_manual(values = rev(c("slategray", "slategray1", "wheat", "wheat4", "pink", "pink4")))+
  theme_classic()+
  theme(legend.position = 'none')+
  ylab('Antibiotic consumption')+
  xlab('Trend in resistance')


## Compare stable vs stabilsing: no clear effect

wilcox.test(amc_median[amc_median$Sector == 'Community'&amc_median$category =='stable',]$Median_AMC,amc_median[amc_median$Sector == 'Community'&amc_median$category =="stabilising",]$Median_AMC)
wilcox.test(amc_median[amc_median$Sector == 'Hospital Sector'&amc_median$category == 'stable',]$Median_AMC,amc_median[amc_median$Sector == 'Hospital Sector'&amc_median$category=='stabilising',]$Median_AMC)

# Correct for 

data.subset = amc_median[amc_median$category == 'stabilising' | amc_median$category == 'stable',] 

summary(glm(Median_AMC ~ category + Antibiotic_class + Country + Pathogen,data = data.subset[data.subset$Sector == 'Hospital Sector',]))
summary(glm(Median_AMC ~ category + Antibiotic_class + Country + Pathogen,data = data.subset[data.subset$Sector == 'Community',]))

## Looking at increasing vs decreasing

amc_median$rough.trend = (amc_median$category == "increasing (s)" | amc_median$category == "increasing (ns)")*1
amc_median[which(amc_median$category == "decreasing (s)" | amc_median$category == "decreasing (ns)"),]$rough.trend = -1

# remove stable and stabalising
trends_to_test = amc_median[amc_median$rough.trend != 0,]

wilcox.test(trends_to_test[trends_to_test$Sector == 'Community'&trends_to_test$rough.trend==-1,]$Median_AMC,trends_to_test[trends_to_test$Sector == 'Community'&trends_to_test$rough.trend==1,]$Median_AMC)
wilcox.test(trends_to_test[trends_to_test$Sector == 'Hospital Sector'&trends_to_test$rough.trend==-1,]$Median_AMC,trends_to_test[trends_to_test$Sector == 'Hospital Sector'&trends_to_test$rough.trend==1,]$Median_AMC)

