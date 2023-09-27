# Exploratory analysis of temporal AMC/AMR correlations.

# Summarise:
# 

rm(list = ls())
library(ggplot2)

# get categoriesed temporal trends
summary_fits = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")

# get the correlations
load(file = "data/AMR_AMC_correlation.RData")

# rename for plotting + get rid of poor fits
summary_fits_toplot <- summary_fits[summary_fits$category!='poor fit',]
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns decreasing',]$category = "decreasing (ns)"
summary_fits_toplot[summary_fits_toplot$category == 'logistic s decreasing',]$category = 'decreasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns increasing',]$category = 'increasing (ns)' 
summary_fits_toplot[summary_fits_toplot$category == 'logistic s increasing',]$category = 'increasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'flat ',]$category = 'stable'
summary_fits_toplot[summary_fits_toplot$category == 'plateauing',]$category = 'stabilising'


# combine two datasets.

temporal_cor_hospital = merge(rho_temporal_HS_AllNew_INPAT_present,summary_fits_toplot, by = c('Country','Pathogen','Antibiotic'))
temporal_cor_community = merge(rho_temporal_Community_AllNew_INPAT_present,summary_fits_toplot, by = c('Country','Pathogen','Antibiotic'))

temporal_cor_hospital_by_trend = ddply(temporal_cor_hospital,.(category),function(x) round(c(mean(x$rho),mean(x$rho)-qnorm(0.975) * sd(x$rho)/sqrt(length(x$rho)),mean(x$rho)+qnorm(0.975) * sd(x$rho)/sqrt(length(x$rho))),2))
temporal_cor_community_by_trend = ddply(temporal_cor_community,.(category),function(x) round(c(mean(x$rho),mean(x$rho)-qnorm(0.975) * sd(x$rho)/sqrt(length(x$rho)),mean(x$rho)+qnorm(0.975) * sd(x$rho)/sqrt(length(x$rho))),2))

colnames(temporal_cor_hospital_by_trend) = c('category','mean.rho','lower.CI','upper.CI')
colnames(temporal_cor_community_by_trend) = c('category','mean.rho','lower.CI','upper.CI')

# not in itself significant - but if we look at sustained trends vs others
m = aov(rho~category,data = temporal_cor_hospital)
summary(m)

# not in itself significant - but if we look at sustained trends vs others
m = aov(rho~category,data = temporal_cor_hospital)
summary(m)

temporal_cor_hospital$sig.trend = temporal_cor_hospital$category == "increasing (s)" | temporal_cor_hospital$category == "decreasing (s)" 

m = aov(rho~sig.trend,data = temporal_cor_hospital)
summary(m)

ggplot(temporal_cor_hospital, aes(x = rho)) +
  geom_histogram(fill = "white", colour = "black") +
  facet_grid(category ~ ., scales = "free")


# hypothesis: correlations reflect concordant trends in antibiotic use and resistance, rather than year-by-year variation.

#### Trends in consumption data ### 

# get consumption data
amc_summary = read.csv(file = "data/summary_AMC_byclass_improved.csv")

# look at correlation between year (not doing linear regression because of large variation in magnitude of DDDs by class)
amc_summary$Year = as.integer(amc_summary$Year)
trends_in_amc= ddply(amc_summary,.(Country,Antibiotic_class,Sector,Antimicrobial.Type), function(x) cor(x$DDD,x$Year,method = 'spearman'))
colnames(trends_in_amc)[5] = 'rho'

# trend by category is indeed what we expect, but not enough power to really say anything conclusive

trends_in_amc = trends_in_amc[trends_in_amc$Antimicrobial.Type == 'AllNew',] 
trends_in_amc = merge(trends_in_amc,summary_fits_toplot,by = c("Country","Antibiotic_class"))
ddply(trends_in_amc,.(Sector,category), function(x) c(mean(x$rho,na.rm = T), mean(x$rho,na.rm = T)-qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho))),mean(x$rho,na.rm = T)+qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho)))))

ggplot(trends_in_amc[trends_in_amc$Sector == 'Hospital Sector',], aes(x=category, y=rho)) + 
  geom_violin() +
  stat_summary(fun.y=median, geom="point", size=2, color="red")

# increase power by lumping concerning together.
trends_in_amc$concerning = trends_in_amc$category == "increasing (s)" | trends_in_amc$category == "increasing (ns)" 
trends_in_amc$rough.trend = (trends_in_amc$category == "increasing (s)" | trends_in_amc$category == "increasing (ns)")*1
trends_in_amc[which(trends_in_amc$category == "decreasing (s)" | trends_in_amc$category == "decreasing (ns)"),]$rough.trend = -1

ddply(trends_in_amc,.(Sector,concerning), function(x) c(mean(x$rho,na.rm = T), mean(x$rho,na.rm = T)-qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho))),mean(x$rho,na.rm = T)+qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho)))))
ddply(trends_in_amc,.(Sector,rough.trend), function(x) c(mean(x$rho,na.rm = T), mean(x$rho,na.rm = T)-qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho))),mean(x$rho,na.rm = T)+qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho)))))

m = aov(rho~rough.trend,data = trends_in_amc[trends_in_amc$Sector == 'Hospital Sector',])
summary(m)

# how about if just how increasing/decreasing? Also not massive difference
ddply(trends_in_amc,.(Sector,category), function(x) mean(x$rho>0,na.rm = T))



      