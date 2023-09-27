# Notes: here, we are not matching between the AMR years and AMC years.
# And there is no threshold for how much AMC data needed for inclusion.

# load AMC and AMR and fits
amc_summary <- read.csv(file = "data/summary_AMC_byclass_filtered.csv")
amr_fit_details <- read.csv(file = "data/results/all_fit_detailsINPAT.csv")
cat.summary = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")

# calculate median consumption per country, sector, type, class
amc_median = ddply(amc_summary,.(Country,Sector,Antimicrobial.Type,Antibiotic_class), summarise, Median_AMC = median(DDD,na.rm = T))

###### functions for calculation ####
calculate_correlations_plateauing = function(antimicrobial.type, amc_median, amr_analyse){
  
  amc_analyse = amc_median[amc_median$Antimicrobial.Type == antimicrobial.type, ]
  to_analyse = merge(amc_analyse, amr_analyse, by = c('Country','Antibiotic_class'))
  testing1 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$m2.plateau,use = "complete.obs",conf.level = 0.95))
  testing1$correlate = 'Plateau.level'
  testing2 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$m2.slope,use = "complete.obs",conf.level = 0.95))
  testing2$correlate = 'Slope'
  testing3 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$median.amr,use = "complete.obs",conf.level = 0.95))
  testing3$correlate = "Median.amr"
  
  # join results into single data frame
  correlations = rbind(rbind(testing1,testing2),testing3)
  return(correlations)
}

calculate_correlations_increasing = function(antimicrobial.type,amc_median,amr_analyse){
  amc_analyse = amc_median[amc_median$Antimicrobial.Type == antimicrobial.type,]
  to_analyse = merge(amc_analyse,amr_analyse,by = c('Country','Antibiotic_class'))
  
  testing2 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$m1.slope,use = "complete.obs",conf.level = 0.95))
  testing2$correlate = 'Slope'
  testing3 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$median.amr,use = "complete.obs",conf.level = 0.95))
  testing3$correlate = "Median.amr"
  
  # join results into single data frame
  correlations = rbind(testing2,testing3)
  return(correlations)
}

calculate_correlations_flat = function(antimicrobial.type,amc_median,amr_analyse){
  
  amc_analyse = amc_median[amc_median$Antimicrobial.Type == antimicrobial.type,]
  to_analyse = merge(amc_analyse,amr_analyse,by = c('Country','Antibiotic_class'))
  testing1 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$m0.plateau,use = "complete.obs",conf.level = 0.95))
  testing1$correlate = 'Plateau.level'
  testing3 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$median.amr,use = "complete.obs",conf.level = 0.95))
  testing3$correlate = "Median.amr"
  
  # join results into single data frame
  correlations = rbind(testing1,testing3)
  return(correlations)
}

#### Calculating the correlations #####

amc.types = c('Oral','Parenteral','All')

amr_analyse = cat.summary[which(cat.summary$category == 'plateauing'),]
amr_analyse = merge(amr_analyse,amr_fit_details[,c('combR',"m2.plateau","m2.slope")])

correlations.plateau = data.frame()
for(i in 1:3){
  temp= calculate_correlations_plateauing(amc.types[i],amc_median,amr_analyse)
  temp$amc.type = amc.types[i]
  correlations.plateau = rbind(correlations.plateau,temp)
}

amr_analyse = cat.summary[which(cat.summary$category == 'logistic s increasing'),]
amr_analyse = merge(amr_analyse, amr_fit_details[,c('combR',"m1.slope")])

correlations.increasing = data.frame()
for(i in 1:3){
  temp= calculate_correlations_increasing(amc.types[i],amc_median,amr_analyse)
  temp$amc.type = amc.types[i]
  correlations.increasing = rbind(correlations.increasing,temp)
}


amr_analyse = cat.summary[which(cat.summary$category == 'flat '),]
amr_analyse = merge(amr_analyse,amr_fit_details[,c('combR',"m0.plateau")])

correlations.flat = data.frame()
for(i in 1:3){
  temp= calculate_correlations_flat(amc.types[i],amc_median,amr_analyse)
  temp$amc.type = amc.types[i]
  correlations.flat = rbind(correlations.flat,temp)
}

correlations.flat$trend = 'stable'
correlations.plateau$trend = 'stabalising'
correlations.increasing$trend = 'increasing (s)'

all.correlations = rbind(rbind(correlations.flat,correlations.plateau),correlations.increasing)

write.csv(all.correlations,paste("data/results/all_correlations_INPAT",".csv",sep=''), row.names = FALSE)


