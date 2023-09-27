# Check the impact of filtering on correlation calculations: calculate correlations using unfiltered version

# TAKE-AWAY - the filtering has minimal impact on the results: using median is already pretty robust.

# FILTERED
amc_summary <- read.csv(file = "data/summary_AMC_byclass_filtered.csv")
amr_logistic_fits <- read.csv(file = "data/results/temporal_trend_INPAT.csv")
# resistance: only logistic fits
amr_analyse = amr_logistic_fits[which(!is.na(amr_logistic_fits$logistic.plateau.K)),]

# calculate median consumption
amc_median = ddply(amc_summary,.(Country,Sector,Antimicrobial.Type,Antibiotic_class),summarise,Median_AMC = median(DDD,na.rm = T))

# ORAL consumption: # caluclation correlation with plateau, slope and median 

antimicrobial.type = 'Oral'
amc_analyse = amc_median[amc_median$Antimicrobial.Type == antimicrobial.type,]
to_analyse = merge(amc_analyse,amr_analyse,by = c('Country','Antibiotic_class'))

testing1 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$logistic.plateau.K,use = "complete.obs",conf.level = 0.95))
testing1$correlate = 'Plateau.level'
testing2 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$slope_plateau,use = "complete.obs",conf.level = 0.95))
testing2$correlate = 'Slope'
testing3 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$Median_AMR,use = "complete.obs",conf.level = 0.95))
testing3$correlate = "Median.amr"

# join results into single data frame
correlations.filtered = rbind(rbind(testing1,testing2),testing3)

# UNFILTERED

# load AMC and AMR
amc_summary <- read.csv(file = "data/summary_AMC_byclass.csv")
amr_logistic_fits <- read.csv(file = "data/results/temporal_trend_INPAT.csv")
# resistance: only logistic fits
amr_analyse = amr_logistic_fits[which(!is.na(amr_logistic_fits$logistic.plateau.K)),]

# calculate median consumption
amc_median = ddply(amc_summary,.(Country,Sector,Antimicrobial.Type,Antibiotic_class),summarise,Median_AMC = median(DDD,na.rm = T))

# ORAL consumption: # caluclation correlation with plateau, slope and median 

antimicrobial.type = 'Oral'
amc_analyse = amc_median[amc_median$Antimicrobial.Type == antimicrobial.type,]
to_analyse = merge(amc_analyse,amr_analyse,by = c('Country','Antibiotic_class'))

testing1 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$logistic.plateau.K,use = "complete.obs",conf.level = 0.95))
testing1$correlate = 'Plateau.level'
testing2 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$slope_plateau,use = "complete.obs",conf.level = 0.95))
testing2$correlate = 'Slope'
testing3 = ddply(to_analyse,.(Pathogen,Antibiotic,Antibiotic_class,Sector),function(x) if(dim(x)[1]>4) SpearmanRho(x$Median_AMC,x$Median_AMR,use = "complete.obs",conf.level = 0.95))
testing3$correlate = "Median.amr"

# join results into single data frame
correlations.unfiltered = rbind(rbind(testing1,testing2),testing3)

