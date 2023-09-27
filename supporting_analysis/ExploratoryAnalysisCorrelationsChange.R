# Exploratory analysis looking at temporal correlation between AMC and *change* in resistance frequency

# This is bascially all noise, no correlations here. 

rm(list = ls())

# get amc data and amr data
amc_summary = read.csv(file = "data/summary_AMC_byclass_improved.csv")
amr_summary <- read.csv("data/summary_AMR.csv")
amr_summary$combC <- paste(amr_summary$Antibiotic_class, amr_summary$Country, amr_summary$patientType,amr_summary$Pathogen, sep = "|")

inoutpat <- "INPAT"

for(mysector in c("Community", "Hospital Sector")){
  for(mytype in c("Oral", "AllNew")){
      
      tmp_amc <- amc_summary[amc_summary$Antimicrobial.Type == mytype & amc_summary$Sector == mysector, ]
      tmp_amr <- amr_summary[amr_summary$patientType==inoutpat, ]

      # compute yearly change (i.e. next year - this year)
      
      # sort by i) year and ii) combc_2
      tmp_amr = tmp_amr[order( tmp_amr[,16], tmp_amr[,1] ),]
      tmp_change = diff(tmp_amr$p)
      
      # remove the ones where different country/class/species
      to_keep = c()
      for(i in 1:length(tmp_change)){
        to_keep = c(to_keep,tmp_amr$combC[i] == tmp_amr$combC[i+1])
      }
      
      tmp_change = tmp_change[to_keep]
      tmp_amr = tmp_amr[c(to_keep,FALSE),]
      tmp_amr = cbind(tmp_amr,tmp_change)
      
      tmp_merged <- merge(tmp_amc, tmp_amr, by = c('Country','Antibiotic_class','Year'), all = F)
      
      #tmp_merged[which(tmp_merged$Country=="France" & tmp_merged$Pathogen=="ESCCOL" & tmp_merged$Antibiotic=="AMP"), ]
      
      rho_temporal <- ddply(tmp_merged, .(Country, Pathogen, Antibiotic),    function(x) c(SpearmanRho(x$DDD, x$tmp_change, use = "complete.obs",conf.level = 0.95), Nyears = nrow(x)))
      rho_temporal <- rho_temporal[rho_temporal$Nyears>4, ]
      if(mysector=="Hospital Sector") mysector2 <- "HS" else mysector2 <- mysector
      assign(x = paste0("rho_temporal_", mysector2, "_", mytype, "_", inoutpat,'_change'), value = rho_temporal)
      
      print(mysector)
      print(mytype)
      print(t.test(rho_temporal$rho))
    }
}
rm(rho_temporal)

rho_temporal_names <- ls()[grepl(pattern = "rho_temporal", ls())] # temporal correlations objects
objects_to_save <- rho_temporal_names
save(list = objects_to_save, file = "data/AMR_AMC_correlation_change.RData")

# compare to correlation between amc and amr frequency (i.e. not change)

load(file = "data/AMR_AMC_correlation.RData")

# Compare with previous AMC vs frequency correlation.

to_compare = merge(rho_temporal_HS_AllNew_INPAT_present,rho_temporal_HS_AllNew_INPAT_change,by = c('Country','Pathogen','Antibiotic'))
t.test(to_compare$rho.x,to_compare$rho.y,paired = TRUE)
cor.test(to_compare$rho.x,to_compare$rho.y)

to_compare = merge(rho_temporal_HS_Oral_INPAT_present,rho_temporal_HS_Oral_INPAT_change,by = c('Country','Pathogen','Antibiotic'))
t.test(to_compare$rho.x,to_compare$rho.y,paired = TRUE)
cor.test(to_compare$rho.x,to_compare$rho.y)

to_compare = merge(rho_temporal_Community_AllNew_INPAT_present,rho_temporal_Community_AllNew_INPAT_change,by = c('Country','Pathogen','Antibiotic'))
t.test(to_compare$rho.x,to_compare$rho.y,paired = TRUE)
cor.test(to_compare$rho.x,to_compare$rho.y)

to_compare = merge(rho_temporal_Community_Oral_INPAT_present,rho_temporal_Community_Oral_INPAT_change,by = c('Country','Pathogen','Antibiotic'))
t.test(to_compare$rho.x,to_compare$rho.y,paired = TRUE)
cor.test(to_compare$rho.x,to_compare$rho.y)

## Look by category

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

temporal_cor_hospital = merge(rho_temporal_HS_AllNew_INPAT_change,summary_fits_toplot, by = c('Country','Pathogen','Antibiotic'))
temporal_cor_community = merge(rho_temporal_Community_AllNew_INPAT_change,summary_fits_toplot, by = c('Country','Pathogen','Antibiotic'))

temporal_cor_hospital_by_trend = ddply(temporal_cor_hospital,.(category),function(x) c(mean(x$rho,na.rm = T), mean(x$rho,na.rm = T)-qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho))),mean(x$rho,na.rm = T)+qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho)))))
temporal_cor_community_by_trend = ddply(temporal_cor_community,.(category),function(x) c(mean(x$rho,na.rm = T), mean(x$rho,na.rm = T)-qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho))),mean(x$rho,na.rm = T)+qnorm(0.975) * sd(x$rho,na.rm = T)/sqrt(sum(!is.na(x$rho)))))

colnames(temporal_cor_hospital_by_trend) = c('category','mean.rho','lower.CI','upper.CI')
colnames(temporal_cor_community_by_trend) = c('category','mean.rho','lower.CI','upper.CI')

# look more specifically at distributions

ggplot(temporal_cor_hospital, aes(x = rho)) +
  geom_histogram(fill = "white", colour = "black") +
  facet_grid(category ~ ., scales = "free")

                   
