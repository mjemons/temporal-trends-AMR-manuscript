rm(list = ls())
##########                                                                CORRELATE AMR AND AMC                                                                ##########

library(minpack.lm)
library(dplyr)

source("code/UsefulFunctions.R")

# get Francois' results
df_coefficients <- read.csv(file = "data/smooth_splines_coefficients.csv")


# it is possible to start here and load AMR and AMC files
# TODO DEAL WITH ANTIMICROBIAL.TYPE
amc_summary <- read.csv(file = "data/2018/csv/summary_AMC_byclass.csv")
amr_summary <- read.csv(file = "data/2018/csv/summary_AMR.csv")
amc_tot <- read.csv(file = "data/2018/csv/total_AMC_byclass.csv")
amr_tot <- read.csv(file = "data/2018/csv/total_AMR.csv")

#amr_summary_ <- amr_summary[which(amr_summary$DateUsedForStatisticsYear < 2020), ]

# slight data cleaning
colnames(amr_summary)[colnames(amr_summary)=="DateUsedForStatisticsYear"] <- "Year"
amr_summary$patientType[amr_summary$patientType=="NULL"] <- "UNK"
amr_summary$patientType[amr_summary$patientType=="O"] <- "OUTPAT" #  Patients that go to the hospital for Dialysis, other Day Hospital Care and to Emergency room should be classified as O
colnames(amr_summary)[colnames(amr_summary)=="ReportingCountry"] <- "Country"
amr_summary <- amr_summary[!is.na(amr_summary$Pathogen),]
amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Antibiotic, amr_summary$patientType, amr_summary$Year, sep = "|")
amc_summary$combC <- paste(amc_summary$class, amc_summary$Sector, amc_summary$Year, sep = "|")

# change country name
amr_summary$Country <- as.character(amr_summary$Country)
amr_summary$Country[amr_summary$Country=="AT"] <- "Austria"
amr_summary$Country[amr_summary$Country=="BE"] <- "Belgium"
amr_summary$Country[amr_summary$Country=="BG"] <- "Bulgaria"
amr_summary$Country[amr_summary$Country=="CY"] <- "Cyprus"
amr_summary$Country[amr_summary$Country=="CZ"] <- "Czech Republic"
amr_summary$Country[amr_summary$Country=="DE"] <- "Germany"
amr_summary$Country[amr_summary$Country=="DK"] <- "Denmark"
amr_summary$Country[amr_summary$Country=="EE"] <- "Estonia"
amr_summary$Country[amr_summary$Country=="EL"] <- "Greece"
amr_summary$Country[amr_summary$Country=="ES"] <- "Spain"
amr_summary$Country[amr_summary$Country=="FI"] <- "Finland"
amr_summary$Country[amr_summary$Country=="FR"] <- "France"
amr_summary$Country[amr_summary$Country=="HR"] <- "Croatia"
amr_summary$Country[amr_summary$Country=="HU"] <- "Hungary"
amr_summary$Country[amr_summary$Country=="IE"] <- "Ireland"
amr_summary$Country[amr_summary$Country=="IT"] <- "Italy"
amr_summary$Country[amr_summary$Country=="LT"] <- "Lithuania"
amr_summary$Country[amr_summary$Country=="LU"] <- "Luxembourg"
amr_summary$Country[amr_summary$Country=="LV"] <- "Latvia"
amr_summary$Country[amr_summary$Country=="MT"] <- "Malta"
amr_summary$Country[amr_summary$Country=="NL"] <- "Netherlands"
amr_summary$Country[amr_summary$Country=="NO"] <- "Norway"
amr_summary$Country[amr_summary$Country=="PL"] <- "Poland"
amr_summary$Country[amr_summary$Country=="PT"] <- "Portugal"
amr_summary$Country[amr_summary$Country=="RO"] <- "Romania"
amr_summary$Country[amr_summary$Country=="SE"] <- "Sweden"
amr_summary$Country[amr_summary$Country=="SI"] <- "Slovenia"
amr_summary$Country[amr_summary$Country=="SK"] <- "Slovakia"
amr_summary$Country[amr_summary$Country=="UK"] <- "United Kingdom"

# create combinations of pathogens, countries, antibiotic:
amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Country, amr_summary$Antibiotic, sep = "|")
all_combR <- unique(amr_summary$combR)
#all_combR <- all_combR[grepl(pattern = "ESCCOL", x = all_combR)]
all_combR <- sort(all_combR)

log1_df <- read.csv('data/first_log_fit.csv')
log2_df <- read.csv('data/second_log_fit.csv')

###### Analysis starts here ###### 

tokeep1 = c() # inpatients for inclusion
tokeep2 = c() # outpatients for inclusion
failed = c() # to keep the failed ones
Slope = c() # to store the initial slopes
for(i in 1:length(all_combR)){
  mycomb <- all_combR[i]
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  tmp2 <- tmp[which(tmp$patientType=="OUTPAT"),]
  #why do we stop ? like this we loose the entire combination instead of just
  #having a missing data point. Alternatively we could just take one of the two
  if(any(duplicated(tmp1$Year))) stop()
  if(any(duplicated(tmp2$Year))) stop()
  # inpatients 
  if((nrow(tmp1) > 5) & (all(tmp1$N > 30))){ # at least 5 years & all points > 100
    #fit = weighted.mean(tmp1$p, tmp1$N)
    #can't calculate AIC of a mean, thus give a linear model with slope zero 
    fit = lm(I(tmp1$p-0*tmp1$Year)~1, weights = tmp1$N)
    keep_temp = c(i,coef(fit), AIC(fit))
    if(is.na(keep_temp[1])){
      failed = rbind(failed,mycomb)
    }
    else{
      Slope = c(Slope,keep_temp[1])
      tokeep1 = rbind(tokeep1,c(i,keep_temp))
    }
  }
  # outpatients 
  if((nrow(tmp2) > 5) & (all(tmp2$N > 30))){ # at least 5 years for one of them
    fit = lm(I(tmp2$p-0*tmp2$Year)~1, weights = tmp2$N)
    keep_temp = c(fit, AIC(fit))
    #this doesn't work, since the next above already breaks the loop, such that not all OUTPAT are counted that fail
    if(is.na(keep_temp[1])){
      failed = rbind(failed,mycomb)
    }
    else{
      tokeep2 = rbind(tokeep2,c(i,keep_temp))
      #if(is.na(keep_temp[2])) failed = rbind(failed,mycomb)
    }
  }
}

#find out which are not present in the logistic fit dataframe
diff_1 <- as.vector(setdiff(tokeep1[,1],log1_df[,1]))
#reduce to those that were successfully fit logistically
tokeep1_reduced_1 <- tokeep1 %>% as.data.frame() %>% filter(!tokeep1[,1] %in% diff_1)

write.csv(tokeep1_reduced_1,"data/first_model_lm.csv", row.names = FALSE)
#logistic_fit = all_combR[tokeep1_reduced[which(tokeep1_reduced[,4]>log1_df[,5])],1]

diff_2 <- as.vector(setdiff(tokeep1[,1],log2_df[,1]))
tokeep1_reduced_2 <- tokeep1 %>% as.data.frame() %>% filter(!tokeep1[,1] %in% diff_2)
write.csv(tokeep1_reduced_2,"data/second_model_lm.csv", row.names = FALSE)