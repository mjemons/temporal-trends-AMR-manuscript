# This script takes to temporal trend data and:
# 1) fits a starndard logistic regression
# 2) logistic regression with a plateau < 1
# 3) categorises temporal trend based on AIC of these fits.

# Q: do you need to save the fits? If not, can get rid.

                                                        ##########

logistic_fits_wrapper = function(patient_type = 'INPAT',time_points = 5,observations= 30){

#### Set up #####

# get the data
amr_summary <- read.csv(file = "data/summary_AMR.csv")
amr_summary <- amr_summary %>% arrange(Year)

# create combinations of pathogens, countries, antibiotic:
amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Country, amr_summary$Antibiotic, sep = "|")
all_combR <- unique(amr_summary$combR)
all_combR <- sort(all_combR)

# create data frame for analysis

df_coefficients = as.data.frame(all_combR)
colnames(df_coefficients)[1] = 'combR'

###### Fitting ###### 

# Fit standard logistic model.
# Only include countries with at least 5 years of data with 30 observations

tokeep1 = c() # inpatients for inclusion
failed = c() # to keep the failed ones
Slope = c() # to store slope parameter
Slope_stderr = c() # to store the std. error of the initial slopes

counter = 0
for(i in 1:length(all_combR)){
  mycomb <- all_combR[i]
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType==patient_type),]
  # inpatients 
  if((nrow(tmp1) > time_points) & (all(tmp1$N > observations))){
    keep_temp = fitting1(tmp1)
    if(is.na(keep_temp[1])){
      counter = counter + 1
      failed = rbind(failed,mycomb)
    }
    if(!is.na(keep_temp[1])){
      Slope = c(Slope,keep_temp[1])
      Slope_stderr = c(Slope_stderr,keep_temp[5])
      tokeep1 = rbind(tokeep1,c(i,keep_temp))
    }
  }
}


trend.increasing =which(tokeep1[,2]>0)
trend.decreasing =which(tokeep1[,2]<0)
stat.increasing = which(tokeep1[,3]>0)
stat.decreasing = which(tokeep1[,4]<0)

all.slopes = which(!is.na(tokeep1[,2]))

stat.unclear = which(tokeep1[,4]>0&tokeep1[,3]<0)

# If regression slope is positive (resistance increasing): fit regression with plateau

increasing = trend.increasing
testing = c()
testing_2 = c()
carrying_K = c()
failed2 = c()
Slope_plateau = c() 
Slope_plateau_stderr = c()
plateau.pseudo.intercept = c()
for(i in 1:length(increasing)){
  x = increasing[i]
  mycomb <- all_combR[tokeep1[x,1]]
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  starval = tokeep1[x,2] 
  temp = fitting2(tmp1,starval)
  aicval = temp[1]
  carrying_K = c(carrying_K, temp[4])
  Slope_plateau = c(Slope_plateau, temp[3])
  Slope_plateau_stderr = c(Slope_plateau_stderr, temp[5])
  plateau.pseudo.intercept = c(plateau.pseudo.intercept,temp[2])
  testing = c(testing,aicval)
  # store index, aic and r^2
  testing_2 = rbind(testing_2, c(tokeep1[x,1], aicval,temp[7],temp[8]))
  if(is.na(aicval))  failed2 = rbind(failed2,mycomb)
}

#choose the model with lowest AIC
plateau = all_combR[tokeep1[increasing[which(tokeep1[increasing,5]>testing)],1]]
index_plateau = tokeep1[increasing[which(tokeep1[increasing,5]>testing)],1]
no.plateau = all_combR[tokeep1[increasing[which(tokeep1[increasing,5]<testing)],1]]
index_no.plateau = tokeep1[increasing[which(tokeep1[increasing,5]<testing)],1]
unclear = all_combR[tokeep1[increasing[which(is.na(testing))],1]]
index_unclear = tokeep1[increasing[which(is.na(testing))],1]

nonsig.decreasing <- setdiff(trend.decreasing,stat.decreasing)

index_nonsig.decrease = tokeep1[nonsig.decreasing,1]
index_sig.decrease = tokeep1[stat.decreasing,1]

# categorise all by trend:
df_coefficients$logistic.trend = 'NA'
df_coefficients$slope = 'NA'
df_coefficients$slope_stderr = 'NA'
df_coefficients$Patient.type = patient_type
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[all.slopes,1]],]$slope = Slope
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[all.slopes,1]],]$slope_stderr = Slope_stderr
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$logistic.trend = 'N.S. increase'
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[trend.decreasing,1]],]$logistic.trend = 'N.S. decrease'
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[stat.decreasing,1]],]$logistic.trend = 'S. decrease'
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[stat.increasing,1]],]$logistic.trend = 'S. increase'

df_coefficients$logistic.plateau = 'NA'
df_coefficients$logistic.plateau.K = 'NA'
df_coefficients$plateau.pseudo.intercept = 'NA'
df_coefficients$slope_plateau = 'NA'
df_coefficients$slope_plateau_stderr = 'NA'
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$slope_plateau = Slope_plateau
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$slope_plateau_stderr = Slope_plateau_stderr
df_coefficients[!df_coefficients$combR %in% plateau,]$slope_plateau = 'NA'
df_coefficients[!df_coefficients$combR %in% plateau,]$slope_plateau_stderr = 'NA'
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$logistic.plateau.K = carrying_K
df_coefficients[!df_coefficients$combR %in% plateau,]$logistic.plateau.K = 'NA'
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$plateau.pseudo.intercept = plateau.pseudo.intercept
df_coefficients[!df_coefficients$combR %in% plateau,]$plateau.pseudo.intercept = 'NA'
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[trend.decreasing,1]],]$logistic.plateau = 'non-significant decrease'
df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[stat.decreasing,1]],]$logistic.plateau = 'significant decrease'
df_coefficients[df_coefficients$combR %in% plateau,]$logistic.plateau = 'plateau'
df_coefficients[df_coefficients$combR %in% unclear,]$logistic.plateau = 'plateau fail'
df_coefficients[df_coefficients$combR %in% no.plateau & df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$logistic.plateau = 'non-significant increase no plateau'
df_coefficients[df_coefficients$combR %in% no.plateau & df_coefficients$combR %in% all_combR[tokeep1[stat.increasing,1]],]$logistic.plateau = 'significant increase no plateau'
df_coefficients[df_coefficients$combR %in% unclear & df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$logistic.plateau = 'non-significant increase no plateau'
df_coefficients[df_coefficients$combR %in% unclear & df_coefficients$combR %in% all_combR[tokeep1[stat.increasing,1]],]$logistic.plateau = 'significant increase no plateau'

# add relevant information to regression coefficient data frame:

# median antibiotic consumption
temp = ddply(amr_summary,.(combR,Country,Pathogen,Antibiotic,Antibiotic_class,patientType),summarise,median(p,na.rm=T))
colnames(temp)[which(colnames(temp) == "..1")] = 'Median_AMR'
colnames(temp)[which(colnames(temp) == "patientType")] = 'Patient.type'

# add data on antibiotic class and median antibiotic consumption
df_full = merge(temp,df_coefficients,by = c('combR','Patient.type'),all.x = FALSE,all.y = TRUE)
# get rid of rows where not enough data for fitting
df_full = df_full[which(df_full$logistic.plateau != "NA"),]

# add r^2 of best fitting model

tokeep1 = as.data.frame(tokeep1)
colnames(tokeep1) = c('index',"slope","slope.ci.low","slope.ci.high","aic","slope.pval","r2",'pseudo.intercept',"best.fit.error")
tokeep1$combR = all_combR[tokeep1[,1]]

testing_2 = as.data.frame(testing_2)
colnames(testing_2) = c('index',"aic","r2","best.fit.error")
testing_2$combR = all_combR[testing_2$index]

# add the pseuo.intersept of standard model and the r2 of best fitting model to df_full

df_full = merge(df_full,tokeep1[,c(7:10)],by = "combR")
is.stable = which(df_full$logistic.plateau == 'plateau')
plateau.r2 = testing_2[which(testing_2$combR %in% df_full[is.stable,]$combR),]$r2
df_full[is.stable,]$r2 = plateau.r2

plateau.error = testing_2[which(testing_2$combR %in% df_full[is.stable,]$combR),]$best.fit.error
df_full[is.stable,]$best.fit.error = plateau.error


# clean up column names
colnames(df_full)[colnames(df_full) == 'r2'] = "r2.best.fit"
colnames(df_full)[colnames(df_full) == 'slope'] = "standard.slope"
colnames(df_full)[colnames(df_full) == 'slope_stderr'] = "standard.slope.stderr"
colnames(df_full)[colnames(df_full) == 'pseudo.intercept'] = "standard.pseudo.intercept"
colnames(df_full)[colnames(df_full) == 'logistic.plateau.K'] = "plateau.k"
colnames(df_full)[colnames(df_full) == 'slope_plateau'] = "plateau.slope"
colnames(df_full)[colnames(df_full) == 'slope_plateau_stderr'] = "plateau.slope.stderr"

# re-order
df_full = df_full[,c(1,2,3,4,5,6,7,8,11,9,10,17,14,15,12,13,16,18)]

# write results to file
#write.csv(tokeep1,paste("data/results/first_log_fit_",patient_type,".csv",sep=''), row.names = FALSE)
#write.csv(testing_2,paste("data/results/second_log_fit_",patient_type,".csv",sep=''), row.names = FALSE)
write.csv(df_full,paste("data/results/temporal_trend_",patient_type,".csv",sep=''), row.names = FALSE)

}