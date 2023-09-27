rm(list = ls())
library(readxl)
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library("ggpubr")
library(sqldf)

if(grepl(pattern="martin", x = getwd())){
  new_data_folder <- "../../data/original/" 
  amc_dir <- "../../data/original/"
} else {
  new_data_folder <- "~/ownCloud/AMR_Sonja_Martin/"
  amc_dir <- "~/Dropbox (Infectious Disease)/ECDCdata/"
}
setwd(amc_dir)

#MDR flag to make the addition of MDR to the existing analysis facultative
MDR <- FALSE

#What we did now was to just get rid off all other combinations and sum oral and parenteral together

# 10/05/2023 FB cleaned the bits of code that seem unnecessary by commenting them out
# maybe can remove it from code if Martin agrees it is unnecessary?

# in case we want to revisit resistance profiles
# resistance_profile <- FALSE

lf <- list.files(amc_dir)
AMR <- fread(paste0(new_data_folder, "AMRTEST.csv"))
file_AMC <- lf[which(grepl(pattern = "AMC_", x = lf))]

##########                                                                CLEAN AND SUMMARISE AMR DATA                                                                ##########

# list of all ABs for which we would have antibiotic resistance data:
all_ABs <- c("AMC","AMK","AMP","AMX","CAZ","CIP","CTX","ETP","GEN","IPM","MEM","TZP","CLR","ERY","LVX","OXA","PEN","FLC","FOX","LNZ","RIF","TEC","VAN","TGC","DAP","COL","CRO","MFX","FEP","TOB","MET","PIP","NAL","OFX","NOR","AZM")

# summarise data by year, country, pathogen, antibiotic, patientType
# and generate a genotype (one line = multiple antibiotics tests) per line

# all_a_wide <- list()
# all_dims <- list()

AMR <- AMR[which(AMR$DateUsedForStatisticsYear < 2020), ]

# clean outpatient / inpatient
AMR$patientType[AMR$patientType=="NULL"] <- "UNK"
AMR$patientType[AMR$patientType=="O"] <- "OUTPAT" #  Patients that go to the hospital for Dialysis, other Day Hospital Care and to Emergency room should be classified as O

#stopifnot(all(a$SIR %in% c("S", "I", "R"))) # NOT FULFILLED FOR NOW (SOME NA)
#amr_summary <- ddply(AMR, .(DateUsedForStatisticsYear, ReportingCountry, Pathogen, Antibiotic, patientType), summarise, N_S = sum(SIR=="S", na.rm = T), N_R = sum(SIR=="R", na.rm = T), N_I = sum(SIR=="I", na.rm = T))
amr_summary <- sqldf("SELECT DateUsedForStatisticsYear, ReportingCountry, Pathogen, Antibiotic, patientType,
            SUM(SIR='S') AS N_S, SUM(SIR='R') AS N_R, SUM(SIR='I') AS N_I
            FROM AMR
            GROUP BY DateUsedForStatisticsYear, ReportingCountry, Pathogen, Antibiotic, patientType") # this gives identical results and is much faster

# FB added this bit 10/05/2023 to investigate MDR
# create isolate-specific combination:
AMR$comb <- paste(AMR$PatientCounter, AMR$IsolateId, AMR$DateUsedForStatisticsYear,  AMR$ReportingCountry, AMR$Pathogen, AMR$patientType, sep = "|")
a <- sqldf("SELECT comb,
            COUNT(SIR) AS N_tested, SUM(SIR='S') AS N_negative, SUM(SIR='R' or SIR='I') AS N_positive, DateUsedForStatisticsYear  AS DateUsedForStatisticsYear, ReportingCountry AS ReportingCountry, Pathogen AS Pathogen, GROUP_CONCAT(Antibiotic) as Antibiotic, patientType AS patientType
            FROM AMR
            GROUP BY comb") # sqldf amazingly faster (x100) than ddply
head(a)

# summarise by year, country, pathogen, patientType
a$N_tested <- as.double(a$N_tested) # need to convert to double for the AVG function to return a double
a2 <- sqldf("SELECT DateUsedForStatisticsYear, ReportingCountry, Pathogen, patientType,
            'MDR' AS Antibiotic, AVG(1.0 * N_tested) AS N_tested, SUM(N_positive > 1) AS N_R, SUM(N_positive < 2) AS N_S, 0 AS N_I
            FROM a
            GROUP BY DateUsedForStatisticsYear, ReportingCountry, Pathogen, patientType")

# amr_summary used to have duplicate entries because some countries had multiple xlsx files
amr_summary$test <- paste(amr_summary$DateUsedForStatisticsYear,
                          amr_summary$ReportingCountry,
                          amr_summary$Pathogen,
                          amr_summary$Antibiotic,
                          amr_summary$patientType, sep = "|")
stopifnot(length(unique(amr_summary$test))==nrow(amr_summary))

# FB adds MDR df to the main df, 10/05/2023
amr_summary$N_tested <- NA
a2$test <- paste(a2$DateUsedForStatisticsYear,
                 a2$ReportingCountry,
                 a2$Pathogen,
                 a2$Antibiotic,
                 a2$patientType, sep = "|")
stopifnot(all(sort(names(amr_summary))==sort(names(a2)))) # check names are the same
good_order <- names(amr_summary)
a2 <- a2[, good_order]
stopifnot(all(names(amr_summary) == names(a2)))
amr_summary2 <- amr_summary
if(MDR == TRUE){
  amr_summary2 <- rbind(amr_summary2, a2) # BIND THE TWO DF TOGETHER
}

# compute frequencies with confidence intervals:
amr_summary2$N <- amr_summary2$N_S + amr_summary2$N_I + amr_summary2$N_R
amr_summary2$p <- (amr_summary2$N_I+amr_summary2$N_R) / amr_summary2$N
amr_summary2$ci_span <- 1.96 * sqrt(amr_summary2$p * (1 - amr_summary2$p) / amr_summary2$N)
amr_summary2$p_min <- amr_summary2$p - amr_summary2$ci_span
amr_summary2$p_max <- amr_summary2$p + amr_summary2$ci_span

dim(amr_summary2)
names(amr_summary2)
# this seems similar to the file that was saved, except the save filed has two extra columns: Antibiotic_long and Antibiotic_class. These colums are added in 'data_cleaning.R'

# WRITE DOWN THE DATA

write.csv(x = amr_summary2, file = paste0(new_data_folder, "temporal_trends_AMR/data/original/summary_AMR.csv"), row.names = F)

################################################################################################################################################################################
##########                                                                CLEAN AND SUMMARISE AMC DATA                                                                ##########
################################################################################################################################################################################


amc <- read_excel(paste0("data/2018/", file_AMC), sheet = 1, col_names = TRUE, col_types = NULL, na = "", skip = 0)
amc <- data.frame(amc)
summary(amc)
#filter negative values, these make no sense as DDD "The DDD is the assumed average maintenance dose per day for a drug used for its main indication in adults."[
amc <- amc %>% filter(amc$DDD.1000.inh.day>=0)

pdf("output/histogram_AMC.pdf", paper = "a4", width = 0, height = 0)
hist(amc$DDD.1000.inh.day, breaks  = 100)
lines(density(amc$DDD.1000.inh.day), col = 'red')
dev.off()
which(amc$DDD.1000.inh.day < 0)

amc <- amc %>% filter(amc$DDD.1000.inh.day>=0)
which(amc$DDD.1000.inh.day < 0)
amc$class <- substr(amc$ATCCode, 1, 4)
table(amc$class)
amc$subclass <- substr(amc$ATCCode, 1, 5)
table(amc$subclass)

# number of countries reporting by class
list_of_countries <- list()
ns <- c()
i <- 1
for(myclass in unique(amc$class)){
  idx <- which(amc$class==myclass)
  which_countries <- unique(amc$Country[idx])
  n_countries <- length(which_countries)
  list_of_countries[[i]] <- which_countries
  ns <- c(ns, n_countries)
  i <- i + 1
}
names(list_of_countries)<-unique(amc$class)

# number of countries reporting by subclass
list_of_countries <- list()
ns <- c()
i <- 1
for(myclass in unique(amc$subclass)){
  idx <- which(amc$subclass==myclass)
  which_countries <- unique(amc$Country[idx])
  n_countries <- length(which_countries)
  list_of_countries[[i]] <- which_countries
  ns <- c(ns, n_countries)
  i <- i + 1
}
names(list_of_countries)<-unique(amc$subclass)
ns
list_of_countries

# J01DA & J01DZ to eliminate
subclass_to_eliminate <- unique(amc$subclass)[which(ns < 5)]

amc_tot_byclass <- ddply(amc, .(class, Sector), summarise, tot_cons = sum(DDD.1000.inh.day, na.rm = T))
amc_summary_byclass <- ddply(amc, .(class, Country, Year, Sector, Antimicrobial.Type), summarise, DDD = sum(DDD.1000.inh.day))
amc_summary_bysubclass <- ddply(amc, .(subclass, Country, Year, Sector, Antimicrobial.Type), summarise, DDD = sum(DDD.1000.inh.day))

# TODO DEAL WITH ANTIMICROBIAL.TYPE SOMETIMES ALL, PERANTERAL, ORAL
table(amc_summary_byclass$Antimicrobial.Type)
combis <- paste(amc_summary_byclass$class, amc_summary_byclass$Country, amc_summary_byclass$Sector, amc_summary_byclass$Year, sep = "|")
table(combis, amc_summary_byclass$Antimicrobial.Type)

amc_summary_byclass$combC <- paste(amc_summary_byclass$class, amc_summary_byclass$Country, amc_summary_byclass$Antimicrobial.Type, amc_summary_byclass$Sector,  sep = "|")
amc_plot <- amc_summary_byclass %>% group_by(combC)
idx <- 1

combC_plot <- amc_summary_byclass['combC'] %>% distinct()
plotlist<-list()
for (i in 1:nrow(combC_plot)){
  print(idx)
  comb = combC_plot[i,'combC']
  amc_plot_sub = amc_plot %>% filter(combC == comb)
  
  #if there are more zeros than >0, delete that combC
  if(sum(amc_plot_sub$DDD>0)<sum(amc_plot_sub['DDD']==0)){
    amc_plot <- amc_plot %>% filter(combC != comb)
    amc_summary_byclass <- amc_summary_byclass %>% filter(combC != comb)
    next
  }
  
  iqr <- IQR(amc_plot_sub$DDD)
  Q1 <- quantile(amc_plot_sub$DDD, .25)
  Q3 <- quantile(amc_plot_sub$DDD, .75)
  #if zero lies not in the IQR, set zero to be the median and interpolate it
  if((Q1 - 1.5*iqr)>0 & sum(amc_plot_sub$DDD==0)>0){
    print(comb)
    amc_plot_sub["DDD"][amc_plot_sub["DDD"] == 0] <- NA_real_ #median(amc_plot_sub$DDD)
    print(amc_summary_byclass["DDD"][amc_summary_byclass["combC"] == comb])
    mask <- amc_summary_byclass["DDD"][amc_summary_byclass["combC"] == comb]==0
    amc_summary_byclass["DDD"][amc_summary_byclass["combC"] == comb][mask]<- NA_real_ #median(amc_plot_sub$DDD)
    print(amc_summary_byclass["DDD"][amc_summary_byclass["combC"] == comb])
  }
  
  plotlist[[idx]] <- ggscatter(amc_plot_sub, x="Year", y="DDD",add="reg.line",add.params = list(size=0.5), conf.int = TRUE,
                               cor.coef = FALSE, cor.method = "pearson", cor.coeff.args = list(size=3),main=(comb))+
    theme_bw(base_size = 7)
  idx <- idx + 1
}
Export <- gridExtra::marrangeGrob(plotlist, nrow = 4, ncol = 3)
# Export to a pdf file
ggsave(filename = "output/AMC_slope.pdf", Export,width = 210, height = 297, units = "mm")

# write summary of AM consumption
write.csv(x = amc_summary_byclass, file = "data/2018/csv/summary_AMC_byclass.csv", row.names = F)
write.csv(x = amc_summary_bysubclass, file = "data/2018/csv/summary_AMC_bysubclass.csv", row.names = F)
write.csv(x = amc_tot_byclass, file = "data/2018/csv/total_AMC_byclass.csv", row.names = F)

# look at consumption over time
allcomb <- unique(amc_summary_byclass[, c("Country", "class", "Sector")])
allcomb <- allcomb[with(allcomb, order(class, Sector, Country)), ]
allcomb$class <- as.character(allcomb$class)
amc_tot_byclass$class <- as.character(amc_tot_byclass$class)
sub_allcomb <- allcomb[allcomb$class %in% amc_tot_byclass$class[amc_tot_byclass$Sector=="Community" & amc_tot_byclass$tot_cons > 500], ] # sub-select largely consumed classes

for(j in 1:nrow(sub_allcomb)){
  mycountry <- sub_allcomb[j, 1]
  myclass <- sub_allcomb[j, 2]
  mysector <- sub_allcomb[j, 3]
  if(mysector == "Community"){
    subtab <- amc_summary_byclass[which(amc_summary_byclass$Country==mycountry & amc_summary_byclass$class == myclass & amc_summary_byclass$Sector ==  mysector),]
    if(nrow(subtab) > 10){
      with(subtab, plot(Year, DDD, type = "o", ylab = "consumption antibiotic", ylim = c(0, 1.5 * max(DDD)), main = paste(mycountry, myclass, mysector, sep = ", ")))
    }
  }
}



