
get_year <- function(nn){
  return(
    sapply(nn, function(x){
      x <- strsplit(x, "_")[[1]]
      return(as.numeric(x[3]))
    })
  )
}
generate_randomised <- function(mytmp){
  stopifnot("N" %in% colnames(mytmp))
  stopifnot("p" %in% colnames(mytmp))
  #binomial sampling
  mytmp$random_p <- rbinom(n = rep(1, nrow(mytmp)), size =  mytmp$N, prob = mytmp$p) / mytmp$N
  return(mytmp)
}
bs_ss <- function(mytmp){ # boostrap CI
  stopifnot("N" %in% colnames(mytmp))
  stopifnot("p" %in% colnames(mytmp))
  mytmp_rand <- generate_randomised(mytmp) # generate randomised dataset
  stopifnot("random_p" %in% colnames(mytmp_rand))
  bs_spl <- smooth.spline(x = mytmp_rand$Year, y = mytmp_rand$random_p, w = 1 / mytmp_rand$N, df = 5)
  bs_pred <- predict(bs_spl)
  bs_pred1 <- predict(bs_spl, deriv = 1)
  bs_pred2 <- predict(bs_spl, deriv = 2)
  return(rbind(bs_pred$x, bs_pred$y, bs_pred1$y, bs_pred2$y))
}

# it is possible to start here and load AMR and AMC files
# TODO DEAL WITH ANTIMICROBIAL.TYPE


amc_summary <- read.csv(file = "data/summary_AMC_byclass.csv")
amr_summary <- read.csv(file = "data/summary_AMR.csv")
#amc_tot <- read.csv(file = "data/2018/csv/total_AMC_byclass.csv")
#amr_tot <- read.csv(file = "data/2018/csv/total_AMR.csv")

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
all_combR <- sort(unique(amr_summary$combR))
