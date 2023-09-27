rm(list = ls())
##########                                                                CORRELATE AMR AND AMC                                                                ##########

# it is possible to start here and load AMR and AMC files
# TODO DEAL WITH ANTIMICROBIAL.TYPE

library(dotenv)
library(dtw) 
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(forestplot)
library(DescTools)
library(confintr)
load_dot_env(file = ".env")

setwd(Sys.getenv(c("HOME_DIR")))

#set flag if you want to do AMR AMC correlation over countries or time
correlation_type = 'Time' # 'Time'
option = 'other' #mean

amc_summary <- read.csv(file = "data/2018/csv/summary_AMC_byclass.csv")
amr_summary <- read.csv(file = "data/2018/csv/summary_AMR.csv")
amc_tot <- read.csv(file = "data/2018/csv/total_AMC_byclass.csv")
amr_tot <- read.csv(file = "data/2018/csv/total_AMR.csv")

# slight data cleaning
colnames(amr_summary)[colnames(amr_summary)=="DateUsedForStatisticsYear"] <- "Year"
colnames(amr_summary)[colnames(amr_summary)=="ReportingCountry"] <- "Country"
amr_summary <- amr_summary[!is.na(amr_summary$Pathogen),]
amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Antibiotic, amr_summary$patientType, amr_summary$Year, sep = "|")
amc_summary$combC <- paste(amc_summary$class, amc_summary$Antimicrobial.Type, amc_summary$Sector, amc_summary$Year, sep = "|")

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

amr_summary$Antibiotic_class[amr_summary$Antibiotic=="AMC"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="AMK"] <- "J01G"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="AMP"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="AMX"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="AZM"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CAZ"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CIP"] <- "J01M"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="COL"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CLO"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CLR"] <- "J01F"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CRO"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CTX"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="DAP"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="DIC"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="DOR"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="ERY"] <- "J01F"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="ETP"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="FEP"] <- "J01R"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="FLC"] <- "J01R"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="FOX"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="GEN"] <- "J01G"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="GEH"] <- "J01G"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="IMP"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="IPM"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="LNZ"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="LVX"] <- "J01M"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="MEM"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="MET"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="MFX"] <- "J01M"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="NAL"] <- "J01M"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="NET"] <- "J01G"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="NOR"] <- "J01R"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="PIP"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="POL"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="OXA"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="OFX"] <- "J01M"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="PEN"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="RIF"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="TEC"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="TGC"] <- "J01A"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="TOB"] <- "J01G"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="TZP"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="VAN"] <- "J01X"


# #clean AMC data such that 'all' is lost and is instead the sum of the subgroups
amc_summary$combC_ <- paste(amc_summary$class, amc_summary$Country, amc_summary$Sector, amc_summary$Year,  sep = "|")
amc_summary_copy <- amc_summary %>% filter(Antimicrobial.Type %in% c('Oral','Parenteral')) %>% group_by(combC_) %>% summarise(sum_DDD = sum(DDD))
amc_summary_copy <- merge(amc_summary, amc_summary_copy, by='combC_')  %>% select(-combC,-Antimicrobial.Type, -DDD) %>% unique()
amc_summary <- amc_summary_copy %>% select(-combC_) %>% rename(DDD = sum_DDD)
amc_summary$uptake <- 'all'
amc_summary$combC <- paste(amc_summary$class, amc_summary$Sector,  sep = "|")

# DEFINE THE RELEVANT ANTIMICROBIAL TYPE TO LOOK AT:
# library(plyr)
# tmp <- ddply(amc_summary, .(class, Country, Sector), summarise,
#              N_types = length(unique(Antimicrobial.Type)), type1 = sort(unique(Antimicrobial.Type))[1],
#              type2 = sort(unique(Antimicrobial.Type))[2], type3 = sort(unique(Antimicrobial.Type))[3],
#              type4 = sort(unique(Antimicrobial.Type))[4], type5 = sort(unique(Antimicrobial.Type))[5],
#              type6 = sort(unique(Antimicrobial.Type))[6])
# tmp$all_types <- do.call(paste, c(tmp[,c("type1", "type2", "type3", "type4", "type5", "type6")], sep="-"))
# sort(
#   table(
#     tmp$all_types
#   )
# )
# table(amc_summary$Antimicrobial.Type)
# amc_summary$relevant_type <- NA
# #for(cc in unique(amc_summary$Country)){
#   for(ss in unique(amc_summary$Sector)){
#     for(cl in unique(amc_summary$class)){
#       idx <- which(
#         #amc_summary$Country==cc &
#           amc_summary$Sector==ss &
#           amc_summary$class==cl
#       )
#       tmp <- amc_summary[idx,]
#       sorted_types <- sort(table(tmp$Antimicrobial.Type), decreasing = T)
#       amc_summary[idx, "relevant_type"] <- names(sorted_types)[1]
#     }
#   }
# #}
# table(amc_summary$relevant_type, useNA = "ifany")
# table(amc_summary$relevant_type==amc_summary$Antimicrobial.Type)
# amc_summary[which(amc_summary$relevant_type=="Oral"),]


# selection of species
# selection of types of resistance
# selection of antibotic class consumed
# selection of sector
# for each year plot the relation between AMC and AMR across countries

if(correlation_type == "Country"){
  
  if(option == 'mean'){
    amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Antibiotic, amr_summary$patientType, amr_summary$Country, sep = "|")
    amc_summary$combC <- paste(amc_summary$class, amc_summary$uptake, amc_summary$Sector, amc_summary$Country,sep = "|")
    #amc_summary_median <- amc_summary %>% group_by(combC) %>% summarise(median_DDD = median(DDD),sd_DDD = sd(DDD), n = n(), uptake=Antimicrobial.Type)
    amc_summary_median <- amc_summary %>% group_by(combC) %>% summarise(median_DDD = median(DDD),sd_DDD = sd(DDD), n = n())
    
    amc_summary <- amc_summary %>% select(-Year,-DDD) %>% distinct() %>% left_join(amc_summary_median, by = 'combC') %>% distinct()
    
    amr_summary_median <- amr_summary %>% group_by(combR) %>% summarise(median_p = median(p),sd_p = sd(p), n = n())
    
    amr_summary <- amr_summary %>% select(-Year,-p, -N_S,-N_I, -N_R, -N, -ci_span, -p_min, -p_max) %>% distinct() %>% left_join(amr_summary_median, by = 'combR') %>% distinct()
    
    amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Antibiotic, amr_summary$patientType, sep = "|")
    amc_summary$combC <- paste(amc_summary$class, amc_summary$uptake, amc_summary$Sector, sep = "|")
    
    # all possible combinations:
    combR <- expand.grid(
      unique(amr_summary$Pathogen),
      unique(amr_summary$Antibiotic),
      unique(amr_summary$patientType)
    )
    names(combR) <- c("Pathogen", "Antibiotic", "patientType")
    combR$comb <- paste(combR$Pathogen, combR$Antibiotic, combR$patientType, sep = "|")
    
    combC <- expand.grid(
      unique(amc_summary$class),
      unique(amc_summary$uptake),
      unique(amc_summary$Sector)
    )
    names(combC) <- c("class", "Type", "Sector")
    combC$comb <- paste(combC$class, combC$Type, combC$Sector, sep = "|")
  
  # all possible combinations:

  # reduce to combinations that do happen often enough
  frequentC <- table(amc_summary$comb)
  frequentC <- frequentC[frequentC > 1]#used to be 15
  frequentC <- names(frequentC)
  combC <- combC[combC$comb %in% frequentC,]
  
  frequentR <- table(amr_summary$combR)
  frequentR <- frequentR[frequentR > 1]#used to be 15
  frequentR <- names(frequentR)
  combR <- combR[combR$comb %in% frequentR,]
  combR <- combR[combR$patientType!="UNK",]
  
  #all_years <- sort(unique(c(amc_summary$Year, amr_summary$Year)))
  # all_comb <- expand.grid(
  #   unique(amr_summary$Pathogen),
  #   unique(amr_summary$Antibiotic),
  #   unique(amr_summary$patientType),
  #   unique(amc_summary$class),
  #   unique(amc_summary$Sector), all_years)
  #names(all_comb) <- c("Pathogen", "Antibiotic", "patientType", "class", "Sector", "Year")
  #all_comb <- with(all_comb, all_comb[order(Pathogen, Antibiotic, patientType, class, Sector, Year), ])
  #for(j in 1:ncol(all_comb)) all_comb[, j] <- as.character(all_comb[, j])
  
  # now compute the coefficient of corrrelation for each combination
  NAvec <- rep(NA, nrow(combR) * nrow(combC))
  all_coeff <- data.frame(Pathogen = NAvec, Antibiotic = NAvec, patientType = NAvec,
                          class = NAvec, abType = NAvec, Sector = NAvec,
                          intercept = NAvec, slope = NAvec, ci_intercept_low = NAvec, ci_slope_low = NAvec, ci_intercept_high = NAvec, ci_slope_high = NAvec,
                          sd_coeff = NAvec, R2 = NAvec, p = NAvec, correlation_coef=NAvec, ci_correlation_coeff_low = NAvec, ci_correlation_coeff_high = NAvec, N_countries = NAvec, N_strains = NAvec)
  idx <- 1
  for(i in 1:nrow(combR)){
    for(j in 1:nrow(combC)){
        if(option == 'mean'){
          
          myspecies <- as.character(combR$Pathogen[i])
          myres <- as.character(combR$Antibiotic[i])
          mytype <- as.character(combR$patientType[i])
          
          myclass <- as.character(combC$class[j])
          mytypeC <- as.character(combC$Type[j])
          mysector <- as.character(combC$Sector[j])
          
          sub_amr <- amr_summary[which(
            amr_summary$Pathogen == myspecies &
              amr_summary$Antibiotic == myres &
              amr_summary$patientType == mytype
          ),]
          
          sub_amc <- amc_summary[which(amc_summary$class == myclass &
                                         amc_summary$uptake == mytypeC &
                                         amc_summary$Sector == mysector),]
        }
          merged <- merge(x = sub_amr, y = sub_amc, by = "Country")
          if(nrow(merged) == 0){
            next
          }
          if(merged$Antibiotic_class == merged$class){
            if(nrow(merged) > 4 & sum(merged$median_DDD, na.rm = TRUE) > 0 & (sum(is.na(merged$median_p)==FALSE)>4)){ # more than 2 countries and some consumption of AB
              cat(i," ", j, "\n")
              #with(merged, plot(DDD, p, xlim = c(0, max(merged$DDD)), ylim = c(0,1), main = paste(myspecies, myyear), pch = 20, ylab = paste("frequency resistance to", myres), xlab = paste("consumption of ", myclass, "in", mysector)))
              #with(merged, text(x = DDD, y = p + 0.01, labels = Country))
              lm0 <- lm(median_p ~ median_DDD, data = merged)
              #abline(lm0$coefficients[1], lm0$coefficients[2])
              sum_lm0 <- summary(lm0)
              ci <- confint(lm0)
              slope <- lm0$coefficients[2]
              correlation <- SpearmanRho(merged$median_DDD, merged$median_p ,use = "complete.obs", conf.level = 0.95)
              
              #if(!is.na(slope)) p_slope <- sum_lm0$coefficients["DDD","Pr(>|t|)"] else p_slope <- NA
              
              if(!is.na(slope)) p_slope <- sum_lm0$coefficients["median_DDD","Pr(>|t|)"] else p_slope <- NA
              all_coeff[idx,] <- c(myspecies, myres, mytype, myclass, mytypeC, mysector, lm0$coefficients[1], slope, c(ci), sum_lm0$sigma, sum_lm0$r.s, p_slope, correlation[1],correlation[2],correlation[3], nrow(merged), sum(merged$N))
              
              idx <- idx + 1
              }
            }
    } # end of loop on consumption combinations
  } # end of loop on resistance combinations
  # eliminate NA
  all_coeff <- all_coeff[-which(apply(all_coeff, 1, function(vec) all(is.na(vec)))), ]
  write.csv(x = all_coeff, file = "data/amr_summaryression_countries_mean.csv", row.names = F)
  }
  else{ # all possible combinations:
    amc_summary$combC <- paste(amc_summary$class,amc_summary$uptake, amc_summary$Sector, amc_summary$Year, sep = "|")
    
    combR <- expand.grid(
      unique(amr_summary$Pathogen),
      unique(amr_summary$Antibiotic),
      unique(amr_summary$patientType),
      unique(amr_summary$Year)
    )
    names(combR) <- c("Pathogen", "Antibiotic", "patientType", "Year")
    combR$comb <- paste(combR$Pathogen, combR$Antibiotic, combR$patientType, combR$Year, sep = "|")
    
    combC <- expand.grid(
      unique(amc_summary$class),
      unique(amc_summary$uptake),
      unique(amc_summary$Sector),
      unique(amc_summary$Year)
    )
    names(combC) <- c("class", "Type", "Sector", "Year")
    combC$comb <- paste(combC$class, combC$Type, combC$Sector, combC$Year, sep = "|")
    
    # reduce to combinations that do happen often enough
    frequentC <- table(amc_summary$comb)
    frequentC <- frequentC[frequentC > 1]#used to be 15
    frequentC <- names(frequentC)
    combC <- combC[combC$comb %in% frequentC,]
    
    frequentR <- table(amr_summary$combR)
    frequentR <- frequentR[frequentR > 1]#used to be 15
    frequentR <- names(frequentR)
    combR <- combR[combR$comb %in% frequentR,]
    combR <- combR[combR$patientType!="UNK",]
    
    #all_years <- sort(unique(c(amc_summary$Year, amr_summary$Year)))
    # all_comb <- expand.grid(
    #   unique(amr_summary$Pathogen),
    #   unique(amr_summary$Antibiotic),
    #   unique(amr_summary$patientType),
    #   unique(amc_summary$class),
    #   unique(amc_summary$Sector), all_years)
    #names(all_comb) <- c("Pathogen", "Antibiotic", "patientType", "class", "Sector", "Year")
    #all_comb <- with(all_comb, all_comb[order(Pathogen, Antibiotic, patientType, class, Sector, Year), ])
    #for(j in 1:ncol(all_comb)) all_comb[, j] <- as.character(all_comb[, j])
    
    # now compute the coefficient of corrrelation for each combination
    NAvec <- rep(NA, nrow(combR) * nrow(combC))
    all_coeff <- data.frame(Pathogen = NAvec, Antibiotic = NAvec, patientType = NAvec,
                            class = NAvec, abType = NAvec, Sector = NAvec, YearR = NAvec, YearC = NAvec,
                            intercept = NAvec, slope = NAvec, ci_intercept_low = NAvec, ci_slope_low = NAvec, ci_intercept_high = NAvec, ci_slope_high = NAvec,
                            sd_coeff = NAvec, R2 = NAvec, p = NAvec, N_countries = NAvec, N_strains = NAvec)
    idx <- 1
    stop = FALSE
    plotlist<-list()
    for(i in 1:nrow(combR)){
      for(j in 1:nrow(combC)){
        myyearR <- as.character(combR$Year[i])
        myyearC <- as.character(combC$Year[j])
        if(myyearR==myyearC){ # for now test correlations for same year only
          
          myspecies <- as.character(combR$Pathogen[i])
          myres <- as.character(combR$Antibiotic[i])
          mytype <- as.character(combR$patientType[i])
          
          myclass <- as.character(combC$class[j])
          mytypeC <- as.character(combC$Type[j])
          mysector <- as.character(combC$Sector[j])
          
          sub_amr <- amr_summary[which(
            amr_summary$Pathogen == myspecies &
              amr_summary$Year == myyearR &
              amr_summary$Antibiotic == myres &
              amr_summary$patientType == mytype
          ),]
          
          sub_amc <- amc_summary[which(amc_summary$class == myclass &
                                         amc_summary$uptake == mytypeC &
                                         amc_summary$Year == myyearC &
                                         amc_summary$Sector == mysector),]
          
          merged <- merge(x = sub_amr, y = sub_amc, by = "Country")
          if(nrow(merged) == 0){
            next
          }
          if(merged$Antibiotic_class == merged$class){
            if(nrow(merged) > 4 & sum(merged$DDD, na.rm = TRUE) > 0){ # more than 2 countries and some consumption of AB
              cat(i," ", j, "\n")
              
              merged_sub <- merged %>% filter(p != 'NA' & DDD != 'NA')
              
              #with(merged, plot(DDD, p, xlim = c(0, max(merged$DDD)), ylim = c(0,1), main = paste(myspecies, myyear), pch = 20, ylab = paste("frequency resistance to", myres), xlab = paste("consumption of ", myclass, "in", mysector)))
              #with(merged, text(x = DDD, y = p + 0.01, labels = Country))
              lm0 <- lm(p ~ DDD, weights = N, data = merged)
              #abline(lm0$coefficients[1], lm0$coefficients[2])
              plotlist[[idx]]<-ggscatter(merged_sub, x = "DDD", y = "p", size=0.5, label = "Country", font.label=c(5),
                                         add = "reg.line", add.params = list(size=0.5), conf.int = TRUE,
                                         cor.coef = FALSE, cor.method = "pearson", cor.coeff.args = list(size=3),
                                         ylab = paste("resistance to", myres), xlab = paste("median consumption of ", myclass, "in", mysector),main = paste(myspecies),ylim=c(0,1),xlim=c(0,max(merged_sub$DDD)))+
                theme_bw(base_size = 7)
                sum_lm0 <- summary(lm0)
              ci <- confint(lm0)
              slope <- lm0$coefficients[2]
              if(!is.na(slope)) p_slope <- sum_lm0$coefficients["DDD","Pr(>|t|)"] else p_slope <- NA
              all_coeff[idx,] <- c(myspecies, myres, mytype, myclass, mytypeC, mysector, myyearR, myyearC, lm0$coefficients[1], slope, c(ci), sum_lm0$sigma, sum_lm0$r.s, p_slope, nrow(merged), sum(merged$N))
              idx <- idx + 1
              
              if(idx == 20){
                stop = TRUE
                break
              }
            }
          }
        }
      } # end of loop on consumption combinations
      if(stop){break}
    } # end of loop on resistance combinations}
    # eliminate NA
    Export <- gridExtra::marrangeGrob(plotlist, nrow = 4, ncol = 3)
    # Export to a pdf file
    ggsave(filename = "output/correlations_AMR_AMC.pdf", Export,width = 210, height = 297, units = "mm")
    
    all_coeff <- all_coeff[-which(apply(all_coeff, 1, function(vec) all(is.na(vec)))), ]
    write.csv(x = all_coeff, file = "data/amr_summaryression_countries.csv", row.names = F)
  }
  
}


if(correlation_type == 'Time'){
  
  amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Antibiotic, amr_summary$patientType, amr_summary$Country, sep = "|")
  amc_summary$combC <- paste(amc_summary$class, amc_summary$uptake, amc_summary$Sector, amc_summary$Country, sep = "|")
  
  # all possible combinations:
  combR <- expand.grid(
    unique(amr_summary$Pathogen),
    unique(amr_summary$Antibiotic),
    unique(amr_summary$patientType),
    unique(amr_summary$Country)
  )
  names(combR) <- c("Pathogen", "Antibiotic", "patientType", "Country")
  combR$comb <- paste(combR$Pathogen, combR$Antibiotic, combR$patientType, combR$Country, sep = "|")
  
  combC <- expand.grid(
    unique(amc_summary$class),
    unique(amc_summary$uptake),
    unique(amc_summary$Sector),
    unique(amc_summary$Country)
  )
  names(combC) <- c("class", "Type", "Sector", "Country")
  combC$comb <- paste(combC$class, combC$Type, combC$Sector, combC$Country, sep = "|")
  
  # reduce to combinations that do happen often enough
  frequentC <- table(amc_summary$comb)
  frequentC <- frequentC[frequentC > 1]#used to be 15
  frequentC <- names(frequentC)
  combC <- combC[combC$comb %in% frequentC,]
  
  frequentR <- table(amr_summary$combR)
  frequentR <- frequentR[frequentR > 1]#used to be 15
  frequentR <- names(frequentR)
  combR <- combR[combR$comb %in% frequentR,]
  combR <- combR[combR$patientType!="UNK",]
  
  
  # now compute the coefficient of corrrelation for each combination
  NAvec <- rep(NA, nrow(combR) * nrow(combC))
  all_coeff_time <- data.frame(Pathogen = NAvec, Antibiotic = NAvec, patientType = NAvec,
                          class = NAvec, abType = NAvec, Sector = NAvec, CountryR = NAvec, CountryC = NAvec, dtwDistance = NAvec, acfLag_0 = NAvec,acfLag_1 = NAvec, acfLag_2 = NAvec,
                          ci_lower_acf = NAvec, ci_higher_acf = NAvec, ci_lower_boot=NAvec, ci_upper_boot=NAvec, N_timepoints = NAvec, N_strains = NAvec)
  idx <- 1
  for(i in 1:nrow(combR)){
    for(j in 1:nrow(combC)){
      mycountryR <- as.character(combR$Country[i])
      mycountryC <- as.character(combC$Country[j])
      if(mycountryR==mycountryC){ # for now test correlations for same year only
        myspecies <- as.character(combR$Pathogen[i])
        myres <- as.character(combR$Antibiotic[i])
        mytype <- as.character(combR$patientType[i])
        
        myclass <- as.character(combC$class[j])
        mytypeC <- as.character(combC$Type[j])
        mysector <- as.character(combC$Sector[j])
        
        sub_amr <- amr_summary[which(
          amr_summary$Pathogen == myspecies &
            amr_summary$Country == mycountryR &
            amr_summary$Antibiotic == myres &
            amr_summary$patientType == mytype
        ),]
        
        sub_amc <- amc_summary[which(amc_summary$class == myclass & 
                                       amc_summary$uptake == mytypeC &
                                       amc_summary$Country == mycountryC & 
                                       amc_summary$Sector == mysector),]
        
        merged <- merge(x = sub_amr, y = sub_amc, by = "Year")
        if(nrow(merged) == 0){
          next
        }
        if(merged$Antibiotic_class == merged$class){
          if(nrow(merged) > 7 & sum(merged$DDD, na.rm = TRUE) > 0 & sum(merged$p, na.rm=TRUE) > 7){ # more than 2 countries and some consumption of AB
            cat(i," ", j, "\n")
            merged <- merged %>% filter(p != 'NA' & DDD != 'NA')
            #with(merged, plot(DDD, p, xlim = c(0, max(merged$DDD)), ylim = c(0,1), main = paste(myspecies, myyear), pch = 20, ylab = paste("frequency resistance to", myres), xlab = paste("consumption of ", myclass, "in", mysector)))
            #with(merged, text(x = DDD, y = p + 0.01, labels = Country))
            #calculate the cross-correlation
            corr <- ccf(merged$DDD, merged$p, lag = 2, ylab = "cross-correlation", pl = FALSE) 
            ci <- ci_cor(merged$DDD, merged$p, type = "bootstrap", R = 999, seed = 1)
            upperCI <- qnorm((1 + 0.95)/2)/sqrt(corr$n.used)
            lowerCI <- -qnorm((1 + 0.95)/2)/sqrt(corr$n.used)
            
            alignement <- dtw(merged$DDD, merged$p, window.type=sakoeChibaWindow,
                              window.size=2)
            distance <- alignement$normalizedDistance
            #if(!is.na(slope)) p_slope <- sum_lm0$coefficients["DDD","Pr(>|t|)"] else p_slope <- NA
            all_coeff_time[idx,] <- c(myspecies, myres, mytype, myclass, mytypeC, mysector, mycountryR, mycountryC, distance, corr$acf[3], corr$acf[4], corr$acf[5],lowerCI,upperCI,ci$interval[1],ci$interval[2], nrow(merged), sum(merged$N))
            idx <- idx + 1
          }
        }
      }
    } # end of loop on consumption combinations
  } # end of loop on resistance combinations
  
  # eliminate NA
  all_coeff_time <- all_coeff_time[-which(apply(all_coeff_time, 1, function(vec) all(is.na(vec)))), ]
  write.csv(x = all_coeff_time, file = "data/amr_summaryression_time.csv", row.names = F)

}
#run here if you don't wish to rerun the correlation again
all_coeff <- read.csv(file = "data/amr_summaryression_countries.csv")
all_coeff$slope <- as.numeric(as.character(all_coeff$slope))
all_coeff$p <- as.numeric(as.character(all_coeff$p))
all_coeff$R2 <- as.numeric(as.character(all_coeff$R2))
all_coeff$N_countries <- as.numeric(as.character(all_coeff$N_countries))
all_coeff$N_strains <- as.numeric(as.character(all_coeff$N_strains))
all_coeff$YearR <- as.numeric(as.character(all_coeff$YearR))
all_coeff$YearC <- as.numeric(as.character(all_coeff$YearC))
hist(all_coeff$R2)

# quantile quantile plot
punif <- (1:nrow(all_coeff)) / nrow(all_coeff)
qqplot(-log10(punif), -log10(all_coeff$p), type = "p", pch = 20, xlim = c(0,10), ylim = c(0,10), xlab = "p-value uniform", ylab = "p-value")
abline(0,1)

ab_class <- data.frame(Antibiotic = c("AMP", "CAZ", "CIP", "CTX", "ERY", "GEH", "GEN", "OXA", "PEN", "RIF", "TOB", "VAN"),
           class = c("J01C", "J01D", "J01M", "J01D", "J01F", "J01G", "J01G", "J01C", "J01C", "J04A", "J01G", "J01X"))

###### determinants of the R2 ######
all_coeff$patientTypeSector <- paste(all_coeff$patientType, all_coeff$Sector, sep= "|")
all_coeff$Antibiotic_class <- ab_class$class[match(all_coeff$Antibiotic, ab_class$Antibiotic)]
all_coeff$MatchingClass <- as.character(all_coeff$Antibiotic_class)==as.character(all_coeff$class)
all_coeff$MatchingSector <- all_coeff$patientTypeSector == "INPAT|Hospital Sector" | all_coeff$patientTypeSector == "OUTPAT|Community"

lm1 <- lm(R2 ~ Pathogen + patientTypeSector + MatchingClass, data = all_coeff)
summary(lm1)
car::Anova(lm1)
lm1$coefficients

View(all_coeff[which(all_coeff$R2 > 0.9 & all_coeff$N_countries > 10), ])

lm1_sp <- lm(R2 ~  patientTypeSector + MatchingClass,
             data = all_coeff[all_coeff$Pathogen=="STRPNE",])
summary(lm1_sp)
lm1_sp$coefficients

lm1_ec <- lm(R2 ~ patientTypeSector + MatchingClass,
             data = all_coeff[all_coeff$Pathogen=="ESCCOL",])
summary(lm1_ec)

lm1_sa <- lm(R2 ~ patientTypeSector + MatchingClass,
             data = all_coeff[all_coeff$Pathogen=="STAAUR",])
summary(lm1_sa)

