rm(list = ls())
##########                                                                CORRELATE PLATEAU AND AMC                                                                ##########

# it is possible to start here and load AMR and AMC files
# TODO DEAL WITH ANTIMICROBIAL.TYPE
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(forestplot)
library(DescTools)

amc_summary <- read.csv(file = "data/2018/csv/summary_AMC_byclass.csv")
amr_summary <- read.csv(file = "data/2018/csv/summary_AMR.csv")
amc_tot <- read.csv(file = "data/2018/csv/total_AMC_byclass.csv")
amr_tot <- read.csv(file = "data/2018/csv/total_AMR.csv")

coeffiecents_df <- read.csv('data/df_coefficients.csv')
coeffiecents_df$combR <- paste(coeffiecents_df$species,coeffiecents_df$country, coeffiecents_df$ab, coeffiecents_df$Patient.type, sep = "|")
# slight data cleaning
colnames(amr_summary)[colnames(amr_summary)=="DateUsedForStatisticsYear"] <- "Year"
colnames(amr_summary)[colnames(amr_summary)=="ReportingCountry"] <- "Country"
amr_summary <- amr_summary[!is.na(amr_summary$Pathogen),]


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

#find out all different antibiotics
amr_unique <- amr_summary['Antibiotic'] %>% distinct()

#make better naming for the plots
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="AMC"] <- "Amoxicillin + clavulanic acid"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="AMK"] <- "Amikacin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="AMP"] <- "Ampicillin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="AMX"] <- "Amoxicillin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="AZM"] <- "Azithromycin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="CAZ"] <- "Ceftazidime"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="CIP"] <- "Ciprofloxacin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="CLO"] <- "Clofazimine"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="CLR"] <- "Clarithromycin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="COL"] <- "Colistin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="CRO"] <- "Ceftriaxone"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="CTX"] <- "Cefotaxime"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="DAP"] <- "Daptomycin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="DIC"] <- "DIC"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="DOR"] <- "Doripenem"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="ERY"] <- "Erythromycin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="ETP"] <- "Ertapenem"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="FEP"] <- "Cefepime"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="FLC"] <- "Fluconazole" #antifungal agent. Why in here?
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="FOX"] <- "Cefoxitin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="GEH"] <- "Gentamicin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="GEN"] <- "Gentamicin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="IMP"] <- "Imipenem"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="LNZ"] <- "Linezolid"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="LVX"] <- "Levofloxacin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="MEM"] <- "Meropenem"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="MET"] <- "Metronidazole"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="MFX"] <- "Moxifloxacin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="NAL"] <- "Nalidixic acid"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="NET"] <- "Netilmicin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="NOR"] <- "Norfloxacine"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="OFX"] <- "Ofloxacin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="OXA"] <- "Oxacillin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="PEN"] <- "Penicillin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="PIP"] <- "Piperacillin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="POL"] <- "POL"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="RIF"] <- "Rifampicin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="TEC"] <- "Teicoplanin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="TGC"] <- "Tigecycline"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="TOB"] <- "Tobramycin"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="TZP"] <- "Piperacillin + Tazobactam"
amr_summary$Antibiotic_long[amr_summary$Antibiotic=="VAN"] <- "Vancomycin"

amr_summary$Antibiotic_class[amr_summary$Antibiotic=="AMC"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="AMK"] <- "J01G"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="AMP"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="AMX"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CAZ"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CIP"] <- "J01M"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="COL"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CRO"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="CTX"] <- "J01D" 
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="ERY"] <- "J01F"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="ETP"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="FEP"] <- "J01R"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="FOX"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="GEN"] <- "J01G"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="GEH"] <- "J01G"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="IMP"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="IPM"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="LNZ"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="LVX"] <- "J01M"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="MEM"] <- "J01D"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="MFX"] <- "J01M"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="NOR"] <- "J01R"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="PIP"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="OXA"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="OFX"] <- "J01M"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="PEN"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="RIF"] <- "J01X"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="TGC"] <- "J01A"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="TOB"] <- "J01G"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="TZP"] <- "J01C"
amr_summary$Antibiotic_class[amr_summary$Antibiotic=="VAN"] <- "J01X"


amc_summary$Antibiotic_class[amc_summary$class=="J01A"] <- "J01A"
amc_summary$Antibiotic_class[amc_summary$class=="J01B"] <- "J01B"
amc_summary$Antibiotic_class[amc_summary$class=="J01C"] <- "J01C"
amc_summary$Antibiotic_class[amc_summary$class=="J01D"] <- "J01D"
amc_summary$Antibiotic_class[amc_summary$class=="J01E"] <- "J01E"
amc_summary$Antibiotic_class[amc_summary$class=="J01F"] <- "J01F"
amc_summary$Antibiotic_class[amc_summary$class=="J01G"] <- "J01G"
amc_summary$Antibiotic_class[amc_summary$class=="J01M"] <- "J01M"
amc_summary$Antibiotic_class[amc_summary$class=="J01R"] <- "J01R"
amc_summary$Antibiotic_class[amc_summary$class=="J01X"] <- "J01X"


#cleaning the AMC data
amc_summary$class[amc_summary$class=="J01A"] <- "J01A Tetracyclines"
amc_summary$class[amc_summary$class=="J01B"] <- "J01B Amphenicols"
amc_summary$class[amc_summary$class=="J01C"] <- "J01C Penicillins"
amc_summary$class[amc_summary$class=="J01D"] <- "J01D Other Betalactams"
amc_summary$class[amc_summary$class=="J01E"] <- "J01E Sulfonamides etc."
amc_summary$class[amc_summary$class=="J01F"] <- "J01F Macrolides etc."
amc_summary$class[amc_summary$class=="J01G"] <- "J01G Aminoglycosides"
amc_summary$class[amc_summary$class=="J01M"] <- "J01M Quinolones"
amc_summary$class[amc_summary$class=="J01R"] <- "J01R Combinations"
amc_summary$class[amc_summary$class=="J01X"] <- "J01X Others"

# renaming of combination annotation
amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Country, amr_summary$Antibiotic, amr_summary$patientType,  sep = "|")
amc_summary$combC <- paste(amc_summary$class, amc_summary$Country, amc_summary$Antimicrobial.Type, amc_summary$Sector,  sep = "|")


# selection of species
# selection of types of resistance
# selection of antibotic class consumed
# selection of sector
# for each year plot the relation between AMC and AMR across countries

# reduce dataset to OUTPAT as we only have fit coefficients for those
amr_summary <- amr_summary %>% filter(patientType=='INPAT')

amr_summary_new <- amr_summary %>% group_by(combR) %>% summarize(sumN = sum(N))

# #plot the AMC consumption over time
# amc_plot <- amc_summary %>% group_by(combC)
# idx <- 1
# 
# combC_plot <- amc_plot['combC'] %>% distinct()
# plotlist<-list()
# for (i in 1:nrow(combC_plot)){
#   print(idx)
#   comb = combC_plot[i,'combC']
#   amc_plot_sub = amc_plot %>% filter(combC == comb)
#   plotlist[[idx]] <- ggscatter(amc_plot_sub, x="Year", y="DDD",add="reg.line",add.params = list(size=0.5), conf.int = TRUE,
#                                cor.coef = FALSE, cor.method = "pearson", cor.coeff.args = list(size=3),main=(comb))+
#     theme_bw(base_size = 7)
#   idx <- idx + 1
# }
# Export <- gridExtra::marrangeGrob(plotlist, nrow = 4, ncol = 3)
# # Export to a pdf file
# ggsave(filename = "output/AMC_slope.pdf", Export,width = 210, height = 297, units = "mm")

# calculate the median of the AMC consumption over the years
amc_summary_median <- amc_summary %>% group_by(combC) %>% summarise(median_DDD = median(DDD),sd_DDD = sd(DDD), n = n(), uptake=Antimicrobial.Type)

# calculate the slope of the AMC over time
amc_summary_slope <- amc_summary %>% group_by(combC) %>% summarise(slope = lm(DDD~Year)$coefficients[2], n = n(), uptake=Antimicrobial.Type, upper=confint(lm(DDD~Year))[2], lower=confint(lm(DDD~Year))[4])
amc_summary_slope <- amc_summary %>% select(-Year,-DDD) %>% distinct() %>% left_join(amc_summary_slope, by = 'combC')
amc_summary_slope <- amc_summary_slope %>% filter(uptake=='Oral') %>% filter(Sector=='Community') %>% distinct()
#collapse the dataframe to be without year and DDD and join the medians from above to it
amc_summary <- amc_summary %>% select(-Year,-DDD) %>% distinct() %>% left_join(amc_summary_median, by = 'combC') %>% distinct()

amc_summary_median <- amc_summary %>% filter(uptake=='Oral')
amc_summary_median <- amc_summary_median %>% filter(Sector=='Community') %>% distinct()
#plot the median DDD and the variance DDD
amc_summary_median[is.na(amc_summary_median)] <- 0
p<- ggplot(amc_summary_median, aes(x=Country, y=median_DDD,fill=class)) + 
  geom_bar(stat = 'identity', position = "dodge2") +
  geom_errorbar(aes(ymin=median_DDD-sd_DDD, ymax=median_DDD+sd_DDD), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  #ylim(-10,50)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("output/medianDDDplot.png", width=15, height=15, dpi=500)

# plot the median DDD and the variance DDD
amc_summary_slope[is.na(amc_summary_slope)] <- 0
p<- ggplot(amc_summary_slope, aes(x=Country, y=slope,fill=class)) + 
  geom_bar(stat = 'identity',position = "dodge2") +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(.9)) +
  #coord_flip(ylim=c(-2,2)) +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("output/DDDSlopeplot.png", width=15, height=15, dpi=500)

# plot the slope as histogram
p<- ggplot(amc_summary_slope, aes(y=slope,fill=class)) + 
  geom_histogram() +
  coord_flip() 
#ylim(-10,50)+
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("output/DDDSlopehistogram.png", width=15, height=15, dpi=500)

write.csv(x = amc_summary_slope, file = "data/amc_consumption_regression.csv", row.names = F)


# add plateau and slope to the amr_summary dataframe
coefficients_summary <- coeffiecents_df %>% select(120:155)
amr_summary <- amr_summary %>% select(-Year,-c(6:13)) %>% distinct() %>% left_join(coefficients_summary, by = 'combR')%>%left_join(amr_summary_new,by='combR')%>%select(-c(8:29)) %>% filter(logistic.plateau.K<=1)#%>%left_join(amr_summary,by='combR')

#perform merging of amr and amc data to visualise the spread of the data or its clustering
amc_summary_median_copy <- amc_summary_median
amr_amc_joined <- left_join(amr_summary, amc_summary_median_copy, by = c('Country', 'Antibiotic_class')) #%>% filter(logistic.plateau.K<=1)

p<- ggplot(amr_amc_joined, aes(x=median_DDD, y=logistic.plateau.K, color=class)) + 
  geom_pointrange(aes(xmin=median_DDD-sd_DDD, xmax=median_DDD+sd_DDD))+
  xlim(-5,20)
  #coord_flip(ylim=c(-2,2)) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("output/DDD_Kplot.png", width=15, height=15, dpi=500)


# renaming of combination annotation
amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Antibiotic, amr_summary$patientType,  sep = "|")
amc_summary$combC <- paste(amc_summary$class,  amc_summary$Sector,  sep = "|")

# #clean AMC data such that 'all' is lost and is instead the sum of the subgroups
amc_summary$combC_ <- paste(amc_summary$class, amc_summary$Country, amc_summary$Sector,  sep = "|")
amc_summary_copy <- amc_summary %>% filter(uptake %in% c('Oral','Parenteral')) %>% group_by(combC_) %>% summarise(median_DDD_all = sum(median_DDD))
amc_summary_copy <- merge(amc_summary, amc_summary_copy, by='combC_')  %>% select(-combC,-Antimicrobial.Type, -median_DDD, -uptake, -sd_DDD, -n) %>% unique()
amc_summary <- amc_summary_copy %>% rename(median_DDD = median_DDD_all)
amc_summary$uptake <- 'all'
amc_summary$combC <- paste(amc_summary$class, amc_summary$uptake, amc_summary$Sector,  sep = "|")

# all possible combinations:
combR <- expand.grid(
  unique(amr_summary$Pathogen),
  #unique(amr_summary$Country),
  unique(amr_summary$Antibiotic),
  unique(amr_summary$patientType)
)
names(combR) <- c("Pathogen",  "Antibiotic", "patientType")
combR$comb <- paste(combR$Pathogen,  combR$Antibiotic, combR$patientType, sep = "|")


combC <- expand.grid(
  unique(amc_summary$class),
  #unique(amc_summary$Country),
  unique(amc_summary$uptake),
  unique(amc_summary$Sector)
)
names(combC) <- c("class", "Type", "Sector")
combC$comb <- paste(combC$class, combC$Type, combC$Sector, sep = "|")

# reduce to combinations that do happen often enough this time without countries, as we will plot over them
frequentC <- table(amc_summary$combC)
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
all_coeff_K <- data.frame(Pathogen = NAvec, Antibiotic = NAvec, patientType = NAvec,
                          class = NAvec,  abType = NAvec, Sector = NAvec,
                          intercept = NAvec, slope = NAvec, ci_intercept_low = NAvec, ci_slope_low = NAvec, ci_intercept_high = NAvec, ci_slope_high = NAvec,
                          sd_coeff = NAvec, R2 = NAvec, p = NAvec, correlation_coef=NAvec, ci_correlation_coeff_low = NAvec, ci_correlation_coeff_high = NAvec, N_countries = NAvec, N_strains = NAvec)
idx <- 1

#correlate the location of the plateau and the AMC consumption
#pdf("output/correlations_plateau_AMC.pdf", paper = "a4", width = 0, height = 0)
#par(mfrow = c(4, 3))
#par(las = 1)
combC <- combC %>% select(-comb) %>% distinct()
combR <- combR %>% select(-comb) %>% distinct()
plotlist<-list()
for(i in 1:nrow(combR)){
  for(j in 1:nrow(combC)){
    myspecies <- as.character(combR$Pathogen[i])
    myres <- as.character(combR$Antibiotic[i])
    mytype <- as.character(combR$patientType[i])
    
    myclass <- as.character(combC$class[j])
    mytypeC <- as.character(combC$Type[j])
    mysector <- as.character(combC$Sector[j])
    
    sub_amr <- amr_summary[which(
      amr_summary$Pathogen == myspecies &
        amr_summary$Antibiotic == myres
    ),]
    
    sub_amc <- amc_summary[which(amc_summary$class == myclass & 
                                   #amc_summary$Antimicrobial.Type == mytypeC &
                                   amc_summary$Sector == mysector),]
    
    merged <- merge(x = sub_amr, y = sub_amc, by = "Country")
    merged <- merged %>% mutate(pairs = median_DDD * logistic.plateau.K)
    
    if(nrow(merged) > 4 & sum(merged$median_DDD, na.rm=TRUE) > 0 & (sum(is.na(merged$logistic.plateau.K)==FALSE)>4) & (sum(is.na(merged$pairs)==FALSE)>4)){ # more than 4 countries and some consumption of AB and a value for the plateau and only "All" as a type
      cat(i," ", j, "\n")
      #print(paste0(myspecies,myres,mytype,myclass,mytypeC,mysector))
      
      merged_sub <- merged %>% filter(logistic.plateau.K != 'NA' & median_DDD != 'NA')
      #with(merged, plot(median_DDD, logistic.plateau.K, xlim = c(0, max(merged$median_DDD)), ylim = c(0,1), main = paste(myspecies), pch = 20, ylab = paste("plateau in resistance to", myres), xlab = paste("median ",mytypeC," consumption of ", myclass, "in", mysector)))
      #with(merged, text(x = median_DDD, y = logistic.plateau.K + 0.01, labels = Country,cex=0.5))
      #for (x in logistic)
      lm0 <- lm(logistic.plateau.K ~ median_DDD, data = merged)
      #correlation <- cor(merged$median_DDD, merged$logistic.plateau.K ,method = "spearman",use = "complete.obs")
      correlation <- SpearmanRho(merged$median_DDD, merged$logistic.plateau.K ,use = "complete.obs", conf.level = 0.95)

      #abline(lm0$coefficients[1], lm0$coefficients[2])
      plotlist[[idx]]<-ggscatter(merged_sub, x = "median_DDD", y = "logistic.plateau.K", size=0.5, label = "Country", font.label=c(5),
                                 add = "reg.line", add.params = list(size=0.5), conf.int = TRUE,
                                 cor.coef = FALSE, cor.method = "pearson", cor.coeff.args = list(size=3),
                                 ylab = paste("plateau in resistance to", myres), xlab = paste("median ",mytypeC," consumption of ", myclass, "in", mysector),main = paste(myspecies),ylim=c(0,1),xlim=c(0,max(merged_sub$median_DDD)))+
        theme_bw(base_size = 7)
      sum_lm0 <- summary(lm0)
      ci <- confint(lm0)
      slope <- lm0$coefficients[2]
      if(!is.na(slope)) p_slope <- sum_lm0$coefficients["median_DDD","Pr(>|t|)"] else p_slope <- NA
      all_coeff_K[idx,] <- c(myspecies, myres, mytype, myclass, mytypeC, mysector, lm0$coefficients[1], slope, c(ci), sum_lm0$sigma, sum_lm0$r.s, p_slope, correlation[1],correlation[2],correlation[3], nrow(merged), sum(merged$N))
      idx <- idx + 1
      
    }
  } # end of loop on consumption combinations
} # end of loop on resistance combinations
Export <- gridExtra::marrangeGrob(plotlist, nrow = 4, ncol = 3)
# Export to a pdf file
ggsave(filename = "output/correlations_plateau_AMC.pdf", Export,width = 210, height = 297, units = "mm")
# eliminate NA
all_coeff_K <- all_coeff_K[-which(apply(all_coeff_K, 1, function(vec) all(is.na(vec)))), ]
write.csv(x = all_coeff_K, file = "data/all_coeff_AMC_K.csv", row.names = F)

#run here if you don't wish to rerun the correlation again
# all_coeff_K <- read.csv(file = "data/all_coeff_AMC_K.csv")
# all_coeff_K$slope <- as.numeric(as.character(all_coeff_K$slope))
# all_coeff_K$p <- as.numeric(as.character(all_coeff_K$p))
# all_coeff_K$R2 <- as.numeric(as.character(all_coeff_K$R2))
# all_coeff_K$N_countries <- as.numeric(as.character(all_coeff_K$N_countries))
# all_coeff_K$N_strains <- as.numeric(as.character(all_coeff_K$N_strains))

#hist(all_coeff_K$R2)

all_coeff <- data.frame(Pathogen = NAvec, Antibiotic = NAvec, patientType = NAvec,
                        class = NAvec, abType = NAvec, Sector = NAvec,
                        intercept = NAvec, slope = NAvec, ci_intercept_low = NAvec, ci_slope_low = NAvec, ci_intercept_high = NAvec, ci_slope_high = NAvec,
                        sd_coeff = NAvec, R2 = NAvec, p = NAvec, correlation_coef=NAvec, ci_correlation_coeff_low = NAvec, ci_correlation_coeff_high = NAvec, N_countries = NAvec, N_strains = NAvec)
idx <- 1
#correlate the slope and the AMC consumption
#pdf("output/correlations_slope_AMC.pdf", paper = "a4", width = 0, height = 0, onefile=TRUE)
plotlist<-list()
for(i in 1:nrow(combR)){
  for(j in 1:nrow(combC)){
    myspecies <- as.character(combR$Pathogen[i])
    myres <- as.character(combR$Antibiotic[i])
    mytype <- as.character(combR$patientType[i])
    
    myclass <- as.character(combC$class[j])
    mytypeC <- as.character(combC$Type[j])
    mysector <- as.character(combC$Sector[j])
    
    sub_amr <- amr_summary[which(
      amr_summary$Pathogen == myspecies &
        amr_summary$Antibiotic == myres
    ),]
    
    sub_amc <- amc_summary[which(amc_summary$class == myclass & 
                                   #amc_summary$Antimicrobial.Type == mytypeC &
                                   amc_summary$Sector == mysector),]
    
    merged <- merge(x = sub_amr, y = sub_amc, by = "Country")
    merged <- merged %>% mutate(pairs = median_DDD * slope)
    if(nrow(merged) > 4 & sum(merged$median_DDD, na.rm=TRUE) > 0 & (sum(is.na(merged$slope)==FALSE)>4) & (sum(is.na(merged$pairs)==FALSE)>4)){ # more than 4 countries and some consumption of AB and a value for the plateau
      cat(i," ", j, "\n")
      
      merged_sub <- merged %>% filter(slope != 'NA' & median_DDD != 'NA')
      
      #with(merged, plot(median_DDD, initial.slope, xlim = c(0, max(merged$median_DDD)), ylim = c(-0.6,0.6), main = paste(myspecies), pch = 20, ylab = paste("initial slope in resistance to", myres), xlab = paste("median ",mytypeC," consumption of ", myclass, "in", mysector)))
      #with(merged, text(x = median_DDD, y = initial.slope + 0.01, labels = Country,cex=0.5))
      lm0 <- lm(slope ~ median_DDD, data = merged)
      #correlation <- cor(merged$median_DDD, merged$slope ,method = "spearman",use = "complete.obs")
      correlation <- SpearmanRho(merged$median_DDD, merged$slope ,use = "complete.obs", conf.level = 0.95)
      #abline(lm0$coefficients[1], lm0$coefficients[2])
      plotlist[[idx]]<-ggscatter(merged_sub, x = "median_DDD", y = "slope", size=0.5, label = "Country", font.label=c(5),
                                 add = "reg.line", add.params = list(size=0.5), conf.int = TRUE,
                                 cor.coef = FALSE, cor.method = "pearson", cor.coeff.args = list(size=3),
                                 xlab = paste("median ",mytypeC," consumption of ", myclass, "in", mysector), ylab = paste("slope of logistic fit in resistance to", myres),main = paste(myspecies),ylim = c(-0.6,0.6),xlim=c(0,max(merged_sub$median_DDD)))+
        theme_bw(base_size = 7)
      
      
      sum_lm0 <- summary(lm0)
      ci <- confint(lm0)
      slope <- lm0$coefficients[2]
      if(!is.na(slope)) p_slope <- sum_lm0$coefficients["median_DDD","Pr(>|t|)"] else p_slope <- NA
      all_coeff[idx,] <- c(myspecies, myres, mytype, myclass, mytypeC, mysector, lm0$coefficients[1], slope, c(ci), sum_lm0$sigma, sum_lm0$r.s, p_slope, correlation[1],correlation[2],correlation[3], nrow(merged), sum(merged$N))
      idx <- idx + 1
      
    }
  } # end of loop on consumption combinations
} # end of loop on resistance combinations
Export <- gridExtra::marrangeGrob(plotlist, nrow = 4, ncol = 3)
# Export to a pdf file
ggsave(filename = "output/correlations_slope_AMC.pdf", Export,width = 210, height = 297, units = "mm")
# eliminate NA
all_coeff <- all_coeff[-which(apply(all_coeff, 1, function(vec) all(is.na(vec)))), ]
write.csv(x = all_coeff, file = "data/all_coeff_AMC_slope.csv", row.names = F)


all_coeff <- data.frame(Pathogen = NAvec, Antibiotic = NAvec, patientType = NAvec,
                        class = NAvec,  abType = NAvec, Sector = NAvec,
                        intercept = NAvec, slope = NAvec, ci_intercept_low = NAvec, ci_slope_low = NAvec, ci_intercept_high = NAvec, ci_slope_high = NAvec,
                        sd_coeff = NAvec, R2 = NAvec, p = NAvec, correlation_coef=NAvec, ci_correlation_coeff_low = NAvec, ci_correlation_coeff_high = NAvec, N_countries = NAvec, N_strains = NAvec)
idx <- 1

#correlate the slope of plateaud cases and the AMC consumption
#pdf("output/correlations_slope_AMC.pdf", paper = "a4", width = 0, height = 0, onefile=TRUE)
plotlist<-list()
for(i in 1:nrow(combR)){
  for(j in 1:nrow(combC)){
    myspecies <- as.character(combR$Pathogen[i])
    myres <- as.character(combR$Antibiotic[i])
    mytype <- as.character(combR$patientType[i])
    
    myclass <- as.character(combC$class[j])
    #mytypeC <- as.character(combC$Type[j])
    mysector <- as.character(combC$Sector[j])
    
    sub_amr <- amr_summary[which(
      amr_summary$Pathogen == myspecies &
        amr_summary$Antibiotic == myres
    ),]
    
    sub_amc <- amc_summary[which(amc_summary$class == myclass & 
                                   #amc_summary$Antimicrobial.Type == mytypeC &
                                   amc_summary$Sector == mysector),]
    
    merged <- merge(x = sub_amr, y = sub_amc, by = "Country")
    merged <- merged %>% mutate(pairs = median_DDD * slope_plateau)
    if(nrow(merged) > 4 & sum(merged$median_DDD, na.rm=TRUE) > 0 & (sum(is.na(merged$slope_plateau)==FALSE)>4) & (sum(is.na(merged$pairs)==FALSE)>4)){ # more than 4 countries and some consumption of AB and a value for the plateau
      cat(i," ", j, "\n")
      
      merged_sub <- merged %>% filter(slope_plateau != 'NA' & median_DDD != 'NA')
      
      #with(merged, plot(median_DDD, initial.slope, xlim = c(0, max(merged$median_DDD)), ylim = c(-0.6,0.6), main = paste(myspecies), pch = 20, ylab = paste("initial slope in resistance to", myres), xlab = paste("median ",mytypeC," consumption of ", myclass, "in", mysector)))
      #with(merged, text(x = median_DDD, y = initial.slope + 0.01, labels = Country,cex=0.5))
      lm0 <- lm(slope_plateau ~ median_DDD, data = merged)
      #correlation <- cor(merged$median_DDD, merged$slope_plateau ,method = "spearman",use = "complete.obs")
      correlation <- SpearmanRho(merged$median_DDD, merged$slope_plateau ,use = "complete.obs", conf.level = 0.95)
      #abline(lm0$coefficients[1], lm0$coefficients[2])
      plotlist[[idx]]<-ggscatter(merged_sub, x = "median_DDD", y = "slope_plateau", size=0.5, label = "Country", font.label=c(5),
                                 add = "reg.line", add.params = list(size=0.5), conf.int = TRUE,
                                 cor.coef = FALSE, cor.method = "pearson", cor.coeff.args = list(size=3),
                                 xlab = paste("median ",mytypeC," consumption of ", myclass, "in", mysector), ylab = paste("slope of logistic fit plateau in resistance to", myres),main = paste(myspecies),ylim = c(-1,2),xlim=c(0,max(merged_sub$median_DDD)))+
        theme_bw(base_size = 7)
      
      
      sum_lm0 <- summary(lm0)
      ci <- confint(lm0)
      slope <- lm0$coefficients[2]
      if(!is.na(slope)) p_slope <- sum_lm0$coefficients["median_DDD","Pr(>|t|)"] else p_slope <- NA
      all_coeff[idx,] <- c(myspecies, myres, mytype, myclass, mytypeC, mysector, lm0$coefficients[1], slope, c(ci), sum_lm0$sigma, sum_lm0$r.s, p_slope, correlation[1],correlation[2],correlation[3], nrow(merged), sum(merged$N))
      idx <- idx + 1
      
    }
  } # end of loop on consumption combinations
} # end of loop on resistance combinations
Export <- gridExtra::marrangeGrob(plotlist, nrow = 4, ncol = 3)
# Export to a pdf file
ggsave(filename = "output/correlations_slope_plateau_AMC.pdf", Export,width = 210, height = 297, units = "mm")
# eliminate NA
all_coeff <- all_coeff[-which(apply(all_coeff, 1, function(vec) all(is.na(vec)))), ]
write.csv(x = all_coeff, file = "data/all_coeff_AMC_slope_plateau.csv", row.names = F)


all_coeff_autocorr <- data.frame(Pathogen = NAvec, Antibiotic = NAvec, patientType = NAvec,
                        class = NAvec, abType = NAvec, Sector = NAvec,
                        intercept = NAvec, slope = NAvec, ci_intercept_low = NAvec, ci_slope_low = NAvec, ci_intercept_high = NAvec, ci_slope_high = NAvec,
                        sd_coeff = NAvec, R2 = NAvec, p = NAvec, correlation_coef=NAvec, ci_correlation_coeff_low = NAvec, ci_correlation_coeff_high = NAvec, N_countries = NAvec, N_strains = NAvec)
idx <- 1

#correlate the slope of logistic plateau vs logistic slope
#pdf("output/correlations_slope_AMC.pdf", paper = "a4", width = 0, height = 0, onefile=TRUE)
plotlist<-list()
for(i in 1:nrow(combR)){
  myspecies <- as.character(combR$Pathogen[i])
  myres <- as.character(combR$Antibiotic[i])
  mytype <- as.character(combR$patientType[i])

  sub_amr <- amr_summary[which(
    amr_summary$Pathogen == myspecies &
      amr_summary$Antibiotic == myres
  ),]

  sub_amc <- amc_summary[which(amc_summary$class == myclass & 
                                 #amc_summary$Antimicrobial.Type == mytypeC &
                                 amc_summary$Sector == mysector),]
  
  merged <- sub_amr
  
  merged <- merged %>% mutate(pairs = logistic.plateau.K * slope)
  print(merged)
  if(nrow(merged) > 4 & sum(merged$logistic.plateau.K, na.rm=TRUE) > 0 & (sum(is.na(merged$slope)==FALSE)>4) & (sum(is.na(merged$pairs)==FALSE)>4)){ # more than 4 countries and some consumption of AB and a value for the plateau
    cat(i," ", j, "\n")
    merged_sub <- merged %>% filter(slope != 'NA' & logistic.plateau.K != 'NA')
    
    #with(merged, plot(median_DDD, initial.slope, xlim = c(0, max(merged$median_DDD)), ylim = c(-0.6,0.6), main = paste(myspecies), pch = 20, ylab = paste("initial slope in resistance to", myres), xlab = paste("median ",mytypeC," consumption of ", myclass, "in", mysector)))
    #with(merged, text(x = median_DDD, y = initial.slope + 0.01, labels = Country,cex=0.5))
    
    lm0 <- lm(logistic.plateau.K ~ slope, data = merged)
    #correlation <- cor(merged$logistic.plateau.K, merged$slope ,method = "spearman",use = "complete.obs")
    correlation <- SpearmanRho(merged$logistic.plateau.K, merged$slope ,use = "complete.obs", conf.level = 0.95)
    #abline(lm0$coefficients[1], lm0$coefficients[2])
    plotlist[[idx]]<-ggscatter(merged_sub, x = "logistic.plateau.K", y = "slope", size=0.5, label = "Country", font.label=c(5),
                               add = "reg.line", add.params = list(size=0.5), conf.int = TRUE,
                               cor.coef = FALSE, cor.method = "pearson", cor.coeff.args = list(size=3),
                               xlab = paste("plateau of logistic fit plateau in resistance to", myres), ylab = paste("slope of logistic fit plateau in resistance to", myres),main = paste(myspecies),ylim = c(-0.6,0.6),xlim=c(0,max(merged_sub$median_DDD)))+
      theme_bw(base_size = 7)
    
    
    sum_lm0 <- summary(lm0)
    ci <- confint(lm0)
    slope <- lm0$coefficients[2]
    if(!is.na(slope)) p_slope <- sum_lm0$coefficients["slope","Pr(>|t|)"] else p_slope <- NA
    all_coeff_autocorr[idx,] <- c(myspecies, myres, mytype, myclass, mytypeC, mysector, lm0$coefficients[1], slope, c(ci), sum_lm0$sigma, sum_lm0$r.s, p_slope, correlation[1],correlation[2],correlation[3], nrow(merged), sum(merged$N))
    idx <- idx + 1
      
  } # end of loop on consumption combinations
} # end of loop on resistance combinations
Export <- gridExtra::marrangeGrob(plotlist, nrow = 4, ncol = 3)
# Export to a pdf file
ggsave(filename = "output/correlations_autocorr.pdf", Export,width = 210, height = 297, units = "mm")
# eliminate NA
all_coeff_autocorr <- all_coeff_autocorr[-which(apply(all_coeff_autocorr, 1, function(vec) all(is.na(vec)))), ]
write.csv(x = all_coeff_autocorr, file = "data/all_coeff_autocorr.csv", row.names = F)

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