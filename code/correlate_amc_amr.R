rm(list = ls())
source("code/Utils.R")
# setwd("~/ownCloud/AMR_Sonja_Martin/temporal_trends_AMR/")

# load AMC and AMR and fits
amc_summary <- read.csv(file = "data/summary_AMC_byclass_improved.csv")


# properties of AMC data
sort(unique(amc_summary$Country))
range(amc_summary$Year)
sort(table(amc_summary$Sector), decreasing = T)
sort(table(amc_summary$Antimicrobial.Type), decreasing = T)


# 2) calculate median consumption per country, sector, type, class (median across years)
amc_median <- ddply(amc_summary,.(Country,Sector,Antimicrobial.Type,Antibiotic_class), summarise, Median_AMC = median(DDD,na.rm = T), mad_AMC = mad(DDD, na.rm = T))
amc_median$cv <- amc_median$mad_AMC / amc_median$Median_AMC

# 3) temporal coefficient of variation for all, community vs. hospital sector

h1 <- hist(amc_median[amc_median$Antimicrobial.Type=="AllNew" & amc_median$Sector=="Community", "cv"], breaks = seq(0,2,0.2), plot = F)
h2 <- hist(amc_median[amc_median$Antimicrobial.Type=="AllNew" & amc_median$Sector=="Hospital Sector", "cv"],  breaks = seq(0,2,0.2), plot = F)
plot(NULL, ylim = c(0, 3), xlim = c(0, 2), xlab = "Coefficient of variation", ylab = "Density", axes = F)
points(h1$mids, h1$density, lwd=  2, col = "blue", type = "l")
points(h2$mids, h2$density, lwd=  2, col = "red", type = "l")
axis(side = 1, seq(0, 2, 0.5)); axis(side = 2, seq(0, 3, 0.5), las = 1)

###############################                                   NOW INTRODUCE AMR FITS                               ###############################

print("Now loading AMR fits...")

for(inoutpat in c("INPAT")){ # do not consider OUTPAT resistance now
  
  # INPAT OR OUTPAT AMR
  amr_fit_details <- read.csv(file = paste0("data/results/all_fit_details", inoutpat, ".csv"))
  cat.summary = read.csv(file = paste0("data/results/temporal_trend_AICc", inoutpat, ".csv"))
  cat.summary$category[cat.summary$category=="flat "] <- "flat" # for some reason there is a space  after the 'flat', correct this
  cat.summary$full.category[cat.summary$full.category=="flat "] <- "flat"
  
  # 3) compute a slope and a plateau when relevant
  stopifnot(nrow(cat.summary) == nrow(amr_fit_details)); stopifnot(all(cat.summary$combR==amr_fit_details$combR))
  cat.summary <- cbind(cat.summary, amr_fit_details) # merge with amr_fit_details
  
  # 3.1) compute plateau for relevant categories (flat, plateauing) NOT ns increase/decrease
  cat.summary$final_plateau <- NA
  cat.summary$final_plateau[cat.summary$category=="flat"] <- cat.summary$m0.plateau[cat.summary$category=="flat"]
  cat.summary$final_plateau[cat.summary$category=="plateauing"] <- cat.summary$m2.plateau[cat.summary$category=="plateauing"]
  table(is.na(cat.summary$final_plateau))
  
  # 3.2) compute slope for relevant categories
  cat.summary$final_slope_v1 <- NA    # derivative at first year
  cat.summary$final_slope_v2 <- NA # slope parameter
  cat.summary$final_slope_v3 <- NA # maximum reachable slope
  cat.summary$final_slope_v4 <- NA # maximum slope reached over the time period
  
  for(myyear in 1999:2019) {cat.summary[, paste0("deriv1", "_", myyear)] <- NA ; cat.summary[, paste0("deriv2", "_", myyear)] <- NA; }
  
  # functions for the derivative of each model evaluated at the first year when data are available
  get_deriv_min_1 <- function(vec) with(vec,              m1.slope * exp(m1.intercept - m1.slope * 1999:2019) / (1 + exp(m1.intercept - m1.slope * 1999:2019))^2)
  get_deriv_min_2 <- function(vec) with(vec, m2.plateau * m2.slope * exp(m2.intercept - m2.slope * 1999:2019) / (1 + exp(m2.intercept - m2.slope * 1999:2019))^2)
  for(i in 1:nrow(cat.summary)) cat.summary[i, paste0("deriv1_", 1999:2019)] <- get_deriv_min_1(cat.summary[i, ])
  for(i in 1:nrow(cat.summary)) cat.summary[i, paste0("deriv2_", 1999:2019)] <- get_deriv_min_2(cat.summary[i, ])
  
  # only logistic s increasing and plateau s increasing
  # derivatives at first year
  cat.summary$final_slope_v1[cat.summary$category=="logistic s increasing"]     <- sapply(which(cat.summary$category=="logistic s increasing"),      function(i) cat.summary[i, paste0("deriv1_", cat.summary$min.year[i])])
  cat.summary$final_slope_v1[cat.summary$full.category=="plateau s increasing"] <- sapply(which(cat.summary$full.category=="plateau s increasing"), function(i) cat.summary[i, paste0("deriv1_", cat.summary$min.year[i])])
  
  # directly the slope parameter
  cat.summary$final_slope_v2[cat.summary$category=="logistic s increasing"] <- cat.summary$m1.slope[cat.summary$category=="logistic s increasing"]
  cat.summary$final_slope_v2[cat.summary$full.category=="plateau s increasing"] <- cat.summary$m2.slope[cat.summary$full.category=="plateau s increasing"]
  
  # maximum rate of increase ( = plateau * slope / 4)
  cat.summary$final_slope_v3[cat.summary$category=="logistic s increasing"] <- 1 * cat.summary$m1.slope[cat.summary$category=="logistic s increasing"] / 4
  cat.summary$final_slope_v3[cat.summary$full.category=="plateau s increasing"] <- cat.summary$m2.plateau[cat.summary$full.category=="plateau s increasing"] * cat.summary$m2.slope[cat.summary$full.category=="plateau s increasing"] / 4
  
  # minimum or maximum rate of increase on the period considered
  cat.summary$final_slope_v4[cat.summary$category=="logistic s increasing"] <- sapply(which(cat.summary$category=="logistic s increasing"),         function(i) max(cat.summary[i, paste0("deriv1_", cat.summary$min.year[i]:cat.summary$max.year[i])], na.rm = T))
  cat.summary$final_slope_v4[cat.summary$full.category=="plateau s increasing"] <- sapply(which(cat.summary$full.category=="plateau s increasing"), function(i) max(cat.summary[i, paste0("deriv2_", cat.summary$min.year[i]:cat.summary$max.year[i])], na.rm = T))
  
  for(mysector in c("Community", "Hospital Sector")){
    if(mysector=="Hospital Sector") mysector2 <- "HS" else mysector2 <- mysector
    for(mytype in c("Oral", "AllNew")){
      
      tmp_amc_median <- amc_median[amc_median$Antimicrobial.Type == mytype & amc_median$Sector == mysector, ]
      tmp_merged <- merge(tmp_amc_median, cat.summary, by = c('Country','Antibiotic_class'), all = F)
      
      # correlate plateau, median, slope
      assign(x = paste0("rho_plateau_", mysector2, "_", mytype, "_", inoutpat), value =  ddply(tmp_merged, .(Pathogen, Antibiotic),    function(x) if(sum(!is.na(x$final_plateau)) > 4)   c(SpearmanRho(x$Median_AMC, x$final_plateau, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$final_plateau)))))
      assign(x = paste0("rho_median_", mysector2, "_", mytype, "_", inoutpat), value = ddply(tmp_merged, .(Pathogen, Antibiotic),     function(x) if(sum(!is.na(x$median.amr)) > 4)       c(SpearmanRho(x$Median_AMC, x$median.amr, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$median.amr)))))
      
      # correlating slopes:
      assign(x = paste0("rho_slope_v1_", mysector2, "_", mytype, "_", inoutpat), value =  ddply(tmp_merged, .(Pathogen, Antibiotic),   function(x) if(sum(!is.na(x$final_slope_v1)) > 4) c(SpearmanRho(x$Median_AMC, x$final_slope_v1, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$final_slope_v1)))))
      assign(x = paste0("rho_slope_v2_", mysector2, "_", mytype, "_", inoutpat), value =  ddply(tmp_merged, .(Pathogen, Antibiotic),   function(x) if(sum(!is.na(x$final_slope_v2)) > 4) c(SpearmanRho(x$Median_AMC, x$final_slope_v2, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$final_slope_v2)))))
      assign(x = paste0("rho_slope_v3_", mysector2, "_", mytype, "_", inoutpat), value =  ddply(tmp_merged, .(Pathogen, Antibiotic),   function(x) if(sum(!is.na(x$final_slope_v3)) > 4) c(SpearmanRho(x$Median_AMC, x$final_slope_v3, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$final_slope_v3)))))
      assign(x = paste0("rho_slope_v4_", mysector2, "_", mytype, "_", inoutpat), value =  ddply(tmp_merged, .(Pathogen, Antibiotic),   function(x) if(sum(!is.na(x$final_slope_v4)) > 4) c(SpearmanRho(x$Median_AMC, x$final_slope_v3, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$final_slope_v4)))))
     
      # slopes only for stabilising traj 
      assign(x = paste0("rho_slope_v1_", mysector2, "_", mytype, "_", inoutpat, "_stabilisingOnly"), value =  ddply(tmp_merged[tmp_merged$full.category == "plateau s increasing",  ], .(Pathogen, Antibiotic),   function(x) if(sum(!is.na(x$final_slope_v1)) > 4) c(SpearmanRho(x$Median_AMC, x$final_slope_v1, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$final_slope_v1)))))
      assign(x = paste0("rho_slope_v2_", mysector2, "_", mytype, "_", inoutpat, "_stabilisingOnly"), value =  ddply(tmp_merged[tmp_merged$full.category == "plateau s increasing",  ], .(Pathogen, Antibiotic),   function(x) if(sum(!is.na(x$final_slope_v2)) > 4) c(SpearmanRho(x$Median_AMC, x$final_slope_v2, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$final_slope_v2)))))
      assign(x = paste0("rho_slope_v3_", mysector2, "_", mytype, "_", inoutpat, "_stabilisingOnly"), value =  ddply(tmp_merged[tmp_merged$full.category == "plateau s increasing",  ], .(Pathogen, Antibiotic),   function(x) if(sum(!is.na(x$final_slope_v3)) > 4) c(SpearmanRho(x$Median_AMC, x$final_slope_v3, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$final_slope_v3)))))
      assign(x = paste0("rho_slope_v4_", mysector2, "_", mytype, "_", inoutpat, "_stabilisingOnly"), value =  ddply(tmp_merged[tmp_merged$full.category == "plateau s increasing",  ], .(Pathogen, Antibiotic),   function(x) if(sum(!is.na(x$final_slope_v4)) > 4) c(SpearmanRho(x$Median_AMC, x$final_slope_v3, use = "complete.obs",conf.level = 0.95), Ncountries = sum(!is.na(x$final_slope_v4)))))
      
      
      if(mysector == "Community" & mytype == "AllNew"){ # example spatial correlations
        
        pdf("output/example_spatial_correlations.pdf", width = 4, height = 6)
        
        sub_tmp <- tmp_merged[which(tmp_merged$Pathogen=="ESCCOL" & tmp_merged$Antibiotic=="CIP"), ]
        
        par(mfrow = c(2,1), mar = c(4,4,1,1))
        plot(NULL, xlim = c(0, 4), ylim = c(0, 0.5), las = 1, xlab = "Median antimicrobial use", ylab = "Plateau resistance")
        lm0 <- lm(final_plateau ~ Median_AMC, data = sub_tmp)
        pp <- predict.lm(object = lm0, newdata = data.frame(Median_AMC=seq(0,4,0.1)), interval = "confidence", level = 0.95)
        polygon(x = c(seq(0,4,0.1), rev(seq(0,4,0.1))), y = c(pp[,"lwr"], rev(pp[, 'upr'])), col = "gray", border = NA)
        points(sub_tmp$Median_AMC, sub_tmp$final_plateau, pch = 20)
        with(lm0, segments(x0=0,y0=coefficients[1],x1=4,y1=coefficients[1]+4*coefficients[2]))
        # SpearmanRho(sub_tmp$Median_AMC, sub_tmp$final_plateau, use = "complete.obs",conf.level = 0.95)
        
        plot(NULL, xlim = c(0, 4), ylim = c(0, 0.05), las = 1, xlab = "Median antimicrobial use", ylab = "Slope of resistance increase")
        lm1 <- lm(final_slope_v4 ~ Median_AMC, data = sub_tmp)
        #pp <- predict.lm(object = lm1, newdata = data.frame(final_slope_v4=seq(0,05,0.001)), interval = "confidence", level = 0.95)
        #polygon(x = c(seq(0,05,0.001), rev(seq(0,05,0.001))), y = c(pp[,"lwr"], rev(pp[, 'upr'])), col = "gray", border = NA)
        points(sub_tmp$Median_AMC, sub_tmp$final_slope_v4, pch = 20)
        with(lm1, segments(x0=0,y0=coefficients[1],x1=4,y1=coefficients[1]+4*coefficients[2]))
        # text(sub_tmp$Median_AMC, sub_tmp$final_slope_v4, sub_tmp$Country)
        # SpearmanRho(sub_tmp$Median_AMC, sub_tmp$final_slope_v4, use = "complete.obs",conf.level = 0.95)
        dev.off()
      }
      
    }
  }
  
}

# the spatial correlations:
rho_spatial_names <- ls()[grepl(pattern = "rho_", ls())]

##################### ############### ATTEMPT AT TEMPORAL CORRELATION #################################### 

# here temporal correlation between AMR and AMC
# (no need to use fit of temporal trends)

amr_summary <- read.csv("data/summary_AMR_filtered.csv")
amr_summary$combC_2 <- paste(amr_summary$Antibiotic_class, amr_summary$Country, amr_summary$patientType, sep = "|")

inoutpat <- "INPAT"
yearshift_vec <- c(-1, 0, 1); names(yearshift_vec) <- c("future", "present", "past")

for(mysector in c("Community", "Hospital Sector")){
  for(mytype in c("Oral", "AllNew")){
    for(yearshift in yearshift_vec){
      
      tmp_amc <- amc_summary[amc_summary$Antimicrobial.Type == mytype & amc_summary$Sector == mysector, ]
      tmp_amr <- amr_summary[amr_summary$patientType==inoutpat, ]
      tmp_amc$Year <- as.character(as.numeric(tmp_amc$Year) + yearshift) # shifted year
      tmp_merged <- merge(tmp_amc, tmp_amr, by = c('Country','Antibiotic_class','Year'), all = F)
      
      #tmp_merged[which(tmp_merged$Country=="France" & tmp_merged$Pathogen=="ESCCOL" & tmp_merged$Antibiotic=="AMP"), ]
      
      rho_temporal <- ddply(tmp_merged, .(Country, Pathogen, Antibiotic),    function(x) c(SpearmanRho(x$DDD, x$p, use = "complete.obs",conf.level = 0.95), Nyears = nrow(x)))
      rho_temporal <- rho_temporal[rho_temporal$Nyears>4, ]
      if(mysector=="Hospital Sector") mysector2 <- "HS" else mysector2 <- mysector
      assign(x = paste0("rho_temporal_", mysector2, "_", mytype, "_", inoutpat, "_", names(yearshift_vec)[yearshift_vec==yearshift]), value = rho_temporal)
      
      print(mysector)
      print(mytype)
      print(t.test(rho_temporal$rho))
      
      if(mysector == "Hospital Sector" & mytype == "AllNew" & yearshift == 0){ # example temporal correlations
        
        pdf("output/example_temporal_correlations.pdf", width = 4, height = 3)
        par(mar = c(4,4,1,1))
        sub_tmp <- tmp_merged[which(tmp_merged$Pathogen=="ESCCOL" & tmp_merged$Antibiotic=="CIP" & tmp_merged$Antibiotic_class=="J01M"), ]
        plot(NULL, xlim = c(0, 0.5), ylim = c(0, 0.5), las = 1, xlab = "Antimicrobial use across years", ylab = "Resistance across years")
        colset <- MetBrewer::met.brewer(name = "Hokusai1", n = 25); names(colset) <-  unique(sub_tmp$Country)
        for(cc in unique(sub_tmp$Country)){
          sub_tmp[which(sub_tmp$Country==cc),] -> truc
          if(nrow(truc)>4){ # more than 4 years as above
            points(truc$DDD, truc$p, pch = 20, col = colset[cc])
            lm2 <- lm(p ~ DDD, data = truc)
            with(lm2, segments(x0=0,y0=coefficients[1],x1=0.5,y1=coefficients[1]+0.5*coefficients[2], col = colset[cc])) 
          }
        }
        lm3 <- lmer(p ~ DDD + (1|Country), data = sub_tmp)
        segments(x0=0,y0=fixef(lm3)[1],x1=0.5,y1=fixef(lm3)[1]+0.5*fixef(lm3)[2], col = "gray", lwd = 4)
        dev.off()
        
      }
      
    }
  }
}
rm(rho_temporal)

rho_temporal_names <- ls()[grepl(pattern = "rho_temporal_", ls())] # temporal correlations objects

#################################### SAVE #################################### 

objects_to_save <- c(rho_spatial_names, rho_temporal_names)
save(list = objects_to_save, file = "data/AMR_AMC_correlation.RData")




