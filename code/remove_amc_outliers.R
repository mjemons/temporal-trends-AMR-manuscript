set.seed(1234)
amc_summary <- read.csv(file = "data/summary_AMC_byclass.csv")

amc_summary$combC =paste(amc_summary$class, amc_summary$Country, amc_summary$Sector,amc_summary$Antimicrobial.Type, sep = "|")
full_combC = unique(amc_summary$combC)

# mark outliers in amc_summary variable
amc_summary$outlier = FALSE

#loop through all AMC combinations and remove all datapoints that are outside of 3* the IQR
for(i in 1:length(full_combC)){
  tmp1 <- amc_summary[which(amc_summary$combC == full_combC[i]),]
  median.tmp1 = median(tmp1$DDD,na.rm = T)
  IQR.tmp1 = IQR(tmp1$DDD,na.rm = T)
  to.remove = (tmp1$DDD > median.tmp1 + 3* IQR.tmp1)|(tmp1$DDD < median.tmp1 - 3* IQR.tmp1)
  amc_summary[which(amc_summary$combC == full_combC[i]),]$outlier = to.remove
}

amc_clean = amc_summary[which(amc_summary$outlier == FALSE),]

#save the flagged data
write.csv(file = "data/summary_AMC_byclass_outliers_flagged.csv",amc_summary)

#save the cleaned data
write.csv(file = "data/summary_AMC_byclass_filtered.csv",amc_clean)

############################################################################################################################################
#################                               CLEAN AMC DATA FURTHER                                                     #################
############################################################################################################################################

# load AMC and AMR and fits
amc_summary <- read.csv(file = "data/summary_AMC_byclass_filtered.csv")
amc_summary$combC_2 <- paste(amc_summary$class, amc_summary$Country, amc_summary$Sector, sep = "|")
amc_summary$X <- NULL
amc_summary$Type <- amc_summary$Antimicrobial.Type

# 0) Improve AMC data by merging oral, parenteral, etc. into "All_new" category when "All" is not present or set to 0
# data.frame with this new category
new_amc_summary <- data.frame(matrix(NA, nrow = 0, ncol = ncol(amc_summary)))
print("Cleaning AMC data...")
for(mycomb in unique(amc_summary$combC_2)){
  idx <- which(amc_summary$combC_2 == mycomb)
  sub_amc_summary <- amc_summary[idx, ]
  years <- unique(sub_amc_summary$Year)
  for(myyear in years){
    
    is_valid_all_DDD <- FALSE
    all_DDD <- sub_amc_summary$DDD[sub_amc_summary$Year==myyear & sub_amc_summary$Antimicrobial.Type=="All"] # value of "All" DDD
    if(length(all_DDD) > 0){
      if(all_DDD > 0 & all(all_DDD > sub_amc_summary$DDD[sub_amc_summary$Year==myyear & sub_amc_summary$Antimicrobial.Type!="All"])){
        is_valid_all_DDD <- TRUE
      }
    }
    
    if(is_valid_all_DDD){ # Case 1: "All" already exists and non-zero, trust it
      myDDD <- all_DDD
      mytype <- "All"
    } else { # Case 2: "All" does not exist or is 0, compute the sum
      types_toconsider <- sub_amc_summary$Year==myyear & sub_amc_summary$Antimicrobial.Type!="All"
      myDDD <- sum(sub_amc_summary$DDD[types_toconsider]) # sum of all antimicrobial types considered
      mytype <- paste(sub_amc_summary$Antimicrobial.Type[types_toconsider], collapse = "|")
    }
    if(is.na(myDDD)) stop("new DDD is NA")
    myvec <- c(
      unique(sub_amc_summary$class),
      unique(sub_amc_summary$Country),
      myyear,
      unique(sub_amc_summary$Sector),
      "AllNew",
      myDDD,
      paste0(unique(sub_amc_summary$combC_2), "|AllNew"),
      unique(sub_amc_summary$Antibiotic_class),
      FALSE,
      unique(sub_amc_summary$combC_2),
      mytype
    )
    stopifnot(length(myvec)==ncol(amc_summary))
    new_amc_summary <- rbind(new_amc_summary, myvec)
  }
}
names(new_amc_summary) <- names(amc_summary)
new_amc_summary$DDD <- as.numeric(new_amc_summary$DDD)
amc_summary <- rbind(amc_summary, new_amc_summary)

# save the new amc

write.csv(amc_summary,"data/summary_AMC_byclass_improved.csv")

############################################################################################################################################
#################                                VISUALISE AMC DATA                                                        #################
############################################################################################################################################

# properties of AMC data

sort(unique(amc_summary$Country))
range(amc_summary$Year)
sort(table(amc_summary$Sector), decreasing = T)
sort(table(amc_summary$Antimicrobial.Type), decreasing = T)

# 1) visualise AMC trajectories
head(amc_summary)

amc.types = c('Oral','Parenteral','All', 'AllNew')
col.types <- RColorBrewer::brewer.pal(n= 4, name = "Set1")
pch.types <- c(20, 20, 20, 1)
lty.types <- c(1,1,1,2)
names(col.types) <- names(pch.types) <- names(lty.types) <- amc.types

plot_amc <- function(mycomb){
  idx <- which(amc_summary$combC_2 == mycomb)
  ymax <- 1.1 *max(amc_summary[idx, "DDD"])
  
  for(mytype in amc.types){
    if(mytype == "Oral"){
      plot_fun <- plot
    } else {
      plot_fun <- points
    }
    sub_idx <- which(amc_summary$combC_2 == mycomb & amc_summary$Antimicrobial.Type == mytype)
    plot_fun(amc_summary[sub_idx, "Year"], amc_summary[sub_idx, "DDD"], type = "o", pch = pch.types[mytype], lty = lty.types[mytype], las = 1, col = col.types[mytype],
             xlim = c(1997, 2016), ylim = c(0, ymax), xlab = "Year", ylab = "DDD", bty = "n", main = mycomb, cex.main = 0.8)
  }
  legend("topleft", legend = amc.types, col = col.types, pch = pch.types, lty = lty.types, bty = "n")
}

file.name = paste("output/supporting/", 'all_consumption_data.pdf',sep = '')
pdf(file.name, paper = "a4", width = 0, height = 0)
par(mfrow = c(4, 3))
for(mycomb in unique(amc_summary$combC_2)) plot_amc(mycomb)
dev.off()
