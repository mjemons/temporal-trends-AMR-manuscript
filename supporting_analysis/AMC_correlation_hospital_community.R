rm(list = ls())

source("code/Utils.R")
set.seed(1234)

# load AMC improved trajectories
amc_summary <- read.csv(file = "data/summary_AMC_byclass_improved.csv")

# check if community/hospital use are correlated across countries
# and if they are correlated within countries, across time.

dim(amc_summary)
head(amc_summary)
table(amc_summary$Antibiotic_class)

mytype <- "AllNew"

tmp <- c()

for(yr in unique(amc_summary$Year)){
  for(myclass in unique(amc_summary$Antibiotic_class)){
  
    crit <- amc_summary$Year==yr & amc_summary$Antibiotic_class == myclass & amc_summary$Antimicrobial.Type ==mytype
    idx1 <- which(crit & amc_summary$Sector == "Community")
    idx2 <- which(crit & amc_summary$Sector == "Hospital Sector")
    ccs <- intersect(amc_summary[idx1, "Country"], amc_summary[idx2, "Country"])
    
    if(length(ccs)>= 5){ # at least five countries
      a1 <- amc_summary[idx1, ][amc_summary[idx1, "Country"] %in% ccs, ]
      a2 <- amc_summary[idx2, ][amc_summary[idx2, "Country"] %in% ccs, ]
      stopifnot(nrow(a1)==nrow(a2))
      
      # add correlation to the table:
      tmp <- rbind(tmp, 
                   c(yr, myclass, SpearmanRho(a1$DDD, a2$DDD, use = "complete.obs",conf.level = 0.95))
      )
    }
  }
} 

tmp <- data.frame(tmp)
names(tmp) <- c("year", "class", "rho", "rho_low", "rho_up")
# transform the rhos in numeric
for(nn in c( "rho", "rho_low", "rho_up")) tmp[, nn] <- as.numeric(as.character(tmp[, nn]))
tmp -> df_spatial
hist(df_spatial$rho)

# now compute correlation across years
tmp <- c()

for(cc in unique(amc_summary$Country)){
  for(myclass in unique(amc_summary$Antibiotic_class)){
    
    crit <- amc_summary$Country==cc & amc_summary$Antibiotic_class == myclass & amc_summary$Antimicrobial.Type ==mytype
    idx1 <- which(crit & amc_summary$Sector == "Community")
    idx2 <- which(crit & amc_summary$Sector == "Hospital Sector")
    yrs <- intersect(amc_summary[idx1, "Year"], amc_summary[idx2, "Year"])
    
    if(length(yrs)>= 5){ # at least five countries
      a1 <- amc_summary[idx1, ][amc_summary[idx1, "Year"] %in% yrs, ]
      a2 <- amc_summary[idx2, ][amc_summary[idx2, "Year"] %in% yrs, ]
      stopifnot(nrow(a1)==nrow(a2))
      
      # add correlation to the table:
      tmp <- rbind(tmp, 
                   c(cc, myclass, SpearmanRho(a1$DDD, a2$DDD, use = "complete.obs",conf.level = 0.95))
      )
    }
  }
} 
tmp <- data.frame(tmp)
names(tmp) <- c("year", "class", "rho", "rho_low", "rho_up")
# transform the rhos in numeric
for(nn in c( "rho", "rho_low", "rho_up")) tmp[, nn] <- as.numeric(as.character(tmp[, nn]))
tmp -> df_temporal

hist(df_temporal$rho)

mean(df_temporal$rho, na.rm = T)
mean(df_spatial$rho, na.rm = T)

median(df_temporal$rho, na.rm = T)
median(df_spatial$rho, na.rm = T)

plot_vec <- function(vec_to_plot, bin_width, xpos, width_scaling_factor = 0.2, plot_median =  T){
  # function to draw a nice violin plot
  med <- median(vec_to_plot, na.rm = T)
  if(plot_median){
    print(med)
    segments(y0 = med, y1 = med, x0 = xpos-width_scaling_factor, x1 = xpos+width_scaling_factor, col = "black", lwd = 3)
  }
  bins <- hist(vec_to_plot, breaks = seq(-1, 1, bin_width), plot = F)
  max_d <- max(bins$density)
  
  for(i in 1:length(bins$mids)){ # for each bin
    sub <- which(vec_to_plot >= bins$mids[i]-bin_width / 2 & vec_to_plot < bins$mids[i]+bin_width/2)
    width <- bins$density[i] / max_d * width_scaling_factor
    points(runif(n = length(sub), min = xpos-width, max = xpos+width), vec_to_plot[sub], pch = 20, col = "gray")
  }
  if(plot_median){
    med <- median(vec_to_plot, na.rm = T)
    print(med)
    segments(y0 = med, y1 = med, x0 = xpos-width_scaling_factor, x1 = xpos+width_scaling_factor, col = "black", lwd = 3)
  }
}
pdf("output/Correlation_use_hospital_community.pdf", width = 6, height = 6)
par(mar = c(1,4,4,1), xpd = T)
plot(NULL, xlab = "", ylab = "Correlation coefficient", axes = F, xlim = c(0,2), ylim = c(-1, 1))
axis(side = 2, at = seq(-1, 1, 0.2), las = 1)
segments(x0 = 0, x1 = 2, y0 = 0, y1 = 0, lty = 2)

plot_vec(vec_to_plot = df_temporal$rho, bin_width = 0.1, xpos = 0.5)
plot_vec(vec_to_plot = df_spatial$rho, bin_width = 0.1, xpos = 1.5)

text(x = 0.5, y = -1, "Correlations across years")
text(x = 1.5, y = -1, "Correlations across countries")
text(x = 1, y = 1.2, "Correlations in antibiotic use\nin hospital vs. community", cex = 1.5)
dev.off()

# Plot legend:
plot_legend <- "
We investigated the correlation between hospital use and community use.
First, we looked at the correlation across years.
For each drug-country combination, we correlated the anbiotic use in the hospital sector with that of the community.
Second, we looked at the correlation across countries.
For each drug-year combination, we correlated the anbiotic use in the hospital sector with that of the community.
We show both correlations across years (left) and correlations across countries (right).
Each point is the correlation coefficient for a drug-country or a drug-year combination.
The horizontal line shows the median.
In all cases, the measure of correlation was Spearman rank correlation coefficient.
"


#vec_to_plot <- df_temporal$rho
#bin_width <- 0.1 # bin width





