rm(list = ls())
library(dplyr)
library(DescTools)
set.seed(1234)

amr_summary <- read.csv(file = "data/summary_AMR_filtered.csv")
amr_summary <- amr_summary %>% arrange(Year)

# focus on INPAT

patient_type = 'INPAT'; time_points = 5; observations = 30; minimum_N_R = 10

# bug drug combinations
amr_summary$bugdrug <- paste(amr_summary$Pathogen, amr_summary$Antibiotic, sep = "_")


# compute SpearmanRho between pairs of AMR trajectories across countries
# across pairs of bug-drug combinations

tmp <- c()

# loop on country (this is a bit long)

for(cc in unique(amr_summary$Country)){ # loop on country
  
  print(cc)
  a1 <- amr_summary[amr_summary$Country==cc & amr_summary$patientType==patient_type, ]
  
  bg <- unique(a1$bugdrug)
  nbg <- length(bg)
  
  # for all pairs of bug-drug combinations
  for(i in 1:(nbg-1)){
    for(j in (i+1):nbg){
      ii <- which(a1$bugdrug==bg[i])
      jj <- which(a1$bugdrug==bg[j])
      yrs <- intersect(a1[ii, "Year"], a1[jj, "Year"]) # years
      
      if(length(yrs)>=time_points){  # five years or more
        
        a1_ii <- a1[ii, ][a1[ii, "Year"] %in% yrs, ]
        a1_jj <- a1[jj, ][a1[jj, "Year"] %in% yrs, ]
        stopifnot(all(a1_ii$Year==a1_jj$Year))
        
        # Apply criterion " where all years had at least 30 isolates tested for drug resistance, and in total at least 10 resistant bacterial isolates"
        if(
          (nrow(a1_ii) >= time_points) & (all(a1_ii$N >= observations)) & (sum(a1_ii$N_I_R) >= minimum_N_R) &
          (nrow(a1_jj) >= time_points) & (all(a1_jj$N >= observations)) & (sum(a1_jj$N_I_R) >= minimum_N_R)
          ){
          rho <- SpearmanRho(a1_ii$p, a1_jj$p, use = "complete.obs",conf.level = 0.95)
          tmp <- rbind(tmp, c(cc, bg[i], bg[j], rho))
        }
      }
    }
  }
} # end of loop on country

tmp <- data.frame(tmp)
names(tmp) <- c("country", "bg1", "bg2", "rho", "rho_low", "rho_up")
head(tmp)

# extract species and antibiotics:
tmp$sp1 <- sapply(tmp$bg1, function(x) strsplit(x, "_")[[1]][1]) 
tmp$sp2 <- sapply(tmp$bg2, function(x) strsplit(x, "_")[[1]][1]) 
tmp$ab1 <- sapply(tmp$bg1, function(x) strsplit(x, "_")[[1]][2]) 
tmp$ab2 <- sapply(tmp$bg2, function(x) strsplit(x, "_")[[1]][2]) 

# create columns same species and same ab
tmp$ss <- tmp$sp1 == tmp$sp2
tmp$sa <- tmp$ab1 == tmp$ab2

# transform the rhos in numeric
for(nn in c( "rho", "rho_low", "rho_up")) tmp[, nn] <- as.numeric(as.character(tmp[, nn]))

dim(tmp)
# 83429    12
table(tmp$country)

table(same_antibio = tmp$sa, same_species = tmp$ss)

# distribution of rho:
# we compare rho across countries, combinations, for same species, different AB:

t.test(tmp$rho[!tmp$sa & !tmp$ss]) # not same AB, different species
hist(tmp$rho[!tmp$sa & !tmp$ss])

t.test(tmp$rho[tmp$ss]) # same species, different AB
hist(tmp$rho[tmp$ss])

t.test(tmp$rho[tmp$sa]) # same AB, different species
hist(tmp$rho[tmp$sa])

summary(lm0 <- lm(rho ~ ss + sa, data = tmp))
car::Anova(lm0)

plot_vec <- function(vec_to_plot, bin_width, xpos, width_scaling_factor = 0.2, point.cex = 1, plot_median =  T){
  # function to draw a nice violin plot
  if(plot_median){
    med <- median(vec_to_plot, na.rm = T)
    print(med)
    segments(y0 = med, y1 = med, x0 = xpos-width_scaling_factor, x1 = xpos+width_scaling_factor, col = "black", lwd = 3)
  }
  bins <- hist(vec_to_plot, breaks = seq(-1, 1, bin_width), plot = F)
  max_d <- max(bins$density)
  
  for(i in 1:length(bins$mids)){ # for each bin
    sub <- which(vec_to_plot >= bins$mids[i]-bin_width / 2 & vec_to_plot < bins$mids[i]+bin_width/2)
    width <- bins$density[i] / max_d * width_scaling_factor
    points(runif(n = length(sub), min = xpos-width, max = xpos+width), vec_to_plot[sub], pch = 20, col = "gray", cex = point.cex)
  }
  segments(y0 = med, y1 = med, x0 = xpos-width_scaling_factor, x1 = xpos+width_scaling_factor, col = "black", lwd = 3)
}

pdf("supporting_analysis/output/Correlation_AMR_trajectories.pdf", width = 6, height = 6)
par(mar = c(1,4,4,1), xpd = T)
plot(NULL, xlab = "", ylab = "Correlation coefficient", axes = F, xlim = c(0,3), ylim = c(-1, 1))
axis(side = 2, at = seq(-1, 1, 0.2), las = 1)
segments(x0 = 0, x1 = 3, y0 = 0, y1 = 0, lty = 2)

# t.test(tmp$rho[!tmp$sa & !tmp$ss]) # not same AB, different species
plot_vec(vec_to_plot = tmp$rho[!tmp$sa & !tmp$ss], bin_width = 0.1, xpos = 0.5, point.cex = 0.2)

# t.test(tmp$rho[tmp$ss]) # same species, different AB
plot_vec(vec_to_plot = tmp$rho[tmp$ss], bin_width = 0.1, xpos = 1.5, point.cex = 0.5)

# t.test(tmp$rho[tmp$sa]) # same AB, different species
plot_vec(vec_to_plot = tmp$rho[tmp$sa], bin_width = 0.1, xpos = 2.5, point.cex = 0.5)

text(x = 0.5, y = -1, "Different bug\nDifferent drug", cex = 0.8)
text(x = 1.5, y = -1, "Same bug\nDifferent drug", cex = 0.8)
text(x = 2.5, y = -1, "Different bug\nSame drug", cex = 0.8)
text(x = 1.4, y = 1.2, "Correlations in resistance trajectories across years", cex = 1.5)
dev.off()

# Plot legend:
plot_legend <- "
We investigated the correlation between pairs of AMR trajectories.
For each country, we considered all pairs of bug-drug combinations.
As for the analysis of resistance trajectories, we focused on combinations where we had at least five overlapping years, all years had at least 30 isolates tested for drug resistance, and in total at least 10 resistant bacterial isolates.
For each pair of trajectories, we correlated resistance across years on the overlapping years.
We show these correlations for pairs of different species, different drugs (left);
for pairs of same species, different drugs (middle);
and for pairs of different species, same drug (right);
Each point is the correlation coefficient for a country and pair of bug-drug combination.
The horizontal line shows the median.
The measure of correlation was Spearman rank correlation coefficient.
"



