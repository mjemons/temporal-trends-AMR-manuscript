
#################################### ANALYSIS CAN START HERE, BY READING THE FILES WITH ALL INFORMATION ####################################
rm(list = ls())

library(dotenv)

load_dot_env(file = ".env")

setwd(Sys.getenv(c("HOME_DIR")))

#source("scripts/source_epidata.R")
source(paste0(Sys.getenv(c("HOME_DIR")),"scripts/source_epidata.R"))
rm(list = ls()[!ls() %in% c("all_combR", "amr_summary")])

df_coefficients <- read.csv(file = "data/smooth_splines_coefficients.csv")
df_ci <- read.csv(file = "data/smooth_splines_coefficients_ci.csv")

# STUDY THE SHAPE OF THE INITIAL INCREASE IN RESISTANCE:

summary(df_coefficients$low_frequency_0d) # ALWAYS < 0.05 by definition of low frequency
summary(df_coefficients$low_frequency_1d)
summary(df_coefficients$low_frequency_2d)
t.test(df_coefficients$low_frequency_1d) # -> at low frequency, increase in frequency 0.0033 per years
t.test(df_coefficients$low_frequency_2d) # at low frequency, significant positive (accelerating) curve


table(df_coefficients$low_frequency_1d_lowerCI > 0, useNA = "ifany") # 45 significant rise
table(df_coefficients$low_frequency_1d_upperCI < 0, useNA = "ifany") # 25 significant decline

table(df_coefficients$low_frequency_2d_lowerCI > 0, useNA = "ifany") # 25 significant acceleration
table(df_coefficients$low_frequency_2d_upperCI < 0, useNA = "ifany") # 24 significant deceleration

# how does initial increase depends on species, country, ab?
summary(lm0 <- lm(low_frequency_1d ~ species + country + ab, data = df_coefficients))
car::Anova(lm0)

summary(lm1 <- lm(low_frequency_2d ~ species + country + ab, data = df_coefficients))
car::Anova(lm1)

# figure first derivative, second derivative by country, by AB, by species
pdf("output/initial_increase_bycountry_species_ab.pdf", width = 6, height = 8)
par(mfrow = c(3,1), mar = c(4,4,1,1), las = 1)
ylims <- c(-0.01, 0.03)
ylabs <- seq(-0.01, 0.03, 0.01)
yposlabels <- 0.025
plot(NULL, xlim = c(0, 17), ylim = ylims, axes = F, xlab = "", ylab = "initial increase in frequency per year")
axis(2, at = ylabs)
abline(h = 0)
idx <- 1
for(cc in unique(df_coefficients$country)){
  all_coeff <- df_coefficients$low_frequency_1d[df_coefficients$country == cc]
  if(sum(!is.na(all_coeff)) > 2){ # at least two ab * species combinations
    y <- mean(all_coeff, na.rm = T)
    ci <- 1.96 * sd(all_coeff, na.rm = T) / sum(!is.na(all_coeff))
    points(x = idx, y = y, pch = 20)
    segments(x0 = idx, y0 = y - ci, x1 = idx, y1 = y + ci)
    text(x = idx, y = yposlabels, cc, srt = 45)
    idx <- idx + 1
  }
}

plot(NULL, xlim = c(0.5, 8), ylim = ylims, axes = F, xlab = "", ylab = "initial increase in frequency per year")
axis(2, at = ylabs)
abline(h = 0)
idx <- 1
for(sp in unique(df_coefficients$species)){
  all_coeff <- df_coefficients$low_frequency_1d[df_coefficients$species == sp]
  if(sum(!is.na(all_coeff)) > 2){ # at least two ab * species combinations
    y <- mean(all_coeff, na.rm = T)
    ci <- 1.96 * sd(all_coeff, na.rm = T) / sum(!is.na(all_coeff))
    points(x = idx, y = y, pch = 20)
    segments(x0 = idx, y0 = y - ci, x1 = idx, y1 = y + ci)
    text(x = idx, y = yposlabels, sp, srt = 45)
    idx <- idx + 1
  }
}
plot(NULL, xlim = c(1, 17), ylim = ylims, axes = F, xlab = "", ylab = "initial increase in frequency per year")
axis(2, at = ylabs)
abline(h = 0)
idx <- 1
for(ab in unique(df_coefficients$ab)){
  all_coeff <- df_coefficients$low_frequency_1d[df_coefficients$ab == ab]
  if(sum(!is.na(all_coeff)) > 2){ # at least two ab * species combinations
    y <- mean(all_coeff, na.rm = T)
    ci <- 1.96 * sd(all_coeff, na.rm = T) / sum(!is.na(all_coeff))
    points(x = idx, y = y, pch = 20)
    segments(x0 = idx, y0 = y - ci, x1 = idx, y1 = y + ci)
    text(x = idx, y = yposlabels, ab, srt = 45)
    idx <- idx + 1
  }
}
dev.off()

# figure second derivative, second derivative by country, by AB, by species
pdf("output/initial_acceleration_bycountry_species_ab.pdf", width = 6, height = 8)
par(mfrow = c(3,1), mar = c(4,4,1,1), las = 1)
ylims <- c(-0.001, 0.012)
ylabs <- seq(0, 0.01, 0.002)
yposlabels <- 0.01
plot(NULL, xlim = c(0, 17), ylim = ylims, axes = F, xlab = "", ylab = "initial acceleration per year")
axis(2, at = ylabs)
abline(h = 0)
idx <- 1
for(cc in unique(df_coefficients$country)){
  all_coeff <- df_coefficients$low_frequency_2d[df_coefficients$country == cc]
  if(sum(!is.na(all_coeff)) > 2){ # at least two ab * species combinations
    y <- mean(all_coeff, na.rm = T)
    ci <- 1.96 * sd(all_coeff, na.rm = T) / sum(!is.na(all_coeff))
    points(x = idx, y = y, pch = 20)
    segments(x0 = idx, y0 = y - ci, x1 = idx, y1 = y + ci)
    text(x = idx, y = yposlabels, cc, srt = 45)
    idx <- idx + 1
  }
}

plot(NULL, xlim = c(0.5, 8), ylim = ylims, axes = F, xlab = "", ylab = "initial acceleration per year")
axis(2, at = ylabs)
abline(h = 0)
idx <- 1
for(sp in unique(df_coefficients$species)){
  all_coeff <- df_coefficients$low_frequency_2d[df_coefficients$species == sp]
  if(sum(!is.na(all_coeff)) > 2){ # at least two ab * species combinations
    y <- mean(all_coeff, na.rm = T)
    ci <- 1.96 * sd(all_coeff, na.rm = T) / sum(!is.na(all_coeff))
    points(x = idx, y = y, pch = 20)
    segments(x0 = idx, y0 = y - ci, x1 = idx, y1 = y + ci)
    text(x = idx, y = yposlabels, sp, srt = 45)
    idx <- idx + 1
  }
}
plot(NULL, xlim = c(1, 17), ylim = ylims, axes = F, xlab = "", ylab = "initial acceleration per year")
axis(2, at = ylabs)
abline(h = 0)
idx <- 1
for(ab in unique(df_coefficients$ab)){
  all_coeff <- df_coefficients$low_frequency_2d[df_coefficients$ab == ab]
  if(sum(!is.na(all_coeff)) > 2){ # at least two ab * species combinations
    y <- mean(all_coeff, na.rm = T)
    ci <- 1.96 * sd(all_coeff, na.rm = T) / sum(!is.na(all_coeff))
    points(x = idx, y = y, pch = 20)
    segments(x0 = idx, y0 = y - ci, x1 = idx, y1 = y + ci)
    text(x = idx, y = yposlabels, ab, srt = 45)
    idx <- idx + 1
  }
}
dev.off()


# STUDY THE PLATEAU:
summary(glm0 <- nnet::multinom(is_plateau ~ species + country + ab, data = df_coefficients))
car::Anova(glm0)


# REDO FIGURE, INDICATING THE INFERRED SECOND-ORDER POLYNOMIAL, AND WHETHER A PLATEAU IS INFERRED OR NOT ON THE FIGURE
all_combR <- sort(unique(amr_summary$combR))
pdf("output/all_temporal_trends.pdf", paper = "a4", width = 0, height = 0)
par(mfrow = c(4, 3))
par(las = 1)
for(i in 1:length(all_combR)){
  mycomb <- all_combR[i]
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  tmp2 <- tmp[which(tmp$patientType=="OUTPAT"),]
  if(any(duplicated(tmp1$Year))) stop()
  if(any(duplicated(tmp2$Year))) stop()
  is_plot <- F
  if((nrow(tmp1) > 5) & (all(tmp1$N > 100))){ # at least 5 years & all points > 100
    print(mycomb)
    plot(tmp1$Year, tmp1$p, col = "red", ylim = c(0,1), xlim = c(1998, 2018), pch = 20, type = "o", ylab = "frequency", xlab = "year", main = mycomb)
    segments(x0 = tmp1$Year, y0 = tmp1$p_min, x1 = tmp1$Year, y1 = tmp1$p_max, col = 'red')
    spl <- smooth.spline(x = tmp1$Year, y = tmp1$p, w = 1 / tmp1$N, df = 5)
    lines(spl$x, spl$y, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5), lwd = 4)
    is_plot <- T
  }
  if(mycomb %in% df_coefficients$combR){
    is_plateau <- df_coefficients$is_plateau[which(df_coefficients$combR == mycomb)]
    if(!is.na(is_plateau)){
      text(x = 2007.5, y = 0.9, labels = is_plateau)
    }
  }
  if((nrow(tmp2) > 5) & (all(tmp2$N > 100))){ # at least 5 years for one of them
    if(is_plot) plot_fun <- points else plot_fun <- plot
    plot_fun(tmp2$Year+0.1, tmp2$p, col = "blue", pch = 20, type = "o", ylim = c(0,1), xlim = c(1998, 2018), ylab = "frequency", xlab = "year", main = mycomb)
    segments(x0 = tmp2$Year+0.1, y0 = tmp2$p_min, x1 = tmp2$Year+0.1, y1 = tmp2$p_max, col = 'blue')
    spl <- smooth.spline(x = tmp2$Year, y = tmp2$p, w = 1 / tmp2$N, df = 5)
    lines(spl$x, spl$y, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5), lwd = 4)
  }
}
dev.off()


