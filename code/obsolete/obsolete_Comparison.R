library(RColorBrewer)
library(lme4)
library(dotenv)

load_dot_env(file = ".env")

setwd(Sys.getenv(c("HOME_DIR")))

all_combR <- sort(unique(amr_summary$combR))

pdf("output/AllCourses.pdf", paper = "a4", width = 0, height = 0)
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
  if((nrow(tmp1) > 5) & (all(tmp1$N > 30))){ # at least 5 years & all points > 100
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
  if(mycomb %in% df_coefficients$combR){
    is_trend <- df_coefficients$logistic.summary[which(df_coefficients$combR == mycomb)]
    if(!is.na(is_trend)){
      text(x = 2007.5, y = 0.8, labels = is_trend)
    }
  }
  if(mycomb %in% df_coefficients$combR){
    is_K <- df_coefficients$logistic.plateau.K[which(df_coefficients$combR == mycomb)]
    if(!is.na(is_trend)){
      text(x = 2007.5, y = 0.7, labels = is_K)
    }
  }
  
  if((nrow(tmp2) > 5) & (all(tmp2$N > 30))){ # at least 5 years for one of them
    if(is_plot) plot_fun <- points else plot_fun <- plot
    plot_fun(tmp2$Year+0.1, tmp2$p, col = "blue", pch = 20, type = "o", ylim = c(0,1), xlim = c(1998, 2023), ylab = "frequency", xlab = "year", main = mycomb)
    segments(x0 = tmp2$Year+0.1, y0 = tmp2$p_min, x1 = tmp2$Year+0.1, y1 = tmp2$p_max, col = 'blue')
    spl <- smooth.spline(x = tmp2$Year, y = tmp2$p, w = 1 / tmp2$N, df = 5)
    lines(spl$x, spl$y, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5), lwd = 4)
  }
}
dev.off()

all_combR <- sort(unique(amr_summary$combR))
Failed = df_coefficients[df_coefficients$logistic.plateau == 'NA',]$combR
pdf("output/testingFitFail.pdf", paper = "a4", width = 0, height = 0)
par(mfrow = c(4, 3))
par(las = 1)

for(i in 1:length(Failed)){
  mycomb <- Failed[i]
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  tmp2 <- tmp[which(tmp$patientType=="OUTPAT"),]
  if(any(duplicated(tmp1$Year))) stop()
  if(any(duplicated(tmp2$Year))) stop()
  is_plot <- F
  if((nrow(tmp1) > 5) & (all(tmp1$N > 30))){ # at least 5 years & all points > 100
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
  if(mycomb %in% df_coefficients$combR){
    is_trend <- df_coefficients$logistic.trend[which(df_coefficients$combR == mycomb)]
    is_plateau <- df_coefficients$logistic.plateau[which(df_coefficients$combR == mycomb)]
    if(!is.na(is_trend)){
      text(x = 2007.5, y = 0.8, labels = paste(is_trend,is_plateau,sep=' '))
    }
  }
  
  if((nrow(tmp2) > 5) & (all(tmp2$N > 30))){ # at least 5 years for one of them
    if(is_plot) plot_fun <- points else plot_fun <- plot
    plot_fun(tmp2$Year+0.1, tmp2$p, col = "blue", pch = 20, type = "o", ylim = c(0,1), xlim = c(1998, 2023), ylab = "frequency", xlab = "year", main = mycomb)
    segments(x0 = tmp2$Year+0.1, y0 = tmp2$p_min, x1 = tmp2$Year+0.1, y1 = tmp2$p_max, col = 'blue')
    spl <- smooth.spline(x = tmp2$Year, y = tmp2$p, w = 1 / tmp2$N, df = 5)
    lines(spl$x, spl$y, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5), lwd = 4)
  }
}
dev.off()


# create combinations of pathogens, antibiotic: (not countries)
amr_summary$combR_ <- paste(amr_summary$Pathogen, amr_summary$Antibiotic, sep = "|")
all_combR <- sort(unique(amr_summary$combR_))

# create color palette for countries
col_countries <- rep(NA, length(unique(amr_summary$Country)))
names(col_countries) <- sort(unique(amr_summary$Country))

col_countries[1:29] <- c(brewer.pal(n = 12, name = "Paired"), brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 5, name = "Set1"))

pdf("output/all_temporal_trends_allcountries.pdf", paper = "a4", width = 0, height = 0)
par(mfrow = c(4, 3))
par(las = 1)
for(i in 1:length(all_combR)){
  mycomb <- all_combR[i]
  print(mycomb)
  tmp <- amr_summary[which(amr_summary$combR_ == mycomb),]
  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  tmp2 <- tmp[which(tmp$patientType=="OUTPAT"),]
  tmp1 <- tmp1[tmp1$N > 10, ]
  tmp2 <- tmp2[tmp2$N > 10, ]
  countries <- unique(tmp1$Country)
  if(nrow(tmp1) > 5){
    idx <- 1
    for(cc in countries){
      if(idx==1) plot_fun <- plot else plot_fun <- points
      plot_fun(tmp1$Year[tmp1$Country==cc], tmp1$p[tmp1$Country==cc], col = col_countries[cc], ylim = c(0,1), xlim = c(1998, 2020), pch = 20, type = "o", ylab = "frequency", xlab = "year", main = mycomb)
      #lines(smooth.spline(x=))
      if(length(tmp1$Year[tmp1$Country==cc])>5){
        spl <- smooth.spline(x = tmp1$Year[tmp1$Country==cc], y =tmp1$p[tmp1$Country==cc], w = 1 / tmp1$N[tmp1$Country==cc], df = 5)
        lines(spl$x, spl$y, col=col_countries[cc], lwd = 3)
        idx <- idx + 1
      }
    }
  }
  if(i == 2){
    for(j in 1:29){
      rect(xleft = 1998, xright = 1999, ybottom = (j-1)/29, ytop = j/29, col = col_countries[j])
      text(x = 2000, y = (j-0.5)/29, names(col_countries)[j], adj = 0, cex = 0.6)
    }
  }
}

dev.off()
