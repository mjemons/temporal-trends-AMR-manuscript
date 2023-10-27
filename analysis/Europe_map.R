rm(list = ls())
# library(dplyr); library(plyr); library(rworldmap); library(ggplot2); setwd("~/ownCloud/AMR_Sonja_Martin/temporal_trends_AMR/")
# Get the world map
worldMap <- getMap()

# Member States of the European Union
Countries <- c("Austria","Belgium","Bulgaria","Croatia",
                   "Cyprus",
                   "Czech Rep.","Denmark","Estonia",
                   "Finland",
                   "France",
                   "Germany","Greece","Hungary","Ireland","Italy","Latvia",
                   "Lithuania",
                   "Luxembourg",
                   "Malta","Netherlands", "Norway", "Poland",
                   "Portugal","Romania","Slovakia",
                   "Slovenia",
                   "Spain",
                   "Sweden",
                   "United Kingdom")
# Select only the index of states member of the E.U.
ind <- which(worldMap$NAME%in%Countries)

# Extract longitude and latitude border's coordinates of members states of E.U. 
europeCoords <- lapply(ind, function(i){
  df <- data.frame(worldMap@polygons[[i]]@Polygons[[1]]@coords)
  df$region =as.character(worldMap$NAME[i])
  colnames(df) <- list("long", "lat", "region")
  return(df)
})

europeCoords <- do.call("rbind", europeCoords)

# Add resistance trends (binomial model) data data for each member
binomial_coeff <- read.csv("data/results/binomial_model.csv")

#get rid off anything that is not country
binomial_coeff <- binomial_coeff %>% filter(!grepl(c("Pathogen"),term))
binomial_coeff <- binomial_coeff %>% filter(!grepl(c("Antibiotic"),term))

# # FB added Slovenia again with the smallest coefficient (arbitrarily, -1.5)
# binomial_coeff$X <- NULL
# #binomial_coeff <- rbind(binomial_coeff, c("CountrySlovenia", 0, 0, 0, 0, 0, 0, "Country", "Slovenia"))
europeCoords$coeff <- binomial_coeff$estimate[match(europeCoords$region, binomial_coeff$local_term)]
europeCoords$coeff[europeCoords$region=="Czech Rep."] <- binomial_coeff$estimate[binomial_coeff$local_term=="Czech Republic"]

# FB added the reference (Austria) at:
europeCoords$coeff[europeCoords$region=="Austria"] <- 0.
europeCoords[europeCoords$coeff < -1.2, "coeff"] <- -1.5 # arbitrarily set at low value for Luxembourg and Slovenia

# Plot the map
P <- ggplot() + geom_polygon(data = europeCoords, aes(x = long, y = lat, group = region, fill = coeff),
                             colour = "black", linewidth = 0.1) +
  coord_map(xlim = c(-13, 35),  ylim = c(32, 71))

P <- P + scale_fill_gradient(name = "Coefficient", low = "#FFFF00FF", high = "#FF0000FF", na.value = "grey50") # change color gradient to yellow to red


P <- P + theme(#panel.grid.minor = element_line(colour = NA), panel.grid.minor = element_line(colour = NA),
  panel.background = element_rect(fill = NA, colour = NA),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  #rect = element_blank(),
  plot.margin = unit(0 * c(-1.5, -1.5, -1.5, -1.5), "lines")) # remove axes

pdf("output/europe_map_binomial.pdf")
plot(P)
dev.off()

# FB ADDED ATTEMPT CORRELATING AB USE IN 2000 WITH BINOMIAL COEFFICIENT
amc_summary <- read.csv("data/summary_AMC_byclass_improved.csv")

ddply(amc_summary[amc_summary$Sector=="Community" & amc_summary$Antimicrobial.Type=="AllNew", ], .(Antibiotic_class), summarise, sum_DDD = sum(DDD, na.rm = T))
.Last.value$Antibiotic_class[.Last.value$sum_DDD > 500]
u_classes <- c("J01A", "J01C", "J01D", "J01F", "J01M") # only keep those that are consumed enough
years <- seq(1998, 2016)

plot_DDD_AMRtrend_correlation <- function(myselection, year, class, country_labels = T, return_tab = F, xlim = NULL, xlab = NULL){
  
  # plot the DDD - AMR trend correlation for a selection of consumption (amc_summary) data
  tmp <- ddply(amc_summary[myselection,], .(Country), summarise, sum_DDD = sum(DDD, na.rm = T))
  tmp <- tmp[tmp$Country=="Austria" | (tmp$Country %in% binomial_coeff$local_term), ]
  tmp <- tmp[!tmp$Country %in% c("Luxembourg", "Slovenia"), ]
  idx <- match(tmp$Country, binomial_coeff$local_term)
  tmp$coeff <- binomial_coeff$estimate[idx]
  tmp$coeff[tmp$Country=="Austria"] <- 0 # this is the reference country
  minDDD <- 0.99 * floor(min(tmp$sum_DDD))
  maxDDD <- 1.01 * ceiling(max(tmp$sum_DDD))
  xDDD <- seq(minDDD, maxDDD, length.out = 100)
  if(is.null(xlim)) xlim <-  c(minDDD, maxDDD)
  if(is.null(xlab)) xlab <-  paste("Consumption of", class, "in", year)
  plot(NULL, pch = 20, las = 1, xlim = xlim, ylim = c(-1.5, 1.5),
       xlab = xlab, ylab = "Coefficient (fraction increasing trend)", axes = F)
  
  ################################
  lm0 <- lm(coeff ~ sum_DDD, data = tmp)
  ################################
  
  ci <- predict(lm0, interval = "confidence", newdata = data.frame(sum_DDD=xDDD))
  print("Summary of linear model:")
  print(summary(lm0))
  print(confint(lm0))
  polygon(x = c(xDDD, rev(xDDD)), y = c(ci[,"lwr"], rev(ci[, 'upr'])), col = "gray", border = NA)
  points(x = xDDD, y = ci[, "fit"], type = "l", lwd= 2)
  points(tmp$sum_DDD, tmp$coeff, pch = 20)
  idx_sel <- match(c("France", "Italy", "Greece"), tmp$Country)
  if(country_labels) text(tmp$sum_DDD[idx_sel], tmp$coeff[idx_sel]+0.02, c("FR", "IT", "GR"), cex = 0.8, adj = 0)
  
  axis(side = 1, at = seq(0, xlim[2], length.out = 5))
  axis(side = 2, at = seq(-1.5, 1.5, 0.5), las = 1)
  
  if(return_tab) return(tmp)
}

# correlation of increasing trends in AMR with level of antibiotic consumption
for(mysector in c("Community", "Hospital Sector")){
  if(mysector=="Community") name <- paste0("output/DDD_increasingTrend_correlation_", "Community.pdf") else name <- paste0("output/DDD_increasingTrend_correlation_", "HS.pdf")
  pdf(name, width = 4*1.5, height = 3*1.5)
  par(mfrow = c(1,1), mar = c(4,4,1,1))
  # AMC data to select:
  selection <- amc_summary$Sector==mysector & amc_summary$Antimicrobial.Type=="AllNew"
  plot_DDD_AMRtrend_correlation(selection, class = "all", year =  "all years")
  dev.off()
}


# plot the DDD-AMR trend correlation by class and by year to see the difference clearly
for(mysector in c("Community", "Hospital Sector")){
  
  if(mysector=="Community") name <- paste0("output/DDD_increasingTrend_correlation_decomposed_", "Community.pdf") else name <- paste0("output/DDD_increasingTrend_correlation_decomposed_", "HS.pdf")
  # figure:
  pdf(name, width = length(u_classes) * 4, height = length(years) * 3)
  par(mfrow = c(length(years), length(u_classes)), mar = c(4,4,1,1))
  
  for(myyear in years){
    for(cc in u_classes){
      
      # AMC data to select:
      selection <- amc_summary$Year==myyear & amc_summary$Sector==mysector & amc_summary$Antimicrobial.Type=="AllNew" & amc_summary$Antibiotic_class==cc
      plot_DDD_AMRtrend_correlation(selection, year = myyear, class = cc)
      
    }
  }
  dev.off()
}



pdf("output/consumption_AMCtrends_correlations.pdf", width = 8, height = 8)

# check temporal trends by countries
par(mfrow = c(2,2))
for(mysector in c("Community", "Hospital Sector")){
  
  matrix_trends <- matrix(NA, nrow = nrow(binomial_coeff), ncol = length(u_classes))
  matrix_intercept2000 <- matrix(NA, nrow = nrow(binomial_coeff), ncol = length(u_classes))
  ff <- function(mat){
    colnames(mat) <- u_classes
    rownames(mat) <- binomial_coeff$local_term
    mat <- data.frame(mat)
    mat$all <- NA
    mat
  }
  matrix_trends <- ff(matrix_trends)
  matrix_intercept2000 <- ff(matrix_intercept2000)
  for(country in binomial_coeff$local_term){
    #print(country)
    for(cc in u_classes){
      idx <- which(amc_summary$Country==country & amc_summary$Sector==mysector & amc_summary$Antimicrobial.Type=="AllNew" & amc_summary$Antibiotic_class==cc)
      if(length(idx) > 5){
        #cat(country, cc, "\n")
        lm0 <- lm(DDD ~ Year, data = amc_summary[idx, ])
        matrix_trends[country, cc] <- lm0$coefficients[2]
        matrix_intercept2000[country, cc] <- mean(amc_summary[idx, "DDD"][amc_summary[idx, "Year"] == 2000], na.rm = T) # lm0$coefficients[1] + lm0$coefficients[2] * 2000
      }
      #amc_summary[idx, ]
    }
    ########################################
    # and trends overall in AMC over time:
    ########################################
    idx <- which(amc_summary$Country==country & amc_summary$Sector==mysector & amc_summary$Antimicrobial.Type=="AllNew" & amc_summary$Antibiotic_class %in% u_classes)
    if(length(idx) > 5){
      amc_summary_all <- ddply(amc_summary[idx, ], .(Year), summarise, DDD = sum(DDD), classes = paste(Antibiotic_class, collapse = "_"))
      #print(table(amc_summary_all$classes))
      amc_summary_all <- amc_summary_all[amc_summary_all$classes=="J01A_J01C_J01D_J01F_J01M", ]
      
      ########################################
      lm0 <- lm(DDD ~ Year, data = amc_summary_all)
      ########################################
      matrix_trends[country, "all"] <- lm0$coefficients[2]
      #print(range(amc_summary_all$Year))
      matrix_intercept2000[country, "all"] <- mean(amc_summary_all[, "DDD"][amc_summary_all[, "Year"] == 2000], na.rm = T) # lm0$coefficients[1] + lm0$coefficients[2] * 2000
    }
  }
  
  # COUNTRIES WITH LARGE CONSUMPTION TEND TO BE THOSE WITH MORE DECREASING USE
  ylim <- c(-0.7, 0.7)
  if(mysector=="Community") xlim <- c(0, 40) else xlim <- c(0, 4)
  plot(
    x = matrix_intercept2000$all,
    y = matrix_trends$all,
    pch = 20, ylim = ylim , xlim = xlim,
    ylab = "Trend in AB use overall", xlab = "Consumption in 2000", las = 1, main = mysector
  )
  #abline(h = 2)
  abline(v = 0); abline(h = 0)
  text(x = matrix_intercept2000$all, y = matrix_trends$all, rownames(matrix_trends), cex = 0.5)
  lm1 <- lm(matrix_trends$all ~ matrix_intercept2000$all)
  abline(lm1, col = "gray", lwd = 2)
  
  for(cc in c(u_classes, "all")){
    print(cc)
    print(
      summary(lm(matrix_trends[, cc] ~ matrix_intercept2000[, cc]))
    )
  }
  
  # NO ASSOCIATION BETWEEN OVERALL TRENDS IN RESISTANCE AND TRENDS IN USE
  stopifnot(all(binomial_coeff$local_term==rownames(matrix_trends)))
  print(summary(lm(binomial_coeff$estimate ~ matrix_trends$all)))
  print(summary(lm(binomial_coeff$estimate ~ matrix_intercept2000$all)))
  fisher.test(
    table(binomial_coeff$estimate>0, matrix_trends$all>0)
  )
  plot(x = matrix_trends$all, y = binomial_coeff$estimate, pch = 20, las = 1, ylim = c(-1.5, 1.5), xlim = ylim,
       ylab = "Coefficient (fraction increasing trend)", xlab = "Trend in AB use overall", main = mysector)
  abline(h = 0); abline(v = 0)
  text(matrix_trends$all, binomial_coeff$estimate, binomial_coeff$local_term, cex = 0.5)
  # no correlation
  subset <- !binomial_coeff$local_term %in% c("Luxembourg", "Slovenia")
  lm2 <- lm(binomial_coeff$estimate[subset] ~ matrix_trends$all[subset])
  abline(lm2, col = "gray", lwd = 2)
  
  if(mysector=="Community"){
    matrix_trends_c <- matrix_trends
    matrix_intercept2000_c <- matrix_intercept2000
  }
  if(mysector=="Hospital Sector"){
    matrix_trends_hs <- matrix_trends
    matrix_intercept2000_hs <- matrix_intercept2000
  }
  
  rm(matrix_trends, matrix_intercept2000)
}
dev.off()

######################################## FINAL FIGURE ########################################
# plot together AB use overall vs. consumption 2000; trend in R vs. consumption in 2000
##############################################################################################

pdf("output/Trends_AMR_AMC_consumption_correlations.pdf", width = 8, height = 8)
# check temporal trends by countries
par(mfrow = c(2,2), xpd=TRUE)
for(mysector in c("Community", "Hospital Sector")){
  
  # COUNTRIES WITH LARGE CONSUMPTION IN COMMUNITY TEND TO BE THOSE WITH MORE DECREASING USE IN HOSPITAL NOT COMMUNITY
  # (NOTE that such a correlation exist when considering hospital sector)
  
  if(mysector=="Community"){
    xlim <- c(0, 40); ylim <- c(-0.7, 0.7)
    matrix_trends <- matrix_trends_c
    matrix_intercept2000 <- matrix_intercept2000_c
  } else {
    xlim <- c(0, 40); ylim <- c(-0.1, 0.1)
    matrix_trends <- matrix_trends_hs
    matrix_intercept2000 <- matrix_intercept2000_c # as per Sonja's suggestion
  }
  plot(
   NULL, ylim = ylim , xlim = xlim,
    ylab = "Trend in total antibiotic use", xlab = "Total community consumption in 2000", las = 1, main = mysector
  )
  #text(x = matrix_intercept2000$all, y = matrix_trends$all, rownames(matrix_trends), cex = 0.5)
  y <- matrix_trends$all
  x <- matrix_intercept2000$all
  lm1 <- lm(y ~ x) 
  #segments(x0 = 0, x1 = 40, y0 = lm1$coefficients[1], y1 = lm1$coefficients[1] + 40*lm1$coefficients[2], col = "gray", lwd = 2)
  
  ci1 <- predict(lm1, interval = "confidence", newdata = data.frame(x = seqx <- seq(xlim[1], xlim[2], length.out = 100)))
  polygon(x = c(seqx, rev(seqx)), y = c(ci1[,"lwr"], rev(ci1[, 'upr'])), col = "gray", border = NA)
  points(x = seqx, y = ci1[, "fit"], type = "l", lwd= 2)
  points( x = matrix_intercept2000$all,
          y = matrix_trends$all,
          pch = 20)
  segments(x0 = xlim[1], x1 = xlim[2], y0 = 0, y1 = 0)
  if(mysector == "Community") text(xlim[1], ylim[2] + 0.2 * (ylim[2]-ylim[1]), "A", cex = 4) else text(xlim[1], ylim[2] + 0.2 * (ylim[2]-ylim[1]), "B", cex = 4)
  
  # plot names of certain countries
  idx_sel <- match(c("France", "Italy", "Greece"), rownames(matrix_trends))
  text(x = matrix_intercept2000$all[idx_sel],
       y = matrix_trends$all[idx_sel],
       labels = c("FR", "IT", "GR"))
  
  # print correlation for each class:
  for(cc in c(u_classes, "all")){
    print(cc)
    print(
      summary(lm(matrix_trends[, cc] ~ matrix_intercept2000[, cc]))
    )
    cat("\n\n")
  } # this correlation holds true for:
  # all classes, in HS
  # somme classes in community; not J01M; not all; pattern less strong in community
}

# 2nd part of the plot: correlation between increasing trends and total consumption in 2000
myyear <- 2000
for(mysector in c("Community")){ # only use community as per Sonja's suggestion
  if(mysector=="Community") xlim <- c(0, 40) else xlim <- c(0, 4)
  ylim <- c(-1.5, 1.5)
  amc_selection <- amc_summary$Sector==mysector & # use community consumption as per Sonja's suggestion
    amc_summary$Antimicrobial.Type=="AllNew" & amc_summary$Year==myyear & amc_summary$Antibiotic_class %in% u_classes
  plot_DDD_AMRtrend_correlation(myselection = amc_selection, class = "all", year =  "2000", country_labels = T, xlim = xlim, xlab = "Total community consumption in 2000")
  segments(x0 = xlim[1], x1 = xlim[2], y0 = 0, y1 = 0)
  if(mysector == "Community") {
    text(xlim[1], ylim[2] + 0.2 * (ylim[2]-ylim[1]), "C", cex = 4) 
  } else {
    text(xlim[1], ylim[2] + 0.2 * (ylim[2]-ylim[1]), "D", cex = 4)
  }
}

dev.off()

