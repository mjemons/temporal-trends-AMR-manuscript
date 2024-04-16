set.seed(1234)
source('code/Utils.R')

# looking at the distribution of fits
# plotting the timecourses of AMR over time, by country, pathogen, drug
# both on p scale and logit(p) scale

# coefficients summary
coefficients_summary = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")
# get the fit details
all.fits <- read.csv(file = "data/results/all_fit_detailsINPAT.csv")

mylogit <- function(p) log(p/(1-p))
mylogit10 <- function(p) log10(p/(1-p))

#rename categories in the trend category
coefficients_summary[coefficients_summary$category == 'logistic ns decreasing',]$trend = "decreasing (ns)"
coefficients_summary[coefficients_summary$category == 'logistic s decreasing',]$trend = 'decreasing (s)'
coefficients_summary[coefficients_summary$category == 'logistic ns increasing',]$trend = 'increasing (ns)' 
coefficients_summary[coefficients_summary$category == 'logistic s increasing',]$trend = 'increasing (s)'
coefficients_summary[coefficients_summary$category == 'flat ',]$trend = 'stable'
coefficients_summary[coefficients_summary$category == 'plateauing',]$trend = 'stabilising'
coefficients_summary[coefficients_summary$category=='poor fit',]$trend = 'poor fit'

coefficients_summary$trend = factor(coefficients_summary$trend, levels = c('decreasing (s)',"decreasing (ns)",'stable','stabilising','increasing (ns)','increasing (s)', 'poor fit'))


# get time series
amr_summary <- read.csv(file = "data/summary_AMR_filtered.csv")
amr_summary <- amr_summary %>% arrange(Year)
amr_summary$combR <- paste(amr_summary$Pathogen, amr_summary$Country, amr_summary$Antibiotic, sep = "|")

# get fit parameters

### write function for plotting

plot_time_courses = function(cat){
  all_combR = unique(coefficients_summary[coefficients_summary$category == cat,]$combR)
  file.name = paste("output/supporting/",gsub(" ","",cat),'.pdf',sep = '')
  pdf(file.name, paper = "a4", width = 0, height = 0)
  par(mfrow = c(4, 3))
  par(las = 1)
  for(i in 1:length(all_combR)){
      #get data
      mycomb <- all_combR[i]
      tmp <- amr_summary[which(amr_summary$combR == mycomb),]
      tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
      tmp1.fit = coefficients_summary[coefficients_summary$combR==mycomb,]
      
      fit.data = all.fits[all.fits$combR==mycomb,]
      
      # do plots
      mytitle = mycomb
      plot(tmp1$Year, tmp1$p, col = "red", ylim = c(0,1), xlim = c(1998, 2019), pch = 20, type = "o", ylab = "frequency", xlab = "year", main = mytitle)
      segments(x0 = tmp1$Year, y0 = tmp1$p_min, x1 = tmp1$Year, y1 = tmp1$p_max, col = 'red')
      if(tmp1.fit$best.model == 'flat'){abline(h = fit.data$m0.plateau,col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4)}
      if(tmp1.fit$best.model == 'plateau'){
        curve(fit.data$m2.plateau/(1+exp(fit.data$m2.intercept-fit.data$m2.slope*x)),add=TRUE, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4)}
      if(tmp1.fit$best.model == 'logistic'){
        curve(1/(1+exp(fit.data$m1.intercept-fit.data$m1.slope*x)),from = 1998, to = 2019,add=TRUE, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4)}
      text(2000,0.9,paste("r2 =", round(tmp1.fit$r2.best.fit,2)),pos=4)
      text(2000,0.8,paste("error =", round(tmp1.fit$error.best.fit,2)),pos=4)
      print(mycomb)
    }
  dev.off()
}

plot_time_courses_logit = function(cat){ # FB added plotting time course on logit scale
  
  all_combR = unique(coefficients_summary[coefficients_summary$category == cat,]$combR)
  file.name = paste("output/supporting/", gsub(" ","",cat),'_logit10.pdf',sep = '')
  pdf(file.name, paper = "a4", width = 0, height = 0)
  par(mfrow = c(4, 3))
  par(las = 1)
  for(i in 1:length(all_combR)){
    #get data
    mycomb <- all_combR[i]
    tmp <- amr_summary[which(amr_summary$combR == mycomb),]
    tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
    tmp1.fit = coefficients_summary[coefficients_summary$combR==mycomb,]
    
    fit.data = all.fits[all.fits$combR==mycomb,]
    
    # do plots
    mytitle = mycomb
    plot(tmp1$Year, mylogit10(tmp1$p), col = "red", ylim = c(-3,3), xlim = c(1998, 2019), pch = 20, type = "o", ylab = "logit10(frequency)", xlab = "year", main = mytitle)
    segments(x0 = tmp1$Year, y0 = mylogit10(tmp1$p_min), x1 = tmp1$Year, y1 = mylogit10(tmp1$p_max), col = 'red')
    if(tmp1.fit$best.model == 'flat'){abline(h = mylogit10(fit.data$m0.plateau),col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4)}
    if(tmp1.fit$best.model == 'plateau'){
      curve(mylogit10(fit.data$m2.plateau/(1+exp(fit.data$m2.intercept-fit.data$m2.slope*x))),add=TRUE, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4)}
    if(tmp1.fit$best.model == 'logistic'){
      curve(mylogit10(1/(1+exp(fit.data$m1.intercept-fit.data$m1.slope*x))),from = 1998, to = 2019,add=TRUE, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4)}
    text(2000, mylogit10(0.9),paste("r2 =", round(tmp1.fit$r2.best.fit,2)),pos=4)
    text(2000, mylogit10(0.8),paste("error =", round(tmp1.fit$error.best.fit,2)),pos=4)
    print(mycomb)}
  dev.off()
}

trends = unique(coefficients_summary$category)
for(i in 1:length(trends)){
  cat = trends[i]
  plot_time_courses(cat)
}
# plot on logit scale
for(i in 1:length(trends)){
  cat = trends[i]
  plot_time_courses_logit(cat)
}


######## Goodness of fit FIGURE 3 ########

colset <- RColorBrewer::brewer.pal(n = 7, name = "Set2"); names(colset) <- unique(coefficients_summary$trend)[c(5, 3, 4, 2, 7, 6, 1)]
colset <- c("slategray", "slategray1", "wheat", "wheat4", "pink", "pink4", 'gray'); names(colset) <- unique(coefficients_summary$trend)[c(5, 3, 4, 2, 7, 6, 1)]

# exemple categories
coefficients_summary[which(coefficients_summary$error.best.fit < 0.01 & coefficients_summary$trend=='increasing (s)'), "combR"]

select_cat <- function(x){
  if(x!='poor fit'){
    possible_comb <- coefficients_summary[which(coefficients_summary$error.best.fit < 0.02 & coefficients_summary$trend==x), "combR"]
  } else {
    possible_comb <- coefficients_summary[which(coefficients_summary$trend==x & coefficients_summary$error.best.fit>0.1), "combR"]
  }
  for(mycomb in possible_comb){
    sub <- which(amr_summary$combR == mycomb & amr_summary$patientType=="INPAT")
    if(length(sub) > 10 & all(amr_summary[sub, "p"] > 0.05)) print(mycomb) # 10 time points and large enough freq
  }
}
select_cat('increasing (s)')
select_cat('increasing (ns)')
select_cat('stabilising')
select_cat('stable')
select_cat('decreasing (ns)')
select_cat('decreasing (s)')
select_cat('poor fit')

ylim <- c(0, 1)

pdf("output/SIfigure1.pdf", width = 3, height = 8)

par(mar = c(4,4,1,1), xpd = F)
layout <- layout(matrix(c(1,2,3), 3, 1, byrow = T))

curve(1/(1+exp(1847.528-0.9226887*x)),ylab = "Frequency", xlab = "Year", las = 1, bty = "n", col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4, ylim = ylim,xlim = c(1998, 2019))
text(x = 2008.5, y = ylim[2]-0.06 , 'Standard logistic model', adj = 0.5, cex = 1.5, col = 'black')
text(x = 2012,0.5, expression(f(t) == frac(1,1+e^(-k(t-k[0])))))


plot(NA,ylim = ylim, xlim = c(1998, 2019), las = 1, bty = "n", ylab = "Frequency", xlab = "Year",)
abline(h = 0.65, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4, ylim = ylim, xlim = c(1998, 2019))
text(x = 2008.5, y = ylim[2]-0.06, 'Flat model', adj = 0.5, cex = 1.5, col = 'black')
text(x = 2012,0.5, expression(f(t) == k[2]))

curve(0.7/(1+exp(1847.528-0.9226887*x)),ylab = "Frequency", xlab = "Year", las = 1, bty = "n", col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4, ylim = ylim, xlim = c(1998, 2019))
text(x = 2008.5, y = ylim[2]-0.06, 'Plateauing logistic model', adj = 0.5, cex = 1.5, col = 'black')
text(x = 2012,0.5, expression(f(t) == frac(k[2],1+e^(-k[1](t-k[0])))))

dev.off()

pdf("output/goodnessOfFitFigure.pdf", width = 11, height = 11)

par(mar = c(4,4,1,1), xpd = F)
layout <- layout(matrix(c(1,2,3,4,5,6,7,8,8), 3, 3, byrow = T))

combs <- c("KLEPNE|Netherlands|CIP",  "PSEAER|Hungary|AMK", "ESCCOL|Austria|CIP", "KLEPNE|Germany|CIP", "STAAUR|Spain|OXA", "STAAUR|Italy|RIF", "ESCCOL|Germany|AMX")
names(combs) <- c("increasing (s)", "increasing (ns)", "stabilising", "stable", "decreasing (ns)", "decreasing (s)", "poor fit")      

combs <- combs[c(6, 5, 4, 3, 2, 1, 7)]

species_acro <- c("ACISPP", "ENCFAE", "ENCFAI", "ESCCOL", "KLEPNE", "PSEAER", "STAAUR", "STRPNE")
species_names <- c("Acinetobacter spp.", "E. faecalis", "E. faecium", "E. coli", "K. pneumoniae", "P. aeruginosa", "S. aureus", "S. pneumoniae")

#get data
for(mycomb in combs){
  
  tmp <- amr_summary[which(amr_summary$combR == mycomb),]
  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  tmp1.fit = coefficients_summary[coefficients_summary$combR==mycomb,]
  
  fit.data = all.fits[all.fits$combR==mycomb,]
  
  # do plots
  mytitle <- mycomb
  for(i in 1:8) mytitle <- gsub(pattern = species_acro[i], replacement = species_names[i], mytitle)
  mytitle <- gsub(pattern = "\\|", replacement = ", ", x = mytitle)
  mytitle2 <- names(combs)[combs==mycomb] # name of category
  
  mycol <- colset[mytitle2]
  if(mytitle2=='poor fit') ylim <- c(0, 1) else ylim <- c(0, 0.5)
  plot(tmp1$Year, tmp1$p, col = mycol, ylim = ylim, xlim = c(1998, 2019), pch = 20, type = "o", lwd = 2, ylab = "Frequency", xlab = "Year", las = 1, bty = "n")
  segments(x0 = tmp1$Year, y0 = tmp1$p_min, x1 = tmp1$Year, y1 = tmp1$p_max, col = mycol, lwd = 2)
  if(tmp1.fit$best.model == 'flat'){abline(h = fit.data$m0.plateau,col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4)}
  if(tmp1.fit$best.model == 'plateau'){
    curve(fit.data$m2.plateau/(1+exp(fit.data$m2.intercept-fit.data$m2.slope*x)),add=TRUE, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4)}
  if(tmp1.fit$best.model == 'logistic'){
    curve(1/(1+exp(fit.data$m1.intercept-fit.data$m1.slope*x)),from = 1998, to = 2019,add=TRUE, col=rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 0.5), lwd = 4)}
  
  text(x = 2008.5, y = ylim[2], mytitle2, adj = 0.5, cex = 1.5, col = mycol)
  text(x = 2008.5, y = ylim[2] - 0.1 * (ylim[2]-ylim[1]), mytitle, adj = 0.5, cex = 1.)
}

# error distribution version 2
plot(NULL, xlim = c(0, 0.2), ylim = c(0, 250), xlab = "Error", ylab = "Counts", axes = F)
axis(side = 1, at = seq(0, 0.2, 0.05))
axis(side = 2, at = seq(0, 250, 50), las= 1)
#this is the part that fails between intel and ARM -> made to 0.3 instead
mb <- seq(0, 0.3, 0.01); n <- length(mb)
cumul <- rep(0, n-1)
for(mycat in unique(coefficients_summary$trend)[c(5, 3, 4, 2, 7, 6, 1)]){
  sub <- which(coefficients_summary$trend==mycat)
  errors <- coefficients_summary$error.best.fit[sub]
  htmp <- hist(errors, breaks = mb, plot = F)
  for(i in 1:(n-1)){
    rect(xleft = mb[i], xright = mb[i+1], ybottom = cumul[i], ytop = cumul[i]+htmp$counts[i], col = colset[mycat])
  }
  cumul <- cumul + htmp$counts
}
legend("topright", legend = names(colset)[c(6, 5, 4, 3, 2, 1)], fill = colset[c(6, 5, 4, 3, 2, 1)])

dev.off()

