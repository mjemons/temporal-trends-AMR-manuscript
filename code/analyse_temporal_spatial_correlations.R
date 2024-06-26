rm(list = ls());
#library(plyr); library(broom); library(MetBrewer)
load(file = "data/AMR_AMC_correlation.RData")
amr_summary <- read.csv(file = "data/summary_AMR_filtered.csv");  corr_AB_class <- table(amr_summary$Antibiotic, amr_summary$Antibiotic_class); corr_AB_class <- apply(corr_AB_class, 1, function(x)names(x)[which(x!=0)]); rm(amr_summary)
source('code/Utils.R')
#ls()

# AllNew = cleaning of antibiotic type (All, Oral, etc.)
# four slopes defined, v2-v3-v4 are better
# application of correlation slope / plateau / etc. only when it makes sense

get_overall_ci <- function(mymean, mylower, myupper){
  # a function to get overall confidence intervals of mean of several mean quantities with each their own CI
  stopifnot(length(mymean)==length(mylower))
  stopifnot(length(mymean)==length(myupper))
  mysd <- (myupper - mylower)/(2 *1.96) # get sd if normality of distribution
  stopifnot(length(mymean)==length(mysd))
  sims <- sapply(1:length(mymean), function(i) rnorm(n = 1000, mean = mymean[i], sd = mysd[i]))
  sims <- rowMeans(sims)
  return(
    quantile(sims,  c(0.025, 0.975))
  )
}

species_acro <- c("ACISPP", "ENCFAE", "ENCFAI", "ESCCOL", "KLEPNE", "PSEAER", "STAAUR", "STRPNE")
species_names <- c("Acinetobacter spp.", "E. faecalis", "E. faecium", "E. coli", "K. pneumoniae", "P. aeruginosa", "S. aureus", "S. pneumoniae")
species_cols <- MetBrewer::met.brewer(name = "VanGogh2", n = 8)
names(species_names) <- names(species_cols) <- species_acro

#################### FOREST PLOTS OF AMR-AMC SPATIAL CORRELATION ####################

forest_plot <- function(rho_object, main = "", add_legend = F) {
  
  stopifnot("Pathogen" %in% names(rho_object))
  stopifnot("Antibiotic" %in% names(rho_object))
  stopifnot("rho" %in% names(rho_object))
  stopifnot("lwr.ci" %in% names(rho_object))
  stopifnot("upr.ci" %in% names(rho_object))
  stopifnot("Ncountries" %in% names(rho_object))
  
  n <- length(rho_object$Pathogen)
  oo <- order(rho_object$rho, decreasing = F)
  par(mar = c(4,4,2,1), xpd = TRUE)
  
  plot(NULL, ylim = c(0, n), xlim = c(-1, 1), bty = "n", axes = F, xlab = "Correlation coefficient", ylab = "", main = main)
  
  # vertical lines:
  segments(x0 = 0, x1 = 0, y0 = 0, y1  = n+1, lty = 1, col = "black")
  #segments(x0 = mean(rho_object$rho), x1 = mean(rho_object$rho), y0 = 0, y1  = n+1, lwd = 2, lty = 2)
  axis(side = 1, at = seq(-1, 1, 0.5), cex.axis = 1)
  
  
  
  mycols <- sapply(rho_object$lwr.ci, function(x) ifelse(x>0,"black", "gray")) #species_cols[rho_object$Pathogen[oo]]
  segments(x0 = rho_object$lwr.ci[oo], x1 = rho_object$upr.ci[oo], y0 = 1:n, y1 = 1:n, lwd = 3, col = mycols[oo])
  segments(x0 = rho_object$lwr.ci[oo], x1 = rho_object$lwr.ci[oo], y0 = (1:n)-0.2, y1 = (1:n)+0.2, lwd = 3, col = mycols[oo])
  segments(x0 = rho_object$upr.ci[oo], x1 = rho_object$upr.ci[oo], y0 = (1:n)-0.2, y1 = (1:n)+0.2, lwd = 3, col = mycols[oo])
  
  #points(x = rho_object$rho[oo], y = 1:n, pch = 1, col = "black")
  points(x = rho_object$rho[oo], y = 1:n, pch = 20, col = mycols[oo], cex = 2)
  text(x = -1.01, y = 1:n,
       labels = paste(rho_object$Pathogen[oo], ", ",
                      rho_object$Antibiotic[oo],
                      paste(", N =", rho_object$Ncountries[oo]), sep = ""),
       cex = 0.5, adj = 1)
  col_classes <- RColorBrewer::brewer.pal(n = 8, name = "Set1"); names(col_classes) <- sort(unique(corr_AB_class)) # colors for classes
  rect(xleft = -0.95, xright = -1, ybottom = (1:n)-0.5, ytop = (1:n)+0.5, col = col_classes[corr_AB_class[rho_object$Antibiotic[oo]]], border = NA) # rectangles
  if(add_legend){
    subset <- which(names(col_classes) %in% corr_AB_class[rho_object$Antibiotic[oo]]) #c(2, 3, 5, 6, 7) # subset of colors presents
    nclasses <- length(subset)
    rect(xleft = -0.8, xright = -0.9, ybottom = seq(n, n-(nclasses-1)*2, - 2)-0.9, ytop = seq(n, n-(nclasses-1)*2, -2)+0.9, col = col_classes[subset], border = NA)
    text(x = -0.78, y = seq(n, n-2*(nclasses-1), - 2),
         labels = c("Tetracyclines", "Penicillins", "Other beta-lactams", "Macrolides", "Aminoglycosides", "Quinolones", "Others", "Antimycobacterials")[subset],
         #labels = names(col_classes),
         adj = 0)
  }
  
  # mean and CI at the bottom
  mymean <- mean(rho_object$rho); myci <- t.test(rho_object$rho)$conf.int
  width <- n/80
  polygon(c(mymean, mymean, myci[2]), c(0.5-1/n-width,0.5-1/n+width,0.5-1/n), col = "black")
  polygon(c(mymean, mymean, myci[1]), c(0.5-1/n-width,0.5-1/n+width,0.5-1/n), col = "black")
  
  #legend(-0.9, y = n, legend = species_names, pch = 20, col = species_cols, bty = "n", cex = 0.9)
  #legend(-0.9, y = n, legend = rep("", 8), pch = 1, col = "black", bty = "n", cex = 0.9)
 # text(x = -1.2, y = 1:n, labels = species_names[rho_object$Pathogen[oo]], adj = 0, cex = 0.5, col = mycols)
}

# all figures
inoutpat <- "INPAT"
for(variable in c(c("plateau", "slope_v2", "slope_v3", "slope_v4", "median"))){
  for(mysector in c("Community", "Hospital Sector")){
    if(mysector=="Hospital Sector") mysector2 <- "HS" else mysector2 <- mysector
    for(mytype in c("Oral", "AllNew")){
      # get rho object:
      rho_object <- get(paste0("rho_", variable, "_", mysector2, "_", mytype, "_", inoutpat))
      truc <- paste0("output/forestplot", gsub(pattern = "_", replacement  = "", variable), "", mysector2, "", mytype, "", inoutpat, ".pdf")
      pdf(truc, width = 6, height = max(4, 6/40 * length(rho_object$rho)))
      forest_plot(rho_object)
      dev.off()
    }
  }
}

# added 16/10/2023: figure 5 directly
pdf("output/figure5.pdf", width = 6*2, height = 6)
par(mfrow = c(1,2))
forest_plot(rho_plateau_Community_AllNew_INPAT, add_legend = T)
text(-1.2, 44, "A", cex = 2) 
forest_plot(rho_slope_v4_Community_AllNew_INPAT)
text(-1.2, 16.1, "B", cex = 2) 
dev.off()

# with this measure of slope, actually weak signal of positive correlation
# plateau:
write.table("Plateau_Community", file = "report.txt", append = F, row.names = F)
write.table(x = broom::tidy(t.test(rho_plateau_Community_AllNew_INPAT$rho)), file = "report.txt", sep = ",", append = T, row.names = F)
write.table("Plateau_Hospital", file = "report.txt", append = T, row.names = F)
write.table(x = broom::tidy(t.test(rho_plateau_HS_AllNew_INPAT$rho)), file = "report.txt", sep = ",", append = T, row.names = F)

# slope:
write.table("Slope_Community", file = "report.txt", append = T, row.names = F)
write.table(x = broom::tidy(t.test(rho_slope_v4_Community_AllNew_INPAT$rho)), file = "report.txt", sep = ",", append = T, row.names = F)
write.table("Slope_Hospital", file = "report.txt", append = T, row.names = F)
write.table(x = broom::tidy(t.test(rho_slope_v4_HS_AllNew_INPAT$rho)), file = "report.txt",  sep = ",", append = T, row.names = F)

# stabilising only:
t.test(rho_slope_v4_Community_AllNew_INPAT_stabilisingOnly$rho)
t.test(rho_slope_v4_HS_AllNew_INPAT_stabilisingOnly$rho)


#################### FIGURES RESTRICTED TO STABILISING ONLY ####################

# stabilising only:
rho_object <- rho_slope_v4_Community_AllNew_INPAT_stabilisingOnly
truc <- paste0("output/forestplotslopev4CommunityAllNewINPAT_stabilisingOnly.pdf")
pdf(truc, width = 6, height = max(4, 6/40 * length(rho_object$rho)))
forest_plot(rho_object)
dev.off()

rho_object <- rho_slope_v4_HS_AllNew_INPAT_stabilisingOnly
truc <- paste0("output/forestplotslopev4HSAllNewINPAT_stabilisingOnly.pdf")
pdf(truc, width = 6, height = max(4, 6/40 * length(rho_object$rho)))
forest_plot(rho_object)
dev.off()

#################### TIME-SHIFT PLOTS OF AMR-AMC TEMPORAL CORRELATION ####################


inoutpat<-"INPAT"

for(mysector in c("Community", "Hospital Sector")){
  if(mysector=="Hospital Sector") mysector2 <- "HS" else mysector2 <- mysector
  
  for(mytype  in c("Oral", "AllNew")){
    
    # remove NA correlations and create three 'rho_object' for past, present, future:
    for(ts in c("_present", "_past", "_future")){
      assign(x = "truc", get(paste0("rho_temporal_", mysector2, "_", mytype, "_", inoutpat, ts)))
      truc <- truc[!is.na(truc$rho), ]
      assign(x = paste0("rho_object", ts), truc)
    }
    
    # define drug bug combinations:
    rho_object_present$db <- paste(rho_object_present$Pathogen, rho_object_present$Antibiotic, sep = "|")
    rho_object_past$db <- paste(rho_object_past$Pathogen, rho_object_past$Antibiotic, sep = "|")
    rho_object_future$db <- paste(rho_object_future$Pathogen, rho_object_future$Antibiotic, sep = "|")
    
    # drug-bug combinations present for past, present, future in at least 5 countries:
    tab_present <- table(rho_object_present$db)
    tab_past <- table(rho_object_past$db)
    tab_future <- table(rho_object_future$db)
    
    db_present <- names(tab_present)[tab_present >= 5]
    db_past <- names(tab_past)[tab_past >= 5]
    db_future <- names(tab_future)[tab_future >= 5]
    all_db <- intersect(intersect(db_present, db_past), db_future)
    
    # mean across countries, for each drug-bug combination
    get_mean_ci <- function(tab){
      
      NAvec <- rep(NA, length(all_db))
      # final data-frame
      out <- list(Pathogen = NAvec, Antibiotic = NAvec, db = all_db, rho = NAvec, lwr.ci = NAvec, upr.ci = NAvec, Ncountries = NAvec)
      for(dd in all_db){
        idx <- which(out$db==dd)
        sub <- which(tab$db==dd)
        ncountries <- length(sub)
        if(ncountries >= 5){
          my_mean_rho <- mean(tab$rho[sub], na.rm = T)
          my_ci_rho <- get_overall_ci(mymean = tab$rho[sub], mylower = tab$lwr.ci[sub], myupper = tab$upr.ci[sub])
          
          out$Pathogen[idx] <- strsplit(dd, "\\|")[[1]][1]
          out$Antibiotic[idx] <- strsplit(dd, "\\|")[[1]][2]
          out$db[idx] <- dd
          out$rho[idx] <- my_mean_rho
          out$lwr.ci[idx] <- my_ci_rho[1]
          out$upr.ci[idx] <- my_ci_rho[2]
          out$Ncountries[idx] <- ncountries
        }
      }
     if(any(out$Ncountries < 5)) stop("all combinations should be represented by at least 5 countries")
     out
    }
    
    truc_present <- get_mean_ci(rho_object_present)
    truc_past <-    get_mean_ci(rho_object_past)
    truc_future <-  get_mean_ci(rho_object_future)

    write.table(paste("Temporal_autocorrelation", mysector2, "|", mytype, "|", inoutpat), file = "report.txt", append = T, row.names = F)
    write.table(x = broom::tidy(t.test(truc_present$rho)), file = "report.txt", sep = ",", append = T, row.names = F)
    
    truc_present$positive_ts <- truc_present$rho > truc_past$rho & truc_present$rho > truc_future$rho
    print(paste(x <- sum(truc_present$positive_ts, na.rm = T), "/", n <- length(truc_present[[1]])))
    print("p-value (binomial test compared to 0.25)")
    print(bt <- binom.test(x = x, n = n, p = 0.25, alternative = "greater"))
    
    # forest plot
    pdf(paste0("output/ForestPlotTimeshift", mysector2, "", mytype, "", inoutpat, ".pdf"), width = 4 * 1.4, height = 1.5 * 3 * 1.4)
    forest_plot(rho_object = truc_present)
    dev.off()
    
    # plot
    pdf(paste0("output/timeshift", mysector2, "", mytype, "", inoutpat, ".pdf"), width = 4 * 1.4, height = 3 * 1.4)
    
    par(mfrow = c(1,1), mar = c(4,4,1,1))
    
    # compare the temporal correlation with the corresponding spatial correlation:
    
    # select only bug|drug combinations that are common to both:
    rho_space <- get(paste0("rho_median_", mysector2, "_", mytype, "_", inoutpat))
    rho_space$db <- paste(rho_space$Pathogen, rho_space$Antibiotic, sep = "|")
    all_db_2 <- intersect(all_db, rho_space$db)
    
    keep_common_db <- function(df, db_sel = all_db) {
      if(is.list(df)) df <- as.data.frame(do.call(cbind, df))
      df <- df[df$db %in% db_sel, ]; df <- df[match(sort(db_sel), df$db), ];
      for(mycol in c("rho", "lwr.ci", "upr.ci", "Ncountries")) df[, mycol] <- as.numeric(df[, mycol])
      return(df)
    }
    
    rho_space <- keep_common_db(df = rho_space, db_sel = all_db_2)
    truc_present_2 <- keep_common_db(truc_present, db_sel = all_db_2)
    truc_past_2 <- keep_common_db(truc_past, db_sel = all_db_2)
    truc_future_2 <- keep_common_db(truc_future, db_sel = all_db_2)
    
    stopifnot(all(rho_space$db == truc_present_2$db))
    stopifnot(all(rho_space$db == truc_past_2$db))
    stopifnot(all(rho_space$db == truc_future_2$db))
    
    plot(NULL, xlim = c(-1, 1), ylim = c(-1, 1), xlab = "Spatial correlation", ylab = "Temporal correlation", pch = 20, las = 1, main = "") #paste(mysector2, "|", mytype, "|", inoutpat))
    segments(x0 = rho_space$rho, x1 = rho_space$rho, y0 = truc_present_2$lwr.ci, y1 = truc_present_2$upr.ci, col = "gray")
    segments(x0 = rho_space$lwr.ci, x1 = rho_space$upr.ci, y0 = truc_present_2$rho, y1 = truc_present_2$rho, col = "gray")
    points(rho_space$rho, truc_present_2$rho, pch = 20, col = species_cols[rho_space$Pathogen])
    abline(h = 0, lty = 2); abline(v = 0, lty = 2)
    summary(lm0 <- lm(truc_present_2$rho ~ rho_space$rho))
    
    legend("topleft", legend = species_names, pch = 20, col = species_cols, bty = "n", cex = 0.8)
    legend("topleft", legend = rep("", 8), pch = 1, col = "black", bty = "n", cex = 0.8)
    
    segments(x0 = -1, x1 = +1, y0 = -1, 1, col = "black", lwd = 1)
    segments(x0 = -1, x1 = +1, y0 = lm0$coefficients[1] - lm0$coefficients[2], y1 =lm0$coefficients[1] + lm0$coefficients[2], col = "red")
    #text(rho_space$rho+0.05, truc_present_2$mean.rho+0.05, labels = truc_present_2$db, cex = 0.4)
    
    
    # p-value, considering that both x and y values are noisy
    # true p-value for relationship (given that both x and y are noisy)
    
    print("Computing p-value...")
    sd_space <- (rho_space$upr.ci - rho_space$lwr.ci)/(2*2*1.96)
    sd_time <- (truc_present_2$upr.ci - truc_present_2$lwr.ci)/(2*2*1.96)
    
    nrep <- 1000
    x <- matrix(NA, nrow = nrep, ncol = length(sd_space))
    y <- matrix(NA, nrow = nrep, ncol = length(sd_time))
    stopifnot(length(sd_space) == length(sd_time)); ndb <- length(sd_time)
    for(i in 1:nrep){
      for(j in 1:ndb){
        samp <- sample(x = 1:ndb, replace = T)
        x[i, j] <- rnorm(n = 1, mean = rho_space$rho[j], sd = sd_space[j])
        y[i, j] <- rnorm(n = 1, mean = truc_present_2$rho[samp[j]], sd = sd_time[samp[j]]) # shuffled value
      }
    }
    all_beta <- c()
    for(i in 1:nrep) all_beta <- c(all_beta, lm(y[i, ] ~ x[i,])$coefficients[2]) # what is coefficient of lm for each shuffled
    pvalue <- 1 - which(sort(all_beta) > lm0$coefficients[2])[1] / nrep
    print("Corrected p-value:"); print(pvalue)
    text(1, -0.8, paste("p =", round(pvalue, 3)), adj = 1)
    
    dev.off()
    
    hist(all_beta, main =  paste(mysector2, "|", mytype, "|", inoutpat))            
    abline(v = lm0$coefficients[2], lwd = 3)
    
  }
}


#################### MULTIPANEL FIGURE 6 ####################


#make multipanel figure 6
#violinplot created in consumption_temporal_trend.R
load(file = "output/violinplot.rdata")

library("vioplot")
colset <- RColorBrewer::brewer.pal(n = 7, name = "Set2"); names(colset) <- unique(h$data$category)[c(1,6,7,2,4,3)]
colset <- c("pink4", "pink", "wheat4", "wheat", "slategray1", 'slategray'); names(colset) <- unique(h$data$category)[c(1,6,7,2,4,3)]

# according to https://stackoverflow.com/questions/14124373/combine-base-and-ggplot-graphics-in-r-figure-window
vp.Right <- grid::viewport(height=unit(1, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"), 
                           y=1, x=0.5)
rho_object <- get_mean_ci(rho_object_present)

pdf(paste0("output/figure_6_v2.pdf"), width =6*2, height = max(4, 6/40 * length(rho_object$rho)))
par(mfrow=c(1,2), xpd = T)
forest_plot(rho_object, add_legend = T)
text(-1.3, 42, "A", cex = 2) 
text(1.1, 42, "B", cex = 2)
plot(NULL, xlim = c(0,6), ylim = c(-1.5,1.5), xlab = "", ylab = "Spearman correlation coefficient use vs. year", axes = F)
axis(side = 2, at = seq(-1,1,0.2), las = 1)
all_cat <- names(table(h$data$category))
for(i in 1:6){
  subset <- which(h$data$category==all_cat[i])
  points(rep(i-0.5, length(subset)), h$data$rho[subset], pch = 20, col = colset[i])
  segments(x0 = i-0.8, x1 = i-0.2, y0 = mean(h$data$rho[subset]), y1 = mean(h$data$rho[subset]), col = colset[i], lwd = 3)
}
text((1:6)-0.5, -1.1, all_cat, srt = 45, col = colset)
#vioplot(h$data$rho ~ h$data$category, col = colset, cex.axis = 0.5, xlab = 'Category', ylab = 'Regression Coefficient [DDD/1000 inhabitants/year/year]', frame.plot = FALSE, vp=vp.Right)
dev.off()

