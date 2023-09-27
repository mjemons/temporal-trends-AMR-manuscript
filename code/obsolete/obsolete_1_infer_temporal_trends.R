rm(list = ls())

##########                                                                INVESTIGATE TEMPORAL TRENDS IN AMR                                                                ##########
library(dotenv)

load_dot_env(file = ".env")

setwd(Sys.getenv(c("HOME_DIR")))

source("scripts/source_epidata.R")

#amr_summary <- amr_summary[which(amr_summary$Year < 2020), ]

# create data-frame to store coefficients:
range_year <- sort(unique(amr_summary$Year))
df_coefficients <- as.data.frame(matrix(NA, nrow = length(all_combR), ncol = 6*length(range_year)))
colnames(df_coefficients) <- apply(expand.grid(c("inpatient_0d_", "inpatient_1d_", "inpatient_2d_", "outpatient_0d_", "outpatient_1d_", "outpatient_2d_"), range_year), 1, paste, collapse="")
n_years <- length(range_year)
n_coeffs <- 6*length(range_year)
df_coefficients$combR <- all_combR

# and a df to store CIs:
df_ci <- as.data.frame(matrix(NA, nrow = length(all_combR), ncol = 12*length(range_year)))
colnames(df_ci) <- apply(expand.grid(c("inpatient_0d_", "inpatient_1d_", "inpatient_2d_",
                                       "outpatient_0d_", "outpatient_1d_", "outpatient_2d_"), range_year, c("_lowerCI", "_upperCI")), 1, paste, collapse="")

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
  
  # INPATIENT
  if((nrow(tmp1) > 5) & (all(tmp1$N > 30))){ # at least 5 years & all points > 100
    print(mycomb)
    plot(tmp1$Year, tmp1$p, col = "red", ylim = c(0,1), xlim = c(1998, 2018), pch = 20, type = "o", ylab = "frequency", xlab = "year", main = mycomb)
    segments(x0 = tmp1$Year, y0 = tmp1$p_min, x1 = tmp1$Year, y1 = tmp1$p_max, col = 'red')
    spl <- smooth.spline(x = tmp1$Year, y = tmp1$p, w = 1 / tmp1$N, df = 5)
    lines(spl$x, spl$y, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.5), lwd = 4)
    is_plot <- T
    
    if(compute_ci <- TRUE){
      all_spl <- lapply(1:100, function(i) bs_ss(tmp1)) # splines on random data tmp1
      all_spl <- do.call(rbind, all_spl)
      # confidence intervals on values, first derivatives, second derivatives
      ci_d0 <- apply(all_spl[seq(2, 400, 4),], 2, function(x) quantile(x, c(0.025, 0.975)))
      ci_d1 <- apply(all_spl[seq(3, 400, 4),], 2, function(x) quantile(x, c(0.025, 0.975)))
      ci_d2 <- apply(all_spl[seq(4, 400, 4),], 2, function(x) quantile(x, c(0.025, 0.975)))
    }
    # save values, derivatives in df table
    pred <- predict(spl)
    pred1 <- predict(spl, deriv = 1)
    pred2 <- predict(spl, deriv = 2)
    df_coefficients[i, paste0("inpatient_0d_", pred$x)] <- pred$y
    df_coefficients[i, paste0("inpatient_1d_", pred$x)] <- pred1$y
    df_coefficients[i, paste0("inpatient_2d_", pred$x)] <- pred2$y
    
    # record CI
    df_ci[i, paste0("inpatient_0d_", pred$x, "_lowerCI")] <- ci_d0["2.5%",]
    df_ci[i, paste0("inpatient_0d_", pred$x, "_upperCI")] <- ci_d0["97.5%",]
    df_ci[i, paste0("inpatient_1d_", pred$x, "_lowerCI")] <- ci_d1["2.5%",]
    df_ci[i, paste0("inpatient_1d_", pred$x, "_upperCI")] <- ci_d1["97.5%",]
    df_ci[i, paste0("inpatient_2d_", pred$x, "_lowerCI")] <- ci_d2["2.5%",]
    df_ci[i, paste0("inpatient_2d_", pred$x, "_upperCI")] <- ci_d2["97.5%",]
    
  }
  # OUTPATIENT
  if((nrow(tmp2) > 5) & (all(tmp2$N > 30))){ # at least 5 years for one of them
    if(is_plot) plot_fun <- points else plot_fun <- plot
    plot_fun(tmp2$Year+0.1, tmp2$p, col = "blue", pch = 20, type = "o", ylim = c(0,1), xlim = c(1998, 2018), ylab = "frequency", xlab = "year", main = mycomb)
    segments(x0 = tmp2$Year+0.1, y0 = tmp2$p_min, x1 = tmp2$Year+0.1, y1 = tmp2$p_max, col = 'blue')
    spl <- smooth.spline(x = tmp2$Year, y = tmp2$p, w = 1 / tmp2$N, df = 5)
    lines(spl$x, spl$y, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5), lwd = 4)
    
    if(compute_ci <- TRUE){ # compute confidence intervals
      all_spl <- lapply(1:100, function(i) bs_ss(tmp2)) # splines on random data generate_randomised(tmp2)
      all_spl <- do.call(rbind, all_spl)
      # confidence intervals on values, first derivatives, second derivatives
      ci_d0 <- apply(all_spl[seq(2, 400, 4),], 2, function(x) quantile(x, c(0.025, 0.975)))
      ci_d1 <- apply(all_spl[seq(3, 400, 4),], 2, function(x) quantile(x, c(0.025, 0.975)))
      ci_d2 <- apply(all_spl[seq(4, 400, 4),], 2, function(x) quantile(x, c(0.025, 0.975)))
    }
    # save values, derivatives in df table
    pred <- predict(spl)
    pred1 <- predict(spl, deriv = 1)
    pred2 <- predict(spl, deriv = 2)
    df_coefficients[i, paste0("outpatient_0d_", pred$x)] <- pred$y
    df_coefficients[i, paste0("outpatient_1d_", pred$x)] <- pred1$y
    df_coefficients[i, paste0("outpatient_2d_", pred$x)] <- pred2$y
    
    # record CI
    df_ci[i, paste0("outpatient_0d_", pred$x, "_lowerCI")] <- ci_d0["2.5%",]
    df_ci[i, paste0("outpatient_0d_", pred$x, "_upperCI")] <- ci_d0["97.5%",]
    df_ci[i, paste0("outpatient_1d_", pred$x, "_lowerCI")] <- ci_d1["2.5%",]
    df_ci[i, paste0("outpatient_1d_", pred$x, "_upperCI")] <- ci_d1["97.5%",]
    df_ci[i, paste0("outpatient_2d_", pred$x, "_lowerCI")] <- ci_d2["2.5%",]
    df_ci[i, paste0("outpatient_2d_", pred$x, "_upperCI")] <- ci_d2["97.5%",]
  }
}
dev.off()

# remove all-NA rows
df_coefficients <- df_coefficients[!(apply(df_coefficients, 1, function(ll) all(is.na(ll[1:n_coeffs])))), ]
df_ci <- df_ci[!(apply(df_ci, 1, function(ll) all(is.na(ll)))), ]

#################################### 1. THE SHAPE OF THE INITIAL INCREASE IN FREQUENCY ####################################
# is the initial increase accelerating or not?

# value_idx <- grepl(pattern = "value", colnames(df_coefficients))
# hist(
#   unlist(df_coefficients[value_idx])
# )

# here we select the first year when the frequency of resistance was low (< 0.05) and add a column of first and second derivatives at this year
df_coefficients$low_frequency_0d <- NA
df_coefficients$low_frequency_1d <- NA
df_coefficients$low_frequency_2d <- NA
df_coefficients$low_frequency_0d_lowerCI <- NA
df_coefficients$low_frequency_1d_lowerCI <- NA
df_coefficients$low_frequency_2d_lowerCI <- NA
df_coefficients$low_frequency_0d_upperCI <- NA
df_coefficients$low_frequency_1d_upperCI <- NA
df_coefficients$low_frequency_2d_upperCI <- NA
df_coefficients$last_1d <- NA
inpatient_0d_cols <- grepl(pattern = "0d", x = colnames(df_coefficients)) & grepl(pattern = "inpatient", x = colnames(df_coefficients)) 
for(i in 1:nrow(df_coefficients)){
  # add column for frequency, first derivative and second derivative at small frequencies
  low_frequency <- df_coefficients[i, inpatient_0d_cols] < 0.05
  if(any(low_frequency, na.rm = T)){
    idx_low_frequency <- which(inpatient_0d_cols)[which(low_frequency)][1] # first year when frequency of resistance was low
    colname <- names(df_coefficients)[idx_low_frequency]
    df_coefficients[i, "low_frequency_0d"] <- df_coefficients[i, idx_low_frequency]
    df_coefficients[i, "low_frequency_1d"] <- df_coefficients[i, idx_low_frequency + 1]
    df_coefficients[i, "low_frequency_2d"] <- df_coefficients[i, idx_low_frequency + 2]
    
    colname_0d_ci <- paste0(colname, c("_lowerCI", "_upperCI"))
    colname_1d_ci <- gsub(pattern = "0d", replacement = "1d", x = colname_0d_ci)
    colname_2d_ci <- gsub(pattern = "0d", replacement = "2d", x = colname_0d_ci)
    df_coefficients[i, c("low_frequency_0d_lowerCI", "low_frequency_0d_upperCI") ] <- df_ci[i, colname_0d_ci]
    df_coefficients[i, c("low_frequency_1d_lowerCI", "low_frequency_1d_upperCI") ] <- df_ci[i, colname_1d_ci]
    df_coefficients[i, c("low_frequency_2d_lowerCI", "low_frequency_2d_upperCI") ] <- df_ci[i, colname_2d_ci]

  }
  # add column for the last first derivative
  non_na <- which(!is.na(df_coefficients[i, inpatient_0d_cols]))
  if(any(non_na)){
    last_value_idx <- which(inpatient_0d_cols)[non_na[length(non_na)]]
    df_coefficients$last_1d[i] <- df_coefficients[i, last_value_idx+1] # first derivative (is at column last_value_idx+1)
  }
}

#################################### 2. DOES RESISTANCE FREQUENCY KEEP INCREASING OR PLATEAUS? ####################################

# test for a plateau: the last 3 points all have non-significant 1st order derivative:
inpatient_1d_cols <- grepl(pattern = "1d", x = colnames(df_coefficients)) & grepl(pattern = "inpatient", x = colnames(df_coefficients)) 
df_coefficients$is_plateau <- NA
for(i in 1:nrow(df_coefficients)){
  non_na <- which(!is.na(df_coefficients[i, inpatient_1d_cols]))
  if(length(non_na) > 3){
    idx_to_consider <- which(inpatient_1d_cols)[non_na[(length(non_na)-2):length(non_na)]] # 3 last non-NA values
    years <- get_year(names(df_coefficients[i, idx_to_consider]))
    if(
      all(df_ci[i, paste0("inpatient_1d_", years, "_lowerCI")] < 0) &
      all(df_ci[i, paste0("inpatient_1d_", years, "_upperCI")] > 0) # if for these three points, the confidence interval includes 0
    ){
      df_coefficients$is_plateau[i] <- "PLATEAU" # then say this is a plateau  
    }
    if(
      all(df_ci[i, paste0("inpatient_1d_", years, "_lowerCI")] > 0) # if for these three points, trend is significantly positive
    ){
      df_coefficients$is_plateau[i] <- "RISE" # then say this is a plateau  
    }
    if(
      all(df_ci[i, paste0("inpatient_1d_", years, "_upperCI")] < 0) # if for these three points, trend is significantly negative
    ){
      df_coefficients$is_plateau[i] <- "DECLINE" # then say this is a plateau  
    }
    if(is.na(df_coefficients$is_plateau[i])) df_coefficients$is_plateau[i] <- "OTHER"
  }
}

df_coefficients$species <- unlist(lapply(strsplit(df_coefficients$combR, "\\|"), function(ll) ll[[1]]))
df_coefficients$country <- unlist(lapply(strsplit(df_coefficients$combR, "\\|"), function(ll) ll[[2]]))
df_coefficients$ab <- unlist(lapply(strsplit(df_coefficients$combR, "\\|"), function(ll) ll[[3]]))

# WRITE UP EVERYTHING:
write.csv(x = df_coefficients, file = "data/smooth_splines_coefficients.csv", row.names = F, col.names = T)
write.csv(x = df_ci, file = "data/smooth_splines_coefficients_ci.csv", row.names = F, col.names = T)



