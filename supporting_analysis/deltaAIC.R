library(dplyr)

# load the fit details
coefficients_summary = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")

coefficients_summary[coefficients_summary$category == 'logistic ns decreasing',]$category = "decreasing (ns)"
coefficients_summary[coefficients_summary$category == 'logistic s decreasing',]$category = 'decreasing (s)'
coefficients_summary[coefficients_summary$category == 'logistic ns increasing',]$category = 'increasing (ns)' 
coefficients_summary[coefficients_summary$category == 'logistic s increasing',]$category = 'increasing (s)'
coefficients_summary[coefficients_summary$category == 'flat ',]$category = 'stable'
coefficients_summary[coefficients_summary$category == 'plateauing',]$category = 'stabilising'

# filter out poor fits as we don't consider them anymore in the analysis after this point
coefficients_summary <- coefficients_summary %>% filter(!(category %in% c('poor fit')))

delta_aicc_df <- data.frame(delta_aicc = coefficients_summary$deltaAICc)

p <- ggplot(delta_aicc_df, aes(x=delta_aicc)) + 
  geom_histogram(binwidth = 1, boundary = 1) + 
  theme_light() + 
  xlab('delta AICc') +
  geom_vline(xintercept = 4, col = 'red') +
  geom_vline(xintercept = 2, col = 'blue')
ggsave('output/histogram_delta_aicc.pdf')

# filter out all fits that have a delta AICc < 2
coefficients_summary_delta_aicc2 <- coefficients_summary %>% filter(deltaAICc < 2) %>%
  select(combR, best.model, second.best.model, category, deltaAICc)

# round deltaAICc values to 3 digits
coefficients_summary_delta_aicc2$deltaAICc <- round(coefficients_summary_delta_aicc2$deltaAICc, 3)
 
write.csv(coefficients_summary_delta_aicc2,paste("data/results/deltaAiccListINPAT",".csv",sep=''), row.names = FALSE)
