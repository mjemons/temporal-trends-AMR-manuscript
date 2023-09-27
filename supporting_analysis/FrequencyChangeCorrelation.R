# coefficients summary
coefficients_summary = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")

coefficients_summary[coefficients_summary$category == 'logistic ns decreasing',]$trend = "decreasing (ns)"
coefficients_summary[coefficients_summary$category == 'logistic s decreasing',]$trend = 'decreasing (s)'
coefficients_summary[coefficients_summary$category == 'logistic ns increasing',]$trend = 'increasing (ns)' 
coefficients_summary[coefficients_summary$category == 'logistic s increasing',]$trend = 'increasing (s)'
coefficients_summary[coefficients_summary$category == 'flat ',]$trend = 'stable'
coefficients_summary[coefficients_summary$category == 'plateauing',]$trend = 'stabilising'
coefficients_summary[coefficients_summary$category=='poor fit',]$trend = 'poor fit'

coefficients_summary$trend = factor(coefficients_summary$trend, levels = c('decreasing (s)',"decreasing (ns)",'stable','stabilising','increasing (ns)','increasing (s)', 'poor fit'))

ggplot(coefficients_summary, aes(x=trend, y=median.amr)) + 
  geom_violin() +
  stat_summary(fun.y=median, geom="point", size=2, color="red")

wilcox.test(coefficients_summary[coefficients_summary$trend == "stabilising",]$median.amr,coefficients_summary[coefficients_summary$trend == "increasing (s)",]$median.amr,paired=FALSE)

