source('code/Utils.R')

### generate summary statistics for the fits

# get the fit details
all.fits = read.csv(file = "data/results/all_fit_detailsINPAT.csv")

print(paste('Increasing trend overall',round(mean(all.fits$m1.slope>0),2)))
print(paste('Sig. increasing trend',round(mean(all.fits$m1.slope>0&all.fits$m1.slope.pval<0.05),2)))
print(paste('Sig. decreasing trend',round(mean(all.fits$m1.slope<0&all.fits$m1.slope.pval<0.05),2)))
print(paste('Median decreasing slope',round(median(all.fits[all.fits$m1.slope<0,]$m1.slope),3)))
print(paste('Median increasing slope',round(median(all.fits[all.fits$m1.slope>0,]$m1.slope),3)))

# get summary 
summary_fits = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")
print(round(table(summary_fits$category)/sum(table(summary_fits$category)),2))

# rename for plotting + get rid of poor fits
summary_fits_toplot <- summary_fits[summary_fits$category!='poor fit',]
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns decreasing',]$category = "decreasing (ns)"
summary_fits_toplot[summary_fits_toplot$category == 'logistic s decreasing',]$category = 'decreasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns increasing',]$category = 'increasing (ns)' 
summary_fits_toplot[summary_fits_toplot$category == 'logistic s increasing',]$category = 'increasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'flat ',]$category = 'stable'
summary_fits_toplot[summary_fits_toplot$category == 'plateauing',]$category = 'stabilising'

summary_fits_toplot$category = factor(summary_fits_toplot$category, levels = c('decreasing (s)',"decreasing (ns)",'stable','stabilising','increasing (ns)','increasing (s)'))

# analyse changes:
# increasing categories re-categorised
increasing.regategorised = table(summary_fits_toplot[summary_fits_toplot$combR %in% all.fits[all.fits$m1.slope>0,]$combR,]$category)/sum(table(summary_fits_toplot[summary_fits_toplot$combR %in% all.fits[all.fits$m1.slope>0,]$combR,]$category))
print(paste('Logistic increasing re-categorised as stable or stabilising:',round(increasing.regategorised["stable"] + increasing.regategorised["stabilising"],2)))

# increasing (s) categories re-categorised
increasing.s.regategorised = table(summary_fits_toplot[summary_fits_toplot$combR %in% all.fits[all.fits$m1.slope>0&all.fits$m1.slope.pval<0.05,]$combR,]$category)/sum(table(summary_fits_toplot[summary_fits_toplot$combR %in% all.fits[all.fits$m1.slope>0&all.fits$m1.slope.pval<0.05,]$combR,]$category))
print(paste('Logistic increasing (s) re-categorised as stabilising:',round(increasing.s.regategorised["stabilising"],2)))

# increasing (s) categories re-categorised for E.coli
increasing.s.regategorised = table(summary_fits_toplot[summary_fits_toplot$combR %in% all.fits[all.fits$m1.slope>0&all.fits$m1.slope.pval<0.05&all.fits$Pathogen == "ESCCOL",]$combR,]$category)/sum(table(summary_fits_toplot[summary_fits_toplot$combR %in% all.fits[all.fits$m1.slope>0&all.fits$m1.slope.pval<0.05&all.fits$Pathogen == "ESCCOL",]$combR,]$category))
print(paste('Logistic increasing (s) re-categorised as stabilising E.coli only:',round(increasing.s.regategorised["stabilising"],2)))
