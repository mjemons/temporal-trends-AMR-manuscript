library(ggpubr)

rm(list = ls())

# get category data:
summary_fits = read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")

# rename for plotting + get rid of poor fits
summary_fits_toplot <- summary_fits[summary_fits$category!='poor fit',]
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns decreasing',]$category = "decreasing (ns)"
summary_fits_toplot[summary_fits_toplot$category == 'logistic s decreasing',]$category = 'decreasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns increasing',]$category = 'increasing (ns)' 
summary_fits_toplot[summary_fits_toplot$category == 'logistic s increasing',]$category = 'increasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'flat ',]$category = 'stable'
summary_fits_toplot[summary_fits_toplot$category == 'plateauing',]$category = 'stabilising'

# write function to produce plots

binary_plot = function(cats,ylabel,plot_title){

summary_fits_toplot$binary_class <- ifelse(summary_fits_toplot$category %in% cats, 1, 0)
summary_fits_toplot$Country <- factor(summary_fits_toplot$Country)

# for some countries, only 1s or 0s: find these:

temp = ddply(summary_fits_toplot,.(Country),summarise,mean(binary_class))
countries.to.remove = temp[temp[,2] == 0 | temp[,2] == 1,1]


mdl <- glm(formula = binary_class ~ Pathogen+Country+Antibiotic_class, family = "binomial", data = summary_fits_toplot)

tt <- broom::tidy(mdl,conf.int=TRUE)
tt <- dplyr::filter(tt, term!="(Intercept)")

tt$global_term <- ''
tt$local_term <- ''

for(i in 1:nrow(tt)){
  tmp <- unlist(strsplit(tt$term[i], "(?<=[a-z])(?=[A-Z])", perl = TRUE))
  tt$global_term[i] <- tmp[1]
  tt$local_term[i] <- tmp[2]
}

tt = tt[!tt$local_term %in% countries.to.remove,]

# change antibiotic_class to Antibiotic

tt[tt$global_term == 'Antibiotic_class',]$global_term = 'Antibiotic class'

# add confidence itnervals - error bars make it very confusing
p <- ggplot(tt, aes(y=estimate,x=reorder(local_term,estimate)))+
  geom_point(width=0.5,stat = 'identity',position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin=conf.low,ymax=conf.high),width=0)+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_y_continuous(limits=c(-5,5),oob=rescale_none) +
  facet_wrap(global_term~., nrow=1, ncol=3, drop=TRUE, scales = "free")+
  theme_light()+
  xlab('')+
  ylab(ylabel)+
  labs(color='Logistic Categories')+
  coord_flip()+
  ggtitle(plot_title)
p
}


cats = c('increasing (s)','increasing (ns)')
p1 = binary_plot(cats,'',"Rising")

cats = c('stable','stabalising')
p2 = binary_plot(cats,'', "Equilibrium")

cats = c('decreasing (s)','decreasing (ns)')
p3 = binary_plot(cats,'Coefficient','Declining')

figure <- ggarrange(p1, p2, p3,
                    ncol = 1, nrow = 3)
figure
ggsave(paste("output/binomialmodel_version2.png"), width=15, height=17, dpi=500)

