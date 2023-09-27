source('code/Utils.R')

### load the results of temporal trend fits and clean

summary_fits <- read.csv(file = "data/results/temporal_trend_AICcINPAT.csv")

# FB added summary including poor fits 27/07/2023
print(signif(table(summary_fits$category=="poor fit")/nrow(summary_fits), 2))

# rename for plotting + get rid of poor fits
summary_fits_toplot <- summary_fits[summary_fits$category!='poor fit',]
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns decreasing',]$category = "decreasing (ns)"
summary_fits_toplot[summary_fits_toplot$category == 'logistic s decreasing',]$category = 'decreasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'logistic ns increasing',]$category = 'increasing (ns)' 
summary_fits_toplot[summary_fits_toplot$category == 'logistic s increasing',]$category = 'increasing (s)'
summary_fits_toplot[summary_fits_toplot$category == 'flat ',]$category = 'stable'
summary_fits_toplot[summary_fits_toplot$category == 'plateauing',]$category = 'stabilising'

summary_fits_toplot$category = factor(summary_fits_toplot$category, levels = c('decreasing (s)',"decreasing (ns)",'stable','stabilising','increasing (ns)','increasing (s)'))

# FB added summary 27/07/2023
print(signif(table(summary_fits_toplot$category, useNA = "ifany")/nrow(summary_fits_toplot), 2))

# Will want to change names a little for plotting

freq_plot <- function(variable.of.interest, coefficients_summary){
  
  count_raw <- ddply(coefficients_summary, .(coefficients_summary[,variable.of.interest]), nrow)
  names(count_raw) <- c(variable.of.interest, "Total")
  
  counts <- ddply(coefficients_summary, .(coefficients_summary[,variable.of.interest],
                                          coefficients_summary$category), nrow)
  names(counts) <- c(variable.of.interest, "logistic", "Count")
  
  counts <- left_join(counts,count_raw, by=variable.of.interest)
  
  counts$Freq <- counts["Count"]/counts["Total"]
  
  # sort by frequency of resistance
  print(counts)
  temp= ddply(coefficients_summary,.(coefficients_summary[,variable.of.interest]),summarise,freq.res = median(median.amr))
  temp <- counts[counts[,'logistic'] %in% c('increasing (s)', 'increasing (ns)'),]
  temp <- temp %>% group_by(temp[[variable.of.interest]]) %>% dplyr::summarise(sum_freq = sum(Freq))
  level.order = temp[wrapr::orderv(temp[,'sum_freq']),1]
  level.order <- as.list(level.order)
  level.order <- level.order[[1]]
  
  #hacky way of adding countries where no increasing category is there for ordering
  if(variable.of.interest=='Country'){
    level.order = c('Luxembourg', 'Slovenia',level.order)
  }
  
  counts[,variable.of.interest] = factor(counts[,variable.of.interest], levels = level.order)
  print(level.order)
  p <- ggplot(counts, aes(x = counts[,variable.of.interest], y = Freq$Count , fill= logistic))+
    #, alpha = Shading)) +
    scale_alpha_discrete(range = c(0.4, 0.8)) +
    geom_bar(stat='identity')  +
    #scale_fill_manual(values=c("#D84315", "#C2185B","#0277BD", "#C0CA33", "#52854C"))+
    scale_fill_manual(values=c("pink4", "pink","wheat4", "wheat", "slategray1","slategray"))+
    ylab('Proportion') +
    geom_text(aes(label=counts$Count),size = 3, position = position_stack(vjust = 0.5),colour='white')+
    xlab(variable.of.interest)+
    labs(fill = "Classification of trend") +
    theme_light()+
    theme(axis.title=element_text(size=10,face="bold"),axis.text.x= element_text(angle = 45, vjust = 0.9, hjust=1, size = 10),legend.title=element_text(size=10,face="bold"), 
          legend.text=element_text(size=10),legend.position="none")
  if(variable.of.interest == 'Antibiotic_class'){
    p = p + xlab('Antibiotic Class')
  }
  if(variable.of.interest == 'Pathogen_long'){
    p = p + xlab('Pathogen')
  }
  p
  ggsave(paste0('output/trends_by_',variable.of.interest,'.svg'), width=25, height=15, dpi=500)
  print(p)
}

for(variable.of.interest in c('Pathogen_long', 'Country', 'Antibiotic_class')){
  pdf(paste0('output/trends_by_',variable.of.interest,'.pdf'), width = 8, height = 6)
  freq_plot(variable.of.interest, summary_fits_toplot)
  dev.off()
}

# freq_plot('Pathogen_long',summary_fits_toplot)
# freq_plot('Country',summary_fits_toplot)
# freq_plot('Antibiotic_class',summary_fits_toplot)


##### Analyse predictors of plateau #####

### multinomial analysis ###
summary_fits_toplot$category = factor(summary_fits_toplot$category, levels = c('increasing (s)','increasing (ns)','stabilising','stable',"decreasing (ns)",'decreasing (s)'))
test = multinom(formula = category ~ Pathogen+Country+Antibiotic_class, family = "binomial", data = summary_fits_toplot)
summary.multinom = summary(test)
car::Anova(test)

ttm <- broom::tidy(test,conf.int=TRUE)
ttm <- dplyr::filter(ttm, term!="(Intercept)")

ttm$global_term <- ''
ttm$local_term <- ''

for(i in 1:nrow(ttm)){
  tmp <- unlist(strsplit(ttm$term[i], "(?<=[a-z])(?=[A-Z])", perl = TRUE))
  ttm$global_term[i] <- tmp[1]
  ttm$local_term[i] <- tmp[2]
}

# add confidence intervals - error bars make it very confusing
p <- ggplot(ttm, aes(y=estimate,x=local_term,colour=y.level))+
  geom_point(width=0.4,stat = 'identity',position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin=conf.low,ymax=conf.high),position = position_dodge(width = 0.75),width=0)+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_y_continuous(limits=c(-10,10),oob=rescale_none) +
  facet_wrap(global_term~., nrow=2, ncol=2, drop=TRUE, scales = "free")+
  theme_light()+
  xlab('Category')+
  ylab('Multinomial Regression Coefficient')+
  labs(color='Logistic Categories')+
  coord_flip()
p
ggsave(paste("output/multinomialmodel.png"), width=15, height=18, dpi=500)

# likelihood ratio test for each variable: all are significant predictors of trend
lrtest(test,"Pathogen")
lrtest(test,"Antibiotic_class")
lrtest(test,"Country")
car::Anova(test)

### univariate approach ###

summary_fits_toplot$binary_class <- ifelse(summary_fits_toplot$category %in% c('increasing (s)', 'increasing (ns)'), 1, 0)

summary_fits_toplot$Country <- factor(summary_fits_toplot$Country)
#get rid of slovenia and luxembourg because they are not binarised and make little to no sense in the analysis
#summary_fits_toplot <- dplyr::filter(summary_fits_toplot, Country!="Slovenia")
#summary_fits_toplot <- dplyr::filter(summary_fits_toplot, Country!="Luxembourg")

mdl <- glm(formula = binary_class ~ Pathogen_long+Country+Antibiotic_class, family = "binomial", data = summary_fits_toplot)

print(mdl)
summary(mdl)

anova(mdl, test="Chisq")

tt <- broom::tidy(mdl,conf.int=TRUE)
tt <- dplyr::filter(tt, term!="(Intercept)")

tt$global_term <- ''
tt$local_term <- ''

for(i in 1:nrow(tt)){
  tmp <- unlist(strsplit(tt$term[i], "(?<=[a-z])(?=[A-Z])", perl = TRUE))
  tt$global_term[i] <- tmp[1]
  tt$local_term[i] <- tmp[2]
}

tt[tt$global_term == 'Antibiotic_class',]$global_term = 'Antibiotic class'
tt[tt$global_term == 'Pathogen_long',]$global_term = 'Pathogen'

# add confidence itnervals - error bars make it very confusing
p <- ggplot(tt, aes(y=estimate,x=reorder(local_term,estimate)))+
  geom_point(width=0.5,stat = 'identity',position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymin=conf.low,ymax=conf.high),width=0)+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_y_continuous(limits=c(-5,5),oob=rescale_none) +
  facet_wrap(global_term~., nrow=1, ncol=3, drop=TRUE, scales = "free")+
  theme_light()+
  xlab('Category')+
  ylab('Binomial Regression Coefficient')+
  labs(color='Logistic Categories')+
  coord_flip()
p
ggsave(paste("output/binomialmodel.pdf"), width=15, height=5, dpi=500)
write.csv(tt, "data/results/binomial_model.csv",sep='')

#plot slopes of logistic fits

all.fits.inpat <- read.csv(paste("data/results/all_fit_detailsINPAT.csv"))
all.fits.inpat$combR.long <- paste(all.fits.inpat$Pathogen, all.fits.inpat$Country ,all.fits.inpat$Antibiotic, "OUTPAT", sep = "|")
all.fits.outpat <- read.csv(paste("data/results/all_fit_detailsOUTPAT.csv"))
all.fits.outpat$combR.long <- paste(all.fits.outpat$Pathogen, all.fits.outpat$Country ,all.fits.outpat$Antibiotic, "INPAT", sep = "|")
all.fits <- rbind(all.fits.inpat, all.fits.outpat)

all.fits$trend = ''
all.fits[all.fits$m1.slope>0 & all.fits$m1.slope.pval < 0.05,]$trend = 'increasing'
all.fits[all.fits$m1.slope>0 & all.fits$m1.slope.pval > 0.05,]$trend = 'non-significant'
all.fits[all.fits$m1.slope<0 & all.fits$m1.slope.pval < 0.05,]$trend = 'decreasing'
all.fits[all.fits$m1.slope<0 & all.fits$m1.slope.pval > 0.05,]$trend = 'non-significant'

# FB added summary below 27/07/2023
print(signif(table(all.fits$trend, useNA = "ifany")/nrow(all.fits), 2))
print(signif(table(all.fits$m1.slope>0, useNA = "ifany")/nrow(all.fits), 2))
any(all.fits$m1.slope==0)

all.fits$Pathogenlong <- ''
all.fits[all.fits$Pathogen == 'ACISPP',]$Pathogenlong = 'Acinetobacter spp.'
all.fits[all.fits$Pathogen == 'ENCFAE',]$Pathogenlong = 'E. faecalis'
all.fits[all.fits$Pathogen == 'ENCFAI',]$Pathogenlong = 'E. faecium'
all.fits[all.fits$Pathogen == 'ESCCOL',]$Pathogenlong = 'E. coli'
all.fits[all.fits$Pathogen == 'KLEPNE',]$Pathogenlong = 'K. pneumoniae'
all.fits[all.fits$Pathogen == 'PSEAER',]$Pathogenlong = 'P. aeruginosa'
all.fits[all.fits$Pathogen == 'STAAUR',]$Pathogenlong = 'S. aureus'
all.fits[all.fits$Pathogen == 'STRPNE',]$Pathogenlong = 'S. pneumoniae'

theme_set(theme_gray())
#for(pathogen in pathogens){
#all.fits.sub.AB <- all.fits[all.fits$Pathogen == pathogen,]
all.fits.sub.AB <- all.fits[order(all.fits$m1.slope),]
p <- ggplot(all.fits.sub.AB, aes(x=reorder(combR.long,m1.slope), y=m1.slope, col = trend)) + 
  geom_point(stat = 'identity',alpha = 1) +
  geom_errorbar(aes(ymin=m1.slope.lower, ymax=m1.slope.upper),width=0,alpha = 1) +
  #coord_flip(ylim=c(-2,2)) +
  ylim(-1.75,1.75)+
  #ggtitle(pathogen)+
  theme_classic()+
  ylab('Slope of logistic fit') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),panel.spacing = unit(2, "lines"), legend.position="none")+ 
  scale_colour_manual(values = c('pink4','slategray','wheat'))+
  scale_x_discrete(expand = c(.02, .02))+
  geom_hline(yintercept=0, linetype="dashed")+
  facet_wrap(~Pathogenlong,scales = "free", nrow=4,ncol=3)
pdf(paste0('output/logisticslopesummary.pdf'), width = 10, height = 10)
print(p)
dev.off()  
print(p)
ggsave(paste("output/logisticslopesummary.png"), width=10, height=10, dpi=500)
