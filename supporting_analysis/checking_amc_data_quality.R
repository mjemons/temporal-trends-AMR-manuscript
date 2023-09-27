library(plyr)
library(ggplot2)
library(gridExtra)

#### With original AMC data ##### 

amc_summary <- read.csv(file = "data/summary_AMC_byclass.csv")

# plot all community consumption separately
# add mean and median

amc_sub = amc_summary[amc_summary$Antimicrobial.Type=='Oral'&amc_summary$Antibiotic_class!="A07A"&amc_summary$Sector=='Community',]
amc_sub$combC <- paste(amc_sub$Antibiotic_class, amc_sub$Country, amc_sub$Sector, sep = "|")
all_combC = unique(amc_sub$combC)

# plot between and within country variation

in.country.stats = ddply(amc_sub,.(Country,Antibiotic_class),summarise,mean = mean(DDD),sd = sd(DDD))
between.country.stats = ddply(in.country.stats,.(Antibiotic_class),summarise,sd = mean(sd))

pdf("supporting_analysis/output/all_consumption_community.pdf", paper = "a4", width = 0, height = 0)
par(mfrow = c(4, 3))
par(las = 1)
for(i in 1:length(all_combC)){
  mycomb <- all_combC[i]
  tmp1 <- amc_sub[which(amc_sub$combC == mycomb),]
  print(mycomb)
  
  plot(tmp1$Year, tmp1$DDD, ylim= c(0,max(tmp1$DDD,na.rm=T)), col = "red", pch = 20, type = "o", ylab = "DDD", xlab = "year", main = mycomb)
  abline(h = mean(tmp1$DDD),col = "cadetblue4")
  abline(h = median(tmp1$DDD),col ="cadetblue3")
}

dev.off() 

# plot countries together. Lazy coding, colours are not consistent across subplots

amc_sub = amc_summary[amc_summary$Antimicrobial.Type=='Oral'&amc_summary$Antibiotic_class!="A07A"&amc_summary$Sector=='Community',]
amc_sub$combC <- paste(amc_sub$Antibiotic_class, amc_sub$Sector, sep = "|")
all_combC = unique(amc_sub$combC)

pdf("supporting_analysis/output/all_consumption_community_countries.pdf", paper = "a4", width = 0, height = 0)
par(mfrow = c(4, 3))
par(las = 1)
for(i in 1:length(all_combC)){
  mycomb <- all_combC[i]
  print(mycomb)
  tmp1 <- amc_sub[which(amc_sub$combC == mycomb),]
  countries = unique(tmp1$Country)
  cols = rainbow(length(countries))
  plot(0,0, ylim= c(0,max(tmp1$DDD,na.rm = T)), xlim = c(min(tmp1$Year),max(tmp1$Year)), pch = 20, type = "o", ylab = "DDD", xlab = "year", main = mycomb)
  for (j in 1:length(countries)){
    tmp2 = tmp1[which(tmp1$Country == countries[j]),]
    lines(tmp2$Year,tmp2$DDD,col = cols[j])
  }
  #  tmp1 <- tmp[which(tmp$patientType=="INPAT"),]
  #  tmp2 <- tmp[which(tmp$patientType=="OUTPAT"),]
  #  if(any(duplicated(tmp1$Year))) stop()
  #  if(any(duplicated(tmp2$Year))) stop()
}

dev.off()

# ok, now look at variability in data

amc_sub = amc_summary[amc_summary$Antimicrobial.Type=='Oral'&amc_summary$Antibiotic_class!="A07A"&amc_summary$Sector=='Community',]
amc_sub$combC <- paste(amc_sub$Antibiotic_class, amc_sub$Sector, sep = "|")
all_combC = unique(amc_sub$combC)

p = list()
for(i in 1:length(all_combC)){
  mycomb <- all_combC[i]
  print(mycomb)
  tmp1 <- amc_sub[which(amc_sub$combC == mycomb),]
  p[[i]] = ggplot(tmp1, aes(x=reorder(Country, DDD, FUN = median), y=DDD, fill=Country))+
    geom_boxplot()+
    ggtitle(mycomb)+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank())
}


ggsave(
  filename = "supporting_analysis/output/all_consumption_community_medians.pdf", 
  plot = marrangeGrob(p, nrow=1, ncol=1), 
  width = 10, height = 6
)

##### With cleaned data ######

# repeat plots with cleaned data

amc_summary <- read.csv(file = "data/summary_AMC_byclass_filtered.csv")

amc_sub = amc_clean[amc_clean$Antimicrobial.Type=='Oral'&amc_clean$Antibiotic_class!="A07A"&amc_clean$Sector=='Community',]
amc_sub$combC <- paste(amc_sub$Antibiotic_class, amc_sub$Country, amc_sub$Sector, sep = "|")
all_combC = unique(amc_sub$combC)

# plot between and within country variation

in.country.stats = ddply(amc_sub,.(Country,Antibiotic_class),summarise,mean = mean(DDD),sd = sd(DDD))
between.country.stats = ddply(in.country.stats,.(Antibiotic_class),summarise,sd = mean(sd))

pdf("supporting_analysis/output/all_consumption_community_clean.pdf", paper = "a4", width = 0, height = 0)
par(mfrow = c(4, 3))
par(las = 1)
for(i in 1:length(all_combC)){
  mycomb <- all_combC[i]
  tmp1 <- amc_sub[which(amc_sub$combC == mycomb),]
  print(mycomb)
  
  plot(tmp1$Year, tmp1$DDD, ylim= c(0,max(tmp1$DDD,na.rm = T)), col = "red", pch = 20, type = "o", ylab = "DDD", xlab = "year", main = mycomb)
  abline(h = mean(tmp1$DDD),col = "cadetblue4")
  abline(h = median(tmp1$DDD),col ="cadetblue3")
}

dev.off()

# plot countries together. Lazy coding, colours are not consistent across subplots

amc_sub = amc_clean[amc_clean$Antimicrobial.Type=='Oral'&amc_clean$Antibiotic_class!="A07A"&amc_clean$Sector=='Community',]
amc_sub$combC <- paste(amc_sub$Antibiotic_class, amc_sub$Sector, sep = "|")
all_combC = unique(amc_sub$combC)

pdf("supporting_analysis/output/all_consumption_community_countries_clean.pdf", paper = "a4", width = 0, height = 0)
par(mfrow = c(4, 3))
par(las = 1)
for(i in 1:length(all_combC)){
  mycomb <- all_combC[i]
  print(mycomb)
  tmp1 <- amc_sub[which(amc_sub$combC == mycomb),]
  countries = unique(tmp1$Country)
  cols = rainbow(length(countries))
  plot(0,0, ylim= c(0,max(tmp1$DDD)), xlim = c(min(tmp1$Year),max(tmp1$Year)), pch = 20, type = "o", ylab = "DDD", xlab = "year", main = mycomb)
  for (j in 1:length(countries)){
    tmp2 = tmp1[which(tmp1$Country == countries[j]),]
    lines(tmp2$Year,tmp2$DDD,col = cols[j])
  }
}

dev.off()

# ok, now look at variability in data

amc_sub = amc_clean[amc_clean$Antimicrobial.Type=='Oral'&amc_clean$Antibiotic_class!="A07A"&amc_clean$Sector=='Community',]
amc_sub$combC <- paste(amc_sub$Antibiotic_class, amc_sub$Sector, sep = "|")
all_combC = unique(amc_sub$combC)

p = list()
for(i in 1:length(all_combC)){
  mycomb <- all_combC[i]
  print(mycomb)
  tmp1 <- amc_sub[which(amc_sub$combC == mycomb),]
  p[[i]] = ggplot(tmp1, aes(x=reorder(Country, DDD, FUN = median), y=DDD, fill=Country))+
    geom_boxplot()+
    ggtitle(mycomb)+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank())
}


ggsave(
  filename = "supporting_analysis/output/all_consumption_community_medians_clean.pdf", 
  plot = marrangeGrob(p, nrow=1, ncol=1), 
  width = 10, height = 6
)



