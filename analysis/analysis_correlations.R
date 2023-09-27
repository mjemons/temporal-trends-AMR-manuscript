# rough preliminary analysis of correlations

#all correaltions are not defined?

ggplot(all.correlations, aes(x=rho,y=correlate))+
  geom_point(position=position_dodgev(height=0.5),aes(colour=Pathogen))+
  facet_grid(trend~ .)

ggplot(all.correlations, aes(x=rho,colour=Pathogen,fill = Pathogen))+
  geom_histogram(alpha=0.5, position="identity",binwidth = 0.05)+
  facet_grid(trend~ correlate)

ddply(all.correlations,.(Pathogen,trend,correlate),summarise,mean=mean(rho))