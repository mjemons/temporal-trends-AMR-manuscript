# this script to compares the logistic fit with a linear fit (?)
# it might be oudated
# note: have not checked this runs

first_model_lm_org <- read.csv("data/first_model_lm.csv")
second_model_lm_org <- read.csv("data/second_model_lm.csv")

#make subcategories for barplot shading if lm fits better or worse
all_lm_coef = c()
tokeep1_plateau = tokeep1 %>% as.data.frame() %>% filter(tokeep1[,1] %in% index_plateau)
second_model_lm = second_model_lm_org %>% as.data.frame() %>% filter(second_model_lm_org[,1] %in% index_plateau) %>% mutate(index = index_plateau)
all_lm_coef = rbind(all_lm_coef, second_model_lm)
plateau_log = all_combR[second_model_lm[which(second_model_lm[,4]>tokeep1_plateau[,5]),1]]
plateau_lm = all_combR[second_model_lm[which(second_model_lm[,4]<tokeep1_plateau[,5]),1]]


tokeep1_no.plateau = tokeep1 %>% as.data.frame() %>% filter(tokeep1[,1] %in% index_no.plateau)
second_model_lm = second_model_lm_org %>% as.data.frame() %>% filter(second_model_lm_org[,1] %in% index_no.plateau) %>% mutate(index = index_no.plateau)
all_lm_coef = rbind(all_lm_coef, second_model_lm)
no.plateau_log = all_combR[second_model_lm[which(second_model_lm[,4]>tokeep1_no.plateau[,5]),1]]
no.plateau_lm = all_combR[second_model_lm[which(second_model_lm[,4]<tokeep1_no.plateau[,5]),1]]

tokeep1_nonsig.decrease =  tokeep1 %>% as.data.frame() %>% filter(tokeep1[,1] %in% index_nonsig.decrease)
first_model_lm = first_model_lm_org %>% as.data.frame() %>% filter(first_model_lm_org[,1] %in% index_nonsig.decrease) %>% mutate(index = index_nonsig.decrease)
all_lm_coef = rbind(all_lm_coef, first_model_lm)
nonsig.decrease_log = all_combR[first_model_lm[which(first_model_lm[,4]>tokeep1_nonsig.decrease[,5]),1]]
nonsig.decrease_lm = all_combR[first_model_lm[which(first_model_lm[,4]<tokeep1_nonsig.decrease[,5]),1]]

tokeep1_sig.decrease =  tokeep1 %>% as.data.frame() %>% filter(tokeep1[,1] %in% index_sig.decrease)
first_model_lm = first_model_lm_org %>% as.data.frame() %>% filter(first_model_lm_org[,1] %in% index_sig.decrease) %>% mutate(index = index_sig.decrease)
all_lm_coef = rbind(all_lm_coef, first_model_lm)
sig.decrease_log = all_combR[first_model_lm[which(first_model_lm[,4]>tokeep1_sig.decrease[,5]),1]]
sig.decrease_lm = all_combR[first_model_lm[which(first_model_lm[,4]<tokeep1_sig.decrease[,5]),1]]

tokeep1_unclear =  tokeep1 %>% as.data.frame() %>% filter(tokeep1[,1] %in% index_unclear)
second_model_lm = second_model_lm_org %>% as.data.frame() %>% filter(second_model_lm_org[,1] %in% index_unclear) %>% mutate(index = index_unclear)
all_lm_coef = rbind(all_lm_coef, second_model_lm)
unclear_log = all_combR[second_model_lm[which(second_model_lm[,4]>tokeep1_unclear[,5]),1]]
unclear_lm = all_combR[second_model_lm[which(second_model_lm[,4]<tokeep1_unclear[,5]),1]]

all_lm_coef <- all_lm_coef %>% as.data.frame() %>% arrange((index))
all_lm <- c(plateau_lm,no.plateau_lm,nonsig.decrease_lm,unclear_lm)

df_coefficients$lm_class = 'NA'
#create categories for barplot shading
df_coefficients[df_coefficients$combR %in% plateau_log,]$lm_class = 'Logistic model fits better'
df_coefficients[df_coefficients$combR %in% plateau_lm,]$lm_class = 'Linear model fits better'
df_coefficients[df_coefficients$combR %in% no.plateau_log & df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$lm_class = 'Logistic model fits better'
df_coefficients[df_coefficients$combR %in% no.plateau_lm & df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$lm_class = 'Linear model fits better'
df_coefficients[df_coefficients$combR %in% no.plateau_log & df_coefficients$combR %in% all_combR[tokeep1[stat.increasing,1]],]$lm_class = 'Logistic model fits better'
#this combination doesn't exist
#df_coefficients[df_coefficients$combR %in% no.plateau_lm & df_coefficients$combR %in% all_combR[tokeep1[stat.increasing,1]],]$lm_class = 'Linear model fits better'
df_coefficients[df_coefficients$combR %in% unclear_lm & df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$lm_class = 'Linear model fits better'
df_coefficients[df_coefficients$combR %in% unclear_log & df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$lm_class = 'Logistic model fits better'
#none of that category either
#df_coefficients[df_coefficients$combR %in% unclear_lm & df_coefficients$combR %in% all_combR[tokeep1[stat.increasing,1]],]$lm_class = 'Linear model fits better'
df_coefficients[df_coefficients$combR %in% unclear_log & df_coefficients$combR %in% all_combR[tokeep1[stat.increasing,1]],]$lm_class = 'Logistic model fits better'
df_coefficients[df_coefficients$combR %in% nonsig.decrease_log,]$lm_class = 'Logistic model fits better'
df_coefficients[df_coefficients$combR %in% nonsig.decrease_lm,]$lm_class = 'Linear model fits better'
df_coefficients[df_coefficients$combR %in% sig.decrease_log,]$lm_class = 'Logistic model fits better'
#doesn't exist either - makes sense it is stat decreasing
#df_coefficients[df_coefficients$combR %in% sig.decrease_lm,]$lm_class = 'Linear model fits better'

if(linear_analysis==TRUE && merge==TRUE){
  #here we merge the plateau value K and the y intercept of the linear model to have a less conservative estimate
  #of the number of actually plateaued cases we have in our data
  df_coefficients[df_coefficients$combR %in% plateau_lm,]$logistic.plateau = 'plateau'
  df_coefficients[df_coefficients$combR %in% no.plateau_lm & df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$logistic.plateau = 'plateau'
  #this combination doesn't exist
  #df_coefficients[df_coefficients$combR %in% no.plateau_lm & df_coefficients$combR %in% all_combR[tokeep1[stat.increasing,1]],]$lm_class = 'Linear model fits better'
  df_coefficients[df_coefficients$combR %in% unclear_lm & df_coefficients$combR %in% all_combR[tokeep1[trend.increasing,1]],]$logistic.plateau = 'plateau'
  #none of that category either
  #df_coefficients[df_coefficients$combR %in% unclear_lm & df_coefficients$combR %in% all_combR[tokeep1[stat.increasing,1]],]$lm_class = 'Linear model fits better'
  df_coefficients[df_coefficients$combR %in% nonsig.decrease_lm,]$logistic.plateau = 'plateau'
  
  df_coefficients$lm_coeff = 'NA'
  
  #set the variable K from the linear model to be the k from the plateau
  df_coefficients[df_coefficients$combR %in% all_combR[tokeep1[all.slopes,1]],]$lm_coeff = all_lm_coef[,3]
  df_coefficients[!df_coefficients$combR %in% all_lm,]$lm_coeff = NA
  
  #df_coefficients$mycol <- ifelse(df_coefficients$logistic.plateau.K, df_coefficients$lm_coeff)
  df_coefficients$z <- dplyr::coalesce(ifelse(df_coefficients$lm_coeff == "", NA, df_coefficients$lm_coeff), df_coefficients$logistic.plateau.K)
  df_coefficients$logistic.plateau.K <- df_coefficients$z
}