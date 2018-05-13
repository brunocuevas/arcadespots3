library('flexsurv')
setwd('../Dropbox/TFG//ArcadeSpots3/')
early_sen = read.csv('early_senescence_events.csv', header = TRUE)
mean_days = aggregate(early_sen$DAYS, by=list(early_sen$PATHO), FUN=mean)
early_sen['days_corrected'] = 100 * early_sen['DAYS'] / as.numeric(mean_days[mean_days['Group.1'] == 'C',]['x'])
early_sen['event'] = 1

#rrr = flexsurvreg(Surv(early_sen_control$DAYS, early_sen_control$event) ~ 1, dist='lnorm')

rrr = flexsurvreg(Surv(early_sen_control$DAYS, early_sen_control$event) ~ 1, dist='lnorm')
early_sen_control =  early_sen[early_sen$PATHO == 'C',]
rrr = flexsurvreg(Surv(early_sen$days_corrected, early_sen$event) ~ early_sen$PATHO, dist='lnorm')

plot(rrr, col = c(1,2,3,4), xlab='days', ylab='alive')
legend(
  x = 'topright',  
  legend = c('CONTROL', 'P0', 'P12', 'P123'), 
  col = c(1,2,3,4), lwd=3)


