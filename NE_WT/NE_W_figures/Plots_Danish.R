setwd('~/Documents/MSc_thesis/NE_WT/NE_W_figures/')


library(dplyr)
library(data.table)
library(ggplot2)
library(fitdistrplus)
library(reshape2)
source('../NE_W_loglik.R',chdir = TRUE)
source('../NE_W_PDF.R',chdir = TRUE)

data(danishuni,package = 'fitdistrplus')
# The univariate dataset was collected at Copenhagen Reinsurance.
# Comprise 2167 fire losses over the period 1980 to 1990.
# They have been adjusted for inflation to reflect 1985 values.
# Expressed in millions of Danish Krone(DKK).

danishuni = data.table(danishuni)

mle_wei = fitdist(danishuni$Loss, 'weibull')
mle_NE_W = NE_W_MLE(danishuni$Loss, init_params = c(1,1,1.5))


New_data = density(danishuni$Loss)
df = data.table(Loss = New_data$x,
                Data = New_data$y,
                Wei = dweibull(New_data$x,
                               shape = mle_wei$estimate['shape'],
                               scale = mle_wei$estimate['scale']),
                NE_W = NE_W_PDF(New_data$x,
                                alpha = mle_NE_W$par[1],
                                gam = mle_NE_W$par[2],
                                beta = mle_NE_W$par[3])
                )
df = melt.data.table(df,value.name = 'probability', variable.name = 'distr',
          id.vars = 'Loss') 
df$distr <- factor(df$distr, levels = c("Wei", "NE_W", "Data"))

ggplot(df)+
  geom_line(aes(x = Loss, y = probability,
                group = distr, colour = distr))+
  scale_y_continuous(trans = 'sqrt')+
  xlab('Loss')+ ylab('Square root Probability')+
  theme(legend.position = c(0.75,0.75))+
  scale_color_manual(name = 'Distribution',
                       breaks = c('Data', 'NE_W','Wei'),
                       labels = c('Data', 'NE_W','Weibull'),
                       values = c('#F8766D','#00BA38','#619CFF') 
                                  #red,green,blue for ggplot
  )+
  xlim(0,150)
  

ggsave('NE_W_Danish.pdf', scale = 0.8, device = 'pdf')




exp_labels <- function(breaks){
  labels <- sprintf("log(%f)", signif(exp(breaks), 7) ) 
  return(labels)
}

ggplot(df)+
  geom_line(data = df[df$distr != 'Data' & abs(exp(df$probability)- 1) > 2e-8 ],
            aes(x = Loss, y = probability,
                group = distr, colour = distr))+
  geom_point(data = df[df$distr == 'Data' & abs(exp(df$probability)- 1) > 2e-8 ],
             aes(x = Loss, y = probability,
                 group = distr, colour = distr))+
  scale_y_continuous(trans = 'log',labels = exp_labels)+
  scale_x_continuous(trans = 'log')+
  xlab('Log(Loss)')+ ylab('Probability')+
  theme(legend.position = c(0.25,0.25))+
  scale_color_discrete(name = 'Distribution',
                       labels = c('Data','NE_W', 'Weibull')
  )
ggsave('NE_W_Danish_log.pdf', scale = 0.8, device = 'pdf')


