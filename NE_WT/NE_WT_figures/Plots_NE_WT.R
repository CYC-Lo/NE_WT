setwd('~/Documents/MSc_thesis/NE_WT/NE_WT_figures/')
source('../NE_STD_sampler.R',chdir = TRUE)
source('../NE_STD_loglik.R',chdir = TRUE)

source('../NE_STD_CDF.R',chdir = TRUE)
source('../NE_STD_PDF.R',chdir = TRUE)

library(ggplot2)
library(latex2exp)
set.seed(1337)
x = NE_STD_sampler(3000000 + 30000 +300,
                   alpha1 = 1,gam1 = 2,beta1 = 3,
                   alpha2 = 1.3,gam2 = 1.8,beta2 = 2)
y = NE_STD_sampler(300,
                   alpha1 = 1,gam1 = 2,beta1 = 3,
                   alpha2 = 1.3,gam2 = 1.8,beta2 = 2)
true_x = seq(min(x[,1]),max(x[,1]),len = 100000)

ggplot()+
  geom_line(aes(x = x[30001:3000000,1], y = x[30001:3000000,2],
                color = 'Samples: 3000000'),size = 1.2)+
  geom_line(aes(x = x[301:30000,1], y = x[301:30000,2],
                color = 'Samples: 30000'),size = 1.2)+
  geom_line(aes(x = x[1:300,1], y = x[1:300,2],
                color = 'Samples: 300'),size = 1.2)+
  geom_line(aes(x = true_x, 
                y = NE_STD_CDF(true_x,
                               alpha1 = 1,gam1 = 2,beta1 = 3,
                               alpha2 = 1.3,gam2 = 1.8,beta2 = 2),
                color = 'Theoretical'))+
  ylab(TeX('$P(X \\leq x)$') )+
  xlab(TeX('$x$') )+
  theme(legend.position = c(0.75,0.25))+
  scale_colour_manual(name = 'Data Type',
                      values = c('#F8766D','#00BA38','#619CFF','black'))
ggsave('NE_WT_sampler_cdf.pdf', scale = 0.8, device = 'pdf')





ggplot()+
  geom_density(aes(x = x[30001:3000000,1],
                   color = 'Samples: 3000000'),size = 1)+
  geom_density(aes(x = x[301:30000,1],
                   color = 'Samples: 30000'),size = 1)+
  geom_density(aes(x = x[1:301,1],
                   color = 'Samples: 300'),size = 1)+
  geom_line(aes(x = true_x, 
                y = NE_STD_PDF(true_x,
                               alpha1 = 1,gam1 = 2,beta1 = 3,
                               alpha2 = 1.3,gam2 = 1.8,beta2 = 2),
                color = 'Theoretical'))+
  ylab(TeX('$P(X = x)$') )+
  xlab(TeX('$x$') )+
  coord_cartesian(xlim = c(-8,8))+
  theme(legend.position = c(0.75,0.75))+
  scale_colour_manual(name = 'Data Type',
                      values = c('#F8766D','#00BA38','#619CFF','black'))


ggsave('NE_WT_sampler_pdf.pdf', scale = 0.8, device = 'pdf')



exp_labels <- function(breaks){
  labels <- sprintf("log(%f)", signif(exp(breaks), 7) ) 
  return(labels)
}
ggplot()+
  geom_density(aes(x = x[30001:3000000,1],
                   color = 'Samples: 3000000'),size = 1)+
  geom_density(aes(x = x[301:30000,1],
                   color = 'Samples: 30000'),size = 1)+
  geom_density(aes(x = x[1:300,1],
                   color = 'Samples: 300'),size = 1)+
  geom_line(aes(x = true_x, 
                y = NE_STD_PDF(true_x,
                               alpha1 = 1,gam1 = 2,beta1 = 3,
                               alpha2 = 1.3,gam2 = 1.8,beta2 = 2),
                color = 'Theoretical'))+
  ylab(TeX('$P(X = x)$') )+
  scale_y_continuous(trans='log', labels = exp_labels)+
  xlab(TeX('$x$') )+
  coord_cartesian(xlim = c(-5,7), ylim = c(1e-8,1))+
  theme(legend.position = c(0.8,0.8))+
  scale_colour_manual(name = 'Data Type',
                      values = c('#F8766D','#00BA38','#619CFF','black'))


ggsave('NE_WT_sampler_logpdf.pdf', scale = 0.8, device = 'pdf')


set.seed(200)
params = list(c(1.6,3,2.8, 1.1,2.3,4.1),
              c(1.2,1.3,5.2, 2.5,2,1.3),
              c(3.5,3,6, 2,3,6))
x_list = list()
mle_params = list()
for (i in 1:3) {
  params_i = params[[i]]
  x = NE_STD_sampler(30000,
                     alpha1 = params_i[1],gam1 = params_i[2],
                     beta1 = params_i[3],alpha2 = params_i[4],
                     gam2 = params_i[5],beta2 = params_i[6])
  x_list[[i]] = x
  mle_params[[i]] = Results = optim(par = c(2,2,2,2,2,3), fn = NE_STD_logLik,
                                    gr = NE_STD_partial_loglik, x = x[,1],
                                    control = list(fnscale = -1)
                                    ) 
}

true_x = seq(-10,10, len= 10000)
#Plot first param set
ggplot()+
  geom_line(aes(x = true_x, 
                y = NE_STD_PDF(true_x,
                               alpha1 = mle_params[[1]]$par[1],
                               gam1 = mle_params[[1]]$par[2],
                               beta1 = mle_params[[1]]$par[3],
                               alpha2 = mle_params[[1]]$par[4],
                               gam2 = mle_params[[1]]$par[5],
                               beta2 = mle_params[[1]]$par[6]),
                color = 'MLE'),size = 1)+
  geom_density(aes(x = x_list[[1]][,1], colour = 'Samples'))+
  ylab(TeX('$P(X= x)$'))+ xlab(TeX('$x$'))+
  scale_color_manual(name = 'Type', values = c('#F8766D','black'))
ggsave('NE_WT_sampler_mle1.pdf', scale = 0.4, device = 'pdf')

#Plot second param set
ggplot()+
  geom_line(aes(x = true_x, 
                y = NE_STD_PDF(true_x,
                               alpha1 = mle_params[[2]]$par[1],
                               gam1 = mle_params[[2]]$par[2],
                               beta1 = mle_params[[2]]$par[3],
                               alpha2 = mle_params[[2]]$par[4],
                               gam2 = mle_params[[2]]$par[5],
                               beta2 = mle_params[[2]]$par[6]),
                color = 'MLE'),size = 1)+
  geom_density(aes(x = x_list[[2]][,1], colour = 'Samples'))+
  ylab(TeX('$P(X= x)$'))+ xlab(TeX('$x$'))+
  scale_color_manual(name = 'Type', values = c('#F8766D','black'))
ggsave('NE_WT_sampler_mle2.pdf', scale = 0.4, device = 'pdf')

#Plot third param set
ggplot()+
  geom_line(aes(x = true_x, 
                y = NE_STD_PDF(true_x,
                               alpha1 = mle_params[[3]]$par[1],
                               gam1 = mle_params[[3]]$par[2],
                               beta1 = mle_params[[3]]$par[3],
                               alpha2 = mle_params[[3]]$par[4],
                               gam2 = mle_params[[3]]$par[5],
                               beta2 = mle_params[[3]]$par[6]),
                color = 'MLE'),size = 1)+
  geom_density(aes(x = x_list[[3]][,1], colour = 'Samples'))+
  ylab(TeX('$P(X= x)$'))+ xlab(TeX('$x$'))+
  scale_color_manual(name = 'Type', values = c('#F8766D','black'))
ggsave('NE_WT_sampler_mle3.pdf', scale = 0.4, device = 'pdf')



