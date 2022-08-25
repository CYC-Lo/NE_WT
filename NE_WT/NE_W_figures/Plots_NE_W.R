setwd('~/Documents/MSc_thesis/NE_WT/NE_W_figures/')


source('../NE_W_PDF.R',chdir = TRUE)
# sets of parameters:
#  alpha,gam,beta
baseline = c(1,2,1.5)
alpha_change = c(0.5,baseline[1],2)
gam_change = c(1,baseline[2],3)
beta_change = c(0.5,baseline[3],2.5)

x = (1:4000)/1000


library(dplyr)
library(reshape2)


library(latex2exp)
library(ggplot2)

#================================
# ALPHA
#================================

alpha_df = data.frame(a = NE_W_PDF(x,alpha = 0.5, gam = 2, beta = 1.5),
                      b= NE_W_PDF(x,alpha = 1, gam = 2, beta = 1.5),
                      c=NE_W_PDF(x,alpha = 2, gam = 2, beta = 1.5),
                      x = x)
alpha_df = melt(alpha_df,id.vars= 'x') 
colnames(alpha_df) <- c('x','alpha','pdf')
levels(alpha_df$alpha) <- c(0.5,1,2)



ggplot(alpha_df)+
  geom_line(aes(x = x, y =pdf, group = alpha, colour = alpha ))+
  ylim(0,3)+
  scale_color_discrete(name = TeX('$\\alpha$'))+
  ylab(TeX('$m(x|\\alpha, \\gamma, \\beta)$'))+
  theme(legend.position=c(0.70,0.70))

ggsave('NE_W_alpha.pdf', scale = 0.4, device = 'pdf')


#================================
# GAMMA
#================================

gamma_df = data.frame(a = NE_W_PDF(x,alpha = 1, gam = 1, beta = 1.5),
                      b= NE_W_PDF(x,alpha = 1, gam = 2, beta = 1.5),
                      c=NE_W_PDF(x,alpha = 1, gam = 3, beta = 1.5),
                      x = x)
gamma_df = melt(gamma_df,id.vars= 'x') 
colnames(gamma_df) <- c('x','gamma','pdf')
levels(gamma_df$gamma) <- c(1,2,3)



ggplot(gamma_df)+
  geom_line(aes(x = x, y =pdf, group = gamma, colour = gamma ))+
  ylim(0,3)+
  scale_color_discrete(name = TeX('$\\gamma$'))+
  ylab(TeX('$m(x|\\alpha, \\gamma, \\beta)$'))+
  theme(legend.position=c(0.70,0.70))


ggsave('NE_W_gamma.pdf', scale = 0.4, device = 'pdf')




#================================
# BETA
#================================
beta_df = data.frame(a = NE_W_PDF(x,alpha = 1, gam = 2, beta = 0.5),
                      b= NE_W_PDF(x,alpha = 1, gam = 2, beta = 1.5),
                      c=NE_W_PDF(x,alpha = 1, gam = 2, beta = 2.5),
                      x = x)
beta_df = melt(beta_df,id.vars= 'x') 
colnames(beta_df) <- c('x','beta','pdf')
levels(beta_df$beta) <- c(1,2,3)



ggplot(beta_df)+
  geom_line(aes(x = x, y =pdf, group =  beta, colour = beta ))+
  ylim(0,3)+
  scale_color_discrete(name = TeX('$\\beta$'))+
  ylab(TeX('$m(x|\\alpha, \\gamma, \\beta)$'))+
  theme(legend.position=c(0.70,0.70))


ggsave('NE_W_beta.pdf', scale = 0.4, device = 'pdf')
# scale_y_continuous(trans='log2')


