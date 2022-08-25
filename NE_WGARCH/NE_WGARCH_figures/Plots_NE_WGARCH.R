setwd('~/Documents/MSc_thesis/NE_WGARCH/NE_WGARCH_figures/')

library(ggplot2)
library(latex2exp)
library(data.table)
# Sample from distribution
source('../NE_WGARCH_sampler.R', chdir = TRUE)
set.seed(1337)
sim = NE_WGARCH_sampler(t = 4000, alpha1 = 1.1, gam1 = 1.5, beta1 = 1.4,
                        alpha2 = 1.7, gam2 = 2.1, beta2 = 2,
                        w = 0.3,a = 0.6, b= 0.2, sigma0 = 0.5)

sim_df = data.frame(matrix(c(sim$return, sim$t), ncol  = 2))
colnames(sim_df) = c('return','t')
sim_df$return = sim_df$return - mean(sim_df$return)

ggplot(data = sim_df)+
  geom_line(aes(x = t, y = return))+
  xlab(TeX('$t$') )+
  ylab(TeX('Return ($r_t$ )'))
ggsave('NE_WGARCH_sample_series.pdf', scale = 0.8, device = 'pdf')



# Fit MLE

#good at 10000
# c(2,2,2,2,2,2,
#   0.5,0.5,0.5,0.1)

#good at 4000
# c(2,2,2,3,3,3,
  # 0.5,0.2,0.7,0.1)
# c(1.1,1.1,1.1,1.1,1.1,1.1,
#   0.5,0.2,0.7,0.1)
source('../NE_WGARCH_loglik.R',chdir = TRUE)
result1 = NE_WGARCH_MLE(data = sim_df, init_params = c(1.1,1.1,1.1,1.1,1.1,1.1,
                                                      0.5,0.2,0.7,1),
                       maxit = 5000)
result1 
set.seed(1337)
result2 = NE_WGARCH_MLE(data = sim_df, init_params = c(1.1,1.1,1.1,1.1,1.1,1.1,
                                                       0.5,0.2,0.7,1),
                        maxit = 10000, method = 'SANN')
result2
result3 = NE_WGARCH_MLE(data = sim_df, init_params = c(1.1,1.1,1.1,1.1,1.1,1.1,
                                                       0.5,0.2,0.7,1),
                        maxit = 10000, method = 'BFGS')
result3


result1.5 =  NE_WGARCH_MLE(data = sim_df, init_params = result2$par,
                           maxit = 5000)
result1.5
result2.5 =  NE_WGARCH_MLE(data = sim_df, init_params = result2$par,
                           maxit = 10000, method = 'SANN')
result2.5
result3.5 = NE_WGARCH_MLE(data = sim_df, init_params = result2$par,
                          maxit = 10000, method = 'BFGS')
result3.5



# Plot results

source('../NE_WGARCH_diag.R',chdir = TRUE)
x_est_1 = sim_df$return[-1]/sqrt(NE_WGARCH_sig(result1,sim_df)[-1])
x_est_2 = sim_df$return[-1]/sqrt(NE_WGARCH_sig(result2,sim_df)[-1])
x_est_3 = sim_df$return[-1]/sqrt(NE_WGARCH_sig(result3,sim_df)[-1])
x_est_1.5 = sim_df$return[-1]/sqrt(NE_WGARCH_sig(result1.5,sim_df)[-1])
x_est_2.5 = sim_df$return[-1]/sqrt(NE_WGARCH_sig(result2.5,sim_df)[-1])
x_est_3.5 = sim_df$return[-1]/sqrt(NE_WGARCH_sig(result3.5,sim_df)[-1])


df_res = data.table(Neld = x_est_1,
                    SANN = x_est_2,
                    BFGS = x_est_3,
                    SANN_Neld = x_est_1.5,
                    SANN_SANN = x_est_2.5,
                    SANN_BFGS = x_est_3.5,
                    type = 'residuals')
df_res = melt.data.table(df_res, id.vars = 'type')
df_res
x_11 = seq(-5,6, length = 5000)
df_MLE = data.table(Neld = NE_STD_PDF(x_11,
                                      alpha1 = result1$par[1],gam1 = result1$par[2],
                                      beta1 = result1$par[3],alpha2 = result1$par[4],
                                      gam2 = result1$par[5],beta2 = result1$par[6]),
                    SANN = NE_STD_PDF(x_11,
                                      alpha1 = result2$par[1],gam1 = result2$par[2],
                                      beta1 = result2$par[3],alpha2 = result2$par[4],
                                      gam2 = result2$par[5],beta2 = result2$par[6]),
                    BFGS = NE_STD_PDF(x_11,
                                      alpha1 = result3$par[1],gam1 = result3$par[2],
                                      beta1 = result3$par[3],alpha2 = result3$par[4],
                                      gam2 = result3$par[5],beta2 = result3$par[6]),
                    SANN_Neld = NE_STD_PDF(x_11,
                                      alpha1 = result1.5$par[1],gam1 = result1.5$par[2],
                                      beta1 = result1.5$par[3],alpha2 = result1.5$par[4],
                                      gam2 = result1.5$par[5],beta2 = result1.5$par[6]),
                    SANN_SANN = NE_STD_PDF(x_11,
                                           alpha1 = result2.5$par[1],gam1 = result2.5$par[2],
                                           beta1 = result2.5$par[3],alpha2 = result2.5$par[4],
                                           gam2 = result2.5$par[5],beta2 = result2.5$par[6]),
                    SANN_BFGS = NE_STD_PDF(x_11,
                                           alpha1 = result3.5$par[1],gam1 = result3.5$par[2],
                                           beta1 = result3.5$par[3],alpha2 = result3.5$par[4],
                                           gam2 = result3.5$par[5],beta2 = result3.5$par[6]),
                    x = x_11,
                    type = 'MLE'
                    )
df_MLE = melt.data.table(df_MLE, id.vars = c('x','type') )
df_MLE
ggplot()+
  geom_density(data = df_res,aes(x = value,colour = type))+
  geom_line(data = df_MLE, aes(x = x, y = value, colour = type))+
  scale_color_manual(name = 'Data Type',values = c('#F8766D','black'))+
  scale_y_continuous(name = TeX('$P(X = x)$'))+
  xlab('x')+
  facet_wrap(~variable)
ggsave('NE_WGARCH_sample_methods.pdf', scale = 0.8, device = 'pdf')

exp_labels <- function(breaks){
  labels <- sprintf("log(%f)", signif(exp(breaks), 7) ) 
  return(labels)
}
ggplot()+
  geom_density(data = df_res,aes(x = value,colour = type))+
  geom_line(data = df_MLE, aes(x = x, y = value, colour = type))+
  scale_color_manual(name = 'Data Type',values = c('#F8766D','black'))+
  scale_y_continuous(name = TeX('$P(X = x)$'),
                     trans = 'log', labels = exp_labels )+
  xlab('x')+
  facet_wrap(~variable)+
  coord_fixed(ylim = c(1e-2,1), xlim = c(-2.15,2.15))
ggsave('NE_WGARCH_sample_methods_log.pdf', scale = 0.8, device = 'pdf')


df_res_known = data.table(Neld = sim$x,
                          SANN = sim$x,
                          BFGS = sim$x,
                          SANN_Neld = sim$x,
                          SANN_SANN = sim$x,
                          SANN_BFGS = sim$x,
                          type = 'residuals')
df_res_known = melt.data.table(df_res_known, id.vars = 'type')
ggplot()+
  geom_density(data = df_res_known, aes(x = value, colour = type))+
  geom_line(data = df_MLE, aes(x = x, y = value, colour = type))+
  scale_color_manual(name = 'Data Type',values = c('#F8766D','black'))+
  scale_y_continuous(name = TeX('$P(X = x)$'))+
  xlab('x')+
  facet_wrap(~variable)
ggsave('NE_WGARCH_sample_methods_known.pdf', scale = 0.8, device = 'pdf')



#=========================
temp = 1
acf(x_est_1,plot = FALSE)$acf %>% head
acf(x_est_2,plot = FALSE)$lag %>% tail
df_acf = data.table(lag = c(acf(x_est_1,plot = FALSE)$lag) ,
                    Neld =c(acf(x_est_1,plot = FALSE)$acf), 
                    SANN =c(acf(x_est_2,plot = FALSE)$acf), 
                    BFGS =c(acf(x_est_3,plot = FALSE)$acf), 
                    SANN_Neld =c(acf(x_est_1.5,plot = FALSE)$acf), 
                    SANN_SANN =c(acf(x_est_2.5,plot = FALSE)$acf),
                    SANN_BFGS =c(acf(x_est_3.5,plot = FALSE)$acf))

df_acf = melt.data.table(df_acf, id.vars = 'lag')
df_acf



ggplot(df_acf)+
  geom_bar(aes(x = lag, y = value, fill = variable),stat = 'identity')+
  facet_wrap(~variable,scale = 'free')+
  ylab('ACF')+
  guides(fill = 'none')
ggsave('NE_WGARCH_sample_ACF.pdf', scale = 0.8, device = 'pdf')


ggplot(df_acf[lag != 0])+
  geom_bar(aes(x = lag, y = value, fill = variable),stat = 'identity')+
  facet_wrap(~variable,scale = 'free')+
  ylab('ACF')+
  guides(fill = 'none')
ggsave('NE_WGARCH_sample_ACF0.pdf', scale = 0.8, device = 'pdf')





#Unidentifiable plot
ggplot(df_res)+
  geom_line(aes(x = rep(1:length(x_est_1),nlevels(variable)),
                y = value,
                group = variable,
                colour = variable))+
  xlim(3800,4000)+
  facet_wrap(~variable)


#=======================
source('../NE_WGARCH_loglik_pen.R',chdir = TRUE)
result1_pen = NE_WGARCH_MLE_pen(data = sim_df, init_params = c(1.1,1.1,1.1,1.1,1.1,1.1,
                                                       0.5,0.2,0.7,1),
                        maxit = 5000)
result1_pen

set.seed(1337)
result2_pen = NE_WGARCH_MLE_pen(data = sim_df, init_params = c(1.1,1.1,1.1,1.1,1.1,1.1,
                                                       0.5,0.2,0.7,1),
                        maxit = 10000, method = 'SANN')
result2_pen
result3_pen = NE_WGARCH_MLE_pen(data = sim_df, init_params = c(1.1,1.1,1.1,1.1,1.1,1.1,
                                                       0.5,0.2,0.7,1),
                        maxit = 10000, method = 'BFGS')
result3_pen






