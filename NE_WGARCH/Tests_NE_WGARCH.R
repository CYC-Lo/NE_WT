setwd('~/Documents/MSc_thesis/NE_WGARCH/')

#==================================
# NE_WGARCH_sampler
#==================================
source('NE_WGARCH_sampler.R')
set.seed(1337)
sim = NE_WGARCH_sampler(t = 10000, alpha1 = 1.1, gam1 = 1.5, beta1 = 1.4,
                        alpha2 = 1.7, gam2 = 2.1, beta2 = 2,
                        w = 0.3,a = 0.6, b= 0.2, sigma0 = 0.5)

sim_df = data.frame(matrix(c(sim$return, sim$t), ncol  = 2))
colnames(sim_df) = c('return','t')



#==================================
# NE_WGARCH_loglik
#==================================
source('NE_WGARCH_loglik.R')


# alpha1,gam1,beta1,alpha2,gam2,beta2,
#                    pre_w,pre_a,pre_b,sigma0
result1 = optim(par = c(1,1,1,1,1,1,
                 1,1,2,1),
              fn = NE_WGARCH_loglik,gr = NE_WGARCH_partial_loglik,
              r0 = sim$return[1], r = sim$return[-1],control = list(fnscale = -1))


result = NE_WGARCH_MLE(data = sim_df)
library(ggplot2)
library(latex2exp)
x_11 = seq(-5,5, length = 4000)
# ggplot()+
#   geom_density(aes(x = sim$x),colour = 'red')+
#   geom_line(aes(x = x_11,
#                 y = NE_STD_PDF(x_11, 
#                                alpha1 = result$par[1],gam1 = result$par[2],
#                                beta1 = result$par[3],alpha2 = result$par[4],
#                                gam2 = result$par[5],beta2 = result$par[6])))


#==================================
# NE_WGARCH_diag
#==================================
source('NE_WGARCH_diag.R')
sigma_2 = NE_WGARCH_sig(result, sim_df)

ggplot()+
  geom_ribbon(aes(x = sim_df$t[1:2000],
                  ymin = -sqrt(sigma_2[1:2000]), ymax = sqrt(sigma_2[1:2000])
                  ), colour = 'red',fill = 'red')+
  geom_point(aes(x = sim_df$t[1:2000], y = sim_df$return[1:2000]))

# ggplot()+
#   geom_ribbon(aes(x = sim_df$t[99500:100000],
#                   ymin = -sqrt(sigma_2[99500:100000]), ymax = sqrt(sigma_2[99500:100000])
#   ), colour = 'red',fill = 'red', alpha = 0.8)+
#   geom_line(aes(x = sim_df$t[99500:100000], y = sim_df$return[99500:100000]))


x_est = sim_df$return[-1]/sqrt(sigma_2[-1])
ggplot()+
  geom_line(aes(x = x_11,
                y = NE_STD_PDF(x_11,
                               alpha1 = result$par[1],gam1 = result$par[2],
                               beta1 = result$par[3],alpha2 = result$par[4],
                               gam2 = result$par[5],beta2 = result$par[6]),
                colour = 'MLE'))+
  geom_density(aes(x = sim$x, colour = 'Residuals'))#+
  # xlab('x')+
  # scale_y_continuous(name = TeX('$P(X = x)$'))+
  # scale_color_manual()

