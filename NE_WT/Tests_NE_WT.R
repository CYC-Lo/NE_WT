#Test for packages

setwd('~/Documents/MSc_thesis/NE_WT/')


#=============================
# NE_W_CDF
#=============================
source('NE_W_CDF.R')
library(ggplot2)
x = 1:10/10
alpha = 0.1
beta = 0.7
gam = 2
M = NE_W_CDF(alpha,gam,beta,x)
ggplot()+
  geom_line(mapping = aes(x = x, y = M))


#=============================
# NE_W_PDF
#=============================
source('NE_W_PDF.R')
library(ggplot2)
x = 1:10/10
alpha = 0.1
beta = 0.7
gam = 2
m = NE_W_PDF(alpha,gam,beta,x)
ggplot()+
  geom_line(mapping = aes(x = x, y = m))


#=============================
# CubicHerm_interp
#=============================
source('CubicHerm_interp.R')
#Test using unif(0,1):
# punif is CDF, dunif is PDF

#check construct_node
n_l = construct_node(0,punif,dunif)
n_r = construct_node(1,punif,dunif)

#check error_split
vals = error_split(n_l,n_r,CDF = punif,PDF = dunif)
#check error_split(null case)
vals_2 = error_split(vals[[1]],vals[[2]],punif,dunif)

#check loop_split(error)
vals = loop_split(c(list(n_l),list(n_r)),error_split,punif,dunif)

#check mono_split(null case)
vals_2 = mono_split(n_l,n_r,punif,dunif)


#check loop_split(mono)
vals_2 = loop_split(vals, mono_split,CDF = punif , dunif)

#check node_lst2mat
vals_2 = node_lst2mat(vals)

#==========================================
#Test using norm(0,1) at range [-1,1]
n_l = construct_node(-10,pnorm,dnorm)
n_r = construct_node(10,pnorm,dnorm)
vals = loop_split(c(list(n_l),list(n_r)),error_split,pnorm,dnorm)


vals = loop_split(vals, mono_split,CDF = pnorm, PDF = dnorm)
vals = node_lst2mat(vals)
vals = vals[,1:92]


spl = splinefunH(vals[2,],vals[1,],vals[3,]^-1)
library(ggplot2)
ggplot()+
  geom_line(aes(x = spl( (1:1000)/1000 ), y =  (1:1000)/1000 ), color = 'red', size = 1.2)+
  geom_line(aes(x = vals[1,],y = vals[2,]), color = 'blue')


#==========================================
#Test with NE_W
source('NE_W_CDF.R')
source('NE_W_PDF.R')
source('CubicHerm_interp.R')

n_l = construct_node(0,NE_W_CDF,NE_W_PDF, alpha = 1,gam = 2,beta = 3)
n_r = construct_node(20,NE_W_CDF,NE_W_PDF, alpha = 1,gam = 2,beta = 3)
vals = loop_split(c(list(n_l),list(n_r)),error_split,NE_W_CDF,NE_W_PDF,alpha = 1,gam = 2,beta = 3)


vals = loop_split(vals, mono_split,NE_W_CDF,NE_W_PDF,alpha = 1, gam = 2,beta = 3)
vals = node_lst2mat(vals)
ggplot()+
  geom_line(aes(x = vals[1,], y = vals[2,]), color= 'red',size = 1.2)+
  geom_line(aes(x = vals[1,], y = NE_W_CDF(vals[1,],alpha = 1, gam = 2, beta = 3)), color = 'green', size = 0.8)+
  geom_line(aes(x = vals[1,], y = construct_node(vals[1,],NE_W_CDF,NE_W_PDF,alpha = 1, gam = 2, beta = 3)[2,]), color = 'black', size = 0.5)


# cut_off_r = which(vals[2,] >= 1-1e-6)[1]
# cut_off_l = rev(which(vals[2,]== 0))[1]
# vals_2 = vals[,cut_off_l:cut_off_r]

# spl = splinefunH(vals_2[2,],vals_2[1,],vals_2[3,]^-1)
# library(ggplot2)
# ggplot()+
#   geom_line(aes(y = (0:200)/10, x = spl( (0:200)/10  ) ), color = 'red', size = 1.2)+
#   geom_line(aes(x = vals[2,],y = vals[1,]), color = 'blue')
#correct set up to interpolate p values wrt to CDF samples


#=============================
# NE_W_sampler
#=============================
source('NE_W_sampler.R')

set.seed(1337)
x = NE_W_sampler(100,alpha = 1,gam = 2,beta = 3)
library(ggplot2)
ggplot()+
  geom_line(aes(x = x[,1], y = x[,2]), color = 'red',size = 1.2)+
  geom_line(aes(x = x[,1], y = NE_W_CDF(x[,1],alpha = 1,gam = 2,beta = 3) ))



#=============================
# NE_WT_CDF
#=============================
source('NE_WT_CDF.R')
ggplot()+
  geom_line(aes(x = (-100:100)/10,
                y = NE_WT_CDF((-100:100)/10,0.2, 2,2,4, 5,1,3  )
  )
  )


#=============================
# NE_WT_PDF
#=============================
source('NE_WT_PDF.R')
ggplot()+
  geom_line(aes(x = (-100:100)/10,
                y = NE_WT_PDF((-100:100)/10,0.4, 2,3,1, 3,2,1/3  )
  )
  )


#=============================
# NE_WT_sampler
#=============================
source('NE_WT_sampler.R')
source('NE_WT_PDF.R')
set.seed(1337)

x = NE_WT_sampler(3000000,theta = 0.4,
                  alpha1 = 1,gam1 = 2,beta1 = 3,
                  alpha2 = 2,gam2 = 2.5,beta2 = 1.5)
library(ggplot2)
ggplot()+
  geom_density(aes(x = x[,1]), color = 'red', size = 1)+
  geom_line(aes(x = x[,1], y = NE_WT_PDF(x[,1],
                                         theta = 0.4,
                                         alpha1 = 1,gam1 = 2,beta1 = 3,
                                         alpha2 = 2,gam2 = 2.5,beta2 = 1.5)
                )
            )+
  scale_y_continuous(trans = 'log10')+
  xlim(-2.5,5)

#=============================
# NE_W_loglik
#=============================
source('NE_W_loglik.R')
source('NE_W_sampler.R')

#sample from NE_W(1,2,3)
set.seed(1337)
x1 = NE_W_sampler(1000000,alpha = 1,gam = 2,beta = 3)


# Verify samples are good
library(ggplot2)
ggplot()+
  geom_line(aes(x = x1[,1], y = x1[,2]), color = 'red',size = 1.2)+
  geom_line(aes(x = x1[,1], y = NE_W_CDF(x1[,1],alpha = 1,gam = 2,beta = 3) ))

ggplot()+
  geom_density(aes(x = x1[,1]), colour = 'red',size = 1)+
  geom_line(aes(x = x1[,1], y = NE_W_PDF(x1[,1],alpha = 1, gam = 2, beta = 3)))

# Test NE_W_MLE
#Manually get results
Results = optim(par = c(2,2,2), fn = NE_W_logLik,
                gr = NE_W_partial_loglik, x = x1[,1],
                control = list(fnscale = -1)
                ) 
#Verify Manual results are good
ggplot()+
  geom_density(aes(x = x1[,1]), colour = 'red',size = 1)+
  geom_line(aes(x = x1[,1], y = NE_W_PDF(x1[,1],alpha = Results$par[1], gam = Results$par[2], beta = Results$par[3])))

Result_MLE = NE_W_MLE(data = x1[,1])
ggplot()+
  geom_density(aes(x = x[,1]), colour = 'red',size = 1.3)+
  geom_line(aes(x = x[,1], y = NE_W_PDF(x[,1],alpha = Result_MLE$par[1], gam = Result_MLE$par[2], beta = Result_MLE$par[3])),
            colour = 'green', size = 1)+
  geom_line(aes(x = x[,1], y = NE_W_PDF(x[,1],alpha = Results$par[1], gam = Results$par[2], beta = Results$par[3])))

#=============================
# NE_W_moments
#=============================
source('NE_W_moments.R')
a1 = integrate(NE_W_PDF, 0, 1.96,alpha = 1, gam = 2, beta = 3)
a2 = moment_func(0,NE_W_PDF, c(0,1.96), alpha = 1, gam = 2, beta = 3)
a1[[1]]== a2[[1]]

#mean
x = moment_func(1,NE_W_PDF, c(0,Inf), alpha = 1, gam = 2, beta = 3)
x1 = NE_W_moment(1,c(0,Inf), alpha = 1, gam = 2, beta = 3)
x1x[[1]] == x1[[1]]

#2nd moment
x = moment_func(2,NE_W_PDF, c(0,Inf), alpha = 1, gam = 2, beta = 3)
x1 = NE_W_moment(2,c(0,Inf), alpha = 1, gam = 2, beta = 3)
x[[1]]==x1[[1]]
x


#=============================
# NE_WT_moments
#=============================
source('NE_WT_PDF.R')
source('NE_WT_moments.R')
#Verify from CDF values
x1 = moment_func(0, pdf = NE_WT_PDF, theta = 0.2,
                 alpha1 = 2,gam1 = 2,beta1 = 1.5,
                 alpha2 = 1,gam2 = 1.8,beta2 = 1.3)[[1]]
x1
# Indeed gets theta, verify by changing theta.
# Theta is the prob of being negative, ie domain =  c(-Inf,0)



#verify E[X^r] = (1-theta)*E[X_1^r] + theta*E[X_2^r]
#r = 1:
x1 = moment_func(1, pdf = NE_WT_PDF, theta = 0.2,
                alpha1 = 2,gam1 = 2,beta1 = 1.5,
                alpha2 = 1,gam2 = 1.8,beta2 = 1.3)[[1]]
x1

theta = 0.2
x2 = (1-theta)*NE_W_moment(1,c(0,Inf), alpha = 2, gam = 2, beta = 1.5)[[1]]+
  (-theta)*NE_W_moment(1,c(0,Inf), alpha = 1, gam = 1.8, beta = 1.3)[[1]]
x2

x3 = NE_WT_moment(1,theta = 0.2,
                  alpha1 = 2,gam1 = 2,beta1 = 1.5,
                  alpha2 = 1,gam2 = 1.8,beta2 = 1.3)[[1]]
x3

#r = 2:
x1 = moment_func(2, pdf = NE_WT_PDF, theta = 0.2,
                 alpha1 = 2,gam1 = 2,beta1 = 1.5,
                 alpha2 = 1,gam2 = 1.8,beta2 = 1.3)[[1]]
x1

theta = 0.2
x2 = (1-theta)*NE_W_moment(2,c(0,Inf), alpha = 2, gam = 2, beta = 1.5)[[1]]+
  (theta)*(-1)^2*NE_W_moment(2,c(0,Inf), alpha = 1, gam = 1.8, beta = 1.3)[[1]]
x2


abs(x2 - x1) < 1e-10

x3 = NE_WT_moment(2,theta = 0.2,
                  alpha1 = 2,gam1 = 2,beta1 = 1.5,
                  alpha2 = 1,gam2 = 1.8,beta2 = 1.3)[[1]]
abs(x3 - x1) < 1e-10



#=============================
# NE_STD_sampler, _CDF, _PDF
#=============================
#NE_STD is the standardised NE_WT
# ie theta = E[X_1]/(E[X_1] - E[X_2]), X_2 is negative side
source('NE_STD_sampler.R')
source('NE_STD_CDF.R')
source('NE_STD_PDF.R')
set.seed(1337)
x = NE_STD_sampler(1000000,
                   alpha1 = 1,gam1 = 2,beta1 = 3,
                   alpha2 = 1.3,gam2 = 1.8,beta2 = 2)

#TESTING
# xvar = NE_STD_var1_sampler(1000000,
#                    alpha1 = 1,gam1 = 2,beta1 = 3,
#                    alpha2 = 1.3,gam2 = 1.8,beta2 = 2)
# Verify samples are good
library(ggplot2)
ggplot()+
  geom_line(aes(x = x[,1], y = x[,2]), color = 'red',size = 1.2)+
  geom_line(aes(x = x[,1], 
                y = NE_STD_CDF(x[,1],
                              alpha1 = 1,gam1 = 2,beta1 = 3,
                              alpha2 = 1.3,gam2 = 1.8,beta2 = 2) ))

ggplot()+
  geom_density(aes(x = x[,1]), colour = 'red',size = 1)+
  geom_line(aes(x = x[,1], 
                y = NE_STD_PDF(x[,1],
                              alpha1 = 1,gam1 = 2,beta1 = 3,
                              alpha2 = 1.3,gam2 = 1.8,beta2 = 2) ))





#=============================
# NE_STD_loglik
#=============================
source('NE_STD_loglik.R')

Results = optim(par = c(2,2,2,2,2,2), fn = NE_STD_logLik,
                gr = NE_STD_partial_loglik, x = x[,1],
                control = list(fnscale = -1)
) 
Results


ggplot()+
  geom_density(aes(x = x[,1]), colour = 'red',size = 1)+
  geom_line(aes(x = x[,1], 
                y = NE_STD_PDF(x[,1],
                               alpha1 = Results$par[1], gam1 = Results$par[2],
                               beta1 = Results$par[3], alpha2 = Results$par[4],
                               gam2 = Results$par[5], beta2 = Results$par[6])))

#Verify the mean is zero
moment_func(r= 1, pdf = NE_STD_PDF, domain = c(-Inf,Inf),
            alpha1 = Results$par[1], gam1 = Results$par[2],
            beta1 = Results$par[3], alpha2 = Results$par[4],
            gam2 = Results$par[5], beta2 = Results$par[6])

Results_MLE = NE_STD_MLE(x = x[,1])



