#load data
setwd("/home/leon/Documents/GitHub/stat/Stats/Data")
data<-read.csv('stackedgroceries.csv')


#get modules
library(bayestestR)
library(rstanarm)
library(ggplot2)
library(tidyverse)
library(mombf)

X= model.matrix(~ X2+X3+X4, data=data)
y= data$values
n= nrow(X)

gseq= exp(seq(log(.01),log(1),length=20))
gseq

#prior elicitation
library(mvtnorm)
V= diag(ncol(X))
beta= rmvnorm(1000, sigma= V)
sse= colSums((X %*% t(beta))^2) / n

r2= double(length(gseq))
for (i in 1:length(gseq)) {
  r2[i]= mean(sse * gseq[i] / (1 + sse * gseq[i]))
}

par(mar=c(4,5,.1,.1), cex.lab=1.3, cex.axis=1.3)
plot(gseq, r2, type='o', xlab='g', ylab=expression(paste('Theoretical ',R^2)))

#comparison with MLE
fit.mle= lm(y ~ X[,-1]) #1st column in x is the intercept, already added by lm
b.mle= coef(fit.mle)
summary(fit.mle)

fit.bayes <- stan_glm(y ~ X[,-1], family = gaussian(link = "identity"), algorithm='sampling', refresh=0)
b.bayes= coef(fit.bayes)

data.frame(mle= b.mle[-1], bayes= b.bayes[-1]) %>% 
  ggplot(aes(x=mle,y=bayes)) + 
  geom_point(shape = "O",size=2) +
  geom_abline(slope=1, intercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  xlab('MLE OLS') +
  ylab('MCMC Bayesian regression') +
  coord_cartesian(xlim=c(1,10),ylim=c(1,10)) +
  theme_classic()

fit.bayesreg <- modelSelection(y=y,x=X, priorCoef=zellnerprior(taustd=1), priorDelta=modelbbprior(1,1))


head(postProb(fit.bayesreg),10)

ci.bayesreg <- coef(fit.bayesreg)[-c(1,nrow(coef(fit.bayesreg))),]
sel.bayesreg <- ci.bayesreg[,4] > 0.5
ci.bayesreg[,1:3]= round(ci.bayesreg[,1:3], 3)  
ci.bayesreg[,4]= round(ci.bayesreg[,4], 4)      
head(ci.bayesreg)

plot(NA, ylim=1.25*range(ci.bayesreg[,1:3]), xlim=c(0,nrow(ci.bayesreg)), ylab='95% CI', xlab='', main='Bayesian Model Selection')
cols= ifelse(beta < ci.bayesreg[ , 1] | beta > ci.bayesreg[, 2], 2, 1)
segments(y0 = ci.bayesreg[, 2], y1 = ci.bayesreg[, 3], x0 = 1:nrow(ci.bayesreg), col = cols)
points(1:4, beta, pch = 16)