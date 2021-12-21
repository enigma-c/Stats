#load data
setwd("/home/leon/Documents/GitHub/stat/Stats/Data")
data<-read.csv('Social_distancing_model_data.csv')


#get modules
library(bayestestR)
library(rstanarm)
library(ggplot2)
library(tidyverse)
library(mombf)

X= model.matrix(~ X1+X2+Covid_topic_prop+covid_cases+strictness_measure, data=data)
y= data$observed
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
  coord_cartesian(xlim=c(0.01,0.4),ylim=c(0.01,0.4)) +
  theme_classic()
vars<-c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)
fit.bayesreg <- modelSelection(y=y,x=X,includevars=vars, priorCoef=zellnerprior(taustd=1), priorDelta=modelbbprior(1,1))

prob<-postProb(fit.bayesreg)
head(postProb(fit.bayesreg),10)

ci.bayesreg <- coef(fit.bayesreg)[1:6,]
sel.bayesreg <- ci.bayesreg[,4] > 0.5
ci.bayesreg[,1:3]= round(ci.bayesreg[,1:3], 3)  
ci.bayesreg[,4]= round(ci.bayesreg[,4], 4)      
head(ci.bayesreg)

plot(NA, ylim=1.25*range(ci.bayesreg[,1:3]), xlim=c(0,nrow(ci.bayesreg)), ylab='95% CI', xlab='', main='Bayesian Model Selection')
cols= ifelse(beta < ci.bayesreg[ , 1] | beta > ci.bayesreg[, 2], 2, 1)
segments(y0 = ci.bayesreg[, 2], y1 = ci.bayesreg[, 3], x0 = 1:nrow(ci.bayesreg), col = cols)
points(1:5, ci.bayesreg[,2], pch = 16)
library(xtable)
xtable(ci.bayesreg)
xtable(prob[1:5,],digits=6)