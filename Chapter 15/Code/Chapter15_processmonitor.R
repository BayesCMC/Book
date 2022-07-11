dat = read.csv("dat.csv", header=TRUE)

library(lme4)
datCurrent=dat[dat$Group==1,] #data for current batches
datCurrent$Time = as.factor(datCurrent$Time)
datCurrent$Batch = as.factor(datCurrent$Batch)
fitReml = lmer(value ~ Time + (1|Batch) + (1|Batch:Time), data=datCurrent)
summary(fitReml) 

library(lmerTest)
anova(fitReml)

## MCMC in R via rstan

library(rstan)

X <- model.matrix(~ as.factor(Time), data = dat) ## design matrix for fixed effects
Z1 <- cbind(X[,1]-apply(X[,2:7], 1, sum), X[,2:7]) ## design matrix for random effects gamma

# rstan program
# power prior a=0.5 with Jeffreys' initial prior

stanDat <- list(Batch=as.integer(dat$Batch), Group=as.integer(dat$Group), y=dat$value, N=nrow(dat), I=10, T=7, K=10, X=X, Z1=Z1, a=0.5) ## 10 batches (5 current, 5 historical), 7 time points, and 10 replicates for each time point and each batch

scode <- "

data {

int<lower=1> N;
int<lower=1> I;
int<lower=1> T;
int<lower=1> K;
vector[N] y;
matrix[N,T] X;
matrix[N,T] Z1;
int<lower=1,upper=2> Group[N];
int<lower=1,upper=I> Batch[N];
real<lower=0,upper=1> a;

}

parameters {

vector[T] eta;/* fixed effects */
vector[I] beta; /* random effect added by batch */
vector[I] gamma1; /* random effect added by batch*(time point 1) */
vector[I] gamma2; /* random effect added by batch*(time point 2) */
vector[I] gamma3; /* random effect added by batch*(time point 3) */
vector[I] gamma4; /* random effect added by batch*(time point 4) */
vector[I] gamma5; /* random effect added by batch*(time point 5) */
vector[I] gamma6; /* random effect added by batch*(time point 6) */
vector[I] gamma7; /* random effect added by batch*(time point 7) */

real<lower=0> sigma_e2; /* residual variance */
real<lower=0> sigma_beta2; /* variance for random effect added by batch */
real<lower=0> sigma_gamma2; /* variance for random effect added by batch*time */

}

transformed parameters {

real<lower=0> sigma_e;
real<lower=0> sigma_beta;
real<lower=0> sigma_gamma;

sigma_e = sqrt(sigma_e2);
sigma_beta = sqrt(sigma_beta2);
sigma_gamma = sqrt(sigma_gamma2);

}

model {

real mu;   /* Conditional Expectation of Y given beta and gamma */

beta ~ normal(0, sigma_beta);
gamma1 ~ normal(0, sigma_gamma);
gamma2 ~ normal(0, sigma_gamma);
gamma3 ~ normal(0, sigma_gamma);
gamma4 ~ normal(0, sigma_gamma);
gamma5 ~ normal(0, sigma_gamma);
gamma6 ~ normal(0, sigma_gamma);
gamma7 ~ normal(0, sigma_gamma);

for (i in 1:N){

mu = X[i,]*eta + beta[Batch[i]] + gamma1[Batch[i]]*Z1[i,1]+gamma2[Batch[i]]*Z1[i,2]+gamma3[Batch[i]]*Z1[i,3] +gamma4[Batch[i]]*Z1[i,4]+gamma5[Batch[i]]*Z1[i,5]+gamma6[Batch[i]]*Z1[i,6]+ gamma7[Batch[i]]*Z1[i,7];

if (Group[i] == 2) 
target += a*normal_lpdf(y[i]|mu, sigma_e);

if (Group[i] == 1)
target += normal_lpdf(y[i]|mu, sigma_e);

}

target += -log(sigma_e2)-log(K*sigma_gamma2+sigma_e2)-log(T*K*sigma_beta2+K*sigma_gamma2+sigma_e2); /* add Jeffreys' initial prior for variance components; assume constant initial prior for fixed effects */  

}
"

fitppj <- stan(model_code=scode, data=stanDat, iter=10000) 

# Jeffreys' prior (a=0)

stanDat <- list(Batch=as.integer(dat$Batch), Group=as.integer(dat$Group), y=dat$value, N=nrow(dat), I=10, T=7, K=10, X=X, Z1=Z1, a=0) ## 10 batches (5 current, 5 historical), 7 time points, and 10 replicates for each time point and each batch

scode <- "

data {

int<lower=1> N;
int<lower=1> I;
int<lower=1> T;
int<lower=1> K;
vector[N] y;
matrix[N,T] X;
matrix[N,T] Z1;
int<lower=1,upper=2> Group[N];
int<lower=1,upper=I> Batch[N];
real<lower=0,upper=1> a;

}

parameters {

vector[T] eta;/* fixed effects */
vector[I] beta; /* random effect added by batch */
vector[I] gamma1; /* random effect added by batch*(time point 1) */
vector[I] gamma2; /* random effect added by batch*(time point 2) */
vector[I] gamma3; /* random effect added by batch*(time point 3) */
vector[I] gamma4; /* random effect added by batch*(time point 4) */
vector[I] gamma5; /* random effect added by batch*(time point 5) */
vector[I] gamma6; /* random effect added by batch*(time point 6) */
vector[I] gamma7; /* random effect added by batch*(time point 7) */

real<lower=0> sigma_e2; /* residual variance */
real<lower=0> sigma_beta2; /* variance for random effect added by batch */
real<lower=0> sigma_gamma2; /* variance for random effect added by batch*time */

}

transformed parameters {

real<lower=0> sigma_e;
real<lower=0> sigma_beta;
real<lower=0> sigma_gamma;

sigma_e = sqrt(sigma_e2);
sigma_beta = sqrt(sigma_beta2);
sigma_gamma = sqrt(sigma_gamma2);

}

model {

real mu;   /* Conditional Expectation of Y given beta and gamma */

beta ~ normal(0, sigma_beta);
gamma1 ~ normal(0, sigma_gamma);
gamma2 ~ normal(0, sigma_gamma);
gamma3 ~ normal(0, sigma_gamma);
gamma4 ~ normal(0, sigma_gamma);
gamma5 ~ normal(0, sigma_gamma);
gamma6 ~ normal(0, sigma_gamma);
gamma7 ~ normal(0, sigma_gamma);

for (i in 1:N){

mu = X[i,]*eta + beta[Batch[i]] + gamma1[Batch[i]]*Z1[i,1]+gamma2[Batch[i]]*Z1[i,2]+gamma3[Batch[i]]*Z1[i,3] +gamma4[Batch[i]]*Z1[i,4]+gamma5[Batch[i]]*Z1[i,5]+gamma6[Batch[i]]*Z1[i,6]+ gamma7[Batch[i]]*Z1[i,7];

if (Group[i] == 2) 
target += a*normal_lpdf(y[i]|mu, sigma_e);

if (Group[i] == 1)
target += normal_lpdf(y[i]|mu, sigma_e);

}

target += -log(sigma_e2)-log(K*sigma_gamma2+sigma_e2)-log(T*K*sigma_beta2+K*sigma_gamma2+sigma_e2); /* add Jeffreys' initial prior for variance components; assume constant initial prior for fixed effects */  

}
"

fitnpj <- stan(model_code=scode, data=stanDat, iter=10000) 

# Power prior a=0.5 with constant initial prior

stanDat <- list(Batch=as.integer(dat$Batch), Group=as.integer(dat$Group), y=dat$value, N=nrow(dat), I=10, T=7, K=10, X=X, Z1=Z1, a=0.5) ## 10 batches (5 current, 5 historical), 7 time points, and 10 replicates for each time point and each batch

scode <- "

data {

int<lower=1> N;
int<lower=1> I;
int<lower=1> T;
int<lower=1> K;
vector[N] y;
matrix[N,T] X;
matrix[N,T] Z1;
int<lower=1,upper=2> Group[N];
int<lower=1,upper=I> Batch[N];
real<lower=0,upper=1> a;

}

parameters {

vector[T] eta;/* fixed effects */
vector[I] beta; /* random effect added by batch */
vector[I] gamma1; /* random effect added by batch*(time point 1) */
vector[I] gamma2; /* random effect added by batch*(time point 2) */
vector[I] gamma3; /* random effect added by batch*(time point 3) */
vector[I] gamma4; /* random effect added by batch*(time point 4) */
vector[I] gamma5; /* random effect added by batch*(time point 5) */
vector[I] gamma6; /* random effect added by batch*(time point 6) */
vector[I] gamma7; /* random effect added by batch*(time point 7) */

real<lower=0> sigma_e; /* residual standard deviation */
real<lower=0> sigma_beta; /* standard deviation for random effect added by batch */
real<lower=0> sigma_gamma; /* standard deviation for random effect added by batch*time */

}

model {

real mu;   /* Conditional Expectation of Y given beta and gamma */

beta ~ normal(0, sigma_beta);
gamma1 ~ normal(0, sigma_gamma);
gamma2 ~ normal(0, sigma_gamma);
gamma3 ~ normal(0, sigma_gamma);
gamma4 ~ normal(0, sigma_gamma);
gamma5 ~ normal(0, sigma_gamma);
gamma6 ~ normal(0, sigma_gamma);
gamma7 ~ normal(0, sigma_gamma);

for (i in 1:N){

mu = X[i,]*eta + beta[Batch[i]] + gamma1[Batch[i]]*Z1[i,1]+gamma2[Batch[i]]*Z1[i,2]+gamma3[Batch[i]]*Z1[i,3] +gamma4[Batch[i]]*Z1[i,4]+gamma5[Batch[i]]*Z1[i,5]+gamma6[Batch[i]]*Z1[i,6]+ gamma7[Batch[i]]*Z1[i,7];

if (Group[i] == 2) 
target += a*normal_lpdf(y[i]|mu, sigma_e);

if (Group[i] == 1)
target += normal_lpdf(y[i]|mu, sigma_e);

}

}
"

fitpp <- stan(model_code=scode, data=stanDat, iter=10000) 

# constant prior (a=0)

stanDat <- list(Batch=as.integer(dat$Batch), Group=as.integer(dat$Group), y=dat$value, N=nrow(dat), I=10, T=7, K=10, X=X, Z1=Z1, a=0) ## 10 batches (5 current, 5 historical), 7 time points, and 10 replicates for each time point and each batch

scode <- "

data {

int<lower=1> N;
int<lower=1> I;
int<lower=1> T;
int<lower=1> K;
vector[N] y;
matrix[N,T] X;
matrix[N,T] Z1;
int<lower=1,upper=2> Group[N];
int<lower=1,upper=I> Batch[N];
real<lower=0,upper=1> a;

}

parameters {

vector[T] eta;/* fixed effects */
vector[I] beta; /* random effect added by batch */
vector[I] gamma1; /* random effect added by batch*(time point 1) */
vector[I] gamma2; /* random effect added by batch*(time point 2) */
vector[I] gamma3; /* random effect added by batch*(time point 3) */
vector[I] gamma4; /* random effect added by batch*(time point 4) */
vector[I] gamma5; /* random effect added by batch*(time point 5) */
vector[I] gamma6; /* random effect added by batch*(time point 6) */
vector[I] gamma7; /* random effect added by batch*(time point 7) */

real<lower=0> sigma_e; /* residual standard deviation */
real<lower=0> sigma_beta; /* standard deviation for random effect added by batch */
real<lower=0> sigma_gamma; /* standard deviation for random effect added by batch*time */

}

model {

real mu;   /* Conditional Expectation of Y given beta and gamma */

beta ~ normal(0, sigma_beta);
gamma1 ~ normal(0, sigma_gamma);
gamma2 ~ normal(0, sigma_gamma);
gamma3 ~ normal(0, sigma_gamma);
gamma4 ~ normal(0, sigma_gamma);
gamma5 ~ normal(0, sigma_gamma);
gamma6 ~ normal(0, sigma_gamma);
gamma7 ~ normal(0, sigma_gamma);

for (i in 1:N){

mu = X[i,]*eta + beta[Batch[i]] + gamma1[Batch[i]]*Z1[i,1]+gamma2[Batch[i]]*Z1[i,2]+gamma3[Batch[i]]*Z1[i,3] +gamma4[Batch[i]]*Z1[i,4]+gamma5[Batch[i]]*Z1[i,5]+gamma6[Batch[i]]*Z1[i,6]+ gamma7[Batch[i]]*Z1[i,7];

if (Group[i] == 2) 
target += a*normal_lpdf(y[i]|mu, sigma_e);

if (Group[i] == 1)
target += normal_lpdf(y[i]|mu, sigma_e);

}

}
"

fitnp <- stan(model_code=scode, data=stanDat, iter=10000) 

# MCMC posterior output
print(fitppj, pars = c("eta", "beta", "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "gamma6", "gamma7", "sigma_e", "sigma_beta", "sigma_gamma"), probs = c(0.025, 0.5, 0.975), digits = 3)
print(fitnpj, pars = c("eta", "beta", "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "gamma6", "gamma7", "sigma_e", "sigma_beta", "sigma_gamma"), probs = c(0.025, 0.5, 0.975), digits = 3)
print(fitpp, pars = c("eta", "beta", "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "gamma6", "gamma7", "sigma_e", "sigma_beta", "sigma_gamma"), probs = c(0.025, 0.5, 0.975), digits = 3)
print(fitnp, pars = c("eta", "beta", "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", "gamma6", "gamma7", "sigma_e", "sigma_beta", "sigma_gamma"), probs = c(0.025, 0.5, 0.975), digits = 3)

############
# Figure 1 #
############

# Trace plot
out = as.array(fitppj)
plot.ts(out[,,"eta[1]"], main="alpha", xlab="Iteractions", ylab="", plot.type="single", col=grey.colors(8))
legend("topleft", legend=c("Chain 1", "Chain 2", "Chain 3", "Chain 4"), col=grey.colors(8), lty=1)

# acf plot
acf(c(extract(fitppj, pars="eta[1]", permuted=FALSE)), main="alpha")

############
# Figure 2 #
############

library(ggpubr)

out = as.matrix(fitppj) # power prior with Jeffreys’ initial prior
chain = rep(1:4, each=5000) 
outplot = data.frame(cbind(out[,c(1,91,92,93)], chain))
outplot$alpha=outplot$eta.1.
outplot$chain=as.factor(outplot$chain)

out1 = as.matrix(fitnpj) # Noninformative prior – Jeffreys’ prior 
outplot1 = data.frame(cbind(out1[,c(1,91,92,93)], chain))
outplot1$alpha=outplot1$eta.1.
outplot1$chain=as.factor(outplot1$chain)

a=ggplot(outplot, aes(x=alpha, lty=chain))+geom_density()+xlim(c(445,455))+theme(axis.title=element_text(size=16))
b=ggplot(outplot, aes(x=sigma_beta, lty=chain))+geom_density()+xlim(c(0,10))+theme(axis.title=element_text(size=16))
c=ggplot(outplot, aes(x=sigma_gamma, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))
d=ggplot(outplot, aes(x=sigma_e, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))

e=ggplot(outplot1, aes(x=alpha, lty=chain))+geom_density()+xlim(c(445,455))+theme(axis.title=element_text(size=16))
f=ggplot(outplot1, aes(x=sigma_beta, lty=chain))+geom_density()+xlim(c(0,10))+theme(axis.title=element_text(size=16))
g=ggplot(outplot1, aes(x=sigma_gamma, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))
h=ggplot(outplot1, aes(x=sigma_e, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))

ggarrange(a,b,c,d,e,f,g,h, labels=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"),ncol=4,nrow=2)


############
# Figure 3 #
############

out = as.matrix(fitppj) # power prior with Jeffreys’ initial prior
chain = rep(1:4, each=5000) 
outplot = data.frame(cbind(out[,c(1,91,92,93)], chain))
outplot$alpha=outplot$eta.1.
outplot$chain=as.factor(outplot$chain)

out1 = as.matrix(fitpp) # power prior with constant initial prior 
outplot1 = data.frame(cbind(out1[,c(1,88,89,90)], chain))
outplot1$alpha=outplot1$eta.1.
outplot1$chain=as.factor(outplot1$chain)

a=ggplot(outplot, aes(x=alpha, lty=chain))+geom_density()+xlim(c(445,455))+theme(axis.title=element_text(size=16))
b=ggplot(outplot, aes(x=sigma_beta, lty=chain))+geom_density()+xlim(c(0,10))+theme(axis.title=element_text(size=16))
c=ggplot(outplot, aes(x=sigma_gamma, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))
d=ggplot(outplot, aes(x=sigma_e, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))

e=ggplot(outplot1, aes(x=alpha, lty=chain))+geom_density()+xlim(c(445,455))+theme(axis.title=element_text(size=16))
f=ggplot(outplot1, aes(x=sigma_beta, lty=chain))+geom_density()+xlim(c(0,10))+theme(axis.title=element_text(size=16))
g=ggplot(outplot1, aes(x=sigma_gamma, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))
h=ggplot(outplot1, aes(x=sigma_e, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))

ggarrange(a,b,c,d,e,f,g,h, labels=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"),ncol=4,nrow=2)

############
# Figure 4 #
############

out = as.matrix(fitnpj) # Jeffreys’ prior
chain = rep(1:4, each=5000) 
outplot = data.frame(cbind(out[,c(1,91,92,93)], chain))
outplot$alpha=outplot$eta.1.
outplot$chain=as.factor(outplot$chain)

out1 = as.matrix(fitnp) # Constant prior 
outplot1 = data.frame(cbind(out1[,c(1,88,89,90)], chain))
outplot1$alpha=outplot1$eta.1.
outplot1$chain=as.factor(outplot1$chain)

a=ggplot(outplot, aes(x=alpha, lty=chain))+geom_density()+xlim(c(445,455))+theme(axis.title=element_text(size=16))
b=ggplot(outplot, aes(x=sigma_beta, lty=chain))+geom_density()+xlim(c(0,10))+theme(axis.title=element_text(size=16))
c=ggplot(outplot, aes(x=sigma_gamma, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))
d=ggplot(outplot, aes(x=sigma_e, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))

e=ggplot(outplot1, aes(x=alpha, lty=chain))+geom_density()+xlim(c(445,455))+theme(axis.title=element_text(size=16))
f=ggplot(outplot1, aes(x=sigma_beta, lty=chain))+geom_density()+xlim(c(0,10))+theme(axis.title=element_text(size=16))
g=ggplot(outplot1, aes(x=sigma_gamma, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))
h=ggplot(outplot1, aes(x=sigma_e, lty=chain))+geom_density()+xlim(c(0,5))+theme(axis.title=element_text(size=16))

ggarrange(a,b,c,d,e,f,g,h, labels=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"),ncol=4,nrow=2)


####################################################################################
# Figure 5 and control limits for power prior (a=0.5) with Jeffreys' initial prior #
####################################################################################

library(ggplot2)

out = as.matrix(fitppj)

y<-list()

for (i in 1:20000){
 
  y_int <- rep(rnorm(5*7,sd= out[i,93]),each=10)
  y_batch <- rnorm(5,sd=out[i,92])
  y1<-as.vector(rnorm(70,mean=X[1:70,]%*%out[i,1:7],sd=out[i,91])) +y_batch[1]+y_int[1:70]
  y2<-as.vector(rnorm(70,mean=X[1:70,]%*%out[i,1:7],sd=out[i,91])) +y_batch[2]+y_int[71:140]
  y3<-as.vector(rnorm(70,mean=X[1:70,]%*%out[i,1:7],sd=out[i,91])) +y_batch[3]+y_int[141:210]
  y4<-as.vector(rnorm(70,mean=X[1:70,]%*%out[i,1:7],sd=out[i,91])) +y_batch[4]+y_int[211:280]
  y5<-as.vector(rnorm(70,mean=X[1:70,]%*%out[i,1:7],sd=out[i,91])) +y_batch[5]+y_int[281:350]  
  y[[i]]<-as.data.frame(c(y1, y2, y3, y4, y5))
  
  y[[i]]$Batch<-rep(1:5, each=70)
  time<-rep(1:7, each=10)
  y[[i]]$Time<-rep(time, 5)
  
  colnames(y[[i]])<-c("value", "Batch", "Time")
  
}

q.min<-unlist(lapply(y, function(x) min(x[, "value"])))
q.max<-unlist(lapply(y, function(x) max(x[, "value"])))

qvec = data.frame(cbind(q.min, q.max))
ggplot(qvec, aes(x=q.max, y=q.min)) + 
geom_point(color="grey") + 
stat_ellipse(level=0.95, color=1, linetype=1, lwd=1.2) +
theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), panel.border=element_rect(color="black", fill=NA, size=1))

lwr<-quantile(q.min, prob=1-pnorm(3))
uppr<-quantile(q.max, prob=pnorm(3))

##################
# Control limits #
##################

## Method (1)

mean(datCurrent$value) - 3*sd(datCurrent$value)
mean(datCurrent$value) + 3*sd(datCurrent$value)

## Method (2)

mean(datCurrent$value)*0.95
mean(datCurrent$value)*1.05

## Method (3)

varvector = as.data.frame(VarCorr(fitReml))
sigma_beta = varvector$sdcor[2]
sigma_gamma = varvector$sdcor[1]
sigma_e = varvector$sdcor[3]
eta = fitReml@beta

y<-list()

for (i in 1:20000){
 
  y_int <- rep(rnorm(5*7,sd=sigma_gamma),each=10)
  y_batch <- rnorm(5,sd=sigma_beta)
  y1<-as.vector(rnorm(70,mean=X[1:70,]%*%eta, sd=sigma_e))+y_batch[1] +y_int[1:70]
  y2<-as.vector(rnorm(70,mean=X[1:70,]%*%eta, sd=sigma_e))+y_batch[2] +y_int[71:140]
  y3<-as.vector(rnorm(70,mean=X[1:70,]%*%eta, sd=sigma_e))+y_batch[3] +y_int[141:210]
  y4<-as.vector(rnorm(70,mean=X[1:70,]%*%eta, sd=sigma_e))+y_batch[4] +y_int[211:280]
  y5<-as.vector(rnorm(70,mean=X[1:70,]%*%eta, sd=sigma_e))+y_batch[5] +y_int[281:350]  
  y[[i]]<-as.data.frame(c(y1, y2, y3, y4, y5))
  
  y[[i]]$Batch<-rep(1:5, each=70)
  time<-rep(1:7, each=10)
  y[[i]]$Time<-rep(time, 5)
  
  colnames(y[[i]])<-c("value", "Batch", "Time")
  
}

q.min<-unlist(lapply(y, function(x) min(x[, "value"])))
q.max<-unlist(lapply(y, function(x) max(x[, "value"])))

lwr<-quantile(q.min, prob=1-pnorm(3))
uppr<-quantile(q.max, prob=pnorm(3))

############
# Figure 6 #
############

library(reshape); library(lattice); library(gridExtra); library(grid);
CL<-read.csv(file="Control Limits.csv", header=TRUE) # Limits developed were saved in the file "Control Limits.csv". There are 3 columns in the file: Approach (5 approaches), LCL (lower control limits), and UCL (upper control limits)
ORDER<-c("3*Sigma", "REML", "Power prior", "Jeffrey prior", "5%*Mean")
LABEL<-rev(c(paste("3", "\u03C3", sep=""), "Monte Carlo with\n REML estimates", "Power Prior with\n Jeffreys' Initial Prior", "Jeffreys' Prior", paste("\u00B1", "5% of Mean", sep="")))
CL.Ordered <- CL[match(ORDER, CL$Approach), ] 
CL.Ordered_m<-melt(CL.Ordered, id.vars=c("Approach"), measure.vars=c("LCL", "UCL"))
trellis.par.set(box.umbrella = list(lty = 1, lwd=3, col="black"), box.rectangle=list(lwd=3, col="black"))

mypanel<-function(x,y,...){
     panel.bwplot(x,y,coef=0,box.ratio=0, col="black", pch="")
     panel.text(CL.Ordered$LCL[1], x=CL.Ordered$LCL[1], y=5.2, cex=1)
     panel.text(CL.Ordered$UCL[1], x=CL.Ordered$UCL[1], y=5.2, cex=1)
     panel.text(CL.Ordered$LCL[2], x=CL.Ordered$LCL[2], y=4.2, cex=1)
     panel.text(CL.Ordered$UCL[2], x=CL.Ordered$UCL[2], y=4.2, cex=1)
     panel.text(CL.Ordered$LCL[3], x=CL.Ordered$LCL[3], y=3.2, cex=1)
     panel.text(CL.Ordered$UCL[3], x=CL.Ordered$UCL[3], y=3.2, cex=1)
     panel.text(CL.Ordered$LCL[4], x=CL.Ordered$LCL[4], y=2.2, cex=1)
     panel.text(CL.Ordered$UCL[4], x=CL.Ordered$UCL[4], y=2.2, cex=1)
     panel.text(CL.Ordered$LCL[5], x=CL.Ordered$LCL[5], y=1.2, cex=1)
     panel.text(CL.Ordered$UCL[5], x=CL.Ordered$UCL[5], y=1.2, cex=1)
}
bwplot(factor(Approach, levels=rev(ORDER), labels=c(LABEL)) ~ value , CL.Ordered_m, ylab="Approach", xlab="Control Limits", col="black", pch="", panel=mypanel)

############
# Figure 7 #
############

library(lattice); library(latticeExtra)
Data=dat
Data$BatchDesc<-ifelse(Data$Batch %in% 1:5, "Current", "Historical")
Data$BatchNo<-ifelse(Data$Batch %in% 6:10, Data$Batch-5, Data$Batch )
LABEL1<-c(paste("\u00B1", " 3", "\u03C3", sep=""), "Monte Carlo with REML estimates", "Power Prior with Jeffreys' Initial Prior", "Jeffreys' Prior", paste("\u00B1", " 5% of Mean", sep=""))

mypanel<-function(x, y, ...){
     panel.xyplot(x, y, ...)
     panel.abline(h=c(CL.Ordered$LCL[1], mean(datCurrent$value), CL.Ordered$UCL[1]), col=c("grey80", "black", "grey80"), lty=c(2, 1, 2), lwd=2)
     panel.abline(h=c(CL.Ordered$LCL[2], CL.Ordered$UCL[2]), col=c("grey75", "grey75"), lty=c(3, 3), lwd=2)
     panel.abline(h=c(CL.Ordered$LCL[3], CL.Ordered$UCL[3]), col=c("grey45", "grey45"), lty=c(4, 4), lwd=2)
     panel.abline(h=c(CL.Ordered$LCL[4], CL.Ordered$UCL[4]), col=c("grey20", "grey20"), lty=c(5, 5), lwd=2)
     panel.abline(h=c(CL.Ordered$LCL[5], CL.Ordered$UCL[5]), col=c("black", "black"), lty=c(6, 6), lwd=2)
     panel.linejoin(x, y, horizontal = FALSE, lwd=2, col=c("grey"))
}
useOuterStrips(xyplot(value ~ Time | factor(BatchNo)*factor(BatchDesc), as.table=TRUE, par.strip.text=list(cex=1), layout=c(5,2,1), data=Data, pch=16, type=c("p"), col="grey", panel=mypanel, ylim=c(420, 475),  ylab="Tablet Weight", key=list(space="bottom", lines=list(col=c("grey80", "grey75", "grey45", "grey20", "black"), lty=c(2,3,4,5,6), lwd=2), text=list(c(LABEL1), cex=1.5))))



