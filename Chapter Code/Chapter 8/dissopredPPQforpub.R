# Program Name: dissopredPPQ.R
# Author: Katherine Giacoletti
# Date: 18Apr2022
#
# Purpose: model disso from 4 PPQ batches to assess process capability as it relates
#
# Data: dissoPPQ.csv
#
# Updates: 
###########################################################
options(repos = 'https://cran.rstudio.com/')
require(MCMCglmm); require(MCMCpack); require(lme4); require(car); require(boot)
require(ggplot2); require(tolerance)

dPPQ = read.csv("dissoPPQ.csv")
dPPQ$Diss = dPPQ$Pct.Dissolved
dPPQ$Batch = as.factor(dPPQ$Batch)

plot(dPPQ2$Batch,dPPQ2$Diss, xlab = "Batch", ylab = "% Dissolved",ylim=c(84,95))
abline(h=85,lty=2)

## no location/no replicate samples per batch
tapply(d$Diss, d$Batch, normtol.int)

fmPPQ = lmer( Diss ~ 1 + ( 1 | Batch), data = d)
summary(fmPPQ)
vcPPQ=as.data.frame(VarCorr(fmPPQ))

bmodPPQ = MCMCglmm( Diss ~ 1, random = ~ Batch, data=d, 
                    nitt = 2505000, burnin = 5000, thin = 100,
                    prior = list(G = list(G1 = list(V = max(vcPPQ$vcov[1],1e-10), nu = 1)))
)

summary(bmodPPQ)
save(bmodPPQ,file="BayesPPQ.RData")
load("BayesPPQ.RData")
plot(bmodPPQ)
traceplot(bmodPPQ$Sol[,1])
par(mfrow=c(2,1))
traceplot(bmod$VCV[,1:2])
apply(sqrt(bmodPPQ$VCV),2,median)
HPDinterval(sqrt(bmodPPQ$VCV[,1]), prob = .95)
HPDinterval(sqrt(bmodPPQ$VCV[,2]), prob = .95)
ecdf(bmodPPQ$Sol[,1])(84.5)
par(mfrow=c(1,1))

Nsim=length(bmodPPQ$Sol[,1])
s1samp<-list()
s2samp<-list()
s3samp<-list()
fs1<-logical(length=Nsim) #T/F fail S1
fs2<-logical(length=Nsim) #T/F fail S2
fs3<-logical(length=Nsim) #T/F fail S3
pfail.s1=0 # prob fail stage 1
pfail.s2=0 # prob fail stage 2
pfail.s3=0 # prob fail stage 3
for (m in 1:Nsim){
     s3samp[[m]]=rnorm(24,bmodPPQ$Sol[m,1],sqrt(sum(bmodPPQ$VCV[m,1:2])))
     s2samp[[m]]=s3samp[[m]][1:12]
     s1samp[[m]]=s3samp[[m]][1:6]
     
     fs1[m]=isTRUE(min(s1samp[[m]])<89.5) # Q =  85, Q+5=90
     
     fs2[m]=isTRUE((fs1[m]=="TRUE")&((min(s2samp[[m]])<69.5)|(mean(s2samp[[m]])<84.5)))
     
     fs3[m]=isTRUE((fs2[m]=="TRUE")&( (min(s3samp[[m]])<59.5)|(mean(s3samp[[m]])<84.5)|(s3samp[[m]][[2]]<69.5) ))
     
}

nfs1=length(fs1[fs1==TRUE]) # number failing Stage 1
nfs2=length(fs2[fs2==TRUE]) # number failing Stage 2
nfs3=length(fs3[fs3==TRUE]) # number failing Stage 3
pfail.s1=nfs1/Nsim
#pfail.s2=ifelse(nfs1==0,NA,nfs2/nfs1)
pfail.s3=ifelse((nfs2==0),NA,nfs3/nfs2)
pfail.USP711 = nfs3/Nsim

