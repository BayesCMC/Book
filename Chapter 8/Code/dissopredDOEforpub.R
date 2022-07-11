# Program Name: dissopredDOE.R
# Author: Katherine Giacoletti
# Date: 4jan2022
#
# Purpose: model disso from DOE to assess intra-batch capability as it relates
#         to 2 factors
#
# Data: disso ex.csv
#
# Updates: 
###########################################################
options(repos = 'https://cran.rstudio.com/')
require(MCMCglmm); require(MCMCpack); require(lme4); require(car); require(boot)

d1 = read.csv("disso ex.csv")
d1$Vessel=as.factor(d1$Vessel)
d1$y=d1$X.Dissolved
d1$x1 = d1$Tablet.Thickness.coded
d1$x2 = d1$Coating.Level.coded

fmod = lmer(y ~ x1 + x2 + I(x1^2) + I(x1*x2) + (1|Location), data=d1)
summary(fmod)
vc=as.data.frame(VarCorr(fmod))

bmod = MCMCglmm(y ~ x1 + x2 + I(x1^2) + I(x1*x2), data = d1, 
                random = ~Location, 
                nitt = 2503000, burnin = 3000, thin = 100,
                prior = list(G = list(G1 = list(V = max(vc$vcov[1],0.000001), nu = 1)))
)
save(bmod,file="DOE disso model.RData")
load("DOE disso model.RData")
summary(bmod)
plot(bmod)
par(mfrow=c(3,2))
traceplot(bmod$Sol[,1:5])
par(mfrow=c(2,1))
traceplot(bmod$VCV[,1:2])
quantile(sqrt(bmod$VCV[,1]),probs=c(0.025,0.5,0.975))
quantile(sqrt(bmod$VCV[,2]),probs=c(0.025,0.5,0.975))
par(mfrow=c(1,1))

## Generate, for each posterior mean & SD, a Stage 1 sample (n=6), Stage 2 sample (n=12), and Stage 3 sample (n=24)
## For each sample, compare to acceptance criteria for that stage of testing
## Count how many pass at each stage and calculate proportion out of total posterior samples = predictive probability of passing at each stage
## Do this over a fine range of x2 & x2 
new.x1 = seq(-1,1,by=0.1)
new.x2 = seq(-1,1,by=0.1)
stage = seq(1,3)

Risk1 = matrix(nrow = length(new.x1), ncol = length(new.x2))
Risk2 = matrix(nrow = length(new.x1), ncol = length(new.x2))
Risk3 = matrix(nrow = length(new.x1), ncol = length(new.x2))
for (i in 1:length(new.x1)){
  for(j in 1:length(new.x2)){
    Nsim=length(bmod$Sol[,1])
    s1samp<-list()
    s2samp<-list()
    s3samp<-list()
    fs1<-logical(length=Nsim) #T/F fail S1
    fs2<-logical(length=Nsim) #T/F fail S2
    fs3<-logical(length=Nsim) #T/F fail S3
    pfail.s1=0 # prob fail stage 1
    pfail.s2=0 # prob fail stage 2
    pfail.s3=0 # prob fail stage 3
    for (m in 1:length(Nsim)){
      s3samp[[m]]=rnorm(24,bmod$Sol[m,1]+bmod$Sol[m,2]*new.x1[i]+bmod$Sol[m,3]*new.x2[j]+bmod$Sol[m,4]*new.x1[j]^2++bmod$Sol[m,5]*new.x1[j]*new.x2[j],sqrt(sum(bmod$VCV[m,])))
      s2samp[[m]]=s3samp[[m]][1:12]
      s1samp[[m]]=s3samp[[m]][1:6]
      
      fs1[m]=isTRUE(min(s1samp[[m]])<84.5)
      
      fs2[m]=isTRUE((fs1[m]=="TRUE")&((min(s2samp[[m]])<64.5)|(mean(s2samp[[m]])<79.5)))
      
      fs3[m]=isTRUE((fs2[m]=="TRUE")&( (min(s3samp[[m]])<54.5)|(mean(s3samp[[m]])<79.5)|(s3samp[[m]][[2]]<64.5) ))
      
    }
    
    nfs1=length(fs1[fs1==TRUE]) # number failing Stage 1
    nfs2=length(fs2[fs2==TRUE]) # number failing Stage 2
    nfs3=length(fs3[fs3==TRUE]) # number failing Stage 3
    pfail.s1=nfs1/Nsim
    pfail.s2=ifelse(nfs1==0,NA,nfs2/nfs1)
    pfail.s3=ifelse((nfs2==0),NA,nfs3/nfs2)
    pfail.USP711=nfs3/Nsim
      
    Risk1[i,j] = pfail.s1
    Risk2[i,j] = pfail.s2
    Risk3[i,j] = pfail.s3
  }
}

write.csv(Risk1,"Predprob Fail Stage 1 USP711 by X1X2.csv")
write.csv(Risk2,"Predprob Fail Stage 2 USP711 by X1X2.csv")
write.csv(Risk3,"Predprob Fail Stage 3 USP711 by X1X2.csv")
