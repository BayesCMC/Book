## Load R libraries
library(brms)       ## Bayesian library built on Stan
library(MCMCpack)   ## Bayesian library for simple modeling
library(MCMCglmm)   ## Bayesian library for hierarchical linear modeling
library(ggplot2)    ## Graphics library

## Chapter/section 6.4.1 Illustrations with brms

  ## Create data
nR = 20
xR = c(88.7, 92.1, 93.5, 93.7, 94.4, 94.6, 95.1, 95.3, 95.9, 97, 
       97.8, 97.9, 98.2, 98.8, 99, 100.2, 101.5, 103.5, 107.6, 109)
nT = 10
xT = c(89.2, 98.5, 100, 100.2, 101.4, 103.6, 104.2, 105.3, 108.9, 112.8)
d = data.frame(x=c(xR, xT), Source=factor( rep(c("RP", "TP"), c(nR, nT)) ))

## Model M1
  ## Set priors  
bprior = prior(normal(100, 20), class=b) + prior(cauchy(0, 1), class = sigma)
  ## Fit model
fit.m1 = brm( x~Source-1, data=d, prior=bprior, chains=4, warmup=2000, iter=2000+10000, cores=4)
## Extract posterior parameters
post.m1 = as.matrix(fit.m1)
## Test eq 6.4
ratio.m1 = (post.m1[,"b_SourceRP"] -
post.m1[,"b_SourceTP"])/post.m1[,"sigma"]
p.m1 = mean( ratio.m1 > -1.5 & ratio.m1 < 1.5 )

## The following performs the MCMC for Model M1 much faster than brm(), but we must use normal/inverse-gamma priors
fit.m1 = MCMCregress( x~Source-1, data=d, seed=241 )


## Model M2
  ## Set priors
bprior = prior(normal(100, 20), class=Intercept) + prior(cauchy(0, 1), class = sigma)
  ## Fit separate models for RP and TP
fit.m2.R = brm( x~(1), data=subset(d, Source=="RP"), prior=bprior,chains=4, warmup=2000,
             iter=2000+10000, cores=4 )
fit.m2.T = brm( x~(1), data=subset(d, Source=="TP"), prior=bprior,chains=4, warmup=2000,
             iter=2000+10000, cores=4 )
  ## Extract posterior parameters
post.m2.R = as.matrix(fit.m2.R)
post.m2.T = as.matrix(fit.m2.T)

  ## Test eq 6.4
ratio.m2 = (post.m2.R[,"b_Intercept"] -
post.m2.T[,"b_Intercept"]) /post.m2.R[,"sigma"]
p.m2 = mean( ratio.m2 > -1.5 & ratio.m2 < 1.5 )  

## The following performs the MCMC for Model M2 much faster than brm(), but we must use normal/inverse-gamma priors
fit.m2.R = MCMCregress( x~(1), data=subset(d, Source=="RP"), seed=5152 )
fit.m2.T = MCMCregress( x~(1), data=subset(d, Source=="TP"), seed=890 )


## Model M3
set.seed(424)
d.R = data.frame( Lot=factor(rep(1:5, each=4)), x=NA )
sigmaR = 5
rhoR = 0.8
sigmaL = sqrt(rhoR)*sigmaR
sigmaE = sqrt(1-rhoR)*sigmaR
d.R$x = 100 + rnorm(5, mean=0, sd=sigmaL)[unclass(d.R$Lot)] + rnorm(nrow(d.R), mean=0, sd=sigmaE)

ggplot(d.R, aes(Lot, x) ) + geom_jitter(width=0.05) + ylab("Value")

  ## Set priors
bprior.m3 = prior(normal(100, 20), class=Intercept) +
                prior(cauchy(0, 1), class = sd) + ## prior for sigmaL
                prior(cauchy(0, 1), class=sigma) ## prior for sigmaE
  ## Fit model            
fit.m3.R = brm( x~(1|Lot), data=d.R, prior=bprior.m3, chains=4, 
    warmup=10000, iter=10000+25000, cores=4 )
    
  ## Refit fit.m2.T to obtain same number of posterior samples as fit.m3.R
fit.m2.T = brm( x~(1), data=subset(d, Source=="TP"), prior=bprior, chains=4, 
    warmup=10000, iter=10000+25000, cores=4 )
    
    
  ## Extract posterior parameters
post.m3.R = as.matrix(fit.m3.R)
post.m2.T = as.matrix(fit.m2.T)

  ## Get posterior for sigma.R = sqrt( sigmaL^2 + sigmaE^2 )
sigma.R = sqrt( post.m3.R[,"sd_Lot__Intercept"]^2 + post.m3.R[,"sigma"]^2 )
  ## Test eq 4
ratio.m3 = (post.m3.R[,"b_Intercept"] - post.m2.T[,"b_Intercept"])/sigma.R
p.m3 = mean( ratio.m3 > -1.5 & ratio.m3 < 1.5 )



## The following performs the MCMC for Model M3 much faster than brm(), but we must use normal/inverse-chisq priors
## The priors are sigma[L]^2 ~ inverse-chisquare(df=0.001)
## sigma[E]^2 ~ inverse-chisquare(df=0.001)
prior.Rcorr = list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))

fit.m3.R = MCMCglmm(fixed=x~(1), random=~Lot, data=d.R, nitt=50000, thin=1, burnin=10000, prior=prior.Rcorr)

