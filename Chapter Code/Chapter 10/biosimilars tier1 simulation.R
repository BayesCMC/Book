## Biosimilar book chapter/section 6.4.2: Bayesian Approach for Equivalence Testing of Tier 1 QAs

#################################################################
##  Simulated power for the unpooled (no correlation),
##                          pooled (no correlation),  
##                          unpooled (correlation in Reference)
################################################################



## Load R libraries
library(mvtnorm)   ## Library to generate multivariate normal data
library(ggplot2)   ## Library for graphics
library(MCMCpack)  ## Bayesian library for simple modeling
library(MCMCglmm)  ## Bayesian library for linear hierarchical modeling
library(tidyr)     ## Library for data manipulation


## Function to create data set for simulations
create.data = function(params, B=1)
{
 ## params is a list with elements:
   ## N = sample size
   ## mu = true mean
   ## sigma = true standard deviation
   ## rho = true correlation among product measurements
   ## n.corr = number of DP batches made from the same DS lots (i.e., # of correlated units per DS lot)
   ##          Note:  N/n.corr must be an integer
 ## B = # of samples to generate


 if ( params$rho==0 | params$n.corr==1 )  ## Create data in the absence of correlation
   Y = matrix( rnorm(B*params$N, mean=params$mu, sd=params$sigma), B, params$N )
 else    ## Create data in the presence of correlation
 {
   Sigma = diag(1, params$N)

   if ( is.null(params$n.corr) )
     params$n.corr = 2
   ## Check if covariance structure is possible
   if ( params$N %% params$n.corr != 0 )
      stop("Partial correlation must have params: 'N' %% 'n.corr' == 0")
   for ( j in 1:(ncol(Sigma)/params$n.corr) )
   {
     index = params$n.corr*(j-1)+1:params$n.corr
     Sigma[index,index][upper.tri(Sigma[index,index], diag=FALSE)] = Sigma[index,index][lower.tri(Sigma[index,index], diag=FALSE)] = params$rho
   }
   Y = rmvnorm(B, mean=rep(params$mu, params$N), sigma=params$sigma^2*Sigma)
 }

 return( Y )
}

set.seed(820)
d.sim = expand.grid( N.R=c(10, 20), mu.R=0,           sigma.R=1, rho.R=c(0, 0.8),
                      N.T = 10,      mu.T=c(1/8, 1.5), sigma.T=c(0.5, 1, 2), rho.T=0,
                      n.ref.corr = c(1, 2, 5), n.test.corr=1,          ## This is "r" in Yang, Novick, and Burdick
                      K = 1.5,
                      Pass.unpooled = NA, Pass.pooled = NA, Pass.unpooled.corr = NA
)
d.sim = subset(d.sim, !(rho.R==0.8 & n.ref.corr==1))                       
d.sim = subset(d.sim, !(rho.R==0 & n.ref.corr > 1))                       
params.R = d.sim[,c("N.R", "mu.R", "sigma.R", "rho.R", "n.ref.corr")]
params.T = d.sim[,c("N.T", "mu.T", "sigma.T", "rho.T", "n.test.corr")]
names(params.R)=names(params.T) = c("N", "mu", "sigma", "rho", "n.corr")

alpha=0.05  ## Significance limit
B=10000

  ## Prior for MCMCglmm. For R and G1, these are both sigma^2 ~ inverse-chisquare(df=0.001).
prior.Rcorr = list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002))) 

for ( i in 1:nrow(params.R) ){
  
  d.R = create.data( as.list(params.R[i,]), B=B )
  d.T = create.data( as.list(params.T[i,]), B=B )
  p.pooled = p.unpooled = p.unpooled.corr = rep(NA, B)
  for ( b in 1:B )
  {
    d.b = data.frame( Y=c(d.R[b,], d.T[b,]), Type=rep(c("R", "T"), c(d.sim[i,"N.R"], d.sim[i,"N.T"])) ) 
      
    fit.pool = MCMCregress(Y~Type, data=d.b, seed=10*(i-1)+b)
    muDiff.over.sigmaR = fit.pool[,"TypeT"]/sqrt(fit.pool[,"sigma2"])
    p.pooled[b] = mean(abs(muDiff.over.sigmaR) < d.sim$K[i])
      
    fit.R = MCMCregress(Y~(1), data=subset(d.b, Type=="R"), seed=1000*(i-1)+b)
    fit.T = MCMCregress(Y~(1), data=subset(d.b, Type=="T"), seed=100*(i-1)+b)    
    muDiff.over.sigmaR = (fit.T[,"(Intercept)"]-fit.R[,"(Intercept)"])/sqrt(fit.R[,"sigma2"])
    p.unpooled[b] = mean(abs(muDiff.over.sigmaR) < d.sim$K[i])
      
    ## Default priors: B$mu=0 and B$V=1e+10, 
    ## The priors for the variance structures (R and G) are lists
    ##   with the expected (co)variances (V) and degree of belief parameter (nu) for the inverseWishart, 
    ##   and also the mean vector (alpha.mu) and covariance matrix (alpha.V) for the redundant working parameters.
    ##   The defaults are nu=0, V=1, alpha.mu=0, and alpha.V=0
      
    if ( d.sim[i,"n.ref.corr"] > 1 ){
      d.br = data.frame(Y=d.R[b,], Batch=factor(rep(1:(d.sim[i,"N.R"]/d.sim[i,"n.ref.corr"]), each=d.sim[i,"n.ref.corr"])))
      fit.Rcorr = MCMCglmm(fixed=Y~(1), random=~Batch, data=d.br, nitt=50000, thin=1, burnin=10000,
                           prior=prior.Rcorr, verbose=FALSE)                    
      muDiff.over.sigmaR = (fit.T[,"(Intercept)"]-fit.Rcorr$Sol[,"(Intercept)"])/sqrt(fit.Rcorr$VCV[,"Batch"]+fit.Rcorr$VCV[,"units"])
      p.unpooled.corr[b] = mean(abs(muDiff.over.sigmaR) < d.sim$K[i])
    }else{
      p.unpooled.corr[b] = p.unpooled[b]
    }
      
      print(b)
  }
  pass.pooled = mean( p.pooled > 1-alpha )
  pass.unpooled = mean( p.unpooled > 1 - alpha )
  pass.unpooled.corr = mean( p.unpooled.corr > 1 - alpha )
    
  d.sim[i,c("Pass.unpooled", "Pass.pooled", "Pass.unpooled.corr")] = c(pass.pooled, pass.unpooled, pass.unpooled.corr)
}


d.plot = tidyr::pivot_longer(d.sim, cols=c("Pass.unpooled", "Pass.pooled", "Pass.unpooled.corr"))

## Test size (type 1 error)
d.plot = subset(d.plot, !(name=="Pass.unpooled.corr" & n.ref.corr==1))
d.plot$x = sapply( strsplit(as.character(d.plot$name), "\\."), function(z){ paste(z[-1], collapse="+") } )              
d.plot$sigT = d.plot$sigma.T


ggplot( subset(d.plot, mu.T==1.5), aes(x, value, shape=factor(N.R)) ) + geom_point(size=3) + 
  facet_grid(sigma.T~n.ref.corr, labeller = label_bquote(rows=(sigma[R]/sigma[T])==.(1/sigma.T), cols=r==.(n.ref.corr)), scale="free_x") +
  xlab("") + ylab("Type I Error") +
  geom_hline( yintercept=0.05, size=1, linetype="dashed" ) +
  theme_grey(base_size=18) + theme(legend.position="top") +
  scale_shape_discrete(name = expression(N[R]))


## Statistical power.  Subset to those scenarios for which type 1 error <= alpha (allowing for Monte Carlo error)  
test.size.05 = subset(d.plot, mu.T==1.5 & value <= 0.07)
test.size.05$Fac = with(test.size.05, paste(N.R, sigma.T, n.ref.corr, name))
d.plot$Fac = with(d.plot, paste(N.R, sigma.T, n.ref.corr, name))

d.power = subset(d.plot, Fac %in% test.size.05$Fac & mu.T==0.125)
ggplot( d.power, aes(x, value, shape=factor(N.R)) ) + geom_point(size=3) + 
  facet_grid(sigma.T~n.ref.corr, labeller = label_bquote(rows=(sigma[R]/sigma[T])==.(1/sigma.T), cols=r==.(n.ref.corr)), scale="free_x") +
  xlab("") + ylab("Statistical Power") +
  geom_hline( yintercept=0.80, size=1, linetype="dashed" ) +
  theme_grey(base_size=18) + theme(legend.position="top") +
  scale_shape_discrete(name = expression(N[R]))