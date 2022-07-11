## Biosimilar book chapter/section 6.5: Bayesian Approach for Equivalence Testing of Tier 2 QAs

################################################################
##  Effect of corelation on tier 2 testing
##  Simulated power for the unpooled (no correlation),
##                          pooled (no correlation),  
##                          unpooled (correlation in Reference)
################################################################
## H0: Fewer than 100p% of test lots lie between -delta and +delta
## Ha: At least 100p% of test lots lie between -delta and +delta

## Load R libraries
library(mvtnorm)   ## Library to generate multivariate normal data
library(ggplot2)   ## Library for graphics
library(lme4)      ## Library for linear mixed effects modeling
library(MCMCpack)  ## Bayesian library for simple modeling
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

#################################################################
##  Effect of corelation on tier 2 testing
##  Simulated power for the unpooled (no correlation),
##                          pooled (no correlation),  
##                          unpooled (correlation in Reference)
################################################################
## H0: Fewer than 100p% of test lots lie between -delta and +delta
## Ha: At least 100p% of test lots lie between -delta and +delta

set.seed(821)
d.sim = expand.grid(  mu=0, sigma=1, 
                       N.R=16, N.T = c(8, 10),
                       n.ref.corr = c(1, 2), rho.R = seq(0, 0.8, by=0.1),
                       n.test.corr = 1, rho.T = 0,
                       Pass.K2.80 = NA, Pass.K2.90 = NA, Pass.K3.80=NA, Pass.K3.90= NA
                    )
d.sim = subset(d.sim, !(rho.R > 0 & n.ref.corr==1))                       
d.sim = subset(d.sim, !(rho.R == 0 & n.ref.corr > 1))                       
params.R = d.sim[,c("N.R", "mu", "sigma", "rho.R", "n.ref.corr")]
params.T = d.sim[,c("N.T", "mu", "sigma", "rho.T", "n.test.corr")]
names(params.R)=names(params.T) = c("N", "mu", "sigma", "rho", "n.corr")


B = 10000                ## Number of Monte Carlo simulations.
for ( i in 1:nrow(params.R) ){
  d.R = create.data( as.list(params.R[i,]), B=B )
  d.T = create.data( as.list(params.T[i,]), B=B )
  pass.K2 = pass.K3 = matrix(NA, B, 2)
  for ( b in 1:B )
  {
    xbar.R = mean(d.R[b,])
    if ( d.sim[i,"n.ref.corr"]==1 )     ## If there is only one batch, then there is only within-batch variability
    {
      sigma.hat.R = sd(d.R[b,])
    }else{                              ## Otherwise, with multiple batches, there are two variance components: batch and within batch
      d.br = data.frame(Y=d.R[b,], Batch=factor(rep(1:(d.sim[i,"N.R"]/d.sim[i,"n.ref.corr"]), each=d.sim[i,"n.ref.corr"])))
      fit.R = lmer( Y~(1|Batch), data=d.br )
      sigma.hat.R = sqrt( VarCorr(fit.R)[["Batch"]][1,1] + sigma(fit.R)^2 )
    }  
    
    delta.K2 = xbar.R + 2*sigma.hat.R
    delta.K3 = xbar.R + 3*sigma.hat.R    
    
    fit.T = MCMCregress(Y~(1), data=data.frame(Y=d.T[b,]), seed=100*(i-1)+b)
    
    prq.K2 = pnorm( (delta.K2-fit.T[,"(Intercept)"])/sqrt(fit.T[,"sigma2"]) ) - pnorm( (-delta.K2-fit.T[,"(Intercept)"])/sqrt(fit.T[,"sigma2"]) )
    prq.K3 = pnorm( (delta.K3-fit.T[,"(Intercept)"])/sqrt(fit.T[,"sigma2"]) ) - pnorm( (-delta.K3-fit.T[,"(Intercept)"])/sqrt(fit.T[,"sigma2"]) )
    pass.K2[b,] = c( mean( prq.K2 >= 0.8 ), mean( prq.K2 >= 0.9 ) )
    pass.K3[b,] = c( mean( prq.K3 >= 0.8 ), mean( prq.K3 >= 0.9) )
  }
  
  d.sim$Pass.K2.80[i] = mean(pass.K2[,1] >= 0.95)
  d.sim$Pass.K2.90[i] = mean(pass.K2[,2] >= 0.95)
  d.sim$Pass.K3.80[i] = mean(pass.K3[,1] >= 0.95)
  d.sim$Pass.K3.90[i] = mean(pass.K3[,2] >= 0.95)
}


## Figure 6.6
d.plot = tidyr::pivot_longer(d.sim, cols=c("Pass.K2.80", "Pass.K2.90", "Pass.K3.80", "Pass.K3.90"))
d.plot$K = sapply( strsplit(as.character(d.plot$name), "\\."), function(z){ as.numeric(substring(z[2], 2)) })
d.plot$q = sapply( strsplit(as.character(d.plot$name), "\\."), function(z){ paste("q = 0.", z[3], sep="") } )

ggplot( d.plot, aes(rho.R, value, linetype=factor(q)) ) + geom_line(size=1) + 
  facet_grid(K~N.T, labeller=label_bquote(rows=K==.(K), cols=N[T]==.(N.T))) +
  xlab(expression(rho[R])) + ylab("Statistical Power") + 
  theme_grey(base_size=18) + theme(legend.title=element_blank())


