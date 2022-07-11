#######################################################################################
###
### Authors: Tony Pourmohamad and Robert Richardson
### Contacts: pourmoht@gene.com and richardson@stat.byu.edu 
### Code for reproducing the results in chapter XXX 
### of the book: Case Studies in Bayesian Methods for Biopharmaceutical CMC 
###
#######################################################################################


#######################################################################################
### Section 1.2.2   Code for the Bayesian Multivariate Normal Model
#######################################################################################

bmn = function(dis_data, B = 10000){
  
  X <- cbind(2-(dis_data[,1]==dis_data[1,1]),dis_data[-1])
  names(X)[1] <- "Group"
  dat_R <- X[X$Group==1,-1]
  dat_T <- X[X$Group==2,-1]
  nlocs <- ncol(dat_R)      
  nreps <- nrow(dat_R)      
  
  Ybar_R <- apply(dat_R,2,mean)
  Ybar_T <- apply(dat_T,2,mean)
  
  anoms <- as.matrix(t(t(dat_R)-Ybar_R))
  S2_R <- matrix(0,nlocs,nlocs)
  for(i in 1:nreps){
    S2_R <- S2_R + (anoms[i,])%*%t(anoms[i,])
  }
  
  anoms <- as.matrix(t(t(dat_T)-Ybar_T))
  S2_T <- matrix(0,nlocs,nlocs)
  for(i in 1:nreps){
    S2_T <- S2_T + (anoms[i,])%*%t(anoms[i,])
  }
  
  nrun <- B
  f2 <- rep(NA,nrun)
  delta <- rep(NA,nrun)
  
  muR <- matrix(NA, nrow = nlocs, ncol = B)
  muT <- matrix(NA, nrow = nlocs, ncol = B)
  
  for(i in 1:nrun){
    
    V_R <- riwish(nreps-1,S2_R)
    mu_R <- rmnorm(1,Ybar_R,V_R)
    
    V_T <- riwish(nreps-1,S2_T)
    mu_T <- rmnorm(1,Ybar_T,V_T)
    
    delta[i] <- max(abs(mu_R-mu_T))
    f2[i] <- 50*log10(1/sqrt(1+(1/length(mu_R))*sum((mu_R-mu_T)^2))*100)
    
    muR[,i] <- mu_R
    muT[,i] <- mu_T
  }
  
  return(list(delta = delta,f2 = f2, muR = muR, muT = muT))  
  
}


#######################################################################################
### Section 1.2.3   Code for the Hierarchical Gaussian Process Model
#######################################################################################

hgp <- function(dis_data,locs,B=1000,n_interp = 30,
                control=list(),adaptive=FALSE){
  draw_mu = function(dat,n,Sigma){
    Xbar <- apply(dat,2,mean)
    Sstar <- solve(solve(S)+n*solve(Sigma))
    mstar <- Sstar%*%(solve(S)%*%m+n*solve(Sigma)%*%(Xbar))
    rmnorm(1,mstar,Sstar)
  }
  
  draw_phiS <- function(){
    alpha <- ifelse(is.null(control$phiS_alpha),(max(locs)-min(locs))/(nloc/2),control$phiS_alpha)
    beta <- ifelse(is.null(control$phiS_beta),1,control$phiS_beta)
    phiSnew <- exp(rnorm(1,log(phiS),prop_phiS))
    Snew <- tau2S*covf2(D,phiSnew)+diag(nloc)*.0001
    alph <- dmnorm(mu_A,m,Snew,log=TRUE)+dmnorm(mu_B,m,Snew,log=TRUE)+dgamma(phiSnew,alpha,beta,log=TRUE)-
      dmnorm(mu_A,m,S,log=TRUE)-dmnorm(mu_B,m,S,log=TRUE)-dgamma(phiS,alpha,beta,log=TRUE)
    if(log(runif(1)) < alph){
      phiS <- phiSnew
      S <- Snew
      acceptS <- acceptS + 1
    }
    phiS
  }
  
  draw_phi <- function(){
    phinew <- exp(rnorm(1,log(phi),prop_phi))
    Sigma_Anew <- tau2_A*covdelta(D,phinew)
    Sigma_Bnew <- tau2_B*covdelta(D,phinew)
    
    alpha <- ifelse(is.null(control$phi_alpha),(max(locs)-min(locs))/(nloc/2),control$phi_alpha)
    beta <- ifelse(is.null(control$phi_beta),1,control$phi_beta)
    
    alph <- sum(dmnorm(dat_A,mu_A,Sigma_Anew,log=TRUE))+sum(dmnorm(dat_B,mu_B,Sigma_Bnew,log=TRUE))+dgamma(phinew,alpha,beta,log=TRUE)-
      sum(dmnorm(dat_A,mu_A,Sigma_A,log=TRUE)+dmnorm(dat_B,mu_B,Sigma_B,log=TRUE))-dgamma(phi,alpha,beta,log=TRUE)
    if(log(runif(1)) < alph){
      phi <- phinew
      Sigma_A <- Sigma_Anew
      Sigma_B <- Sigma_Bnew
      accept <- accept + 1
    }
    phi
  }
  
  draw_tau2S <- function(){
    alpha <- ifelse(is.null(control$tau2S_alpha),mean(diag(var(X[,-1]))),control$tau2S_alpha)
    beta <- ifelse(is.null(control$tau2S_beta),mean(diag(var(X[,-1]))),control$tau2S_beta)
    rigamma(1,alpha+nloc*2/2,beta+1/2*(mu_A-m)%*%solve(covf2(D,phiS),mu_A-m)+
              1/2*(mu_B-m)%*%solve(covf2(D,phiS),mu_B-m))
  }
  
  draw_tau2<- function(dat,mu,n){
    alpha <- ifelse(is.null(control$tau2_alpha),mean(diag(var(X[,-1]))),control$tau2_alpha)
    beta <- ifelse(is.null(control$tau2_beta),mean(diag(var(X[,-1]))),control$tau2_beta)
    wsse <- 0
    anom <- as.matrix(t(t(dat)-mu))
    wsse <- sum(diag(anom%*%solve(covdelta(D,phi))%*%t(anom)))
    rigamma(1,alpha+(nloc*n)/2,beta+1/2*wsse)
  }
  
  mu_post <- function(mu,tau2){
    S22 <- tau2S*covf2(D22,phiS)
    S12 <- tau2S*covf2(D12,phiS)
    mnew <- 50+t(S12)%*%solve(S,mu-50)
    signew <- S22 - t(S12)%*%solve(S,(S12))
    munew <- rmnorm(1,mnew,signew)
    c(mu,munew)
  }
  
  # Initialize Data
  X <- cbind(2-(dis_data[,1]==dis_data[1,1]),dis_data[-1])
  names(X)[1] <- "Group"
  nloc <- length(locs)
  otherlocs <- seq(min(locs),max(locs),length=n_interp)
  alllocs <- c(locs,otherlocs)
  nreps_A <- sum(X[,1]==1)
  nreps_B <- sum(X[,1]==2)
  D <- abs(outer(locs,locs,"-"))
  D22 <- abs(outer(otherlocs,otherlocs,"-"))
  D12 <- abs(outer(locs,otherlocs,"-"))
  dat_A <- X[X$Group==1,-1]
  dat_B <- X[X$Group==2,-1]
  dat <- X
  
  # Initialize mu, sigma2,tau2,phi,tau2S
  # 1st group is denoted with _A and second group with _B
  mu_A <- rep(0,nloc)
  mu_B <- rep(0,nloc)
  tau2 <- ifelse(is.null(control$tau2_starting),mean(diag(var(X[,-1]))),control$tau2_starting)
  tau2_A <- tau2 # Covariance of 1st group mean
  tau2_B <- tau2 #Covariance of 2nd group mean
  phi <- ifelse(is.null(control$phi_starting),(max(locs)-min(locs))/(nloc/2),control$phi_starting)
  phiS <- phi
  tau2S <- ifelse(is.null(control$tau2S_starting),tau2*5,control$tau2S_starting) # Covariance of mean function
  covdelta <- function(D,phi)matern(D,phi,3/2)
  Sigma_A <- tau2_A*covdelta(D,phi) # 1st group
  Sigma_B <- tau2_B*covdelta(D,phi) # second group
  m <- rep(50,nloc)
  covf2 <- function(D,phiS)exp(-abs(D/phiS))
  S <- tau2S*covf2(D,phiS)
  # Metropolis Hastings proposal parameters
  prop_phiS <- ifelse(is.null(control$phiS_prop),1,control$phiS_prop)
  prop_phi <- ifelse(is.null(control$phi_prop),1,control$phi_prop)
  
  # initialize collection vectors
  accept <- 0
  acceptS <- 0
  nrun <- B
  f2 <- rep(NA,nrun)
  delta <- rep(NA,nrun)
  mucollect_A <- matrix(0,length(alllocs),nrun)
  mucollect_B <- matrix(0,length(alllocs),nrun)
  phivec <- matrix(0,nrun,5)
  
  # Run chain
  for(i in 1:nrun){
    # Draw mu_A|.
    mu_A <- draw_mu(dat_A,nreps_A,Sigma_A)
    # Draw mu_B|.
    mu_B <- draw_mu(dat_B,nreps_B,Sigma_B)
    # Draw tau2S
    tau2S <- draw_tau2S()
    S <- tau2S*covf2(D,phiS)
    
    # Draw phiS
    phiS <- draw_phiS()
    S <- tau2S*covf2(D,phiS)+diag(nloc)*.0001
    
    # Draw tau2_A
    tau2_A <- draw_tau2(dat_A,mu_A,nreps_A)
    Sigma_A <- tau2_A*covdelta(D,phi)
    # Draw tau2_B
    tau2_B <- draw_tau2(dat_B,mu_B,nreps_B)
    Sigma_B <- tau2_B*covdelta(D,phi)
    
    # Draw phi
    phi <- draw_phi()
    phivec[i,] <- c(phi,phiS,sqrt(tau2_A),sqrt(tau2_B),sqrt(tau2S))
    
    # posterior predictive for f2
    mustar_A <- mu_post(mu_A,tau2_A)
    mustar_B <- mu_post(mu_B,tau2_B)
    mucollect_A[,i] <- mustar_A[order(alllocs)]
    mucollect_B[,i] <- mustar_B[order(alllocs)]
    f2[i] <- 50*log10(1/sqrt(1+(1/length(mustar_A))*sum((mustar_A-mustar_B)^2))*100)
    delta[i] <- max(abs(mustar_A-mustar_B))
    
    if(adaptive & i>25){
      prop_phi = sd(phivec[1:i,1])/3
      prop_phiS = sd(phivec[1:i,2])/3
    }
  }
  cov_pars <- data.frame(phivec)
  names(cov_pars) <- c("phi","psi","sigma_R","sigma_T","tau")
  mucollect_A <- data.frame(mucollect_A)
  names(mucollect_A) <- rep("muR", ncol(mucollect_A))
  mucollect_B <- data.frame(mucollect_B)
  names(mucollect_B) <- rep("muT", ncol(mucollect_B))
  mu_pars = cbind(mucollect_A,mucollect_B)
  

  return(list(delta = mean(delta), f2 = mean(f2), mcmc_chains = list(delta=delta,f2=f2,mu_pars = mu_pars,cov_pars = cov_pars)))
}


#######################################################################################
### Necessary packages to load
#######################################################################################

library(geoR)
library(ggplot2)
library(LaplacesDemon)
library(MCMCpack)
library(mnormt)
library(pscl)

#######################################################################################
### Example from section 1.3
#######################################################################################

# Read in the data from Ocana, et. al (2009)
X <- matrix(c(
19.78, 37.61, 48.53, 60.62, 63.34, 68.72, 75.76, 83.42, 
22.05, 36.74, 46.12, 53.90, 61.35, 67.35, 72.85, 77.88,
21.50, 34.58, 44.88, 53.46, 61.09, 66.83, 73.11, 77.28,
20.24, 35.23, 50.27, 66.49, 75.87, 79.84, 83.33, 87.57,
17.66, 33.51, 47.97, 59.29, 67.52, 70.05, 77.04, 80.88,
18.16, 33.75, 49.02, 58.26, 71.63, 74.98, 80.38, 84.33,
20.94, 33.31, 47.40, 59.26, 68.65, 75.82, 79.98, 82.96,
19.43, 31.40, 44.18, 56.20, 67.36, 75.03, 79.70, 83.98,
19.12, 32.93, 47.79, 60.19, 69.71, 76.26, 80.84, 83.95,
19.76, 31.67, 44.92, 56.79, 66.83, 74.41, 78.77, 82.30,
20.64, 32.42, 45.80, 58.11, 69.46, 76.81, 81.57, 85.83,
19.22, 32.14, 46.60, 58.16, 64.18, 75.79, 81.49, 83.94, 
35.70, 48.81, 57.05, 65.01, 71.92, 75.32, 80.64, 83.52,
36.10, 52.88, 61.79, 67.54, 74.47, 78.71, 80.71, 83.23,
41.66, 52.52, 63.22, 70.21, 76.71, 80.00, 83.25, 85.61,
35.49, 50.39, 58.99, 65.68, 72.59, 76.58, 80.52, 83.00,
36.06, 52.09, 63.46, 69.84, 75.23, 78.57, 81.53, 84.37,
35.53, 51.63, 62.59, 68.68, 75.64, 79.43, 82.99, 85.81,
32.70, 46.85, 57.51, 65.00, 71.26, 75.24, 79.18, 81.21,
33.85, 44.89, 53.64, 60.87, 67.34, 71.18, 73.51, 76.42,
32.33, 46.58, 56.46, 64.38, 70.30, 75.17, 77.83, 79.86,
32.80, 46.37, 56.66, 65.37, 70.04, 75.33, 77.11, 80.45,
33.86, 46.72, 57.45, 65.40, 71.05, 75.78, 76.22, 79.38,
32.46, 45.60, 55.04, 62.78, 68.54, 73.06, 75.73, 77.68),
nrow = 24, byrow = TRUE)

X <- data.frame(group = rep(c("Reference", "Test"), each = 12), X)

ybar_R <- apply(X[X$group == "Reference", -1], 2, mean)
ybar_T <- apply(X[X$group == "Test", -1], 2, mean)
f2 <- 50*log10(1/sqrt(1+(1/length(ybar_R))*sum((ybar_R-ybar_T)^2))*100)

#pdf("Dissolution_Data.pdf")
par(ps=15)
matplot(t(X[,-1]), pch = rep(c(19, 17), each = 12), 
        col = rep(c("gray65","black"), each = 12),
        xlab = "Time Points", 
        ylab = "Percentage Dissolved")
legend("bottomright", c("Reference", "Test"), pch = c(19, 17),
       col = c("gray65", "black"), bty = "n")
#dev.off()


#######################################################################################
### Equally spaced time points
#######################################################################################

set.seed(4)

tp <- seq(10, 80, 10)

B <- 10000
burnin <- 0.1 * B

method1 <- bmn(X, B)

method1$delta <- method1$delta[-(1:burnin)]
method1$f2<- method1$f2[-(1:burnin)]
method1$mu_R <- method1$muR[,-(1:burnin)]
method1$mu_T <- method1$muT[,-(1:burnin)]

apply(method1$mu_R, 1, ESS)
apply(method1$mu_T, 1, ESS)

N <- B-burnin # Number of samples
chains <- data.frame(samples = rep(c(1:N, 1:N), each = ncol(X) - 1),
                     group = rep(c("muR", "muT"), each = (ncol(X) - 1)*N), 
                     timepoint = paste("Time Point", rep(1:(ncol(X) - 1), 2 * N)),
                     values = c(c(method1$mu_R), c(method1$mu_T)))

g <- ggplot(chains, aes(samples, values)) + 
  geom_line() + 
  labs(x = "Iterations", y = "Posterior Sample Values") + 
  facet_wrap(group~timepoint) +
  theme(text = element_text(size = 16))
#ggsave("muchains.pdf", width = 17.2, height = 8 )

g <- ggplot(chains, aes(values)) + 
  geom_density() + 
  labs(y = "Posterior Density") + 
  facet_wrap(group~timepoint, scales = "free") +
  theme(text = element_text(size = 16))
#ggsave("muposteriors.pdf", width = 17.2, height = 8)

q <- sum(method1$delta < 15 & method1$f2 > 50) / (B - burnin)

alpha <- 0.05
f2_cred <- c(quantile(method1$f2, alpha / 2),quantile(method1$f2, 1 - alpha / 2))
delta_cred <- c(quantile(method1$delta, alpha / 2),quantile(method1$delta, 1 - alpha / 2))

#pdf("bmnF2.pdf")
par(ps = 15)
plot((density(method1$f2)), xlab = expression(F[2]), ylab = "Posterior Density", 
      main = expression("Posterior Density of "*F[2]))
abline(v = 50, lty = 2)
#dev.off()

#pdf("bmndelta.pdf")
par(ps = 15)
plot((density(method1$delta)), xlab = expression(delta), ylab = "Posterior Density",
     main = expression("Posterior Density of "*delta))
abline(v = 15, lty = 2)
#dev.off()




B <- 100000
burnin <- 0.1 * B
thin <- 10

method2 <- hgp(X, locs = tp, B = B, n_interp = 100)

phi <- method2$mcmc_chains$cov_pars$phi[-c(1:burnin)][seq(1, (B-burnin), thin)]
psi <- method2$mcmc_chains$cov_pars$psi[-c(1:burnin)][seq(1, (B-burnin), thin)]
sigma_R <- method2$mcmc_chains$cov_pars$sigma_R[-c(1:burnin)][seq(1, (B-burnin), thin)]
sigma_T <- method2$mcmc_chains$cov_pars$sigma_T[-c(1:burnin)][seq(1, (B-burnin), thin)]
tau <- method2$mcmc_chains$cov_pars$tau[-c(1:burnin)][seq(1, (B-burnin), thin)]

ESS(phi)
ESS(psi)
ESS(sigma_R)
ESS(sigma_T)
ESS(tau)

grid <- sort(c(tp, seq(min(tp), max(tp), length = 100)))
#pdf("equalHGP.pdf")
par(ps=15)
matplot(tp, t(X[,-1]), pch = rep(c(19, 17), each = 12), 
        col = rep(c("gray65","black"), each = 12),
        xlab = "Time (Minutes)", 
        ylab = "Percentage Dissolved")
legend("bottomright", c("Reference", "Test"), pch = c(19, 17),
       col = c("gray65", "black"), bty = "n")

grid1 <- (1:B)[-(1:burnin)][seq(1, (B-burnin), thin)]
grid2 <- ((B+1):(2*B))[-(1:burnin)][seq(1, (B-burnin), thin)]
lines(grid, apply(method2$mcmc_chains$mu_pars[,grid1], 1, mean), col = "gray65", lwd = 2)
lines(grid, apply(method2$mcmc_chains$mu_pars[,grid2], 1, mean), col = "black", lwd = 2)

lower <- apply(method2$mcmc_chains$mu_pars[,grid1], 1, quantile, prob = 0.025)
upper <- apply(method2$mcmc_chains$mu_pars[,grid1], 1, quantile, prob = 0.975)
polygon(c(grid, rev(grid)), c(lower, rev(upper)), col = alpha("gray65",.2), border = NA)

lower <- apply(method2$mcmc_chains$mu_pars[,grid2], 1, quantile, prob = 0.025)
upper <- apply(method2$mcmc_chains$mu_pars[,grid2], 1, quantile, prob = 0.975)
polygon(c(grid, rev(grid)), c(lower, rev(upper)), col = alpha("black",.2), border = NA)
#dev.off()

chains <- data.frame(samples = rep(1:((B-burnin)/thin), times = 5), 
                     parameter = rep(c("phi", "psi", "sigma_R", "sigma_T", "tau"), each = (B-burnin)/thin),
                     values = c(phi, psi, sigma_R, sigma_T, tau))
chains$parameter <- factor(chains$parameter, 
                           labels = c(expression(phi), expression(psi), expression(sigma[R]),
                                      expression(sigma[T]), expression(tau)))

g <- ggplot(chains, aes(samples, values)) + 
  geom_line() + 
  labs(x = "Iterations", y = "Posterior Sample Values") + 
  facet_wrap(~parameter, scales = "free", labeller = label_parsed) +
  theme(text = element_text(size = 22))
#ggsave("equaltraceplots.pdf", width = 17, height = 9)

g <- ggplot(chains, aes(values)) + 
  geom_density() + 
  labs(x = "Values", y = "Posterior Density") + 
  facet_wrap(~parameter, scales = "free", labeller = label_parsed) +
  theme(text = element_text(size = 22))
#ggsave("equalposteriors.pdf", width = 17, height = 9)

q <- sum(method2$mcmc_chains$delta[grid1] < 15 & method2$mcmc_chains$f2[grid1] > 50) / ((B - burnin)/thin)

alpha <- 0.05
f2_cred <- c(quantile(method2$mcmc_chains$f2[grid1], alpha / 2),quantile(method2$mcmc_chains$f2[grid1], 1 - alpha / 2))
delta_cred <- c(quantile(method2$mcmc_chains$delta[grid1], alpha / 2),quantile(method2$mcmc_chains$delta[grid1], 1 - alpha / 2))

#pdf("equalF2.pdf")
par(ps = 15)
plot(density(method2$mcmc_chains$f2[grid1]), xlab = expression(F[2]), ylab = "Posterior Density", 
     main = expression("Posterior Density of "*F[2]))
abline(v = 50, lty = 2)
#dev.off()

#pdf("equaldelta.pdf")
par(ps = 15)
plot(density(method2$mcmc_chains$delta[grid1]), xlab = expression(delta), ylab = "Posterior Density",
     main = expression("Posterior Density of "*delta))
abline(v = 15, lty = 2)
#dev.off()


#######################################################################################
### Unequally spaced time points
#######################################################################################

set.seed(4)

B <- 100000
burnin <- 0.1 * B
thin <- 10

tp <- c(5, 15, 30, 45, 60, 90, 120, 180)

method2 <- hgp(X, locs = tp, B = B, n_interp = 100, 
               control = list(phiS_alpha = 100, phiS_beta = .1))


phi <- method2$mcmc_chains$cov_pars$phi[-c(1:burnin)][seq(1, (B-burnin), thin)]
psi <- method2$mcmc_chains$cov_pars$psi[-c(1:burnin)][seq(1, (B-burnin), thin)]
sigma_R <- method2$mcmc_chains$cov_pars$sigma_R[-c(1:burnin)][seq(1, (B-burnin), thin)]
sigma_T <- method2$mcmc_chains$cov_pars$sigma_T[-c(1:burnin)][seq(1, (B-burnin), thin)]
tau <- method2$mcmc_chains$cov_pars$tau[-c(1:burnin)][seq(1, (B-burnin), thin)]

ESS(phi)
ESS(psi)
ESS(sigma_R)
ESS(sigma_T)
ESS(tau)

grid <- sort(c(tp, seq(min(tp), max(tp), length = 100)))
#pdf("unequalHGP.pdf")
par(ps=15)
matplot(tp, t(X[,-1]), pch = rep(c(19, 17), each = 12), 
        col = rep(c("grey65","black"), each = 12),
        xlab = "Time (Minutes)", 
        ylab = "Percentage Dissolved")
legend("bottomright", c("Reference", "Test"), pch = c(19, 17),
       col = c("grey65", "black"), bty = "n")

grid1 <- (1:B)[-(1:burnin)][seq(1, (B-burnin), thin)]
grid2 <- ((B+1):(2*B))[-(1:burnin)][seq(1, (B-burnin), thin)]
lines(grid, apply(method2$mcmc_chains$mu_pars[,grid1], 1, mean), col = "grey65", lwd = 2)
lines(grid, apply(method2$mcmc_chains$mu_pars[,grid2], 1, mean), col = "black", lwd = 2)

lower <- apply(method2$mcmc_chains$mu_pars[,grid1], 1, quantile, prob = 0.025)
upper <- apply(method2$mcmc_chains$mu_pars[,grid1], 1, quantile, prob = 0.975)
polygon(c(grid, rev(grid)), c(lower, rev(upper)), col = alpha("grey65",.2), border = NA)

lower <- apply(method2$mcmc_chains$mu_pars[,grid2], 1, quantile, prob = 0.025)
upper <- apply(method2$mcmc_chains$mu_pars[,grid2], 1, quantile, prob = 0.975)
polygon(c(grid, rev(grid)), c(lower, rev(upper)), col = alpha("black",.2), border = NA)
#dev.off()

chains <- data.frame(samples = rep(1:((B-burnin)/thin), times = 5), 
                     parameter = rep(c("phi", "psi", "sigma_R", "sigma_T", "tau"), each = (B-burnin)/thin),
                     values = c(phi, psi, sigma_R, sigma_T, tau))
chains$parameter <- factor(chains$parameter, 
                           labels = c(expression(phi), expression(psi), expression(sigma[R]),
                                      expression(sigma[T]), expression(tau)))

g <- ggplot(chains, aes(samples, values)) + 
  geom_line() + 
  labs(x = "Iterations", y = "Posterior Sample Values") + 
  facet_wrap(~parameter, scales = "free", labeller = label_parsed) +
  theme(text = element_text(size = 22))
#ggsave("unequaltraceplots.pdf", width = 17, height = 9)

g <- ggplot(chains, aes(values)) + 
  geom_density() + 
  labs(x = "Values", y = "Posterior Density") + 
  facet_wrap(~parameter, scales = "free", labeller = label_parsed) +
  theme(text = element_text(size = 22))
#ggsave("unequalposteriors.pdf", width = 17, height = 9)

q <- sum(method2$mcmc_chains$delta[grid1]  < 15 & method2$mcmc_chains$f2[grid1]   > 50) / ((B - burnin)/thin)

c(quantile(method2$mcmc_chains$f2[grid1] , alpha / 2),quantile(method2$mcmc_chains$f2[grid1] , 1 - alpha / 2))
c(quantile(method2$mcmc_chains$delta[grid1] , alpha / 2),quantile(method2$mcmc_chains$delta[grid1] , 1 - alpha / 2))

#pdf("unequalF2.pdf")
par(ps = 15)
plot(density(method2$mcmc_chains$f2[grid1]), xlab = expression(F[2]), ylab = "Posterior Density", 
     main = expression("Posterior Density of "*F[2]))
abline(v = 50, lty = 2)
#dev.off()

#pdf("unequaldelta.pdf")
par(ps = 15)
plot(density(method2$mcmc_chains$delta[grid1]), xlab = expression(delta), ylab = "Posterior Density",
     main = expression("Posterior Density of "*delta))
abline(v = 15, lty = 2)
#dev.off()





