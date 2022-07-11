

library(rstan);library(ggplot2); library(dplyr)


## prior predictive simulation
priorpredcode = "
  data {
    real IP;
    real  multiplier.mu;
    real  multiplier.sigma;
  }
  transformed data {
    real H;
    real musigma;
    H = multiplier.sigma *IP;
    musigma = multiplier.mu *IP;
  }  
  parameters {
    real<lower=1> nu;
    real mu;
    real<lower=0> sigma;
  }
  model {
    nu ~ gamma(2, 0.1);
    mu ~ normal(100, musigma);
    sigma ~ normal(0, H);
   }
  generated quantities{
    real  y_sim;
    y_sim = student_t_rng(nu, mu, sigma);
}
"
IP.sim <- 5   
a1 <- c(3,5, 10,20) #  levels of multiple.mu  
a2 <- c(5,10,20, 50) # levels of multiple.sigma
a.m <- expand.grid(a1, a2)  # a matrix for all 16 prior settings


for (i in 1:nrow(a.m)) {
  datasim.list <- list(
    IP=IP.sim, multiplier.mu=a.m[i,1], multiplier.sigma= a.m[i,2]
  )

  standmodel <- stan_model(model_code=priorpredcode) 
  
  stanfit.sim <- sampling(standmodel,
                          data=datasim.list,
                          algorithm="NUTS",
                          control=list(adapt_delta=0.85),
                          iter=3000,
                          warmup=500,
                          chains=4,
                          cores=2,
                          seed=9182567)
  
  y_sim <- round(as.matrix(stanfit.sim),2)
  ysim <- data.frame(multiplier.mu = rep(a.m[i,1], nrow(y_sim)),multiplier.sigma = rep(a.m[i,2], nrow(y_sim)), y_sim[,"y_sim"])
}


## plot prior predictive distribution

ysim$multiple.mu <- as.factor(ysim$multiple.mu )
ysim$multiple.sigma <- as.factor(ysim$multiple.sigma)
p = ggplot(data= ysim.plot, aes(x=y_sim, color=multiple.sigma)) +
  facet_wrap(~ multiple.mu, nrow=2, labeller=label_both) +
  geom_density()+  xlim(c(0,200))+
  stat_function(fun=dnorm,n=101, args=list(mean=100, sd=25), colour="black", linetype=2,lwd=1, alpha=0.8) +
  xlim(c(0,200))+ylab("")+theme(legend.position="right")
p+scale_color_grey(start = 0.1, end = .8)+  theme_classic()


## Prior sensitivity with simulated examples
set.seed(87621056)
y.all <- rnorm(24, 100, 2*IP.sim)
N.m <- list(3, 6, 12, 24)
y.m <- list(y.all[1:3], y.all[1:6], y.all[1:12], y.all[1:24])

postpredcode = "
  data {
    int<lower=1> n;
    real<lower=0> y[n];
    real IP;
    real  multiplier.mu;
    real  multiplier.sigma;
  }
  transformed data {
    real H;
    real musigma;
    H = multiplier.sigma *IP;
    musigma = multiplier.mu *IP;
  }  
parameters {
    real<lower=1> nu;
    real mu;
    real<lower=0> sigma;
  }
  model {
    nu ~ gamma(2, 0.1);
    mu ~ normal(100, musigma);
    sigma ~ normal(0, H);
     y ~ student_t(nu, mu, sigma);
  }
  generated quantities{
    real  y_pred;
    y_pred = student_t_rng(nu, mu, sigma);
  }
"

            
for (j in 1:4) {   # index for four levels of sample size
  for (i in 1:nrow(a.m)) {  # index for prior settings
    N.sim <- c(N.m[[j]])
    y.sim <- c(y.m[[j]])
    datasim.list <- list(
      n=N.sim, y=y.sim, meanY=100, IP=IP.sim, multiplier.mu = a.m[i,1], multiplier.sigma = a.m[i,2]  )
    standmodel <- stan_model(model_code=postpredcode) 
    stanfit.sim <- sampling(standmodel,
                            data=datasim.list,
                            algorithm="NUTS",
                            control=list(adapt_delta=0.85),
                            iter=1000,
                            warmup=500,
                            chains=4,
                            seed=9182567)
    
    summary.post <- monitor(stanfit.sim, permute=FALSE, inc_warmup=FALSE,digits_summary = 2)
    summary.m <- data.frame(param=rownames(summary.post),size=rep(N.sim,nrow(summary.post)), 
                            mean.y=rep(mean(y.sim), nrow(summary.post)), 
                            sd.y=rep(sd(y.sim), nrow(summary.post)),
                            multiplier.mu = rep(a.mu, nrow(summary.post)),
                            multiplier.sigma = rep(a.y, nrow(summary.post)), summary.post)  
    write.table(summary.m, file="postsummary-potency-mu100.csv", sep = ",", row.names=FALSE, append = T)        
  }
}

##  posterior sample summary statistics
mc.m <-  read.csv("postsummary-potency-mu100.csv", header=TRUE)
mc.m$multiplier.mu <- as.factor(mc.m$multiplier.mu)
mc.m$multiplier.sigma <- as.factor(mc.m$multiplier.sigma)
mc.m$size <- factor(mc.m$size,levels=c("3", "6","12","24"), labels=c("3", "6","12","24")) 
ggplot(mc.m, aes(x=param, y=Rhat)) +
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1)+
  ylim(c(0.995, 1.03))+  xlab("")+  theme_minimal()


library(VCA)
rhat.m <- mc.m %>% filter(param == "y_pred")  # only look at y_pred as an example

varPlot(Rhat~multiplier.mu/multiplier.sigma, rhat.m, 
        VCnam=list(text=c("multiplier.mu", "multiplier.sigma"), cex=0.9, col="black"),
        VarLab = list(cex = 0.9, adj = c(0.5, 0.5)),
        YLabel = list(text ="y_pred: Rhat", side = 2, line =3.5, cex = 1.0),
        mar=c(0.1,6,0.5,1), htab = 0.2,
        Points = list(pch = list(var="size", pch=1:4), cex = 1, col = "black"),
        Mean = NULL,  ylim=c(0.995, 1.02))
legend( x = 0, y = 1.02, pch=1:4, pt.cex=c(1,1,1,1), 
        legend=c("size=3", "size=6", "size=12", "size=24"), cex=0.9)



# A case study
myData0 <- round(c(81.38162,  95.55772,  98.71117,  94.61979, 101.02236,  104.56085, 
                   116.30879, 103.16120,  98.71117 ),2)
process <- c(rep("Process_1", "Process_2"), each=c(6,3))
breaks <- round(c(mean(myData0[1:6]), mean(myData0[1:6])+3*sd(myData0[1:6]),
                  mean(myData0[1:6])-3*sd(myData0[1:6])),1)
example <- data.frame(lot=factor(1:length(myData0)),process=process, Assay=myData0)

example %>%
  ggplot(aes(x=lot, y=Assay, shape=process)) +
  scale_y_continuous(limits = c(70, 130),breaks=breaks )+  geom_point(size=3) + 
  geom_hline(yintercept= breaks[1], lty=2, col="black") +  
  geom_hline(yintercept= breaks[2:3], lty=2, col="darkgray") +
  theme_classic()+  
  theme(legend.position="top", axis.text = element_text(size = 12),axis.title = element_text(size = 12))

IP.sim <- 10  
multiplier.mu <- 2
multiplier.sigma <- 5
nchain <- 4
datasim.list <- list(n=n, y=y, meanY=100, IP=IP.sim, a_mu=a.mu, a_y=a.y)
standmodel <- stan_model(model_code=postpredcode) 
stanfit.sim <- sampling(standmodel,
                        data=datasim.list,
                        algorithm="NUTS",
                        iter=1000, #vary iteration and warmup numbers to change the chain size
                        warmup=500,
                        chains=nchain,
                        cores=2,
                        seed=182567)
y_sim <- round(as.matrix(stanfit.sim),2)
write.csv(y_sim, file="mcsamples_example2000.csv")
summary.post <- round(monitor(stanfit.sim, permute=FALSE, inc_warmup=FALSE),4)
write.csv(summary.post, file=" mcmcdiag_n2000.csv")


library(bayesplot)
lp_cp <- log_posterior(stanfit.sim)
np_cp <- nuts_params(stanfit.sim)

color_scheme_set("gray")
mcmc_trace(stanfit.sim, pars = c("mu","sigma","nu"), 
           facet_args = list(ncol = 1, strip.position = "left"))+  xlab("Post-warmup iteration")
mcmc_acf(stanfit.sim, pars = c("mu","sigma","nu"), lags = 10)
mcmc_pairs((stanfit.sim), pars = c("mu", "sigma", "nu"), np = np_cp, off_diag_args = list(size = 0.75))
mcmc_nuts_divergence(np_cp, lp_cp)

# posterior accuracy 
# load the files with different post buin-in MCMC samples
mcse.n2000 <- read.csv("mcmcdiag_n2000.csv", header=TRUE)
mcse.n6000 <- read.csv("mcmcdiag_n6000.csv", header=TRUE)
mcse.n10000 <- read.csv("mcmcdiag_n10000.csv", header=TRUE)
mcse.n20000 <- read.csv("mcmcdiag_n20000.csv", header=TRUE)
mcse.n30000 <- read.csv("mcmcdiag_n30000.csv", header=TRUE)

mcse.n <- data.frame(size=rep(c(2000, 6000, 10000, 20000, 30000), each=5), 
                     rbind(mcse.n2000,mcse.n6000, mcse.n10000,mcse.n20000, mcse.n30000)) %>%
  filter(X %in% c( "mu", "sigma", "nu", "y_pred"))

mcse.n %>%
  filter(X=="mu") %>%
  ggplot(aes(x=factor(size), y=Q50))+
  geom_point()+  geom_line()+
  geom_errorbar(aes(ymin=Q50-MCSE_Q50,ymax=Q50+MCSE_Q50), width=.3,position="dodge")+
  ylab("mean")+  xlab("iterations")+  theme_bw()

mcse.n %>%
  ggplot(aes(x=factor(size), y=MCSE_Q50/sd, shape=X, group=X))+
  geom_point()+  geom_line()+   ylab("MCSE_Q50/sd")+
  xlab("iterations")+  ylim(c(0,0.05))+  theme_classic2()


## prior to posterior comparison
y_pred <-  read.csv("mcsamples_example20000.csv", header=TRUE)
y_sim <-  read.csv("priorpredmcsamples_example20000.csv", header=TRUE)
vardistr <- data.frame(bayes=rep(c("posterior","prior"), c(length(y_pred[,"mu"]),length(y_sim[,"mu"]))), 
                       mu=c(y_pred[,"mu"],y_sim[,"mu"]), sigma=c(y_pred[,"sigma"],y_sim[,"sigma"]),  nu=c(y_pred[,"nu"],y_sim[,"nu"]), 
                       y=c(y_pred[,"y_pred"],y_sim[,"y_pred"]))

library(ggh4x)
scales <- list(
  scale_x_continuous(minor_breaks = c(80,90, 100, 110,120),  limits=c(70, 130)),
  scale_x_continuous(minor_breaks = c(10, 30, 50, 80),    limits=c(0,100)),
  scale_x_continuous(minor_breaks = c(10, 30, 50, 80),  limits=c(0,100)),
  scale_x_continuous(minor_breaks = c(70, 80,90, 100, 110,120, 130), limits=c(50, 150))
)
vardistr %>%
  pivot_longer(cols=c("mu", "sigma", "nu", "y"),
               names_to="parameters",
               values_to="values") %>%
  dplyr::select("parameters", "values", "bayes") -> plot11

p <- ggplot(data=plot11, aes(x = values, fill=bayes)) +
  geom_density(alpha=0.3) +  xlab("") +
  facet_wrap( ~ factor(parameters), scales = "free", nrow = 2)+
  facetted_pos_scales(x = scales) +
  geom_vline(data=group.mean, aes(xintercept = values, color=bayes), linetype="dashed")+
  theme_classic()+   theme(legend.position="bottom",legend.title=element_blank())+
  theme(strip.text.x=element_text(size=16, colour="black"),
        strip.background = element_rect(fill="white", colour="gray",size=2))
p+scale_color_grey(start = 0.1, end = .5)+ scale_fill_grey(start = 0.1, end = .5)+ theme_classic() + 
  theme(legend.position="bottom",legend.title=element_blank())


bayesboot::plotPost(as.numeric(na.omit(y_sim[,"y_pred "])), col=gray, credMass = 0.95, HDItextPlace =0.95,  xlab="values")






