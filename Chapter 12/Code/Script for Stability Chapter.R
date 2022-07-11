
#' ---
#' title: "For Book Chapter on Bayesian Stability"
#' author: "Ji Young Kim"
#' date: "`r format(Sys.time(), '%d-%b-%Y')`"
#' output:
#'  html_document:
#'    css: style.css
#' ---

options(repos = 'https://cran.rstudio.com/') 
mute <- function(x)suppressWarnings(suppressPackageStartupMessages(x)) 

# loading libraries  ---------------------------------------------------
mute(library(tidyverse))
mute(library(plyr))
mute(library(knitr))
mute(library(pander))
mute(library(readxl))
mute(library(lme4))
mute(library(lmerTest))
mute(library(car))
mute(library(latex2exp))
mute(library(patchwork))
mute(library(gt))
mute(library(emmeans))
mute(library(writexl))
mute(library(MASS))
mute(library(ggsci))
mute(library(RColorBrewer))
mute(library(rstan))
# mute(library(rstanarm))
mute(library(HDInterval))
mute(library(merTools))
# mute(library(brms))
mute(library(merDeriv))
mute(library(pander))
mute(library(geomtextpath))

select <- dplyr::select
rename <- dplyr::rename
fill <- tidyr::fill

R.Version()$version.string  


# create folders
dir.create("./Plots", showWarnings = F)
dir.create("./Tables", showWarnings = F)

# To fix the results for replication

set.seed(123)
seed <- sample.int(.Machine$integer.max, 1)
seed

d <- read.csv("batches.csv") %>%
  rename(Time = x, Y = y) %>%
  mutate(Batch = as.character(Batch)) %>%
  arrange(Batch, Time) %>%
  mutate(Group = case_when(Batch %in% 1:3 ~ "Historical", T ~ "FSS"))

d # Batches 1, 2, 3 are historical data to be used to derive priors

# Spec Limits
LSL <- 90
USL <- 110


# Plot the data
p1 <- ggplot(data = d, aes(x = Time, y = Y)) +
  theme_bw() +
  scale_shape_manual(values = 1:length(unique(d$Batch))) +
  scale_colour_grey(start = 0.1, end = 0.7, na.value = "black", aesthetics = "colour") +
  # scale_color_d3()+
  # scale_color_manual(values = brewer.pal(n = 11, name = "BrBG")) +
  geom_hline(yintercept = LSL, col = 1, linetype = 1) +
  geom_hline(yintercept = USL, col = 1, linetype = 1) +
  geom_point(aes(color = Batch, shape = Batch)) +
  geom_smooth(aes(color = Batch, linetype = Group), method = "lm", formula = y ~ x, se = F) +
  scale_x_continuous(breaks = unique(d$Time)) +
  guides(linetype = "none") +
  labs(x = "Time", y = "Response", title = "Example Data") # + theme(legend.position = "bottom")

p1
ggsave("./Plots/p1.png", p1, height = 4, width = 5, dpi=1200)


# Priors based on info from batches 1~3

lmr <- lmer(Y ~ Time + (1 |Batch), data = d %>% filter(Batch %in% 1:3))
summary(lmr)
ranova(lmr)
confint(lmr)

# Estimates of random effects variances
vars_re <- data.frame(summary(lmr)$varcor); vars_re
sigma2_b <- vars_re$vcov[1]; sigma2_b
sigma2_r <- vars_re$vcov[2]; sigma2_r

# Standard error of random effects variances
se_vars <- vcov(lmr, full = TRUE); se_vars

write_xlsx(data.frame(as.matrix(se_vars)), "./Tables/se_vars.xlsx")

se.sigma2_b <- sqrt(diag(se_vars)[3]); se.sigma2_b
se.sigma2_r <- sqrt(diag(se_vars)[4]); se.sigma2_r

df_b <- 2 * (sigma2_b / se.sigma2_b)^2; df_b
df_r <- 2 * (sigma2_r / se.sigma2_r)^2; df_r

# Based on these, priors for sigma2_b and sigma2_r are

# sigma2_b ~ IG(df_b/2, df_b * sigma2_b^2 / 2)
# sigma2_r ~ IG(df_r/2, df_r * sigma2_r^2 / 2)


df_b/2 # 0.75
df_b * sigma2_b^2 / 2 # 0.95

df_r/2 # 8.50
df_r * sigma2_r^2 / 2 # 12.90

d1 <- d %>% filter(as.numeric(Batch) > 3) # %>% droplevels()
d1

###################################
### ICH approach for shelf life ###
###################################

fit.f1 <- lm(Y ~ Time*Batch, data=d1)
Anova(fit.f1) # P-value for interaction greater than 0.25

fit.f2 <- lm(Y ~ Time + Batch, data=d1)
Anova(fit.f2) # P-value for Batch smaller than 0.25
summary(fit.f2)

# final model

fit.f <- fit.f2 

unique(d1 %>% select(Batch, Time)) # max extrapolated = 24


Newdata = expand.grid(Batch = unique(d1$Batch), Time = seq(0,40, by = 0.01))
Pred <- predict(fit.f, Newdata, interval = "confidence")
Predicted <- data.frame(Newdata, Pred) %>% mutate(within.spec = lwr >= LSL & upr <= USL)

SL.f = Predicted %>% filter(within.spec == F) %>% select(Time) %>% min()
SL.f

p2 <- ggplot(data = Predicted, aes(x = Time, y = fit)) +
  theme_bw() +
  scale_shape_manual(values = 1:length(unique(d1$Batch))) +
  scale_colour_grey(start = 0.1, end = 0.4, na.value = "black", aesthetics = "colour") +
  scale_colour_grey(start = 0.1, end = 0.4, na.value = "black", aesthetics = "fill") +
  scale_linetype_manual(values = c("dotdash", "dashed", "longdash")) +

  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Batch), alpha = 0.3) +
  geom_line(aes(color=Batch, linetype = Batch)) +

  geom_hline(yintercept = LSL, col = 1, linetype = 1) +
  geom_hline(yintercept = USL, col = 1, linetype = 1) +

  scale_x_continuous(breaks = seq(0, 40, by=5)) +
  labs(x = "Time", y = "Response", title = "95% Confidence Bands by Batch") # + theme(legend.position = "bottom")

p2

ggsave("./Plots/p2.png", p2, height = 4, width = 5, dpi=1200)


###########################
## Random Batch Approach ##
###########################


standata = list(

  N = nrow(d1),
  N_Batch = length(unique(d1$Batch)),
  Batch = as.numeric(as.factor(d1$Batch)),
  Y = d1$Y,
  Time = d1$Time
)


stanmodel_no_priors ="

  data{

    int N ;
    int N_Batch;
    vector[N] Y;
    vector[N] Time;
    int Batch[N];
  }

  parameters{

  real alpha;
  real beta;
  real<lower=0> sigma2_batch_int;
  real<lower=0> sigma2_resid;

  // generated values for random effects

  vector[N_Batch] u_batch_int; // random intercept by batch
}

  model {

  vector[N] mu; // mean vector

// noninformative Priors

    alpha ~ normal(0, 1e3);
    beta ~ normal(0, 1e3);
    
    sigma2_resid ~ inv_gamma(0.0001, 0.0001);
    sigma2_batch_int ~ inv_gamma(0.0001, 0.0001);


// Random effects distributions

u_batch_int ~ normal(0, sqrt(sigma2_batch_int));

// Mean vector

for (i in 1:N){

mu[i] = alpha +  u_batch_int[Batch[i]] + beta * Time[i] ;

}

// Distribution of Y

Y ~ normal(mu, sqrt(sigma2_resid)); // vector of length N

}

"


stanmodel_info_priors ="

  data{

    int N ;
    int N_Batch;
    vector[N] Y;
    vector[N] Time;
    int Batch[N];
  }

  parameters{

  real alpha;
  real beta;
  real<lower=0> sigma2_batch_int;
  real<lower=0> sigma2_resid;

  // generated values for random effects

  vector[N_Batch] u_batch_int; // random intercept by batch
}

  model {

  vector[N] mu; // mean vector

// Informative Priors on variances

    alpha ~ normal(0, 1e3);
    beta ~ normal(-0.22134, 0.03093);

    sigma2_batch_int ~ inv_gamma(0.75, 0.95);
    sigma2_resid ~ inv_gamma(8.50, 12.90);


// Random effects distributions

u_batch_int ~ normal(0, sqrt(sigma2_batch_int));

// Mean vector

for (i in 1:N){

mu[i] = alpha +  u_batch_int[Batch[i]] + beta * Time[i] ;

}

// Distribution of Y

Y ~ normal(mu, sqrt(sigma2_resid)); // vector of length N

}

"


fit_no_priors <- stan(model_code = stanmodel_no_priors,
                      data = standata,
                      chains = 4,
                      iter   = 130000,
                      warmup = 30000,
                      cores = 4,
                      thin=10,
                      seed = 5555)

fit_info_priors <- stan(model_code = stanmodel_info_priors,
                        data = standata,
                        chains = 4,
                        iter   = 130000,
                        warmup = 30000,
                        cores = 4,
                        thin=10,
                        seed = 5555)

print(fit_no_priors)
print(fit_info_priors)

stanfits <- list(fit_no_priors, fit_info_priors)
label <- c("Noninformative prior", "Informative Priors")

set.seed(seed)

for (j in 1:2){
  
  print(label[j])
  fit <- stanfits[[j]]
  cat("\n")
  
  print(fit)
  print(stan_ac(fit))
  
  # Traceplot in gray scale for paper
  print(stan_trace(fit)  + scale_colour_grey(start = 0.7, end = 0.3, na.value = "black", aesthetics = "colour"))
    
  post <- data.frame(extract(fit))
  post = post[order(post[,1]),]
  post = post[sample(nrow(post)), ]
  
  print(kable(head(post, 10), digits=3))

  write_xlsx(post,  paste0("./Tables/posterior draws ", j, ".xlsx"))

  # Parameters that were drawn
  colnames(post)

  # Intercept for a future batch
  int.fut <- rnorm(nrow(post), post$alpha, sqrt(post$sigma2_batch_int))

  # Fixed slope
  slope.fut <- post$beta

  # Draw histograms
  p.density.int <- ggplot(data=data.frame(int.fut), aes(x=int.fut)) +
    theme_bw() + xlim(90, 120) +
    labs(x = "Intercept of a Future Batch") +
    geom_density()

  p.density.slope <- ggplot(data=data.frame(slope.fut), aes(x=slope.fut)) +
    theme_bw() + xlim(-0.75, 0.25) +
    labs(x = "Slope of a Future Batch") +
    geom_density()

  print(p.density.int + p.density.slope)
  
  ggsave(paste0("./Plots/int slope ", j, ".png"), height = 4, width = 10, dpi=1200)
  

  # Posterior predictive distribution of the Shelf life a future batch
  SLs <- data.frame(int.fut, slope.fut) %>%
    mutate(sl = case_when(slope.fut < 0 ~ (LSL - int.fut)/slope.fut,
                          slope.fut > 0 ~ (USL - int.fut)/slope.fut,
                          slope.fut == 0 ~ Inf))


  # Burdick's book suggests to take 5th percentile as the shelf life
  SL.fut.batch <- data.frame(Batch = "Future",
                             SL = quantile(SLs$sl, 0.05, na.rm=T))
  print(SL.fut.batch)
  floor(SL.fut.batch$SL)

  # Distributions of raw and trimmed shelf life estimates
  print(summary(SLs$sl))
  
  p.sl <- ggplot(data = SLs, aes(x=sl)) +
    geom_histogram(binwidth = 2500) +
    labs(x = "Shelf Life", title = "Shelf life distribution(Entire posterior distribution)")+
    theme_bw(); p.sl


  p.sl.tr <- ggplot(data = SLs %>% filter(sl < 250 & sl > 0), aes(x=sl)) +
    geom_histogram(binwidth = 2) +
    geom_vline(xintercept = round(SL.fut.batch$SL), col=2, linetype = "dashed") +
    # annotate("text", x=230, y = -100, label="Red line at the 5th percentile", size=3.5, col = 2) +
    labs(x = "Shelf Life", title = "Posterior distribution of shelf life between 0 and 250 months", y = "Response") +
    theme_bw(); p.sl.tr


  print(p.sl.tr)

  ggsave(paste0("./Plots/p.sl.tr ", j, ".png"), p.sl.tr, height = 4, width = 6, dpi=1200)
  


  # IRL

  # Generate K release values within spec range
  release <- seq(LSL, USL, length=10000)

  # The probability that the value at the shelf-life is within spec

  prob_withinspec <- do.call(rbind, lapply(release, function(k) {
    
    # k is each release value
    Y_R <- rnorm(nrow(post), k, sqrt(post$sigma2_resid))
    Y_SR <- Y_R + slope.fut * 36 # Shelf life set at 36 months
    prob_k <- mean(Y_SR >= LSL & Y_SR <= USL & Y_R >= LSL & Y_R <= USL)

    data.frame(k, prob_k)

})) # %>% arrange(prob_k)

  p.irl <- ggplot(data=prob_withinspec, aes(x=prob_k, y=k)) +
    #geom_line() +
    geom_point(size=0.2, alpha = 0.3, color = "gray")+
    
    theme_bw() +
    # labs(x=paste("P(Within Spec at the Shelf Life of", floor(SL.fut.batch$SL), "Months)" ), y = "Value at Release") +
    labs(x=paste("P(Within Spec at the Shelf Life of", 36, "Months)" ), y = "Value at Release") +

    geom_textvline(xintercept = 0.95, label = "Success Prob = 95%", hjust = 0.5, angle = 90,
                   linetype = 2, vjust = 0.5, color = "1") +
    geom_texthline(yintercept = LSL, label = "LSL", hjust = 0.5, angle = 90,
                   linetype = 1, vjust = 0.5, color = "1") +
    geom_texthline(yintercept = USL, label = "USL", hjust = 0.5, angle = 90,
                   linetype = 1, vjust = 0.5, color = "1")

  
  print(p.irl)
  
  ggsave(paste0("./Plots/p.irl ", j, ".png"), p.irl, height = 4, width = 5, dpi=1200)

  IRL <- prob_withinspec %>% filter(prob_k >= 0.95) %>% select(k) %>% range()
  
  print("Results")

  print(data.frame(J=j, SL.fut.batch, LRL = IRL[1], ULR = IRL[2]))

}


# Allen's method:

n <- 1 # number of replicate assays used for batch release
lm_d1 <- lm(Y ~ Batch + Time, data = d1)
summary(lm_d1)

assayvar <- sigma(lm_d1)^2
avg.slope <- summary(lm_d1)$coef["Time", "Estimate"]; avg.slope
se.slope <- summary(lm_d1)$coef["Time", "Std. Error"];se.slope

df.assayvar <- lm_d1$df.residual
df.avg.slope <- df.assayvar # assumed the same as residual df

df.t <- (assayvar + se.slope^2)^2 / ( (assayvar^2/ df.assayvar)+(se.slope^2^2 / df.avg.slope) ); df.t
tvalue <- qt(0.95, df.t); tvalue

# with SL of 36
LSL - avg.slope * 36 + tvalue * sqrt(se.slope^2 * 36^2 + (assayvar/n))
