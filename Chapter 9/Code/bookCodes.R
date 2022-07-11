library(tidyverse)
library(readxl)
library(brms)
library(shinystan)
library(rstan)
library(tidybayes)
library(mvtnorm)

### Load data
dataCS <- read_excel("bookData.xlsx") %>%
  mutate(runK = paste0(k, ".",run))

### Case study 1: only single dissolution time point

## Bayesian model fit

fitCS1 <- brm(mvbind(md30, NIRc) ~ factor1*factor2 + (1|q|runK), data = dataCS,
              iter = 15000, warmup = 1000, chains = 3, thin = 5,
              seed = 1234, cores = 3, 
              control = list(adapt_delta = 0.95, max_treedepth = 10)) 

launch_shinystan(fitCS1)


## Model results
summary(fitCS1)

chainsCS1 <- posterior_samples(fitCS1)[,1:14]

# Table 1
table1 <- matrix(c(round(colMeans(chainsCS1), 2)[c(1,3,4,5,2,6:14)],
                   round(apply(chainsCS1, 2, sd),2)[c(1,3,4,5,2,6:14)],
                   apply(chainsCS1, 2, function(x) paste0("(", 
                                                          paste0(round(quantile(x, c(0.025, 0.975)), 2), 
                                                                 collapse = ", "), ")"))[c(1,3,4,5,2,6:14)]), 
                 ncol = 3, byrow = FALSE)


## Conditional model
# mu_Q30
round(mean(chainsCS1[,"b_md30_Intercept"]), 2)
# theta
round(mean(chainsCS1[,"cor_runK__md30_Intercept__NIRc_Intercept"]*
             chainsCS1[,"sd_runK__md30_Intercept"]/chainsCS1[,"sd_runK__NIRc_Intercept"]), 2)
#mu_NIR
round(mean(chainsCS1[,"b_NIRc_Intercept"]), 2)

# model itself
predictCS1target <- function(NIRvalue){
  round(mean(chainsCS1[,"b_md30_Intercept"]), 2)  + 
    round(mean(chainsCS1[,"cor_runK__md30_Intercept__NIRc_Intercept"]*
                 chainsCS1[,"sd_runK__md30_Intercept"]/chainsCS1[,"sd_runK__NIRc_Intercept"]), 2)*
                (NIRvalue - round(mean(chainsCS1[,"b_NIRc_Intercept"]), 2))
}

## Simulation of performance

resultsMean <- matrix(numeric(nrow(chainsCS1)*2), ncol = 2)
colnames(resultsMean) <- c("BatchTrueQ30", "MeanPredQ30ind")

set.seed(12345)
for (i in 1:nrow(chainsCS1)){
  sampleBatch <- matrix(rep(
    rmvnorm(1, 
            mean = unlist(chainsCS1[i, c("b_md30_Intercept", "b_NIRc_Intercept")]), 
            sigma = matrix( c(
              chainsCS1[i, "sd_runK__md30_Intercept"]^2,
              rep(chainsCS1[i, "sd_runK__md30_Intercept"]*
                    chainsCS1[i, "sd_runK__NIRc_Intercept"]*
                    chainsCS1[i, "cor_runK__md30_Intercept__NIRc_Intercept"], 2),
              chainsCS1[i, "sd_runK__NIRc_Intercept"]^2
            ), byrow = TRUE, nrow=2)), 30), byrow = TRUE, nrow = 30) 
  sampleError <- rmvnorm(30, mean = c(0,0), 
                         sigma = matrix( c(
                           chainsCS1[i, "sigma_md30"]^2,
                           rep(chainsCS1[i, "sigma_md30"]*chainsCS1[i, "sigma_NIRc"]*chainsCS1[i, "rescor__md30__NIRc"], 2),
                           chainsCS1[i, "sigma_NIRc"]^2
                         ), byrow = TRUE, nrow=2))
  sample <- sampleBatch + sampleError
  dissoPred <- predictCS1target(sample[,2])
  resultsMean[i,] <- c(colMeans(sampleBatch)[1], mean(dissoPred)) # last two are same for linear predictor, but not if curvature
}

resultsMean <- as.data.frame(resultsMean)

# Figure 1 
resultsMean %>%
  select(BatchTrueQ30, MeanPredQ30ind) %>%
  pivot_longer(cols = BatchTrueQ30:MeanPredQ30ind,  names_to = "Variable", values_to = "Q30" ) %>%
  ggplot(aes(x = Q30, fill = Variable)) +
  geom_density(alpha = 0.4) + 
  scale_fill_manual(labels = c("Batch mean", "Predicted value"), values = c("dodgerblue", "midnightblue")) +
  labs(x = "Q30",
       title = "Dissolution at 30 minutes") + theme_bw()

# Table 2
round(apply(resultsMean, 2, median) , 2)
apply(resultsMean, 2, function(x) 
  paste0("(", paste0(round(quantile(x, c(0.025, 0.975)),2), collapse = ", "), ")"))
apply(resultsMean, 2, function(x) 
  paste0("(", paste0(round(quantile(x, c(0.005, 0.995)),2), collapse = ", "), ")"))

# Figure 3
resultsMeanDiff <- resultsMean %>%
  mutate(Difference = BatchTrueQ30 - MeanPredQ30ind)

resultsMeanDiff %>%
  ggplot(aes(x = Difference)) +
  geom_density(alpha = 0.4, fill="lightskyblue", col = "royalblue") + 
  labs(x = "Difference in Q30",
       title = "Difference in prediction and batch mean at 30 minutes") + theme_bw()

# Table 3 
round(median(resultsMeanDiff$Difference),2)
paste0("(", paste0(round(quantile(resultsMeanDiff$Difference, c(0.025, 0.975)),2), 
                   collapse = ", "), ")")
paste0("(", paste0(round(quantile(resultsMeanDiff$Difference, c(0.005, 0.995)),2), 
                   collapse = ", "), ")")





### Case study 2: full profile prediction

## Transform data in long format
dataCSLong <- gather(data = dataCS, DissoTime, Value, md05:md60)
dataCSLong$DissoTime <- as.numeric(unlist(str_extract_all(dataCSLong$DissoTime,"\\(?[0-9,.]+\\)?")))
dataCSLong$Value <- as.numeric(as.character(dataCSLong$Value))

## Run Weibull model per tablet
set.seed(1234)
startingValue <- list(theta1 = 100, theta2 = 7, theta3 = 0.8) # set starting values
coefsNew <- NULL

for(i in 1:nrow(dataCS)){
  # get data from a single tablet
  tablet_i <- dataCSLong[which(dataCSLong$tabx == i),]
  
  # fit Weibull
  fixedWeibullTablet_i <- nls(Value ~ theta1*(1 - exp(-exp(theta3 * (log(DissoTime) - log(theta2))))),
                              data = tablet_i,
                              start = startingValue)
  coefsNew <- rbind(coefsNew, coef(fixedWeibullTablet_i))
} 
coefsNew <- as.data.frame(coefsNew)

# save the Weibull parameters' estimates (to original data set, not long format one)
dataCS$WUB <- coefsNew$theta1
dataCS$Wlambda <- coefsNew$theta2
dataCS$Wk <- coefsNew$theta3


## Fit Bayesian model
fitCS2 <- brm(mvbind(WUB, Wlambda, Wk, NIRc) ~ factor1*factor2 + (1|q|runK), data = dataCS,
              iter= 30000, warmup = 10000, chains = 3, thin = 5,
              seed = 1234, cores = 3, 
              control = list(adapt_delta = 0.95, max_treedepth = 10)) 
launch_shinystan(fitCS2)


## Model results
summary(fitCS2)
chainsCS2 <- posterior_samples(fitCS2)[,1:36]
# Table 4
table4 <- matrix(c(round(colMeans(chainsCS2), 2)[c(1,5,6,7,2,8,9,10,3,11,12,13,4,14,15,16,17:20, 24:30, 34:36)],
                   round(apply(chainsCS2, 2, sd),2)[c(1,5,6,7,2,8,9,10,3,11,12,13,4,14,15,16,17:20, 24:30, 34:36)],
                   apply(chainsCS2, 2, function(x) paste0("(", paste0(round(quantile(x, c(0.025, 0.975)), 2), 
                                                                      collapse = ", "), ")"))[c(1,5,6,7,2,8,9,10,3,11,12,13,4,14,15,16,17:20, 24:30, 34:36)]), ncol = 3, byrow = FALSE)

## Conditional model
# mu_Weibull
round(mean(chainsCS2[,"b_WUB_Intercept"]), 2)
round(mean(chainsCS2[,"b_Wlambda_Intercept"]), 2)
round(mean(chainsCS2[,"b_Wk_Intercept"]), 2)
# theta
round(mean(chainsCS2[,"cor_runK__WUB_Intercept__NIRc_Intercept"]*
             chainsCS2[,"sd_runK__WUB_Intercept"]/chainsCS2[,"sd_runK__NIRc_Intercept"]), 2)
round(mean(chainsCS2[,"cor_runK__Wlambda_Intercept__NIRc_Intercept"]*
             chainsCS2[,"sd_runK__Wlambda_Intercept"]/
             chainsCS2[,"sd_runK__NIRc_Intercept"]), 2)
round(mean(chainsCS2[,"cor_runK__Wk_Intercept__NIRc_Intercept"]*
             chainsCS2[,"sd_runK__Wk_Intercept"]/chainsCS2[,"sd_runK__NIRc_Intercept"]), 2)
# mu_NIR
round(mean(chainsCS2[,"b_NIRc_Intercept"]), 2)

## Weibull function
weibullF <- function(wub, wl, wk, t){
  wub*(1-exp(-exp( wk*(log(t)-log(wl)))))
}

## Conditional model itself
predictCS2target <- function(NIRvalue){
  wub = round(mean(chainsCS2[,"b_WUB_Intercept"]), 2) + 
    round(mean(chainsCS2[,"cor_runK__WUB_Intercept__NIRc_Intercept"]*
                 chainsCS2[,"sd_runK__WUB_Intercept"]/chainsCS2[,"sd_runK__NIRc_Intercept"]), 2)*
                (NIRvalue - round(mean(chainsCS2[,"b_NIRc_Intercept"]), 2))
  wl = round(mean(chainsCS2[,"b_Wlambda_Intercept"]), 2) + 
    round(mean(chainsCS2[,"cor_runK__Wlambda_Intercept__NIRc_Intercept"]*
                 chainsCS2[,"sd_runK__Wlambda_Intercept"]/chainsCS2[,"sd_runK__NIRc_Intercept"]), 2)*
                (NIRvalue - round(mean(chainsCS2[,"b_NIRc_Intercept"]), 2))
  wk = round(mean(chainsCS2[,"b_Wk_Intercept"]), 2) + 
    round(mean(chainsCS2[,"cor_runK__Wk_Intercept__NIRc_Intercept"]*
                 chainsCS2[,"sd_runK__Wk_Intercept"]/chainsCS2[,"sd_runK__NIRc_Intercept"]), 2)*
                (NIRvalue - round(mean(chainsCS2[,"b_NIRc_Intercept"]), 2))
  
  t30pred <- weibullF(wub, wl, wk, 30) # we compare at dissolution 30 minutes
  
  return(t30pred)
}


## Simulation of performance
resultsMeanCS2 <- matrix(numeric(nrow(chainsCS2)*2), ncol = 2)
colnames(resultsMeanCS2) <- c("BatchTrueQ30", "MeanPredQ30ind")

set.seed(12345)
for (i in 1:nrow(chainsCS2)){
  sampleBatch <- matrix(rep(rmvnorm(1, 
                                    mean = unlist(chainsCS2[i, c("b_WUB_Intercept", "b_Wlambda_Intercept", 
                                                                 "b_Wk_Intercept", "b_NIRc_Intercept")]), 
                                    sigma = matrix( c(
                                      chainsCS2[i, "sd_runK__WUB_Intercept"]^2,
                                      chainsCS2[i, "sd_runK__WUB_Intercept"]*
                                        chainsCS2[i, "sd_runK__Wlambda_Intercept"]*
                                        chainsCS2[i, "cor_runK__WUB_Intercept__Wlambda_Intercept"],
                                      chainsCS2[i, "sd_runK__WUB_Intercept"]*
                                        chainsCS2[i, "sd_runK__Wk_Intercept"]*
                                        chainsCS2[i, "cor_runK__WUB_Intercept__Wk_Intercept"],
                                      chainsCS2[i, "sd_runK__WUB_Intercept"]*
                                        chainsCS2[i, "sd_runK__NIRc_Intercept"]*
                                        chainsCS2[i, "cor_runK__WUB_Intercept__NIRc_Intercept"],
                                      chainsCS2[i, "sd_runK__WUB_Intercept"]*
                                        chainsCS2[i, "sd_runK__Wlambda_Intercept"]*
                                        chainsCS2[i, "cor_runK__WUB_Intercept__Wlambda_Intercept"],
                                      chainsCS2[i, "sd_runK__Wlambda_Intercept"]^2,
                                      chainsCS2[i, "sd_runK__Wlambda_Intercept"]*
                                        chainsCS2[i, "sd_runK__Wk_Intercept"]*
                                        chainsCS2[i, "cor_runK__Wlambda_Intercept__Wk_Intercept"],
                                      chainsCS2[i, "sd_runK__Wlambda_Intercept"]*
                                        chainsCS2[i, "sd_runK__NIRc_Intercept"]*
                                        chainsCS2[i, "cor_runK__Wlambda_Intercept__NIRc_Intercept"],
                                      chainsCS2[i, "sd_runK__WUB_Intercept"]*
                                        chainsCS2[i, "sd_runK__Wk_Intercept"]*
                                        chainsCS2[i, "cor_runK__WUB_Intercept__Wk_Intercept"],
                                      chainsCS2[i, "sd_runK__Wlambda_Intercept"]*
                                        chainsCS2[i, "sd_runK__Wk_Intercept"]*
                                        chainsCS2[i, "cor_runK__Wlambda_Intercept__Wk_Intercept"],
                                      chainsCS2[i, "sd_runK__Wk_Intercept"]^2,
                                      chainsCS2[i, "sd_runK__Wk_Intercept"]*
                                        chainsCS2[i, "sd_runK__NIRc_Intercept"]*
                                        chainsCS2[i, "cor_runK__Wk_Intercept__NIRc_Intercept"],    
                                      chainsCS2[i, "sd_runK__WUB_Intercept"]*
                                        chainsCS2[i, "sd_runK__NIRc_Intercept"]*
                                        chainsCS2[i, "cor_runK__WUB_Intercept__NIRc_Intercept"],
                                      chainsCS2[i, "sd_runK__Wlambda_Intercept"]*
                                        chainsCS2[i, "sd_runK__NIRc_Intercept"]*
                                        chainsCS2[i, "cor_runK__Wlambda_Intercept__NIRc_Intercept"],
                                      chainsCS2[i, "sd_runK__Wk_Intercept"]*
                                        chainsCS2[i, "sd_runK__NIRc_Intercept"]*
                                        chainsCS2[i, "cor_runK__Wk_Intercept__NIRc_Intercept"],  
                                      chainsCS2[i, "sd_runK__NIRc_Intercept"]^2
                                    ), byrow = TRUE, nrow=4)), 30), byrow = TRUE, nrow = 30) 
  sampleError <- rmvnorm(30, mean = c(0,0,0,0), 
                         sigma = matrix( c(
                           chainsCS2[i, "sigma_WUB"]^2,
                           chainsCS2[i, "sigma_WUB"]*
                             chainsCS2[i, "sigma_Wlambda"]*
                             chainsCS2[i, "rescor__WUB__Wlambda"],
                           chainsCS2[i, "sigma_WUB"]*
                             chainsCS2[i, "sigma_Wk"]*
                             chainsCS2[i, "rescor__WUB__Wk"],
                           chainsCS2[i, "sigma_WUB"]*
                             chainsCS2[i, "sigma_NIRc"]*
                             chainsCS2[i, "rescor__WUB__NIRc"],
                           chainsCS2[i, "sigma_WUB"]*
                             chainsCS2[i, "sigma_Wlambda"]*
                             chainsCS2[i, "rescor__WUB__Wlambda"],
                           chainsCS2[i, "sigma_Wlambda"]^2,
                           chainsCS2[i, "sigma_Wlambda"]*
                             chainsCS2[i, "sigma_Wk"]*
                             chainsCS2[i, "rescor__Wlambda__Wk"],
                           chainsCS2[i, "sigma_Wlambda"]*
                             chainsCS2[i, "sigma_NIRc"]*
                             chainsCS2[i, "rescor__Wlambda__NIRc"],
                           chainsCS2[i, "sigma_WUB"]*
                             chainsCS2[i, "sigma_Wk"]*
                             chainsCS2[i, "rescor__WUB__Wk"],
                           chainsCS2[i, "sigma_Wlambda"]*
                             chainsCS2[i, "sigma_Wk"]*
                             chainsCS2[i, "rescor__Wlambda__Wk"],
                           chainsCS2[i, "sigma_Wk"]^2,
                           chainsCS2[i, "sigma_Wk"]*
                             chainsCS2[i, "sigma_NIRc"]*
                             chainsCS2[i, "rescor__Wk__NIRc"],  
                           chainsCS2[i, "sigma_WUB"]*
                             chainsCS2[i, "sigma_NIRc"]*
                             chainsCS2[i, "rescor__WUB__NIRc"],
                           chainsCS2[i, "sigma_Wlambda"]*
                             chainsCS2[i, "sigma_NIRc"]*
                             chainsCS2[i, "rescor__Wlambda__NIRc"],
                           chainsCS2[i, "sigma_Wk"]*
                             chainsCS2[i, "sigma_NIRc"]*
                             chainsCS2[i, "rescor__Wk__NIRc"],  
                           chainsCS2[i, "sigma_NIRc"]^2
                         ), byrow = TRUE, nrow=4))
  sample <- sampleBatch + sampleError
  dissoPred <- predictCS2target(sample[,4])
  resultsMeanCS2[i,] <- c(weibullF(sampleBatch[1,1], sampleBatch[1,2], sampleBatch[1,3], 30), 
                          mean(dissoPred)) # 
}

resultsMeanCS2 <- as.data.frame(resultsMeanCS2)

# Figure 3 
resultsMeanCS2 %>%
  select(BatchTrueQ30, MeanPredQ30ind) %>%
  pivot_longer(cols = BatchTrueQ30:MeanPredQ30ind,  names_to = "Variable", values_to = "Q30" ) %>%
  ggplot(aes(x = Q30, fill = Variable)) +
  geom_density(alpha = 0.4) + 
  scale_fill_manual(labels = c("Batch mean", "Predicted value"), values = c("dodgerblue", "midnightblue")) +
  labs(x = "Q30",
       title = "Dissolution at 30 minutes") + theme_bw()

# Table 4
round(apply(resultsMeanCS2, 2, median) , 2)
apply(resultsMeanCS2, 2, function(x) paste0("(", paste0(round(quantile(x, c(0.025, 0.975)),2), collapse = ", "), ")"))
apply(resultsMeanCS2, 2, function(x) paste0("(", paste0(round(quantile(x, c(0.005, 0.995)),2), collapse = ", "), ")"))

# Figure 4
resultsMeanCS2Diff <- resultsMeanCS2 %>%
  mutate(Difference = BatchTrueQ30 - MeanPredQ30ind)

resultsMeanCS2Diff %>%
  ggplot(aes(x = Difference)) +
  geom_density(alpha = 0.4, fill="lightskyblue", col = "royalblue") + 
  labs(x = "Difference in Q30",
       title = "Difference in prediction and batch mean at 30 minutes") + theme_bw()

# Table 3 
round(median(resultsMeanCS2Diff$Difference),2)
paste0("(", paste0(round(quantile(resultsMeanCS2Diff$Difference, c(0.025, 0.975)),2), collapse = ", "), ")")
paste0("(", paste0(round(quantile(resultsMeanCS2Diff$Difference, c(0.005, 0.995)),2), collapse = ", "), ")")