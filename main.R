#main model：GFI
library(rstan)

NGF_model <- stan_model(file = 'NGF.stan')
GFI_model <- stan_model(file = 'GF2.stan')

# an example：
source('TPWD.R')
source('GenData.R')

set.seed(2414214)

Ntest <- c(25,20,15,10)
Rtest <- c(20,16,12,8)
Stest <- c(40,60,80,100)

lambdatest <- 2
rhotest <- 5
deltatest <- lambdatest + log(rhotest)
btest <- 70
mtest <- 2

Datatest <- Gendata(Ntest,Rtest,Stest,mtest,lambdatest,deltatest,btest)
Datatest

stan_data <- list(
  n = nrow(Datatest[,-c(2,3)]),
  TS = Datatest[,-c(2,3)], # Include T and S
  n_sample = nrow(Datatest[,-3]),
  SAMPLE = Datatest[,-3] # Include T, status, and S
) 
#Actually this code isn't well enough, a reconstruction by yourself is recommended.

fitNGF <- sampling(NGF_model, data = stan_data,
                 iter = 5000, warmup = 3000, chains = 1)
print(fitNGF)

fitGFI <- sampling(GFI_model, data = stan_data,
                 iter = 5000, warmup = 3000, chains = 1)
print(fitGFI)

