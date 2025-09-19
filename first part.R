# Start with a model where survival probability and detection probability depend on time and sex

library(tidyverse)
library(nimble)
set.seed(1604)
set.seed <- c(1,2,3,4)

data <- read.csv("./data/input.csv")
str(data)
det <- data %>%
  mutate("1997" = x97, "1998" = x98, "1999" = x99, "2000" = x00, "2001" = x01, "2002" = x02, "2003" = x03, "2004" = x04) %>%
  select(`1997`:`2004`) %>%
  as.matrix()

SEX <- as.integer(data$sex)
K <- max(SEX)

y <- det
y
T <- ncol(y)      # nb d’occasions
N <- nrow(y)      # nb d’individus

# instant de première capture
first <- apply(y, 1, function(x) min(which(x != 0)))

# mettre NA avant première capture
for(i in 1:N){
  if(first[i] > 1) y[i, 1:(first[i]-1)] <- NA
}

y

my.data <- list(y = y+1)

zinits <- function(y) {
  T <- length(y)
  first <- min(which(y != 0))
  last <- max(which(y == 1))
  z <- rep(NA, T)
  if(last == T) {
    z[first:T] <- 1
  } else {
    if(last == (T - 1)) {
      z[first:(T - 1)] <- 1
      z[T] <- sample(1:2, size = 1)
    } else{
      death <- sample((last + 1):T, size = 1)
      z[first:(death - 1)] <- 1
      z[death:T] <- 2
    }
  }
  return(z)
}

#### Model 1 : p depends on time and sex (categorical) ; phi depends on time and sex (categorical) ####

hmm.p.time.sex_phi.time.sex <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  for(t in 1:(T-1)){
    for(k in 1:K){
      p[t,k] ~ dunif(0,1)
      phi[t,k] ~ dunif(0,1)
    }
  }
  
  for (i in 1:N) {
    for (t in 1:(T-1)) {
      omega[1,1,t,i] <- 1 - p[t, SEX[i]]
      omega[1,2,t,i] <- p[t, SEX[i]]
      omega[2,1,t,i] <- 1
      omega[2,2,t,i] <- 0
      
      gamma[1,1,t,i] <- phi[t, SEX[i]]
      gamma[1,2,t,i] <- 1 - phi[t, SEX[i]]
      gamma[2,1,t,i] <- 0
      gamma[2,2,t,i] <- 1
    }
  }
  
  for(i in 1:N){
    z[i, first[i]] ~ dcat(delta[1:2])
    for(j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,j-1,i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2,j-1,i])
    }
  }
})


my.constants <- list(N = N,
                     T = T,
                     first = first,
                     SEX = SEX,
                     K = K)

initial.values1 <- list(
  p   = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values2 <- list(
  p   = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values3 <- list(
  p   = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values4 <- list(
  p   = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)
initial.values <- c(initial.values1, initial.values2, initial.values3, initial.values4)
params <- c("phi", "p")        

mcmc_p.time.sex_phi.time.sex <- nimbleMCMC(code = hmm.p.time.sex_phi.time.sex,       
                                           constants = my.constants,
                                           data = my.data,
                                           inits = initial.values,
                                           monitors = params,
                                           niter = 20000,
                                           nburnin = 5000,
                                           nchains = 4,
                                           summary = TRUE,
                                           setSeed = set.seed,
                                           WAIC = TRUE)

mcmc_p.time.sex_phi.time.sex$summary
mcmc_p.time.sex_phi.time.sex$WAIC

#### Model 2 : p depends on sex (categorical) and time (random effect) ; phi depends on time and sex (categorical) ####

hmm.p.timeRE.sex_phi.timeCAT.sex <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  
  for(t in 1:(T-1)){
    for(k in 1:K){
      phi[t,k] ~ dunif(0,1)
    }
  }
  for (t in 1:(T-1)){
    eps_p[t] ~ dnorm(0, sd = sdeps_p)
  }
  
  alpha_sex_p[1] <- 0
  for (k in 2:K){
    alpha_sex_p[k] ~ dnorm(0, sd = 1.5)
  }
  
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i,t]) <- mup + alpha_sex_p[SEX[i]] + eps_p[t]
      omega[1,1,i,t] <- 1 - p[i,t]
      omega[1,2,i,t] <- p[i,t]
      omega[2,1,i,t] <- 1
      omega[2,2,i,t] <- 0
      
      gamma[1,1,i,t] <- phi[t, SEX[i]]
      gamma[1,2,i,t] <- 1 - phi[t, SEX[i]]
      gamma[2,1,i,t] <- 0
      gamma[2,2,i,t] <- 1
    }
  }
  
  # priors
  mup ~ dnorm(0, sd= 1.5) # prior intercept on the logit scale
  sdeps_p ~ dunif(0,10) # prior standard deviation for the random effect
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,i,j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2,i, j-1])
    }
  }
})

zinits <- function(y) {
  T <- length(y)
  first <- min(which(y != 0))
  last <- max(which(y == 1))
  z <- rep(NA, T)
  if(last == T) {
    z[first:T] <- 1
  } else {
    if(last == (T - 1)) {
      z[first:(T - 1)] <- 1
      z[T] <- sample(1:2, size = 1)
    } else{
      death <- sample((last + 1):T, size = 1)
      z[first:(death - 1)] <- 1
      z[death:T] <- 2
    }
  }
  return(z)
}

my.constants <- list(N = N,    
                     T = T,     
                     first = first,
                     SEX = SEX,
                     K = K)

initial.values1 <- list(mup = rnorm(1,0,1),
                        sdeps_p = runif(1,0,3),
                        alpha_sex_p = c(0,rnorm(1,0,1.5)),
                        phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
                        z= t(apply(y, 1, zinits)))

initial.values2 <- list(mup = rnorm(1,0,1),
                        sdeps_p = runif(1,0,3),
                        alpha_sex_p = c(0,rnorm(1,0,1.5)),
                        phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
                        z= t(apply(y, 1, zinits)))

initial.values3 <- list(mup = rnorm(1,0,1),
                        sdeps_p = runif(1,0,3),
                        alpha_sex_p = c(0,rnorm(1,0,1.5)),
                        phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
                        z= t(apply(y, 1, zinits)))

initial.values4 <- list(mup = rnorm(1,0,1),
                        sdeps_p = runif(1,0,3),
                        alpha_sex_p = c(0,rnorm(1,0,1.5)),
                        phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
                        z= t(apply(y, 1, zinits)))

initial.values <- c(initial.values1, initial.values2, initial.values3, initial.values4)

params <- c("mup","phi","alpha_sex_p","sdeps_p")

mcmc_p.timeRE.sex_phi.timeCAT.sex <- nimbleMCMC(code = hmm.p.timeRE.sex_phi.timeCAT.sex,       
                                                constants = my.constants,
                                                data = my.data,
                                                inits = initial.values,
                                                monitors = params,
                                                niter = 20000,
                                                nburnin = 5000,
                                                nchains = 4,
                                                summary = TRUE,
                                                setSeed = set.seed,
                                                WAIC = TRUE)

mcmc_p.timeRE.sex_phi.timeCAT.sex$WAIC

#### Model 3 : p depends on sex (categorical) ; phi depends on time and sex (categorical) ####

hmm.p.sex_phi.timeCAT.sex <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  for(t in 1:(T-1)){
    for(k in 1:K){
      phi[t,k] ~ dunif(0,1)
    }
  }
  
  for(k in 1:K){
    p[k] ~ dunif(0,1)
  }
  
  for(i in 1:N){
    omega[1,1,i] <- 1 - p[SEX[i]]
    omega[1,2,i] <- p[SEX[i]]
    omega[2,1,i] <- 1
    omega[2,2,i] <- 0
  }
  
  
  for (i in 1:N) {
    for (t in 1:(T-1)) {
      
      gamma[1,1,t,i] <- phi[t, SEX[i]]
      gamma[1,2,t,i] <- 1 - phi[t, SEX[i]]
      gamma[2,1,t,i] <- 0
      gamma[2,2,t,i] <- 1
    }
  }
  
  for(i in 1:N){
    z[i, first[i]] ~ dcat(delta[1:2])
    for(j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,j-1,i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2,i])
    }
  }
})


my.constants <- list(N = N,
                     T = T,
                     first = first,
                     SEX = SEX,
                     K = K)

initial.values1 <- list(
  p   = runif(K,0,1),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values2 <- list(
  p   = runif(K,0,1),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values3 <- list(
  p   = runif(K,0,1),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values4 <- list(
  p   = runif(K,0,1),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values <- c(initial.values1, initial.values2, initial.values3, initial.values4)
params <- c("phi", "p")        

mcmc_p.sex_phi.timeCAT.sex <- nimbleMCMC(code = hmm.p.sex_phi.timeCAT.sex,       
                                         constants = my.constants,
                                         data = my.data,
                                         inits = initial.values,
                                         monitors = params,
                                         niter = 20000,
                                         nburnin = 5000,
                                         nchains = 4,
                                         summary = TRUE,
                                         setSeed = set.seed,
                                         WAIC = TRUE)

mcmc_p.sex_phi.timeCAT.sex$summary
mcmc_p.sex_phi.timeCAT.sex$WAIC

#### Model 4 : p constant ; phi depends on time and sex (categorical) ####

hmm.p._phi.timeCAT.sex <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  for(t in 1:(T-1)){
    for(k in 1:K){
      phi[t,k] ~ dunif(0,1)
    }
  }
  
  p ~ dunif(0,1)
  omega[1,1] <- 1 - p
  omega[1,2] <- p
  omega[2,1] <- 1
  omega[2,2] <- 0
  
  
  for (i in 1:N) {
    for (t in 1:(T-1)) {
      
      gamma[1,1,t,i] <- phi[t, SEX[i]]
      gamma[1,2,t,i] <- 1 - phi[t, SEX[i]]
      gamma[2,1,t,i] <- 0
      gamma[2,2,t,i] <- 1
    }
  }
  
  for(i in 1:N){
    z[i, first[i]] ~ dcat(delta[1:2])
    for(j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,j-1,i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})


my.constants <- list(N = N,
                     T = T,
                     first = first,
                     SEX = SEX,
                     K = K)

initial.values1 <- list(
  p   = runif(1,0,1),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values2 <- list(
  p   = runif(1,0,1),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values3 <- list(
  p   = runif(1,0,1),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values4 <- list(
  p   = runif(1,0,1),
  phi = matrix(runif((T-1)*K, 0, 1), nrow = T-1, ncol = K),
  z   = t(apply(y, 1, zinits))
)

initial.values <- c(initial.values1, initial.values2, initial.values3, initial.values4)
params <- c("phi", "p")        

mcmc_p._phi.timeCAT.sex <- nimbleMCMC(code = hmm.p._phi.timeCAT.sex,       
                                      constants = my.constants,
                                      data = my.data,
                                      inits = initial.values,
                                      monitors = params,
                                      niter = 20000,
                                      nburnin = 5000,
                                      nchains = 4,
                                      summary = TRUE,
                                      setSeed = set.seed,
                                      WAIC = TRUE)

mcmc_p._phi.timeCAT.sex$summary
mcmc_p._phi.timeCAT.sex$WAIC

#### Model 5 : p depends on sex ; phi depends on time(random effect) and sex (categorical) ####

hmm.p.sex_phi.timeRE.sex <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  for (k in 1:K){
    alpha_sex_phi[k] ~ dnorm(0, sd = 1.5)
  }
  for (t in 1:(T-1)){
    eps_phi[t] ~ dnorm(0, sd = sdeps_phi)
  }
  
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(phi[i,t]) <- muphi + alpha_sex_phi[SEX[i]] + eps_phi[t]
      gamma[1,1,i,t] <- phi[i,t]
      gamma[1,2,i,t] <- 1-phi[i,t]
      gamma[2,1,i,t] <- 0
      gamma[2,2,i,t] <- 1
    }
  }
  
  for(k in 1:K){
    p[k] ~ dunif(0,1)
  }
  
  for(i in 1:N){
    omega[1,1,i] <- 1 - p[SEX[i]]
    omega[1,2,i] <- p[SEX[i]]
    omega[2,1,i] <- 1
    omega[2,2,i] <- 0
  }
  
  muphi ~ dnorm(0, sd= 1.5) # prior intercept on the logit scale
  sdeps_phi ~ dunif(0,10) # prior standard deviation for the random effect
  
  for(i in 1:N){
    z[i, first[i]] ~ dcat(delta[1:2])
    for(j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,i,j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2,i])
    }
  }
})


my.constants <- list(N = N,
                     T = T,
                     first = first,
                     SEX = SEX,
                     K = K)

initial.values1 <- list(
  p   = runif(K,0,1),
  muphi = rnorm(1,0,1),
  sdeps_phi = runif(1,0,3),
  alpha_sex_phi = rnorm(2,0,1.5),
  z   = t(apply(y, 1, zinits))
)

initial.values2 <- list(
  p   = runif(K,0,1),
  muphi = rnorm(1,0,1),
  sdeps_phi = runif(1,0,3),
  alpha_sex_phi = rnorm(2,0,1.5),
  z   = t(apply(y, 1, zinits))
)

initial.values3 <- list(
  p   = runif(K,0,1),
  muphi = rnorm(1,0,1),
  sdeps_phi = runif(1,0,3),
  alpha_sex_phi = rnorm(2,0,1.5),
  z   = t(apply(y, 1, zinits))
)

initial.values4 <- list(
  p   = runif(K,0,1),
  muphi = rnorm(1,0,1),
  sdeps_phi = runif(1,0,3),
  alpha_sex_phi = rnorm(2,0,1.5),
  z   = t(apply(y, 1, zinits))
)

initial.values <- c(initial.values1, initial.values2, initial.values3, initial.values4)
params <- c("p", "muphi","alpha_sex_phi","sdeps_phi")

mcmc_p.sex_phi.timeRE.sex <- nimbleMCMC(code = hmm.p.sex_phi.timeRE.sex,       
                                        constants = my.constants,
                                        data = my.data,
                                        inits = initial.values,
                                        monitors = params,
                                        niter = 20000,
                                        nburnin = 5000,
                                        nchains = 4,
                                        summary = TRUE,
                                        setSeed = set.seed,
                                        WAIC = TRUE)

mcmc_p.sex_phi.timeRE.sex$summary
mcmc_p.sex_phi.timeRE.sex$WAIC
MCMCsummary(mcmc_p.sex_phi.timeRE.sex$samples)
#### Model 6 : p depends on sex ; phi depends on time(continuous) and sex (categorical) ####

TIME <- c(1,2,3,4,5,6,7,8)
TIME <- as.numeric(scale(TIME))

hmm.p.sex_phi.timeCONT.sex <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  
  a_phi[1] <- 0
  for (s in 2:K){
    a_phi[s] ~ dnorm(0, sd = 1.5)
  }
  for (i in 1:N){
    for (t in 1:(T-1)){
      logit(phi[i,t]) <- a0 + a_phi[SEX[i]] + beta_time*TIME[t]
      gamma[1,1,i,t] <- phi[i,t]      # alive -> alive
      gamma[1,2,i,t] <- 1 - phi[i,t]  # alive -> dead
      gamma[2,1,i,t] <- 0             # dead  -> alive
      gamma[2,2,i,t] <- 1             # dead  -> dead
    }
  }
  
  for(k in 1:K){
    p[k] ~ dunif(0,1)
  }
  
  for(i in 1:N){
    omega[1,1,i] <- 1 - p[SEX[i]]
    omega[1,2,i] <- p[SEX[i]]
    omega[2,1,i] <- 1
    omega[2,2,i] <- 0
  }
  
  a0 ~ dnorm(0, sd = 1.5)     # survival intercept
  beta_time ~ dnorm(0, sd = 1.5) # prior intercept
  mup ~ dnorm(0, sd= 1.5) # prior intercept on the logit scale
  sdeps_p ~ dunif(0,10) # prior standard deviation for the random effect
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,i,j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i])
    }
  }
})

initial.values1 <- list(
  p = runif(K,0,1),
  a0 = rnorm(1,0,1),
  z = t(apply(y, 1, zinits)),
  beta_time = rnorm(1,0,1),
  a_phi = c(0, rnorm(1, 0, 0.5))
)

initial.values2 <- list(
  p = runif(K,0,1),
  a0 = rnorm(1,0,1),
  z = t(apply(y, 1, zinits)),
  beta_time = rnorm(1,0,1),
  a_phi = c(0, rnorm(1, 0, 0.5))
)

initial.values3 <- list(
  p = runif(K,0,1),
  a0 = rnorm(1,0,1),
  z = t(apply(y, 1, zinits)),
  beta_time = rnorm(1,0,1),
  a_phi = c(0, rnorm(1, 0, 0.5))
)

initial.values4 <- list(
  p = runif(K,0,1),
  a0 = rnorm(1,0,1),
  z = t(apply(y, 1, zinits)),
  beta_time = rnorm(1,0,1),
  a_phi = c(0, rnorm(1, 0, 0.5))
)



my.constants <- list(N=N, T=T, first=first, SEX=SEX, K=K, TIME=TIME)  # ou TIME=TIME

initial.values <- c(initial.values1, initial.values2, initial.values3, initial.values4)

params <- c("p","a0","beta_time","a_phi","p","phi")

mcmc_p.sex_phi.timeCONT.sex <- nimbleMCMC(code = hmm.p.sex_phi.timeCONT.sex,       
                                          constants = my.constants,
                                          data = my.data,
                                          inits = initial.values,
                                          monitors = params,
                                          niter = 20000,
                                          nburnin = 5000,
                                          nchains = 4,
                                          summary = TRUE,
                                          setSeed = set.seed,
                                          WAIC = TRUE)
mcmc_p.sex_phi.timeCONT.sex$WAIC

#### Model 7 : p depends on sex ; phi depends on sex (categorical) ####

hmm.p.sex_phi.sex <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  
  for(k in 1:K){
    phi[k] ~ dunif(0,1)
    p[k] ~ dunif(0,1)
  }
  
  for(i in 1:N){
    gamma[1,1,i] <- phi[SEX[i]]      # alive -> alive
    gamma[1,2,i] <- 1 - phi[SEX[i]]  # alive -> dead
    gamma[2,1,i] <- 0             # dead  -> alive
    gamma[2,2,i] <- 1 
    
    omega[1,1,i] <- 1 - p[SEX[i]]
    omega[1,2,i] <- p[SEX[i]]
    omega[2,1,i] <- 1
    omega[2,2,i] <- 0
  }
  
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,i])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, i])
    }
  }
})

initial.values1 <- list(
  p = runif(K,0,1),
  phi = runif(K,0,1),
  z = t(apply(y, 1, zinits))
)

initial.values2 <- list(
  p = runif(K,0,1),
  phi = runif(K,0,1),
  z = t(apply(y, 1, zinits))
)

initial.values3 <- list(
  p = runif(K,0,1),
  phi = runif(K,0,1),
  z = t(apply(y, 1, zinits))
)

initial.values4 <- list(
  p = runif(K,0,1),
  phi = runif(K,0,1),
  z = t(apply(y, 1, zinits))
)
my.constants <- list(N=N, T=T, first=first, SEX=SEX, K=K)  

initial.values <- c(initial.values1, initial.values2, initial.values3, initial.values4)

params <- c("p","phi")

mcmc_p.sex_phi.sex <- nimbleMCMC(code = hmm.p.sex_phi.sex,       
                                 constants = my.constants,
                                 data = my.data,
                                 inits = initial.values,
                                 monitors = params,
                                 niter = 20000,
                                 nburnin = 5000,
                                 nchains = 4,
                                 summary = TRUE,
                                 setSeed = set.seed,
                                 WAIC = TRUE)
mcmc_p.sex_phi.sex$WAIC

#### Model 8 : p depends on sex ; phi depends on time (random effect) ####

hmm.p.sex_phi.timeRE <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  for (t in 1:(T-1)){
    eps_phi[t] ~ dnorm(0, sd = sdeps_phi)
  }
  
  for(t in 1:(T-1)){
    logit(phi[t]) <- muphi + eps_phi[t]
    gamma[1,1,t] <- phi[t]
    gamma[1,2,t] <- 1-phi[t]
    gamma[2,1,t] <- 0
    gamma[2,2,t] <- 1
  }
  
  for(k in 1:K){
    p[k] ~ dunif(0,1)
  }
  
  for(i in 1:N){
    omega[1,1,i] <- 1 - p[SEX[i]]
    omega[1,2,i] <- p[SEX[i]]
    omega[2,1,i] <- 1
    omega[2,2,i] <- 0
  }
  
  muphi ~ dnorm(0, sd= 1.5) # prior intercept on the logit scale
  sdeps_phi ~ dunif(0,10) # prior standard deviation for the random effect
  
  for(i in 1:N){
    z[i, first[i]] ~ dcat(delta[1:2])
    for(j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,j-1])
      y[i,j] ~ dcat(omega[z[i,j], 1:2,i])
    }
  }
})


my.constants <- list(N = N,
                     T = T,
                     first = first,
                     SEX = SEX,
                     K = K)

initial.values1 <- list(
  p   = runif(K,0,1),
  muphi = rnorm(1,0,1),
  sdeps_phi = runif(1,0,3),
  z   = t(apply(y, 1, zinits))
)

initial.values2 <- list(
  p   = runif(K,0,1),
  muphi = rnorm(1,0,1),
  sdeps_phi = runif(1,0,3),
  z   = t(apply(y, 1, zinits))
)

initial.values3 <- list(
  p   = runif(K,0,1),
  muphi = rnorm(1,0,1),
  sdeps_phi = runif(1,0,3),
  z   = t(apply(y, 1, zinits))
)

initial.values4 <- list(
  p   = runif(K,0,1),
  muphi = rnorm(1,0,1),
  sdeps_phi = runif(1,0,3),
  z   = t(apply(y, 1, zinits))
)

initial.values <- c(initial.values1, initial.values2, initial.values3, initial.values4)
params <- c("p", "muphi","sdeps_phi")

mcmc_p.sex_phi.timeRE <- nimbleMCMC(code = hmm.p.sex_phi.timeRE,       
                                    constants = my.constants,
                                    data = my.data,
                                    inits = initial.values,
                                    monitors = params,
                                    niter = 20000,
                                    nburnin = 5000,
                                    nchains = 4,
                                    summary = TRUE,
                                    setSeed = set.seed,
                                    WAIC = TRUE)

MCMCsummary(mcmc_p.sex_phi.timeRE$samples)
mcmc_p.sex_phi.timeRE$WAIC

#### Model 9 : everything constant ####

hmm.p._phi. <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  phi~ dunif(0,1)
  p ~ dunif(0,1)
  
  omega[1,1] <- 1 - p
  omega[1,2] <- p
  omega[2,1] <- 1
  omega[2,2] <- 0
  gamma[1,1] <- phi
  gamma[1,2] <- 1 - phi
  gamma[2,1] <- 0
  gamma[2,2] <- 1
  
  for(i in 1:N){
    z[i, first[i]] ~ dcat(delta[1:2])
    for(j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})


my.constants <- list(N = N,
                     T = T,
                     first = first)

initial.values1 <- list(
  p   = runif(1,0,1),
  phi = runif(1,0,1),
  z   = t(apply(y,1,zinits))
)

initial.values2 <- list(
  p   = runif(1,0,1),
  phi = runif(1,0,1),
  z   = t(apply(y,1,zinits))
)

initial.values3 <- list(
  p   = runif(1,0,1),
  phi = runif(1,0,1),
  z   = t(apply(y,1,zinits))
)

initial.values4 <- list(
  p   = runif(1,0,1),
  phi = runif(1,0,1),
  z   = t(apply(y,1,zinits))
)

initial.values <- c(initial.values1, initial.values2, initial.values3, initial.values4)
params <- c("phi", "p")        

mcmc_p._phi. <- nimbleMCMC(code = hmm.p._phi.,       
                           constants = my.constants,
                           data = my.data,
                           inits = initial.values,
                           monitors = params,
                           niter = 20000,
                           nburnin = 5000,
                           nchains = 4,
                           summary = TRUE,
                           setSeed = set.seed,
                           WAIC = TRUE)

mcmc_p._phi.$summary
mcmc_p._phi.$WAIC

save(list = c("mcmc_p.time.sex_phi.time.sex","mcmc_p.timeRE.sex_phi.timeCAT.sex",
              "mcmc_p.sex_phi.timeCAT.sex","mcmc_p._phi.timeCAT.sex","mcmc_p.sex_phi.timeRE.sex",
              "mcmc_p.sex_phi.timeCONT.sex","mcmc_p.sex_phi.sex","mcmc_p.sex_phi.timeRE",
              "mcmc_p._phi."
),
file = "./data/Models_decroissant_seed1604_VF.RData"
)


