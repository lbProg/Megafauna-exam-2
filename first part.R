# Start with a model where survival probability and detection probability depend on time and sex

library(tidyverse)
library(nimble)

set.seed(1604)
chain_seeds <- c(1,2,3,4)

data <- read.csv("./data/input.csv")
y <- data |>
  select(x97:x04) |>
  as.matrix()

T <- ncol(y)      
N <- nrow(y)    

SEX <- as.integer(data$sex)
K <- length(unique(SEX))

first <- apply(y, 1, function(x) min(which(x != 0)))

for(i in 1:N){
  if(first[i] > 1) y[i, 1:(first[i]-1)] <- NA
}

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

# alpha stands for intercepts
# beta stands for coefficients 

#### Model 1 : p depends on time and sex (categorical) ; phi depends on time and sex (categorical) ####

hmm.p.time.sex_phi.time.sex <- nimbleCode({
  
  delta[1] <- 1
  delta[2] <- 0
  
  # Intercepts
  alpha_p0   ~ dnorm(0, sd = 1.5)
  alpha_phi0 ~ dnorm(0, sd = 1.5)
  
  beta_p_time[1]   <- 0
  beta_phi_time[1] <- 0
  for(t in 2:(T-1)){
    beta_p_time[t]   ~ dnorm(0, sd = 1.5)
    beta_phi_time[t] ~ dnorm(0, sd = 1.5)
  }
  
  beta_p_sex[1]   <- 0
  beta_phi_sex[1] <- 0
  for(k in 2:K){
    beta_p_sex[k]   ~ dnorm(0, sd = 1.5)
    beta_phi_sex[k] ~ dnorm(0, sd = 1.5)
  }
  
  # Likelihoods
  
    for(t in 1:(T-1)){
    for(k in 1:K){
      logit(p[t,k])   <- alpha_p0   + beta_p_time[t]   + beta_p_sex[k]
      logit(phi[t,k]) <- alpha_phi0 + beta_phi_time[t] + beta_phi_sex[k]
    }
    }
  
  for (i in 1:N) {
    for (t in 1:(T-1)) {
      # Detection
      omega[1,1,t,i] <- 1 - p[t, SEX[i]]
      omega[1,2,t,i] <- p[t, SEX[i]]
      omega[2,1,t,i] <- 1
      omega[2,2,t,i] <- 0
      
      # Transition 
      gamma[1,1,t,i] <- phi[t, SEX[i]]
      gamma[1,2,t,i] <- 1 - phi[t, SEX[i]]
      gamma[2,1,t,i] <- 0
      gamma[2,2,t,i] <- 1
    }
  }
  
  for(i in 1:N){
    z[i, first[i]] ~ dcat(delta[1:2])
    for(j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2, j-1, i])
      y[i,j] ~ dcat(omega[z[i,j],   1:2, j-1, i])
    }
  }
})



my.constants <- list(N = N,
                     T = T,
                     first = first,
                     SEX = SEX,
                     K = K)

makeinits <- function() {
  list(
    alpha_p0   = rnorm(1, 0, sd = 1.5),
    alpha_phi0 = rnorm(1, 0, sd = 1.5),
    
    beta_p_time   = c(0, rnorm(T-2, 0, sd = 1.5)),
    beta_phi_time = c(0, rnorm(T-2, 0, sd = 1.5)),
    
    beta_p_sex   = c(0, rnorm(K-1, 0, sd = 1.5)),
    beta_phi_sex = c(0, rnorm(K-1, 0, sd = 1.5)),
    
    z = t(apply(y, 1, zinits))
  )
}
initial.values <- replicate(4, makeinits(), simplify = FALSE)
params <- c(
  "alpha_p0","alpha_phi0",
  "beta_p_time","beta_phi_time",
  "beta_p_sex","beta_phi_sex",
  "p","phi"
)    

mcmc_p.time.sex_phi.time.sex <- nimbleMCMC(code = hmm.p.time.sex_phi.time.sex,       
                                           constants = my.constants,
                                           data = my.data,
                                           inits = initial.values,
                                           monitors = params,
                                           niter = 20000,
                                           nburnin = 5000,
                                           nchains = 4,
                                           summary = TRUE,
                                           setSeed = chain_seeds,
                                           WAIC = TRUE)

mcmc_p.time.sex_phi.time.sex$summary
mcmc_p.time.sex_phi.time.sex$WAIC

#### Model 2 : p depends on sex (categorical) and time (random effect) ; phi depends on time and sex (categorical) ####

hmm.p.timeRE.sex_phi.timeCAT.sex <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  
  for (t in 1:(T-1)){
    eps_p[t] ~ dnorm(0, sd = sdeps_p)
  }
  
  beta_p_sex[1] <- 0
  for (k in 2:K){
    beta_p_sex[k] ~ dnorm(0, sd = 1.5)
  }
  
  alpha_phi0 ~ dnorm(0, sd = 1.5)
  beta_phi_time[1] <- 0
  for(t in 2:(T-1)) beta_phi_time[t] ~ dnorm(0, sd = 1.5)
  beta_phi_sex[1] <- 0
  for(k in 2:K)    beta_phi_sex[k]  ~ dnorm(0, sd = 1.5)
  
  for(t in 1:(T-1)){
    for(k in 1:K){
      logit(phi[t,k]) <- alpha_phi0 + beta_phi_time[t] + beta_phi_sex[k]
    }
  }
  
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(p[i,t]) <- mup + beta_p_sex[SEX[i]] + eps_p[t]
      omega[1,1,i,t] <- 1 - p[i,t]
      omega[1,2,i,t] <- p[i,t]
      omega[2,1,i,t] <- 1
      omega[2,2,i,t] <- 0
      
      gamma[1,1,t,i] <- phi[t, SEX[i]]
      gamma[1,2,t,i] <- 1 - phi[t, SEX[i]]
      gamma[2,1,t,i] <- 0
      gamma[2,2,t,i] <- 1
    }
  }
  
  # priors
  mup ~ dnorm(0, sd= 1.5) # prior intercept on the logit scale
  sdeps_p ~ dunif(0,10) # prior standard deviation for the random effect
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2,j-1,i])
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

makeinits <- function() {
  list(
    mup = rnorm(1,0,1),
    sdeps_p = runif(1,0,3),
    beta_p_sex = runif(c(0, rnorm(K-1,0,1.5))),
    alpha_phi0 = rnorm(1, 0, sd = 1.5),
    beta_phi_time = c(0, rnorm(T-2, 0, sd = 1.5)),
    beta_phi_sex = c(0, rnorm(K-1, 0, sd = 1.5)),
    z = t(apply(y, 1, zinits))
  )
}


initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c(
  "mup","sdeps_p","beta_p_sex",
  "alpha_phi0","beta_phi_time","beta_phi_sex")

mcmc_p.timeRE.sex_phi.timeCAT.sex <- nimbleMCMC(code = hmm.p.timeRE.sex_phi.timeCAT.sex,       
                                                constants = my.constants,
                                                data = my.data,
                                                inits = initial.values,
                                                monitors = params,
                                                niter = 20000,
                                                nburnin = 5000,
                                                nchains = 4,
                                                summary = TRUE,
                                                setSeed = chain_seeds,
                                                WAIC = TRUE)

mcmc_p.timeRE.sex_phi.timeCAT.sex$WAIC

#### Model 3 : p depends on sex (categorical) ; phi depends on time and sex (categorical) ####

hmm.p.sex_phi.timeCAT.sex <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  for(k in 1:K){
    p[k] ~ dunif(0,1)
  }
  
  for(i in 1:N){
    omega[1,1,i] <- 1 - p[SEX[i]]
    omega[1,2,i] <- p[SEX[i]]
    omega[2,1,i] <- 1
    omega[2,2,i] <- 0
  }

  alpha_phi0 ~ dnorm(0, sd = 1.5)
  beta_phi_time[1] <- 0
  for(t in 2:(T-1)) beta_phi_time[t] ~ dnorm(0, sd = 1.5)
  beta_phi_sex[1] <- 0
  for(k in 2:K)    beta_phi_sex[k]  ~ dnorm(0, sd = 1.5)
  
  for(t in 1:(T-1)){
    for(k in 1:K){
      logit(phi[t,k]) <- alpha_phi0 + beta_phi_time[t] + beta_phi_sex[k]
    }
  }
  
  for(i in 1:N){
    for(t in 1:(T-1)){
      gamma[1,1,t,i] <-     phi[t, SEX[i]]
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

makeinits <- function() {
  list(
    p = runif(K,0,1),
    alpha_phi0 = rnorm(1, 0, sd = 1.5),
    beta_phi_time = c(0, rnorm(T-2, 0, sd = 1.5)),
    beta_phi_sex = c(0, rnorm(K-1, 0, sd = 1.5)),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c(
  "p",
  "alpha_phi0","beta_phi_time","beta_phi_sex")

mcmc_p.sex_phi.timeCAT.sex <- nimbleMCMC(code = hmm.p.sex_phi.timeCAT.sex,       
                                         constants = my.constants,
                                         data = my.data,
                                         inits = initial.values,
                                         monitors = params,
                                         niter = 20000,
                                         nburnin = 5000,
                                         nchains = 4,
                                         summary = TRUE,
                                         setSeed = chain_seeds,
                                         WAIC = TRUE)

mcmc_p.sex_phi.timeCAT.sex$summary
mcmc_p.sex_phi.timeCAT.sex$WAIC

#### Model 4 : p constant ; phi depends on time and sex (categorical) ####

hmm.p._phi.timeCAT.sex <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  p ~ dunif(0,1)
  omega[1,1] <- 1 - p
  omega[1,2] <- p
  omega[2,1] <- 1
  omega[2,2] <- 0
  
  alpha_phi0 ~ dnorm(0, sd = 1.5)
  beta_phi_time[1] <- 0
  for(t in 2:(T-1)) beta_phi_time[t] ~ dnorm(0, sd = 1.5)
  beta_phi_sex[1] <- 0
  for(k in 2:K)    beta_phi_sex[k]  ~ dnorm(0, sd = 1.5)
  
  for(t in 1:(T-1)){
    for(k in 1:K){
      logit(phi[t,k]) <- alpha_phi0 + beta_phi_time[t] + beta_phi_sex[k]
    }
  }
  
  for(i in 1:N){
    for(t in 1:(T-1)){
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
makeinits <- function() {
  list(
    p = runif(1,0,1),
    alpha_phi0 = rnorm(1, 0, sd = 1.5),
    beta_phi_time = c(0, rnorm(T-2, 0, sd = 1.5)),
    beta_phi_sex = c(0, rnorm(K-1, 0, sd = 1.5)),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c(
  "p",
  "alpha_phi0","beta_phi_time","beta_phi_sex")

mcmc_p._phi.timeCAT.sex <- nimbleMCMC(code = hmm.p._phi.timeCAT.sex,       
                                      constants = my.constants,
                                      data = my.data,
                                      inits = initial.values,
                                      monitors = params,
                                      niter = 20000,
                                      nburnin = 5000,
                                      nchains = 4,
                                      summary = TRUE,
                                      setSeed = chain_seeds,
                                      WAIC = TRUE)

mcmc_p._phi.timeCAT.sex$summary
mcmc_p._phi.timeCAT.sex$WAIC

#### Model 5 : p depends on sex ; phi depends on time(random effect) and sex (categorical) ####

hmm.p.sex_phi.timeRE.sex <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  beta_sex_phi[1] <- 0
  for (k in 2:K){
    beta_sex_phi[k] ~ dnorm(0, sd = 1.5)
  }
  for (t in 1:(T-1)){
    eps_phi[t] ~ dnorm(0, sd = sdeps_phi)
  }
  
  for(i in 1:N){
    for(t in 1:(T-1)){
      logit(phi[i,t]) <- muphi + beta_sex_phi[SEX[i]] + eps_phi[t]
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

makeinits <- function() {
  list(
    p   = runif(K,0,1),
    muphi = rnorm(1,0,1),
    sdeps_phi = runif(1,0,3),
    beta_sex_phi = rnorm(2,0,1.5),
    z   = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c("p", "muphi","beta_sex_phi","sdeps_phi", "eps_phi")

mcmc_p.sex_phi.timeRE.sex <- nimbleMCMC(code = hmm.p.sex_phi.timeRE.sex,       
                                        constants = my.constants,
                                        data = my.data,
                                        inits = initial.values,
                                        monitors = params,
                                        niter = 20000,
                                        nburnin = 5000,
                                        nchains = 4,
                                        summary = TRUE,
                                        setSeed = chain_seeds,
                                        WAIC = TRUE)

mcmc_p.sex_phi.timeRE.sex$summary
mcmc_p.sex_phi.timeRE.sex$WAIC

#### Model 6 : p depends on sex ; phi depends on time(continuous) and sex (categorical) ####

TIME <- c(1,2,3,4,5,6,7,8)
TIME <- as.numeric(scale(TIME))

hmm.p.sex_phi.timeCONT.sex <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  
  beta_phi_sexe[1] <- 0
  for (s in 2:K){
    beta_phi_sexe[s] ~ dnorm(0, sd = 1.5)
  }
  for (i in 1:N){
    for (t in 1:(T-1)){
      logit(phi[i,t]) <- alpha_phi0 + beta_phi_sexe[SEX[i]] + beta_time*TIME[t]
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
  
  alpha_phi0 ~ dnorm(0, sd = 1.5)     # survival intercept
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


makeinits <- function() {
  list(
    p = runif(K,0,1),
    alpha_phi0 = rnorm(1,0,1),
    z = t(apply(y, 1, zinits)),
    beta_time = rnorm(1,0,1),
    beta_phi_sexe = c(0, rnorm(K-1, 0, 0.5))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

my.constants <- list(N=N, 
                     T=T, 
                     first=first, 
                     SEX=SEX, 
                     K=K, 
                     TIME=TIME)  


params <- c("p","alpha_phi0","beta_time","beta_phi_sexe","p","phi")

mcmc_p.sex_phi.timeCONT.sex <- nimbleMCMC(code = hmm.p.sex_phi.timeCONT.sex,       
                                          constants = my.constants,
                                          data = my.data,
                                          inits = initial.values,
                                          monitors = params,
                                          niter = 20000,
                                          nburnin = 5000,
                                          nchains = 4,
                                          summary = TRUE,
                                          setSeed = chain_seeds,
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


makeinits <- function() {
  list(
    p = runif(K,0,1),
    phi = runif(K,0,1),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

my.constants <- list(N=N, T=T, first=first, SEX=SEX, K=K)  

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
                                 setSeed = chain_seeds,
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

makeinits <- function() {
  list(
    p   = runif(K,0,1),
    muphi = rnorm(1,0,1),
    sdeps_phi = runif(1,0,3),
    z   = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

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
                                    setSeed = chain_seeds,
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

makeinits <- function() {
  list(
    p   = runif(1,0,1),
    phi = runif(1,0,1),
    z   = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

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
                           setSeed = chain_seeds,
                           WAIC = TRUE)

mcmc_p._phi.$summary
mcmc_p._phi.$WAIC

save(list = c("mcmc_p.time.sex_phi.time.sex","mcmc_p.timeRE.sex_phi.timeCAT.sex",
              "mcmc_p.sex_phi.timeCAT.sex","mcmc_p._phi.timeCAT.sex","mcmc_p.sex_phi.timeRE.sex",
              "mcmc_p.sex_phi.timeCONT.sex","mcmc_p.sex_phi.sex","mcmc_p.sex_phi.timeRE",
              "mcmc_p._phi."
),
file = "./data/Models_part1.RData"
)


