#### Model 1 : everything constant ####

hmm.p._phi. <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  phi~ dunif(0,1)
  omega[1,1] <- 1 - p
  omega[1,2] <- p
  omega[2,1] <- 1
  omega[2,2] <- 0
  
  p ~ dunif(0,1)
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


my.constants <- list(N = N, T = T, first = first)

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

#### Model 2 : p depends on time (random effect) ####

hmm.p.timeRE_phi. <- nimbleCode({
  delta[1] <- 1          # Pr(alive t = 1) = 1
  delta[2] <- 0          # Pr(dead t = 1) = 0
  
  for (t in 1:(T-1)){
    logit(p[t]) <- mup + eps_p[t] # eps is the random effect (epsilon)
    eps_p[t] ~ dnorm(0, sd = sdeps_p)
    omega[1,1,t] <- 1 - p[t]    # Pr(alive t -> non-detected t)
    omega[1,2,t] <- p[t]        # Pr(alive t -> detected t)
    omega[2,1,t] <- 1           # Pr(dead t -> non-detected t)
    omega[2,2,t] <- 0           # Pr(dead t -> detected t)
  }

  phi ~ dunif(0, 1)          # prior detection
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0           # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1           # Pr(dead t -> dead t+1)
  
  
  # priors
  mup ~ dnorm(0, sd= 1.5) # prior intercept on the logit scale
  sdeps_p ~ dunif(0,10) # prior standard deviation for the random effect
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2, j-1])
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
                     first = first)

makeinits <- function() {
  list(
    mup = rnorm(1,0,1),
    sdeps_p = runif(1,0,3),
    phi = runif(1,0,1),
    z = t(apply(y, 1, zinits))
  )
}


initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c(
  "mup","sdeps_p","phi")

mcmc_p.timeRE_phi. <- nimbleMCMC(code = hmm.p.timeRE_phi.,       
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

mcmc_p.timeRE_phi.

#### Model 3 : p depends on time (categorical) ####

hmm.p.timeCAT_phi. <- nimbleCode({
  delta[1] <- 1          
  delta[2] <- 0 
  
  for(t in 1:(T-1)){
    p[t] ~ dunif(0,1)
    omega[1,1,t] <- 1 - p[t]
    omega[1,2,t] <- p[t]
    omega[2,1,t] <- 1
    omega[2,2,t] <- 0
  }
  
  phi ~ dunif(0, 1) 
  gamma[1,1] <- phi      
  gamma[1,2] <- 1 - phi  
  gamma[2,1] <- 0           
  gamma[2,2] <- 1           
  
  
  for(i in 1:N){
    z[i, first[i]] ~ dcat(delta[1:2])
    for(j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2,j-1])
    }
  }
})


my.constants <- list(N = N,
                     T = T,
                     first = first
                     )

makeinits <- function() {
  list(
    p = runif((T-1),0,1),
    phi = runif(1,0,1),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c(
  "p","phi")

mcmc_p.timeCAT_phi. <- nimbleMCMC(code = hmm.p.timeCAT_phi.,       
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

mcmc_p.timeCAT_phi.$WAIC

#### Model 4 : p depends on sex (categorical) ####

hmm.p.sex_phi. <- nimbleCode({
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

  phi ~ dunif(0, 1)          # prior detection
  gamma[1,1] <- phi      # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi  # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0           # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1           # Pr(dead t -> dead t+1)
  
  for(i in 1:N){
    z[i, first[i]] ~ dcat(delta[1:2])
    for(j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
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
    phi = runif(1,0,1),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c(
  "p","phi")

mcmc_p.sex_phi. <- nimbleMCMC(code = hmm.p.sex_phi.,       
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

mcmc_p.sex_phi.$WAIC


#### Model 5 : p depends on sex (categorical) and phi depends on time (RE) ####

hmm.p.sex_phi.timeRE <- nimbleCode({
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
  
  for (t in 1:(T-1)){
    logit(phi[t]) <- muphi + eps_phi[t] # eps is the random effect (epsilon)
    eps_phi[t] ~ dnorm(0, sd = sdeps_phi)
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
  }
  
  muphi ~ dnorm(0, sd = 1.5)
  sdeps_phi ~ dunif(0, 10)
  
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
    p = runif(K,0,1),
    muphi = rnorm(1,0,1),
    sdeps_phi = runif(1,0,3),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c(
  "muphi","sdeps_phi","p")

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

mcmc_p.sex_phi.timeRE$WAIC

#### Model 6 : p depends on sex (categorical) and phi depends on time (categorical) ####

hmm.p.sex_phi.timeCAT <- nimbleCode({
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
  
  for(t in 1:(T-1)){
    phi[t] ~ dunif(0,1)   # un p par année
    gamma[1,1,t] <- phi[t]
    gamma[1,2,t] <- 1 - phi[t]
    gamma[2,1,t] <- 0
    gamma[2,2,t] <- 1
  }

  
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
    p = runif(K,0,1),
    phi = runif(T-1,0,1),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c(
  "phi","p")

mcmc_p.sex_phi.timeCAT <- nimbleMCMC(code = hmm.p.sex_phi.timeCAT,       
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

mcmc_p.sex_phi.timeCAT$WAIC

#### Model 7 : p depends on sex (categorical) and phi depends on time (continuous) ####
TIME = c(1,2,3,4,5,6,7,8)
TIME = as.numeric(TIME)

hmm.p.sex_phi.timeCONT <- nimbleCode({
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
  
  for (t in 1:(T-1)){
    logit(phi[t]) <- alpha_phi0 + beta_time*(TIME[t]-1)
    gamma[1,1,t] <- phi[t]     # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
  }
  alpha_phi0 ~ dnorm(0, sd = 1.5) # prior intercept
  beta_time ~ dnorm(0, sd = 1.5) # prior slope
  
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
                     TIME = TIME,
                     first = first,
                     SEX = SEX,
                     K = K)

makeinits <- function() {
  list(
    p = runif(K,0,1),
    alpha_phi0 = rnorm(1,0,1),
    beta_time = rnorm(1,0,1),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c(
  "alpha_phi0","beta_time","p")

mcmc_p.sex_phi.timeCONT <- nimbleMCMC(code = hmm.p.sex_phi.timeCONT,       
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

mcmc_p.sex_phi.timeCONT$WAIC

#### Model 8 : p depends on sex ; phi depends on time (random effect) and sex (categorical) ####

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
    beta_sex_phi = c(0, rnorm(K-1,0,1.5)),
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

#### OXY ####
#### OXY CONTINUE ####

hmm.p.sex_phi.timeRE.oxyCONT <- nimbleCode({
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
  
  beta_time[1] <- 0
  for(t in 2:(T-1)){
    beta_time[t] ~ dnorm(0, sd = 1.5)
  }
  
  for (t in 1:(T-1)){
    logit(phi[t]) <- alpha_phi0 + beta_oxy * OXY[t] + beta_time[t]
    gamma[1,1,t] <- phi[t]      # Pr(alive t -> alive t+1)
    gamma[1,2,t] <- 1 - phi[t]  # Pr(alive t -> dead t+1)
    gamma[2,1,t] <- 0           # Pr(dead t -> alive t+1)
    gamma[2,2,t] <- 1           # Pr(dead t -> dead t+1)
  }
  
  alpha_phi0 ~ dnorm(0, sd = 1.5)
  beta_oxy   ~ dnorm(0, sd = 1.5) 
  
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
                     K = K,
                     OXY = log_oxy_sc)

makeinits <- function() {
  list(
    p = runif(K,0,1),
    alpha_phi0 = rnorm(1,0,1),
    beta_oxy = rnorm(1,0,1),
    beta_time = c(0, rnorm(T-2,0,1.5)),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c("p", "alpha_phi0", "beta_oxy", "beta_time")

mcmc_p.sex_phi.timeRE.oxyCONT <- nimbleMCMC(code = hmm.p.sex_phi.timeRE.oxyCONT,       
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

mcmc_p.sex_phi.timeRE.oxyCONT$WAIC
mcmc_p.sex_phi.timeCAT.oxyCONT <- mcmc_p.sex_phi.timeRE.oxyCONT

#### OXY CATEGORIELLE (10 déciles par sex) ####

hmm.p.sex_phi.timeRE.oxyCAT <- nimbleCode({
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
  
  beta_time[1] <- 0
  for(t in 2:(T-1)){
    beta_time[t] ~ dnorm(0, sd = 1.5)
  }
  
  beta_oxy[1] <- 0
  for(g in 2:G){
    beta_oxy[g] ~ dnorm(0, sd = 1.5)
  }
  
  for(i in 1:N){
    for (t in 1:(T-1)){
      logit(phi[i,t]) <- alpha_phi0 + beta_oxy[oxyd[i]] + beta_time[t]
      gamma[1,1,i,t] <- phi[i,t]      # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t]  # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0           # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1           # Pr(dead t -> dead t+1)
    }
  }
  
  alpha_phi0 ~ dnorm(0, sd = 1.5)
  
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
                     K = K,
                     oxyd = oxyd,
                     G = G)

makeinits <- function() {
  list(
    p = runif(K,0,1),
    alpha_phi0 = rnorm(1,0,1),
    beta_oxy = c(0, rnorm(G-1,0,1)),
    beta_time = c(0, rnorm(T-2,0,1.5)),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c("p", "alpha_phi0", "beta_oxy", "beta_time")

mcmc_p.sex_phi.timeRE.oxyCAT <- nimbleMCMC(code = hmm.p.sex_phi.timeRE.oxyCAT,       
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

mcmc_p.sex_phi.timeRE.oxyCAT$WAIC
mcmc_p.sex_phi.timeCAT.oxyCAT <- mcmc_p.sex_phi.timeRE.oxyCAT
#### 2 derniers déciles pas poolés pour les 2 sexes ####

hmm.p.sex_phi.timeRE.oxyCAT2D <- nimbleCode({
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
  
  beta_time[1] <- 0
  for(t in 2:(T-1)){
    beta_time[t] ~ dnorm(0, sd = 1.5)
  }
  
  beta_oxy[1] <- 0
  for(g in 2:G){
    beta_oxy[g] ~ dnorm(0, sd = 1.5)
  }
  
  for(i in 1:N){
    for (t in 1:(T-1)){
      logit(phi[i,t]) <- alpha_phi0 + beta_oxy[oxyd[i]] + beta_time[t]
      gamma[1,1,i,t] <- phi[i,t]      # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t]  # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0           # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1           # Pr(dead t -> dead t+1)
    }
  }
  
  alpha_phi0 ~ dnorm(0, sd = 1.5)
  
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
                     K = K,
                     oxyd = oxyd_2D,
                     G = G)

makeinits <- function() {
  list(
    p = runif(K,0,1),
    alpha_phi0 = rnorm(1,0,1),
    beta_oxy = c(0, rnorm(G-1,0,1)),
    beta_time = c(0, rnorm(T-2,0,1.5)),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c("p", "alpha_phi0", "beta_oxy", "beta_time")

mcmc_p.sex_phi.timeRE.oxyCAT2D <- nimbleMCMC(code = hmm.p.sex_phi.timeRE.oxyCAT2D,       
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

mcmc_p.sex_phi.timeRE.oxyCAT2D$WAIC
mcmc_p.sex_phi.timeCAT.oxyCAT2D <- mcmc_p.sex_phi.timeRE.oxyCAT2D

#### 3 derniers déciles pas poolés pour les 2 sexes ####

hmm.p.sex_phi.timeRE.oxyCAT3D <- nimbleCode({
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
  
  beta_time[1] <- 0
  for(t in 2:(T-1)){
    beta_time[t] ~ dnorm(0, sd = 1.5)
  }
  
  beta_oxy[1] <- 0
  for(g in 2:G){
    beta_oxy[g] ~ dnorm(0, sd = 1.5)
  }
  
  for(i in 1:N){
    for (t in 1:(T-1)){
      logit(phi[i,t]) <- alpha_phi0 + beta_oxy[oxyd[i]] + beta_time[t]
      gamma[1,1,i,t] <- phi[i,t]      # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t]  # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0           # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1           # Pr(dead t -> dead t+1)
    }
  }
  
  alpha_phi0 ~ dnorm(0, sd = 1.5)
  
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
                     K = K,
                     oxyd = oxyd_3D,
                     G = G)

makeinits <- function() {
  list(
    p = runif(K,0,1),
    alpha_phi0 = rnorm(1,0,1),
    beta_oxy = c(0, rnorm(G-1,0,1)),
    beta_time = c(0, rnorm(T-2,0,1.5)),
    z = t(apply(y, 1, zinits))
  )
}

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c("p", "alpha_phi0", "beta_oxy", "beta_time")

mcmc_p.sex_phi.timeRE.oxyCAT3D <- nimbleMCMC(code = hmm.p.sex_phi.timeRE.oxyCAT3D,       
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

mcmc_p.sex_phi.timeRE.oxyCAT3D$WAIC
mcmc_p.sex_phi.timeCAT.oxyCAT3D <- mcmc_p.sex_phi.timeRE.oxyCAT3D
#### 2 derniers déciles pour jsplu et 3 derniers déciles pour jsplu ####

hmm.p.sex_phi.timeRE.oxyCAT_2DF3DM <- nimbleCode({
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
  
  beta_time[1] <- 0
  for(t in 2:(T-1)){
    beta_time[t] ~ dnorm(0, sd = 1.5)
  }
  
  beta_oxy[1,1] <- 0
  for(g in 2:G){
    beta_oxy[1,g] ~ dnorm(0, sd = 1.5)
  }
  
  for(g in 1:G){
    beta_oxy[2,g] ~ dnorm(0, sd = 1.5)
  }
  
  for(i in 1:N){
    for (t in 1:(T-1)){
      logit(phi[i,t]) <- alpha_phi0 + beta_oxy[SEX[i], oxyd[i]] + beta_time[t]
      gamma[1,1,i,t] <- phi[i,t]      # Pr(alive t -> alive t+1)
      gamma[1,2,i,t] <- 1 - phi[i,t]  # Pr(alive t -> dead t+1)
      gamma[2,1,i,t] <- 0           # Pr(dead t -> alive t+1)
      gamma[2,2,i,t] <- 1           # Pr(dead t -> dead t+1)
    }
  }
  
  alpha_phi0 ~ dnorm(0, sd = 1.5)
  
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
                     K = K,
                     oxyd = oxyd_mf,
                     G = G)

makeinits <- function() {
  B <- matrix(NA, nrow=K, ncol=G)
  if (G > 1) B[1, 2:G] <- rnorm(G-1, 0, 1)
  if (K > 1) B[2:K, 1:G] <- rnorm((K-1)*G, 0, 1)
  B[1,1] <- 0
  
  list(
    p = runif(K, 0, 1),
    alpha_phi0 = rnorm(1, 0, 1),
    beta_time = c(0, rnorm(T-2, 0, 1.5)),
    beta_oxy = B,
    z = t(apply(y, 1, zinits))
  )
}
params <- c("p", "alpha_phi0", "beta_time", "beta_oxy")

initial.values <- replicate(4, makeinits(), simplify = FALSE)

params <- c("p", "alpha_phi0", "beta_oxy", "beta_time")

mcmc_p.sex_phi.timeRE.oxyCAT_2DF3DM <- nimbleMCMC(code = hmm.p.sex_phi.timeRE.oxyCAT_2DF3DM,       
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

mcmc_p.sex_phi.timeRE.oxyCAT_2DF3DM$WAIC
mcmc_p.sex_phi.timeCAT.oxyCAT_2DF3DM <- mcmc_p.sex_phi.timeRE.oxyCAT_2DF3DM

save(list = c("mcmc_p._phi.","mcmc_p.timeRE_phi.","mcmc_p.timeCAT_phi.","mcmc_p.sex_phi.",
              "mcmc_p.sex_phi.timeRE","mcmc_p.sex_phi.timeCAT","mcmc_p.sex_phi.timeCONT",
              "mcmc_p.sex_phi.timeRE.sex","mcmc_p.sex_phi.timeCAT.oxyCONT",
              "mcmc_p.sex_phi.timeCAT.oxyCAT","mcmc_p.sex_phi.timeCAT.oxyCAT2D",
              "mcmc_p.sex_phi.timeCAT.oxyCAT3D","mcmc_p.sex_phi.timeCAT.oxyCAT_2DF3DM"
),
file = "./data/Models_croissant_seed1604_VF.RData"
)

mcmc_p.sex_phi.timeCAT.oxyCAT_2DF3DM$WAIC

### En réalité : il faudrait regrouper les 9 premiers pour les femelles et les 8 premiers pour les mâles, puis le 10 tout seul pour les femelles et 9 et 10 ensemble pour les mâles

oxyd_3D
