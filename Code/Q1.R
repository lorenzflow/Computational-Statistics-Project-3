# plot density
f <- function(x,y) exp(-x^2/100 - (y+3/100*x^2-3)^2)

x_vec <- seq(-15,15,0.1)
y_vec <- seq(-5,10,0.1)
grid <- expand.grid(x=x_vec,y=y_vec)
f_points <- f(grid$x, grid$y)

dat <- data.frame(x = c(grid$x), y=c(grid$y), f=c(f_points))

fig <- plot_ly(x = dat$x, y=dat$y, z=dat$f, type = 'mesh3d') 
fig
# simple metropolis hastings
MetrHastw <- function(x0, sigmapropx, sigmapropy, nsteps, f){
  X <- matrix(NA,nrow=2, ncol=nsteps+1)
  X[,1] <- x0
  
  for (i in 2:(nsteps+1)){
    x <- rnorm(1,mean=X[1,(i-1)], sd=sigmapropx)
    y <- rnorm(1,mean=X[2,(i-1)], sd=sigmapropy)
    if (runif(1)<=min(exp(log(f(x,y))-log(f(X[1,i-1], X[2,i-1]))), 1)) X[,i] <- c(x,y)
    else X[,i] <- X[,i-1]
  }
  X[,-1]
}
# algorithm taking into account scale and dependency 
# by specifying the covariance matrix of the proposal
MetrHastw_cov <- function(x0, sigmapropx, sigmapropy, cova, nsteps, f){
  X <- matrix(NA,nrow=2, ncol=nsteps+1)
  X[,1] <- x0
  
  for (i in 2:(nsteps+1)){
    prop <- rmvnorm(1,mean=X[,i-1], sigma=matrix(c(sigmapropx^2,cova,cova,sigmapropy^2), ncol=2))
    x <- prop[1]
    y <- prop[2]
    if (runif(1)<=min(exp(log(f(x,y))-log(f(X[1,i-1], X[2,i-1]))), 1)) X[,i] <- c(x,y)
    else X[,i] <- X[,i-1]
  }
  X[,-1]
}


# tuning parameters
X_mh_adj <- MetrHastw_cov(c(0,3), 6.884244,2.122618, 0.2264226,10^5,f)
cov <- cov(X_mh_adj[1,], X_mh_adj[2,])
sd_x <- sd(X_mh[1,])
sd_y <- sd(X_mh[2,])
cov <- 0
sd_x <- 1
sd_y <- 1

for (i in 1:5){
  
  print(c(sd_x, sd_y,cov))
  Dummy <- MetrHastw_cov(c(0,3), sd_x,sd_y, cov,10^5,f)
  # update parameters
  cov <- cov(Dummy[1,10^4:10^5], Dummy[2,10^4:10^5])
  sd_x <- sd(Dummy[1,10^4:10^5])
  sd_y <- sd(Dummy[2,10^4:10^5])
  
}
# Then trying several values around this estimate
# multiple chain convergence diagnostics
set.seed(0)
N <- 5*10^4 # chain length
X0s <- replicate(6, rnorm(2,c(0,3), 8))
samplesX <- matrix(rep(0,dim(X0s)[2]*N), nrow=dim(X0s)[2])
samplesY <- matrix(rep(0,dim(X0s)[2]*N), nrow=dim(X0s)[2])
for (i in 1:dim(X0s)[2]){
  sample <- MetrHastw_cov(X0s[,i], 7 ,2.1, 0.22,N,f)
  samplesX[i,] <- sample[1,]
  samplesY[i,] <- sample[2,]
}


chainsX <- lapply(1:dim(X0s)[2], function(i) mcmc(samplesX[i,]))
chainsY <- lapply(1:dim(X0s)[2], function(i) mcmc(samplesY[i,]))
X0s
traceplot(chainsX, main='Traceplots for X')
traceplot(chainsY, main='Traceplots for Y')
gelman.diag(mcmc.list(chainsX))
gelman.diag(mcmc.list(chainsY))
gelman.plot(mcmc.list(chainsX), main='Gelman Plot for X')
gelman.plot(mcmc.list(chainsY), main='Gelman Plot for Y')

# sample
set.seed(0)
X_mh_final <- MetrHastw_cov(c(0,3), 7 ,2.1, 0.22,10^5,f)
plot(X_mh_final[1,], type='l')
plot(X_mh_final[2,], type='l')

plot(X_mh_final[1,(500:10^5)], X_mh_final[2,(500:10^5)], col=rgb(red=0, green=0, blue=1, alpha=0.01), xlab='x', ylab='y', main='Joint Distribution f(x,y)', xlim=c(-25,25), ylim=c(-20,20))
acf(X_mh_final[1,10^4:10^5], lag.max=200, main='Acf of X')
acf(X_mh_final[2,10^4:10^5], lag.max=200, main='Acf of Y')

#mean estimates
mean(X_mh_final[1,(10^4+1):10^5])
mean(X_mh_final[2,(10^4+1):10^5])

#numerical integration for finding means
IntY <- function(x) sapply(x, function(b) integrate(function(y) f(b,y),-Inf,Inf)$value)
IntX <- function(y) sapply(y, function(b) integrate(function(x) f(x,b),-Inf,Inf)$value)
k <- 1/integrate(function(x) IntY(x), -Inf, Inf)$value

k*integrate(function(x) x*IntY(x), -Inf, Inf)$value
k*integrate(function(y) y*IntX(y), -Inf, Inf)$value

# adaptive proposal algorithm
c_2 <- 2.4/sqrt(2)
AP <- function(X0=c(0,3), H, U, nsteps, f){
  
  X <- matrix(NA,nrow=2, ncol=nsteps+H)
  X[,1] <- X0
  X[,2:(H)] <- MetrHastw_cov(X0, 7, 2, 0.22, H-1, f)
  
  Kt <- X[,1:H]
  Kt_centered <- Kt - apply(Kt, 1, mean)
  
  for (i in (H+1):(nsteps+H)){
    
    Y <- X[,i-1] +  c_2/sqrt(H-1)* Kt_centered  %*% rnorm(H)
    
    if (runif(1)<=min(f(Y[1], Y[2])/f(X[1,i-1], X[2,i-1]), 1)) X[,i] <- Y
    else X[,i] <- X[,i-1]
    
    if ((i%%U)==0){
      Kt <- X[,(i-H+1):(i)]
      Kt_centered <- Kt - apply(Kt, 1, mean)
    }
  }
  X[,(1+H):(H+nsteps)]
}

# ap multiplc chain conv diagnostics
set.seed(0)
N <- 5*10^4 # chain length
X0s <- replicate(6, rnorm(2,c(0,3), 8))
samplesX_ap <- matrix(rep(0,dim(X0s)[2]*N), nrow=dim(X0s)[2])
samplesY_ap <- matrix(rep(0,dim(X0s)[2]*N), nrow=dim(X0s)[2])
for (i in 1:dim(X0s)[2]){
  sample_ap <- AP(X0s[,i],H=200, U=200, nsteps=N, f)
  samplesX_ap[i,] <- sample_ap[1,]
  samplesY_ap[i,] <- sample_ap[2,]
}


chainsX_ap <- lapply(1:dim(X0s)[2], function(i) mcmc(samplesX_ap[i,]))
chainsY_ap <- lapply(1:dim(X0s)[2], function(i) mcmc(samplesY_ap[i,]))
X0s
traceplot(chainsX_ap, main='Traceplots for X')
traceplot(chainsY_ap, main='Traceplots for Y')
gelman.diag(mcmc.list(chainsX_ap))
gelman.diag(mcmc.list(chainsY_ap))
gelman.plot(mcmc.list(chainsX_ap), main='Gelman Plot for X')
gelman.plot(mcmc.list(chainsY_ap), main='Gelman Plot for Y')

# sample
set.seed(0)
X_AP2 <- AP(c(0,3),H=200, U=200, nsteps=10^5, f)
plot(X_AP2[1,(10^4:10^5)], X_AP2[2,(10^4:10^5)], col=rgb(red=0, green=0, blue=1, alpha=0.01), xlab='x', ylab='y', main='Joint Distribution f(x,y)', xlim=c(-25,25), ylim=c(-20,20))
acf(X_AP2[1,(10^4:10^5)], lag.max=200, main='Acf for X')
acf(X_AP2[2,(10^4:10^5)], lag.max=200, main='Acf for Y')
plot(X_AP2[1,(10^4:10^5)], type='l')
plot(X_AP2[2,(10^4:10^5)], type='l')



#effective sample size and estimates
mean(X_AP2[1,10^4:10^5])
mean(X_AP2[2,10^4:10^5])

ESS_ap2 <- c(effectiveSize(mcmc(X_AP2[1,10^4:10^5])), effectiveSize(mcmc(X_AP2[2,10^4:10^5])))
ESS_ap2
ESS_mh <- c(effectiveSize(mcmc(X_mh_final[1,10^4:10^5])), effectiveSize(mcmc(X_mh_final[2,10^4:10^5])))
ESS_mh
# time algorithms
system.time(replicate(10, AP(c(0,3),H=200, U=200, nsteps=5*10^4, f)))
system.time(replicate(10,MetrHastw_cov(c(0,3), 7 ,2.1, 0.22,5*10^4,f)))


