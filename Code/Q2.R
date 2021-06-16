a <- 8
g <- function(theta) 4/10*exp(-(a-theta)^2) + 6/10 * exp(-(30-theta)^2)

x <- seq(0,38,0.01)
plot(x,g(x), type='l', main='g(.) up to proportionality')

MetrHastw1 <- function(X0, sigmaprop, nsteps, f){
  X <- numeric(nsteps+1)
  X[1] <- X0
  
  for (i in 2:(nsteps+1)){
    Y <- rnorm(1,mean=X[i-1], sd=sigmaprop)
    if (runif(1)<=min(f(Y)/f(X[i-1]), 1)) X[i] <- Y
    else X[i] <- X[i-1]
  }
  X[-1]
}
# not exploring entire state space
X0s <- seq(0,35,3)
chains <- lapply(1:length(X0s), function(i) mcmc(MetrHastw1(X0s[i], 1, 10^3, g)))
traceplot(chains, ylim=c(0,36), main='Traceplots for proposal SD=1')

X0s <- c(1,17,30,40)
chains <- lapply(1:length(X0s), function(i) mcmc(MetrHastw1(X0s[i], 20, 10^3, g)))
traceplot(chains, main='Traceplot for proposal SD=20')

# parallel tempering algorithm
g_temp <- function(theta,T) g(theta)^(1/T)
H <- function(x) -log(g(x))
# swapping probabilities
q_B <- function(l,m,M) {
  if(2<=l & l<M & (m==(l-1) | m==l+1)) return(0.5)
  if((l==1 & m==2) | (l==M & m==M-1)) return(1)
  else return(0)
}
rq_B <- function(l,M){
  if(l==1) return(l+1)
  if(l==M) return(l-1)
  if(2<=l & l<M) return(ifelse(runif(1)<1/2, l-1, l+1))
}

tempered <- function(Ts,x0s, sigmaprops, nsteps){
  M <- length(Ts)
  X <- matrix(rep(0,M*(nsteps+1)), nrow=M, ncol=nsteps+1)
  X[,1] <- x0s
  accepted_mh <- 0
  accepted_sw <- 0
  # prob accept swap
  for (i in 2:nsteps+1){
    #proposal
    prop <- X[,i-1] + rnorm(M,mean=rep(0,M), sd=sigmaprops)
    U <- runif(M)
    ind <- c(U<= sapply(1:M, function(j) exp(log(g_temp(prop[j],Ts[j]))- log(g_temp(X[j,i-1],Ts[j])))))
    # mh update for each density
    X[ind,i] <- prop[ind]
    X[!ind,i] <- X[!ind, i-1]
    accepted_mh <- accepted_mh + sum(ind)
    # between moves
    for(l in sample(1:M)){
      m <- rq_B(l,M)
      if (runif(1)<=min(1, exp((H(X[l,i])-H(X[m,i]))*(1/Ts[l]-1/Ts[m])))) {
        X[c(l,m),i] <- X[c(m,l),i]
        accepted_sw <- accepted_sw + 1}
    }
  }
  print(accepted_mh/(nsteps*M))
  print(accepted_sw/(nsteps*M))
  return(X[,-1])
}

#plot tempered distributions
Ts <- c(1,2,5,10,20, 50, 100, 200)
x <- seq(-10,50,0.01)
k <- 1/integrate(g,lower=-Inf,upper=Inf)$value
plot(x,g(x)*k, type='l', main='Tempered Distributions', xlim=c(-10,50), ylim=c(0,.4), ylab=expression(g[T](x)))
for(i in 1:length(Ts)) lines(x, g_temp(x,Ts[i])/integrate(function(y) g_temp(y,Ts[i]), lower=-Inf, upper=Inf)$value, col=i)
legend('topright', legend=c("T=1", "T=2","T=5","T=10","T=20","T=50","T=100","T=200"),
       col=1:(length(Ts)), lty=1, cex=0.8)

# multiple chains convergence diagnostics
X0s <- replicate(6, rnorm(length(Tsq), 0 ,10))
chains <- lapply(1:6, function(i) mcmc(tempered(Tsq, X0s[,i], c(0.3, 0.6, 1, 1.6, 4, 9,14), 10^4)[1,]))
gelman.diag(mcmc.list(chains))
traceplot(chains, main='Traceplots')
gelman.plot(mcmc.list(chains), main='Gelman Plot')
# sample
set.seed(0)
Tsq <- c(1,2, 5, 10,20, 50, 100)
X_tempq <- tempered(Tsq, rep(0,length(Tsq)), c(0.3, 0.6, 1, 1.6, 4, 9,14), 10^5)
hist(X_tempq[1,10^3:10^5], breaks=150, main='Histogram of the sample', xlab='x')
plot(X_tempq[1,10^3:10^5], type='l', main='Traceplot', ylab=expression(X[t]))
acf(X_tempq[1,10^3:10^5], lag.max=100, main='ACF')

#plot true density against estimated density
x <- seq(0,35,0.1)
plot(density(X_tempq[1,10^3:10^5], adjust=0.1 ),ylim=c(0,1), main='True density and estimated density')
lines(x, g(x)*k, col='red')
legend('topright', legend=c("Estimated", "True"),
       col=c('black', 'red'), lty=1, cex=0.8)

# estimate and numerical integration
mean(X_tempq[1,10^4:10^5]>20)
integrate(function(x) k*g(x), lower=20, upper=Inf)


