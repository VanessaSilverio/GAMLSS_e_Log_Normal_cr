LOGNOcr <- function (mu.link="identity", sigma.link="log",nu.link="logit") {
  mstats <- checklink("mu.link", "Log Normal cr", substitute(mu.link), c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "Log Normal cr", substitute(sigma.link), c("inverse", "log", "identity", "own"))
  vstats <- checklink("nu.link", "Log Normal cr", substitute(nu.link),  c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  
  structure(
    list(family = c("LOGNOcr", "Log Normal cure rate"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE),
         nopar = 3, 
         type = "Continuous",
         
         nu.link = as.character(substitute(nu.link)), 
         nu.linkfun = vstats$linkfun,
         nu.linkinv = vstats$linkinv,
         nu.dr = vstats$mu.eta,
         
         mu.link = as.character(substitute(mu.link)), 
         sigma.link = as.character(substitute(sigma.link)), 
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         
         
         dldm = function(y,mu,sigma) {dldm <- (log(y)-mu)/sigma^2;  dldm},
         d2ldm2 = function(sigma) -1/sigma^2,
         dldd = function(y,mu,sigma)  {dldd <- (1/(sigma^3))*((log(y)-mu)^2-sigma^2); dldd},
         d2ldd2 = function(sigma) -2/sigma^2,
         d2ldmdd = function(y)  rep(0,length(y)),
         
         dldv = function(y,nu) -1/(1-nu),
         d2ldv2 = function(y,nu) -1/(1-nu)^2,
         d2ldmdv = function(y,mu,sigma,nu){
           nd = gamlss:::numeric.deriv(dLOGNOcr(y, mu, sigma, nu, log = TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd, "gradient"))
           d2ldmdv = -dldm * (-1/(1-nu))
           d2ldmdv}, 
         d2ldddv = function(y,mu,sigma,nu){
           nd = gamlss:::numeric.deriv(dLOGNOcr(y, mu, sigma,nu, log = TRUE), "sigma", delta = 1e-04)
           dldm = as.vector(attr(nd, "gradient"))
           d2ldmdv = -dldm * (-1/(1-nu))
           d2ldmdv}, 
         
         
         
         
         G.dev.incr  = function(y,mu,sigma,nu,...) -2*dLOGNOcr(x=y, mu=mu,  sigma = sigma, nu=nu,log = TRUE), #
         rqres = expression(  rqres(pfun="pLOGNOcr", type="Continuous", y=y, mu=mu, sigma=sigma,nu=nu)),
         mu.initial = expression({    mu <- (log(y)+mean(log(y)))/2  }),
         sigma.initial = expression({sigma <- rep(sd(log(y)),length(y)) }),  
         mu.valid = function(mu) all(mu > 0), 
         sigma.valid = function(sigma)  all(sigma > 0),
         nu.initial = expression(nu <-rep(0.5, length(y))),
         nu.valid = function(nu) all(nu >= 0.001) && all(nu < .99),
         
         y.valid = function(y)  all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------

dLOGNOcr <- function(x, mu=0, sigma=1,nu=0,  log = FALSE)
{
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
  if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
  fy1 <- (1-nu)*dlnorm(x=x, meanlog = mu, sdlog = sigma) 
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy}

#----------------------------------------------------------------------------------------
pLOGNOcr <- function(q, mu=0, sigma=1,nu=0,  lower.tail = TRUE, log.p = FALSE)
{  
  #      if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
  if (any(q < 0) )  stop(paste("y must be >=0", "\n", ""))  
  cdf1 <- (1-nu)*plnorm(q=q, meanlog = mu, sdlog = sigma)
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
  
}
#----------------------------------------------------------------------------------------
qLOGNOcr <- function(p, mu=0, sigma=1,nu=0,  lower.tail = TRUE, log.p = FALSE )
{ 
  #if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
  if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (any(p < nu))  stop(paste("p must be less than nu", "\n", ""))
  
  q <-qlnorm(p=p/(1-nu), meanlog = mu, sdlog = sigma, lower.tail = lower.tail, log.p = log.p)
  q
}
#----------------------------------------------------------------------------------------
rLOGNOcr <- function(n, mu=0, sigma=1,nu=0)
{
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
  
  n2=round(nu*n,0) #cure
  n1=n-n2 #failure
  
  uni1<- runif(n = n1,0,1) # failure
  r1 <- rlnorm(uni1,mu,sigma) #failure
  
  dp<-sd(r1);
  r2 <- runif(n2,max(r1)+0.5*dp,max(r1)+2*dp) #cure
  
  Response<<-c(r1,r2)
  
  delta1<-rep(1,n1); #failure
  delta2<-rep(0,n2) #censure
  Delta<<-c(delta1,delta2)
  print("Varibles Response and Delta(1=failure, 0=censored) has been successfully generated")
}
#----------------------------------------------------------------------------------------
