ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(c("MASS"))

dqpois <- Vectorize(function(x, lambda, alpha = 0, log = FALSE){
  if(x <= 150){
    logP <- -log(factorial(x)) +
      x*log(lambda/(1 + alpha*lambda)) +
      (x - 1)*log(1 + alpha*x) - 
      lambda*(1 + alpha*x)/(1 + alpha*lambda)
  }else{
    logP <- -0.5*log(2*pi*x) - x*log(x) + x + 
      x*log(lambda/(1 + alpha*lambda)) +
      (x - 1)*log(1 + alpha*x) - 
      lambda*(1 + alpha*x)/(1 + alpha*lambda)
  }
  if(log){
    return(logP)
  }else{
    return(exp(logP))
  }
}, c("x", "lambda", "alpha"))
pqpois <- Vectorize(function(q, lambda, alpha = 0, lower.tail = TRUE, log.p = FALSE){
  prob <- sum(dqpois(0:q, lambda, alpha, FALSE))
  if(!lower.tail) prob <- 1 - prob
  if(log.p) prob <- log(prob)
  return(prob)
}, c("q", "lambda", "alpha"))
rqpois <- Vectorize(function(n, lambda, alpha = 0, max.samp = 5000){
  sample(x = 0:max.samp, size = n, replace = T, prob = dqpois(0:max.samp, lambda, alpha))
}, "lambda", "alpha")

rgpm <- function(n ,mu, theta) rnbinom(n, mu/(theta-1), 1/theta)
dgpm <- function(n ,mu, theta) dnbinom(n, mu/(theta-1), 1/theta)
pgpm <- function(n ,mu, theta) pnbinom(n, mu/(theta-1), 1/theta)
qgpm <- function(n ,mu, theta) qnbinom(n, mu/(theta-1), 1/theta)