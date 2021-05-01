source("lib/R/base_functions.R")

family <- quasipoisson(link = "log")

minusloglik.qpois <- function(y, delta, X, beta, alpha = 0, family){
  lambda <- family$linkinv(X %*% beta)
  if(alpha != 0){
    - sum(dqpois(y[delta], lambda[delta], alpha, log = TRUE)) -
      sum(pqpois(y[!delta], lambda[!delta], alpha, lower.tail = FALSE, log.p = TRUE))
  }else{
    - sum(dpois(y[delta], lambda[delta], log = TRUE)) -
      sum(ppois(y[!delta], lambda[!delta], lower.tail = FALSE, log.p = TRUE))
  }
}

y <- tmp$y[,1]
lambda <- exp(tmp$X %*% tmp$beta)
alpha <- tmp$alpha
delta <- tmp$y[,2]
X <- tmp$X
beta <- tmp$beta

deriv.lambda <- function(y, lambda, alpha){
  (y - lambda)/(lambda*(1 + alpha*lambda)^2)
}

deriv.alpha <- function(y, lambda, alpha){
  y*(y - 1)/(1 + alpha*y) -
    y*lambda/(1 + alpha*lambda) -
    lambda*(y - lambda)/(1 + alpha*lambda)^2
}

deriv.censored <- Vectorize(function(y, lambda, alpha, deriv.by){
  -sum(deriv.by(0:y, lambda, alpha))/pqpois(y, lambda, alpha, lower.tail = FALSE)
}, c("y", "lambda", "alpha"))

gradloglik.qpois <- function(y, delta, X, beta, alpha, family){
  lambda <- family$linkinv(X %*% beta)
  c(
  t(c(deriv.lambda(y[delta], lambda[delta], alpha),
    deriv.censored(y[!delta], lambda[!delta], alpha, deriv.lambda))*
    family$mu.eta(X %*% beta)) %*% X,
  sum(deriv.alpha(y[delta], lambda[delta], alpha)) +
        sum(deriv.censored(y[!delta], lambda[!delta], alpha, deriv.alpha))
  )
}

minusloglik.qpois(y, delta, X, beta, alpha, quasipoisson())
optim(c(beta, alpha)+0.2,
      fn = function(pars) minusloglik.qpois(y, delta, X,
                                       pars[1:(length(pars) - 1)], pars[length(pars)],
                                       quasipoisson(link = "log")),
      gr = function(pars) gradloglik.qpois(y, delta, X,
                                      pars[1:(length(pars) - 1)], pars[length(pars)],
                                      quasipoisson(link = "log")),
      method = "BFGS")
alpha
