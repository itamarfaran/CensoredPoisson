source("lib/R/base_functions.R")

generate.cqpois.samp <- function(n = NULL, p = NULL, X = NULL, beta = NULL,
                                 alpha = NULL, cens.prob = 0, cens.range = c(0,1)){
  if(is.null(c(n, p, X))) stop("Must specify n,p or X")
  if(is.null(X)){
    if(is.null(c(n,p))) stop("Must specify both n,p")
    X <- cbind(rep(1, n), matrix(rnorm(n*(p - 1), 1, 0.6), n, p - 1))
  }else{
    n <- nrow(X)
    if(any(X[,1] != 1)) X <- cbind(rep(1, n), X)
    p <- ncol(X)
  }
  if(is.null(beta)){
    beta <- rnorm(p, 0.5, 0.005)
  }else{
    if(length(beta) != p) stop("X & beta not of compatible lengths")
  }
  if(is.null(alpha)) alpha <- runif(1, 0, 1)
  
  if(alpha > 0){
    yclean <- rqpois(1, exp(X %*% beta), alpha)
  }else{
    yclean <- rpois(n, exp(X %*% beta))
  }
  delta <- sample(0:1, nrow(X), T, c(cens.prob, 1 - cens.prob))
  y <- yclean
  y[delta == 0] <- floor(y[delta == 0]*runif(sum(1 - delta), cens.range[1], cens.range[2]))
  
  Y <- cbind(y, delta)
  colnames(Y) <- c("y", "event")
  
  return(list(y = Y, yclean = yclean, X = X, beta = beta, alpha = alpha))
}

tmp <- generate.cqpois.samp(100, 9, cens.prob = 0.2, alpha = 0.1)
