
rstar.stage1 <- function(log.output,
                         inflation,
                         b2.constraint=NA) {

  stage <- 1
  
  T <- length(log.output) - 4

  x.og <- cbind(rep(1,T+4), 1:(T+4))
  y.og <- log.output
  output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)) * 100

  log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
  g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  xi.00 <- c(100*g.pot[3:1])

  y.is <- output.gap[5:(T+4)]
  x.is <- cbind(output.gap[4:(T+3)], output.gap[3:(T+2)])
  b.is <- solve(t(x.is) %*% x.is, t(x.is) %*% y.is)
  r.is <- y.is - x.is %*% b.is
  s.is <- sqrt(sum(r.is^2) / (length(r.is)-(dim(x.is)[2])))

  y.ph <- inflation[5:(T+4)]
  x.ph <- cbind(inflation[4:(T+3)],
                (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3,
                output.gap[4:(T+3)])                
  b.ph <- solve(t(x.ph) %*% x.ph, t(x.ph) %*% y.ph)
  r.ph <- y.ph - x.ph %*% b.ph
  s.ph <- sqrt(sum(r.ph^2) / (length(r.ph)-(dim(x.ph)[2])))
  
  y.data <- cbind(100 * log.output[5:(T+4)],
                  inflation[5:(T+4)])
  x.data <- cbind(100 * log.output[4:(T+3)],
                  100 * log.output[3:(T+2)],
                  inflation[4:(T+3)],
                  (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3)

  initial.parameters <- c(b.is, b.ph[1], b.ph[3], 0.85, s.is, s.ph, 0.5)

  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep(Inf,length(initial.parameters)))
  
  if (!is.na(b2.constraint)) {
      if (initial.parameters[4] < b2.constraint) {
          initial.parameters[4] <- b2.constraint
      }
      theta.lb[4] <- b2.constraint
  }

  P.00 <- calculate.covariance(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, NA, NA, xi.00)
  
  f <- function(theta) {return(-log.likelihood.wrapper(theta, y.data, x.data, stage, NA, NA, xi.00, P.00)$ll.cum)}
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8))
  theta <- nloptr.out$solution

  log.likelihood <- log.likelihood.wrapper(theta, y.data, x.data, stage, NA, NA, xi.00, P.00)$ll.cum

  states <- kalman.states.wrapper(theta, y.data, x.data, stage, NA, NA, xi.00, P.00)

  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)
  
  potential.smoothed  <- as.vector(states$smoothed$xi.tT[,1])/100
  output.gap.smoothed <- y.data[,1] - (potential.smoothed * 100)

  return.list                <- list()
  return.list$theta          <- theta
  return.list$log.likelihood <- log.likelihood
  return.list$states         <- states
  return.list$xi.00          <- xi.00
  return.list$P.00           <- P.00
  return.list$potential.filtered  <- potential.filtered
  return.list$output.gap.filtered <- output.gap.filtered
  return.list$potential.smoothed  <- potential.smoothed
  return.list$output.gap.smoothed <- output.gap.smoothed
  return(return.list)
}
