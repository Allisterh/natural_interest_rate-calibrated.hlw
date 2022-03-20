
unpack.parameters.stage1 <- function(parameters, y.data, x.data, xi.00, P.00) {
  A         <- matrix(0, 4, 2)
  A[1:2, 1] <- parameters[1:2]
  A[1, 2]   <- parameters[4] 
  A[3, 2]   <- parameters[3]  
  A[4, 2]   <- 1-parameters[3]  
  
  H         <- matrix(0, 3, 2)  
  H[1, 1]   <- 1
  H[2:3, 1] <- -parameters[1:2] 
  H[2, 2]   <- -parameters[4]  

  R         <- diag(c(parameters[6]^2, parameters[7]^2)) 
  Q         <- matrix(0, 3, 3)
  Q[1, 1]   <- parameters[8]^2 

  F <- matrix(0, 3, 3)
  F[1, 1] <- F[2, 1] <- F[3, 2] <- 1

  y.data[, 1] <- y.data[, 1] - 1:dim(y.data)[1] * parameters[5] # g
  x.data[, 1] <- x.data[, 1] - 0:(dim(x.data)[1]-1) * parameters[5]
  x.data[, 2] <- x.data[, 2] - -1:(dim(x.data)[1]-2) * parameters[5]
      
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}
