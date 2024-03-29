
unpack.parameters.stage2 <- function(parameters, y.data, x.data, lambda.g, xi.00=NA, P.00=NA) {
  A         <- matrix(0, 2, 7)
  A[1, 1:2] <- parameters[1:2]  
  A[1, 3:4] <- parameters[3]/2  
  A[1, 7]   <- parameters[4]     
  A[2, 1]   <- parameters[7]     
  A[2, 5]   <- parameters[6]     
  A[2, 6]   <- 1 - parameters[6] 
  A         <- t(A)
  
  H         <- matrix(0, 2, 4)
  H[1, 1  ] <- 1
  H[1, 2:3] <- -parameters[1:2] 
  H[1, 4  ] <- parameters[5]   
  H[2, 2]   <- -parameters[7]   
  H         <- t(H)

  R         <- diag(c(parameters[8]^2, parameters[9]^2)) 
  Q         <- matrix(0, 4, 4)
  Q[1, 1]   <- parameters[10]^2              
  Q[4, 4]   <- (lambda.g * parameters[10])^2 

  F <- matrix(0, 4, 4)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2] <- F[4,4] <- 1
  
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}
