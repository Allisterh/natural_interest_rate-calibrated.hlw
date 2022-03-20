
unpack.parameters.stage3 <- function(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00) {
  A         <- matrix(0, 2, 6)   
  A[1, 1:2] <- parameters[1:2]   
  A[1, 3:4] <- parameters[3]/2 
  A[2, 1]   <- parameters[5]    
  A[2, 5]   <- parameters[4]   
  A[2, 6]   <- 1 - parameters[4]
  A         <- t(A)

  
  H         <- matrix(0, 2, 7)
  H[1, 1]   <- 1
  H[1, 2:3] <- -parameters[1:2]      
  H[1, 4:5] <- -parameters[3] * 2    
  H[1, 6:7] <- -parameters[3]/2     
  H[2, 2]   <- -parameters[5]       
  H         <- t(H)

  R         <- diag(c(parameters[6]^2, parameters[7]^2)) 

  Q         <- matrix(0, 7, 7)
  Q[1, 1]   <- (1+lambda.g^2)*parameters[8]^2                  
  Q[1, 4]   <- Q[4, 1] <- Q[4, 4]<- (lambda.g*parameters[8])^2 
  Q[6, 6]   <- (lambda.z*parameters[6]/parameters[3])^2        

  F <- matrix(0, 7, 7)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2] <- F[4,4] <- F[5,4]<- F[6,6] <- F[7,6] <- 1
 
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}



