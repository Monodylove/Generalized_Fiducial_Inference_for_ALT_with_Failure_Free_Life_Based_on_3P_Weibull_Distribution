phi <- function(x){
  1/x #Arrhenuis
  #x #log-linear
}

term1 <- function(ts,delta,b){
  t<- ts[1]
  S<- ts[2]
  t-exp(delta+b*phi(S))
}

term2 <- function(ts,lambda,b){
  t<- ts[1]
  S<- ts[2]
  exp(lambda+b*phi(S))
}

f <- function(ts,m,delta,lambda,b){
  temp1 <- term1(ts,delta,b)
  temp2 <- term2(ts,lambda,b)
  value1 <- m*((temp1^(m-1))/(temp2^m))*exp(-(temp1/temp2)^m)#(m-1可能包含非负约束)
  if (is.na(value1)) {
    return(exp(-100))
  } 
  else{return(value1)}
}

R <- function(ts,m,delta,lambda,b){
  temp1 <- term1(ts,delta,b)
  temp2 <- term2(ts,lambda,b)
  value2 <- exp(-(temp1/temp2)^m)
  if (is.na(value2)) {
    return(exp(-100))
    # return(1)    
  }
  else{return(value2)}
}

lnf <- function(ts, m, delta, lambda, b) {
  value3 <- f(ts, m, delta, lambda, b)
  return(log(value3))
}

lnR <- function(ts, m, delta, lambda, b) {
  value4 <- R(ts, m, delta, lambda, b)
  return(log(value4))
}

JL2 <- function(SAMPLE, m, delta, lambda, b) {
  n <- nrow(SAMPLE)
  A2 <- matrix(0, n, n)
  B <- matrix(0, n, 4)
  C <- matrix(0, n, 4)
  
  for (i in 1:n) {
    ts <- SAMPLE[i, c(1,4)]
    status <- SAMPLE[i, 2]
    
    if (status == 1) {
      A2[i, i] <- 1 / ( f(ts,m,delta,lambda,b)^2 )
    } else {
      A2[i, i] <- 1 / ( R(ts,m,delta,lambda,b)^2 )
    }

    t_val <- SAMPLE[i, 1]  
    S_val <- SAMPLE[i, 4] 
    R_val <- R(ts,m,delta,lambda,b)
    B0 <- m * R_val * log(R_val)
    B[i, ] <- rep(B0, 4)
    
    C[i, 1] <- - log(term1(ts,delta,b) / term2(ts,lambda,b)) / m
    C[i, 2] <- (t_val - term1(ts,delta,b)) / term1(ts,delta,b)
    C[i, 3] <- 1
    C[i, 4] <- - t_val * phi(S_val) / term1(ts,delta,b)
  }

  JM <- (t(B) * t(C)) %*% A2 %*% (B * C)
  L2 <- sqrt(det(JM))
  return(L2)
}

lnL <- function(SAMPLE,params){
  Fail_S <- SAMPLE[which(SAMPLE[,2]==1),c(1,4)]
  Sum_S <- SAMPLE[,c(1,4)]
  m <- params[1]
  delta <- params[2]
  lambda <- params[3]
  b <- params[4]
  fail_part <- sum(apply(Fail_S, 1, lnf, m=m, delta=delta, lambda=lambda, b=b))
  sum_part <- sum(apply(Sum_S, 1, lnR, m=m, delta=delta, lambda=lambda, b=b))
  J <- JL2(SAMPLE, m, delta, lambda, b)
  return(fail_part + sum_part + log(J))
}
