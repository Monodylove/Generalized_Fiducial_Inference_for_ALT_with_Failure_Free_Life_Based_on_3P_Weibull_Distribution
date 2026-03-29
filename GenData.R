Gendata <- function(N,R,S,m,lambda,delta,b){
  k <- length(N) #baseline set as 0, k stands for stress level
  Sample <- matrix(NA,ncol = 4)
  #rho <- exp(delta-lambda)
  Eta <- numeric(k) #scale parameter
  Gammay <- numeric(k) #position parameter
  #m:shape parameter rho:accelerate factor
  for (i in 1:k) {
    Temp <- matrix(data=NA,nrow=N[i],ncol = 4)
    Eta[i] <- exp(lambda+b*phi(S[i]))
    Gammay[i] <- exp(delta+b*phi(S[i])) 
    data <- rweibull(N[i],shape = m,scale = Eta[i])
    data <- data + Gammay[i]
    Temp[seq(1,N[i]),1] <- sort(data)
    Temp[seq(1,R[i]),2] <- rep(1,R[i])
    Temp[seq(R[i]+1,N[i]),1] <- rep(sort(data)[R[i]],N[i]-R[i]) 
    Temp[seq(R[i]+1,N[i]),2] <- rep(0,N[i]-R[i])  
    Temp[seq(1,N[i]),3] <- rep(i,N[i])
    Temp[seq(1,N[i]),4] <- rep(S[i],N[i])    
    Sample <- rbind(Sample,Temp)
    Sample <- as.matrix(na.omit(Sample))
  }
  na.omit(Sample)
}