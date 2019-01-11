#Algorithm 2.1

regression<- function(X,y,k,sigma,Xs){
  
  #Calculate K(X,X)
  VielX <- matrix(NA,nrow=length(X),ncol=length(X))
  for(i in 1:length(X)){
    Xi<-matrix(X[i],nrow=length(X),ncol=1) 
    VielX[,i] <- Xi 
  }
  K <- matrix(mapply(k,VielX,X),nrow=length(X),byrow=TRUE)
  #K calculated
  
  #Calculate k* (ks)
  XS<-rep(Xs,length(X))
  ks <- mapply(k,X,XS)
  #end
  
  #calculate all other variables directly
  L <- chol(K+sigma*diag(length(X)))
  
  alpha<-solve(L^T,solve(L,y))
  
  fs <- ks^T%*%alpha
  
  v <- solve(L,ks)
  
  Vfs <- k(xs,xs)-v^T%*%v
  
  logp <- -0.5*y^T%*%alpha-sum(log(diag(L)))-length(X)/2*log(2*pi)
  
  return(c(fs,Vfs,logp))
}  
