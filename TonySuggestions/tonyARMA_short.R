ARMApqREMLfunct <- function(par,nX,p,q){  
  
##generates LL for ARMApq model.
##translated from matlab function ARMApqREMLfunct.m by Nic Ziebarth and Tony Ives
##10 March, 2009

##par should have values for both AR and MA components.
##nX is the data, p = # of AR, q = # of MA

##reoccuring issue is what should diagExtend output when size = 0. This comes up
##when p=1 and q=1. I think any 1 x 1 matrix will work. 

## this is standardized so that a(1)=1 and is not estimated
b <- par[1:p];
lengthPar <- length(par);

##this slight change relative to Matlab case. If q=0 in old version, then a =[], which ##was no good
if(q > 0){ 
   a <- c(1,par[(p+1):lengthPar]);
   }
else{
   a <- 1
   }

B <- diagExtend(p,-1);
sizeB <- dim(B);
B[1,1:p] <- b;
eigs <- eigen(B)$values;
eigB <- max(abs(eigs));
 
if (eigB>=1 | max(abs(a))>10){
   LL<-10^10;
   }
    
T <- length(nX);
    
## solve for stationary distribution
    
k <- max(p,q);

##deal with case of no MA terms separately
if (q > 0){
   A <- array(0,dim=c(k,q));
   A[1,1:q] <- a[2:(q+1)];
   if(k > 1){
      B <- diagExtend(k,-1);
      B[1,1:p] <- b;
      }
   else{
      B <- array(1,dim=c(1,1))*b;
      }

   CC <- cbind(B,A);

   otherPart <- rbind(array(0,dim=c(k,q)), diagExtend(q,1));
   CC <- rbind(CC,t(otherPart));
 
   Ve <- array(0,dim=c(k+q,k+q));
   Ve[1,1] <- 1;
   Ve[1,k+1] <- 1;
   Ve[k+1,1] <- 1;
   Ve[k+1,k+1] <- 1;

   sizeVe <- dim(Ve);
   vecVe <- array(Ve,dim=c(sizeVe[1]*sizeVe[2],1));
   C <- solve(diag(1,(k+q)^2) - kronecker(CC,CC))%*%vecVe;

   C <- array(C,dim=c(k+q,k+q));

   Vstationary<-C[1:k,1:k];
   }
else{

   ##this is particular example of issue when p=1
   if(p > 1){
      B <- diagExtend(p,-1);
      sizeB <- dim(B);
      B[1,1:sizeB[2]] <- b;
      }
   else{
      B <- array(1,dim=c(1,1))*b;
      }

   Ve <- array(0,dim=c(p,p));
   Ve[1,1] <- 1;
 	
   sizeVe <- dim(Ve);
   vecVe <- array(Ve,dim=c(sizeVe[1]*sizeVe[2],1));	
   C <- solve(diag(1,(p+q)^2) - kronecker(B,B))%*%vecVe;

   Vstationary <- array(C,dim=c(p+q,p+q));
}
 
## solve for V_T
    
##MA component
    
##This was added to separately deal with case of no MA lags
if(q > 0){
   AA <- diag(1,q+1);
   for (i in 1:q){
      AA <- AA + a[i+1]*diagExtend(q+1,-i);
      }
   AA <- t(AA)%*%AA;
   sizeAA <- dim(AA);

   aa <- AA[1,1:sizeAA[2]];
   AA <- array(0,dim=c(T,T));
   for (i in 1:q){
      AA <- AA + aa[i+1]*diagExtend(T,-i);
      }
 	
   AA <- AA+t(AA)+aa[1]*diag(1,T);
   }
else{
   AA <- diag(1,T);
   }

BB <- diag(1,k);
if(k > 1){ 
   for (i in 1:p){
      BB <- BB - b[i]*diagExtend(k,-i);
      }
   }
A1 <- BB%*%Vstationary%*%t(BB);
AA[1:k,1:k] <- A1;

## AR component
W <- diag(1,T);
for (i in 1:p){
   W <- W - b[i]*diagExtend(T,-i);
}
invW <- solve(W);
 

V <- invW%*%AA%*%t(invW);
    
## compute LL function for extant data
    
pick <- which(is.na(nX)==FALSE);
T <- length(pick);
X <- nX[pick];
V <- V[pick,pick];
    
invV <- solve(V);
    
U <- array(1,dim=c(T,1));
mu <- solve(t(U)%*%invV%*%U)%*%t(U)%*%invV%*%X;
H <- X - mu; 
 
## condensed ML LL
##s2 <- (t(H)%*%invV%*%H)/T;
##LL <- .5*(T*log(2*pi)+T*log(s2)+log(det(V))+T);

## condensed REML LL
s2 <- (t(H)%*%invV%*%H)/(T-1);
LL <- .5*((T-1)*log(2*pi)+(T-1)*log(s2)+log(det(V))+log(det(t(U)%*%invV%*%U))+(T-1));


if(abs(Im(LL)) > 10^-6) {
   LL <- 10^10;
   }
LL
}
## end ARMApqREMLfunct



diagExtend <- function(p,diagInd){

##This is a function with the functionality of diag in Matlab
##In particular, this generates a matrix of size p x p with ones on
##the diagInd diagonal. Restriction is abs(diagInd) <= p. Returns 0 matrix otw.

diagNew <- abs(diagInd);
if(diagNew <= p){
   firstPart <- diag(1,p,p-diagNew);
   secondPart <- array(0,dim=c(p,diagNew));

   A <- cbind(secondPart,firstPart);
   if(diagInd < 0){
      A <- t(A);
      }
   }
   else{
      A <- c(0)
   }
A
}
## end diagExtend
