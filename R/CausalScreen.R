CausalScreen <-
function(X,P1,P2,l,k,y,z=NULL){
  
# Format X to be a matrix
  tempNames <- colnames(X)
  if(!is.null(X)){
    x<-matrix(0,nrow=nrow(X),ncol=ncol(X))
    for(pp in 1:ncol(X)){x[,pp]<-X[,pp]}
  }
  X <- x
  colnames(X) <- tempNames
  
  # Format P1 to be a matrix
  if(!is.null(P1)){
    p1<-matrix(0,nrow=nrow(P1),ncol=ncol(P1))
    for(pp in 1:ncol(P1)){p1[,pp]<-P1[,pp]}
  }
  P1 <- p1
  
  # Format P2 to be a matrix
  if(!is.null(P2)){
    p2<-matrix(0,nrow=nrow(P2),ncol=ncol(P2))
    for(pp in 1:ncol(P2)){p2[,pp]<-P2[,pp]}
  }
  P2 <- p2
  
  
CSErrorCheck(X,P1,P2,l,k,y,z)
  
###############################
# E(X|P)
############################### 
#expectation of childs genotype given parents

nSNP<-ncol(X)
n <- nrow(X)
matR<-matrix(NA,ncol=5,nrow=nSNP)
colnames(matR)<-c("snpOrder","ScreeningStepOrder","ScreenStat","Statistic","pvalue")
names(matR)<-c("snpOrder","ScreeningStepOrder","ScreenStat","Statistic","pvalue")
matR[,"snpOrder"]<-c(1:nSNP)
rownames(matR)<-colnames(X)

  # E(X|P)
  ExPM<-matrix(0,nrow=nrow(X),ncol=ncol(X))
  for(xx in 1:nSNP){
  	p1<-P1[,xx]
  	p2<-P2[,xx]
  ExP<-rep(0,n)  
  ExP[(p1==0&p2==2)|(p1==2&p2==0)]<-1
  ExP[(p1==2&p2==2)]<-2
  ExP[(p1==1&p2==0)|(p1==0&p2==1)]<-0.5
  ExP[(p1==1&p2==2)|(p1==2&p2==1)]<-1.5
  ExP[(p1==1&p2==1)]<-1
  ExPM[,xx]<-ExP
    }

###############################
# screening statistic and test statistic
###############################  

for(nn in 1:nSNP){
# Test Stat without genetic info 
  
# If z is NULL, fit model2 without z
if(is.null(z)){
  model2<-lm(y~k+l+ExPM[,nn])
}else{
  model2<-lm(y~k+l+ExPM[,nn]+z)
}
  
T2SCREEN<-  (ExPM[,nn])*(y-mean(y)-(coef(model2)[2])*(k-mean(k)))
  
# If z is NULL, fit fitK2 without z
if(is.null(z)){
  fitK2<-lm(k~l+ExPM[,nn])
}else{
  fitK2<-lm(k~l+ExPM[,nn]+z)
}

muK2<-fitK2$fitted.values
sigK2<-summary(fitK2)$sigma^2
SIG2SCREEN<- T2SCREEN-mean((ExPM[,nn])*k)*((k-muK2)/sigK2)*model2$residuals 
STATscre<-(sum(T2SCREEN))^2/(n*var(SIG2SCREEN))
#STATscre<-(sum(T2SCREEN))^2/(n*var(T2SCREEN))
matR[nn,"ScreenStat"]<-STATscre

# Test Stat with genetic info

if(is.null(z)){
  model1<-lm(y~k+l+ExPM[,nn])
}else{
  model1<-lm(y~k+l+ExPM[,nn]+z)
}

T1<-  (X[,nn]-ExPM[,nn])*(y-mean(y)-(coef(model1)[2])*(k-mean(k)))

if(is.null(z)){
  fitK1<-lm(k~l+ExPM[,nn])
}else{
  fitK1<-lm(k~l+ExPM[,nn]+z)
}

muK1<-fitK1$fitted.values
sigK1<-summary(fitK1)$sigma^2
SIG1<- T1-mean((X[,nn]-ExPM[,nn])*k)*((k-muK1)/sigK1)*model1$residuals 
STAT1<-sum(T1)^2/(n*var(SIG1))
matR[nn,"Statistic"]<-STAT1
matR[nn,"pvalue"]<-1-pchisq(STAT1,df=1)

}

###############################
# order based on screening stat
############################### 

matR2<-matR[order(matR[,"ScreenStat"],decreasing=TRUE),]
matR2[,"ScreeningStepOrder"]<-c(1:nSNP)

matR2
}
