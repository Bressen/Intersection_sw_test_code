############################################################
#R-script for producing Figure 2 in paper###################
############################################################

#############################################################
#critvals and p-values
#############################################################
library(mvtnorm)

qprob<-function(corr){
  Sigma<-matrix(c(1,corr,corr,1),2,2)
  0.5-pmvnorm(upper=c(0,0),mean=c(0,0),corr=Sigma)[1]
}


prob.fct<-function(c,alpha,corr){
  q<-qprob(corr)
  0.5*pchisq(c,1,lower.tail=F)+q*pchisq(c,2,lower.tail=F)-alpha
}



q.fct<-function(alpha,corr){
  uniroot(prob.fct,alpha=alpha,corr=corr,interval=c(0,10))$root
}


###############################################################
#Calculating test statistics###################################
#Calculating p-values##########################################
###############################################################



testvals<-function(thetahat1,se1,thetahat2,se2,noninf1,noninf2,corr,alpha){
  z1<-(thetahat1-noninf1)/se1
  z2<-(thetahat2-noninf2)/se2
  zmin<-min(z1,z2)
  zmax<-max(z1,z2)
  SignWald.intersect<-ifelse(zmax>=0 & zmin<=(corr*zmax),1,0)*zmax*zmax+
    ifelse(zmax>=0 & zmin>(corr*zmax),(zmax*zmax+zmin*zmin-2*corr*zmax*zmin)/(1-corr*corr),0)
  critval.intersect<-q.fct(alpha,corr)
  pval.intersect<-ifelse(SignWald.intersect>0,prob.fct(SignWald.intersect,alpha,corr)+alpha,1)
}

zmax.helper<-function(c,frac,power,alpha,corr){
  critval.zmax<-qmvnorm(1-alpha,sigma=matrix(c(1,corr,corr,1),2,2))$quantile
  1-pmvnorm(upper=rep(critval.zmax,2),mean=c(frac*c,c),sigma=matrix(c(1,corr,corr,1),2,2))[1]-power
}

zmax.noncent<-function(frac,power,alpha, corr){
  uniroot(zmax.helper,frac=frac,power=power, alpha=alpha,corr=corr,interval=c(0,10))$root
}

powersim<-function(zmax,zmin,corr,Nsim,alpha){
  thetahat.std<-rmvnorm(Nsim,mean=c(zmin,zmax),sigma=matrix(c(1,corr,corr,1),2,2))
  reject<-rep(0,Nsim)
  reject.zmax<-rep(0,Nsim)
  critval.zmax<-qmvnorm(1-alpha,sigma=matrix(c(1,corr,corr,1),2,2))$quantile
  for(i in 1:Nsim){
    reject[i]<-ifelse(testvals(thetahat.std[i,1],1,thetahat.std[i,2],1,0,0,corr,alpha)<alpha,1,0)
    reject.zmax[i]<-ifelse(max(thetahat.std[i,])>critval.zmax,1,0)
  }
  return(c(mean(reject),mean(reject.zmax)))
}


#############################################################################
#Simulation 1################################################################
#for a range of correlations#################################################
#for each correlation zmax chosen to yield a certain power of min-p test#####
#zmin chosen to be less than zmax############################################
#############################################################################

alpha<-0.025
power<-0.9
corr.range<-c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75)
frac.range<-seq(-1,1,by=0.05)
zmax.range<-matrix(0,length(frac.range),length(corr.range))
zmin.range<-matrix(0,length(frac.range),length(corr.range))
for (i in 1:length(corr.range)){
 for (j in 1:length(frac.range)){
   zmax.range[j,i]<-zmax.noncent(frac.range[j],power,alpha,corr.range[i])   

   } 
}

zmin.range<-zmax.range*frac.range

powersim.range<-matrix(0,length(frac.range),length(corr.range))

set.seed(12345)
for(i in 1:length(corr.range)){
  for(j in 1:length(frac.range)){
powersim.range[j,i]<-powersim(zmax.range[j,i],zmin.range[j,i],corr.range[i],100000,alpha)[1]
print(c(i,j))

  }
}

plot(x=frac.range,y=powersim.range[,1],type="l",xlim=c(-1,1),ylim=c(0.85,1),col=1,lwd=2,ylab="Power",xlab=expression(z[min]/z[max]))
abline(h=0.90,lty=2)
points(x=frac.range,y=powersim.range[,2],type="l",col=2,lwd=2)
points(x=frac.range,y=powersim.range[,3],type="l",col=3,lwd=2)
points(x=frac.range,y=powersim.range[,4],type="l",col=4,lwd=2)
points(x=frac.range,y=powersim.range[,5],type="l",col=5,lwd=2)
points(x=frac.range,y=powersim.range[,6],type="l",col=6,lwd=2)
points(x=frac.range,y=powersim.range[,7],type="l",col=7,lwd=2)
legend("topleft",legend=c("Correlation=-0.75", "Correlation=-0.50","Correlation=-0.25","Correlation=0","Correlation=0.25","Correlation=0.50","Correlation=0.75"),text.col=1:7)

write.csv(powersim.range, file.choose())
