##########################################################
#Analysis of FLOW data####################################
##########################################################

library(targeted)


#########################################################
#Estimates of target parameters##########################
#########################################################

thetahat<-c(0.0315,2.681,3.153)

#########################################################
#Estimated variance######################################
#########################################################

Sigmahat<-matrix(c(0.403,  -1.560,  4.480,
                  -1.560, 846.241, 822.862,
                  4.480, 822.862, 890.064),3,3)

vcov.thetahat<-Sigmahat/3532 #3532 out of 3533 randomized subjects included in analysis. 

#########################################################
#P-values: signed wald test##############################
#########################################################

data.frame(pvals123=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=rep(1,3)/3)$p.value,
pvals12=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=rep(1,3)/3,index=c(1,2))$p.value,
pvals13=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=rep(1,3)/3,index=c(1,3))$p.value,
pvals23=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=rep(1,3)/3,index=c(2,3))$p.value,
pvals123.1=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=c(0.2,0.4,0.4))$p.value,
pvals12.1=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=c(0.2,0.4,0.4),index=c(1,2))$p.value,
pvals13.1=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=c(0.2,0.4,0.4),index=c(1,3))$p.value,
pvals1=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=rep(1,3)/3,index=c(1))$p.value,
pvals2=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=rep(1,3)/3,index=c(2))$p.value,
pvals3=targeted::test_intersection_sw(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),weights=rep(1,3)/3,index=c(3))$p.value)

#########################################################
#P-values: Minimum p-value test##########################
#########################################################


data.frame(pvals123=targeted::test_zmax_onesided(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3))$p.value,
           pvals12=targeted::test_zmax_onesided(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),index=c(1,2))$p.value,
           pvals13=targeted::test_zmax_onesided(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),index=c(1,3))$p.value,
           pvals23=targeted::test_zmax_onesided(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),index=c(2,3))$p.value,
           pvals1=targeted::test_zmax_onesided(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),index=c(1))$p.value,
           pvals2=targeted::test_zmax_onesided(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),index=c(2))$p.value,
           pvals3=targeted::test_zmax_onesided(par=thetahat,vcov=vcov.thetahat,noninf=rep(0,3),index=c(3))$p.value)

