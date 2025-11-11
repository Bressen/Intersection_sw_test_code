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

Sigmahat <- matrix(c(0.403,  -1.560,  4.480,
                  -1.560, 846.241, 822.862,
                  4.480, 822.862, 890.064),3,3)

vcov.thetahat <- Sigmahat/3532 #3532 out of 3533 randomized subjects included in analysis. 

#########################################################
#Creating estimate object################################
#########################################################

est_obj <- lava::estimate(coef=thetahat,vcov=vcov.thetahat)

#########################################################
#Signed wald test equal weights #########################
#########################################################

sw_equalweights <- lava::closed_testing(
  est_obj,
  test = targeted::test_intersection_sw,
  noninf = rep(0, 3),
  weights = rep(1, 3)/3,
  nsim.null=100000
)
sw_equalweights$p.value
summary(sw_equalweights)

#########################################################
#Signed wald test upweighted H2,H3 ######################
#########################################################

sw_upweighted <- lava::closed_testing(
  est_obj,
  test = targeted::test_intersection_sw,
  noninf = rep(0, 3),
  weights = c(0.2,0.4,0.4),
  nsim.null=100000
)
sw_upweighted$p.value
summary(sw_upweighted)

#########################################################
#P-values: Minimum p-value test##########################
#########################################################

p_min <- lava::closed_testing(
  est_obj,
  test = targeted::test_zmax_onesided,
  noninf = rep(0, 3)
)
p_min$p.value
summary(p_min)
