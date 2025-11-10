library(lava)
library(mets)
library(targeted)
future::plan("multicore")
progressr::handlers(global=TRUE)

################################
#Function for simulating data###
################################

simdata<-function(n, # sample-size
                  mu, sigma, 
                  lambda, tau, gamma,
                  trteff1, trteff2
) {
  A <- rbinom(n, 1, 0.5) # Randomized treatment
  Times <- rexp(n, rate = lambda + trteff1 * A) # Time to terminal event
  R <- (Times <= tau) # Terminal event before landmark time
  Y <- (1-R) * rnorm(n, mean = mu + trteff2 * A, sd = sigma) # Outcome among those still alive
  Ytilde <- (1-R) * Y + gamma * R # Utility
  return(data.frame(A, R, Y, Ytilde))
}

##############################
#Estimating target parms######
##############################

est <- function(dat) {
  e1 <- lm(I(1-R) ~ A, data=dat) |>
    lava::estimate(keep="A", labels="theta1") # E(R=0|A=1)
  e2 <- lm(Y ~ A * R, data=dat) |>
    lava::estimate(keep="A", labels="theta2")  # E(Y|R=0, A=1)
  e3 <- lm(Ytilde ~ A, data=dat) |>
    lava::estimate(keep="A", labels="theta3") # E(~Y|A=1)
  merge(e1, e2, e3)
}


#########################################
#Setting up one-run sims#################
#########################################

onerun <- function(..., n, mu, lambda, tau,
                   trteff1, trteff2, sigma, gamma) {
  
  tmp<-simdata(n,mu, sigma, lambda, tau, 
               gamma,trteff1, trteff2)
  esttmp<-est(tmp)
  
  pvals.adj.sw.w1<-lava::closed_testing(
    esttmp,
    test = targeted::test_intersection_sw,
    noninf = rep(0, 3),
    weights = rep(1, 3)/3
  )$p.value
  names(pvals.adj.sw.w1)<-c("theta1.sw.w1","theta2.sw.w1","theta3.sw.w1")
  
  pvals.adj.sw.w2<-lava::closed_testing(
    esttmp,
    test = targeted::test_intersection_sw,
    noninf = rep(0, 3),
    weights = c(0.2,0.4,0.4)
  )$p.value
  names(pvals.adj.sw.w2)<-c("theta1.sw.w2","theta2.sw.w2","theta3.sw.w2")
  
  pvals.adj.zmax<-lava::closed_testing(
    esttmp,
    test = targeted::test_zmax_onesided,
    noninf = rep(0, 3)
  )$p.value
  names(pvals.adj.zmax)<-c("theta1.zmax","theta2.zmax","theta3.zmax")
  
  return(c(pvals.adj.sw.w1,pvals.adj.sw.w2,pvals.adj.zmax))
}

##########################################################
#Summarizing simulation results for each scenario#########
##########################################################

summary.simres<-function(alpha,simres,n){
 
 type=c("Equal weights","Up-weighted H2, H3","Minimum p-value")
 Sample_size<-rep(n,3)
 H1<-with(simres,c(mean(theta1.sw.w1<alpha),mean(theta1.sw.w2<alpha),mean(theta1.zmax<alpha)))
 H2<-with(simres,c(mean(theta2.sw.w1<alpha),mean(theta2.sw.w2<alpha),mean(theta2.zmax<alpha)))
 H3<-with(simres,c(mean(theta3.sw.w1<alpha),mean(theta3.sw.w2<alpha),mean(theta3.zmax<alpha)))
 H12<-with(simres,c(mean(theta1.sw.w1<alpha & theta2.sw.w1<alpha),
                    mean(theta1.sw.w2<alpha & theta2.sw.w2<alpha),
                    mean(theta1.zmax<alpha &theta2.zmax<alpha)))
 H13<-with(simres,c(mean(theta1.sw.w1<alpha & theta3.sw.w1<alpha),
                    mean(theta1.sw.w2<alpha & theta3.sw.w2<alpha),
                    mean(theta1.zmax<alpha &theta3.zmax<alpha)))
 H23<-with(simres,c(mean(theta2.sw.w1<alpha & theta3.sw.w1<alpha),
                    mean(theta2.sw.w2<alpha & theta3.sw.w2<alpha),
                    mean(theta2.zmax<alpha &theta3.zmax<alpha)))
 H123<-with(simres,c(mean(theta1.sw.w1<alpha & theta2.sw.w1<alpha & theta3.sw.w1<alpha),
                     mean(theta1.sw.w2<alpha & theta2.sw.w2<alpha & theta3.sw.w2<alpha),
                     mean(theta1.zmax<alpha &theta2.zmax<alpha &theta3.zmax<alpha)))
 
data.frame(Sample_size,type,H1,H2,H3,H12,H13,H23,H123) 
 
}

###############################
#Global parameters#############

set.seed(123456) #seed

Nsim <- 10000 #number of simulated datasets per scenario
alpha <-0.025 #Nominal significance level
mu <- 40 #mean of Y among those still alive in placebo arm
sigma<-15 #standard deviation of Y among those still alive 
lambda <- 0.07 #rate of terminal event in placebo arm
tau <- 2 #landmark time
trteff1 <- -0.018 #treatment effect terminal event
trteff2 <- 2.7 #treatment effect in Y among those still alive
gamma <- 15 #composite value

###############################
#simulating scenarios##########
###############################


sim_args<-list(n=200, mu=mu, lambda=lambda, tau=tau,
trteff1=trteff1, trteff2=trteff2, sigma=sigma, gamma=gamma) #Creating list of scenario specific parameters

#Sample size 200

sim_args$n<-200

res200 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres200<-summary.simres(alpha,res200,sim_args[[1]]) #summarizing simulations

#Sample size 500

sim_args$n<-500

res500 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres500<-summary.simres(alpha,res500,sim_args[[1]]) #summarizing simulations

#Sample size 1000

sim_args$n<-1000

res1000 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres1000<-summary.simres(alpha,res1000,sim_args[[1]]) #summarizing simulations

#Sample size 2000

sim_args$n<-2000

res2000 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres2000<-summary.simres(alpha,res2000,sim_args[[1]]) #summarizing simulations


#Sample size 3500

sim_args$n<-3500

res3500 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres3500<-summary.simres(alpha,res3500,sim_args[[1]]) #summarizing simulations

#Merging results
results<-rbind(sumres200,sumres500,sumres1000,sumres2000,sumres3500)

#Creating *tex table

library(knitr)

knitr::kable(results, digits=4, format="latex", booktabs=TRUE)

#########################################################
#Resulting table#########################################

\begin{tabular}{rlrrrrrrr}
\toprule
Sample\_size & type & H1 & H2 & H3 & H12 & H13 & H23 & H123\\
\midrule
200 & Equal weights & 0.0841 & 0.1692 & 0.1996 & 0.0250 & 0.0618 & 0.1456 & 0.0250\\
200 & Up-weighted H2, H3 & 0.0408 & 0.2033 & 0.2465 & 0.0250 & 0.0407 & 0.1794 & 0.0250\\
200 & Minimum p-value & 0.0678 & 0.1528 & 0.1996 & 0.0216 & 0.0473 & 0.1274 & 0.0216\\
500 & Equal weights & 0.1812 & 0.4071 & 0.4706 & 0.0920 & 0.1592 & 0.3788 & 0.0920\\
500 & Up-weighted H2, H3 & 0.1264 & 0.4502 & 0.5371 & 0.0920 & 0.1264 & 0.4268 & 0.0920\\
\addlinespace
500 & Minimum p-value & 0.1503 & 0.3693 & 0.4705 & 0.0835 & 0.1324 & 0.3414 & 0.0835\\
1000 & Equal weights & 0.3500 & 0.7291 & 0.8055 & 0.2727 & 0.3403 & 0.7152 & 0.2727\\
1000 & Up-weighted H2, H3 & 0.3169 & 0.7529 & 0.8438 & 0.2727 & 0.3169 & 0.7423 & 0.2727\\
1000 & Minimum p-value & 0.3211 & 0.6871 & 0.8027 & 0.2624 & 0.3132 & 0.6735 & 0.2624\\
2000 & Equal weights & 0.5998 & 0.9578 & 0.9807 & 0.5774 & 0.5986 & 0.9556 & 0.5774\\
\addlinespace
2000 & Up-weighted H2, H3 & 0.5934 & 0.9617 & 0.9875 & 0.5774 & 0.5934 & 0.9608 & 0.5774\\
2000 & Minimum p-value & 0.5938 & 0.9499 & 0.9802 & 0.5750 & 0.5929 & 0.9474 & 0.5750\\
3500 & Equal weights & 0.8421 & 0.9991 & 1.0000 & 0.8414 & 0.8421 & 0.9991 & 0.8414\\
3500 & Up-weighted H2, H3 & 0.8420 & 0.9991 & 1.0000 & 0.8414 & 0.8420 & 0.9991 & 0.8414\\
3500 & Minimum p-value & 0.8417 & 0.9986 & 1.0000 & 0.8411 & 0.8417 & 0.9986 & 0.8411\\
\bottomrule
\end{tabular}

