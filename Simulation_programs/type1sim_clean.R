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

  return(c(pvals123=targeted::test_intersection_sw(esttmp,noninf=rep(0,3),weights=rep(1,3)/3)$p.value,
  pvals12=targeted::test_intersection_sw(esttmp,noninf=rep(0,3),weights=rep(1,3)/3,index=c(1,2))$p.value,
  pvals13=targeted::test_intersection_sw(esttmp,noninf=rep(0,3),weights=rep(1,3)/3,index=c(1,3))$p.value,
  pvals23=targeted::test_intersection_sw(esttmp,noninf=rep(0,3),weights=rep(1,3)/3,index=c(2,3))$p.value,
  pvals123.1=targeted::test_intersection_sw(esttmp,noninf=rep(0,3),weights=c(0.2,0.4,0.4))$p.value,
  pvals12.1=targeted::test_intersection_sw(esttmp,noninf=rep(0,3),weights=c(0.2,0.4,0.4),index=c(1,2))$p.value,
  pvals13.1=targeted::test_intersection_sw(esttmp,noninf=rep(0,3),weights=c(0.2,0.4,0.4),index=c(1,3))$p.value
  ))
}

##########################################################
#Summarizing simulation results for each scenario#########
##########################################################

summary.simres<-function(alpha,simres,sim_args){
  
  type1.w1<-rep(0,4)
  type1.w1[1]<-with(simres,mean(pvals123<alpha))
  type1.w1[2]<-with(simres,mean(pvals12<alpha))
  type1.w1[3]<-with(simres,mean(pvals13<alpha))
  type1.w1[4]<-with(simres,mean(pvals23<alpha))
  type1.w2<-rep(0,4)
  type1.w2[1]<-with(simres,mean(pvals123.1<alpha))
  type1.w2[2]<-with(simres,mean(pvals12.1<alpha))
  type1.w2[3]<-with(simres,mean(pvals13.1<alpha))
  type1.w2[4]<-NA
  Hypothesis<-factor(1:4,labels=c("H123","H12","H13","H23"))
  Sample_size<-rep(sim_args$n,4)
  mu<-rep(sim_args$mu,4)
  sumdata<-data.frame(Hypothesis,Sample_size,mu,type1.w1,type1.w2)
  newnames<-paste(c("type1.w1","type1.w2"),".lambda=",sim_args$lambda)
  names(sumdata)<-c("Hypothesis","Sample_size","mu",newnames)
  return(sumdata)
}

###############################
#Global parameters#############

set.seed(12345) #seed

Nsim <- 10000 #number of simulated datasets per scenario
alpha <-0.025 #Nominal significance level
sigma<-15 #standard deviation of Y among those still alive 
tau <- 2 #landmark time
trteff1 <- 0 #treatment effect terminal event
trteff2 <- 0 #treatment effect in Y among those still alive
gamma <- 15 #composite value


###############################
#simulating scenarios##########
###############################


sim_args<-list(n=200, mu=40, lambda=0.05, tau=tau,
               trteff1=trteff1, trteff2=trteff2, sigma=sigma, gamma=gamma) #Creating list of scenario specific parameters
#Sample size 200

sim_args$n<-200

#mu=40 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)

sim_args$mu<-40
sim_args$lambda<-0.05

res200.mu1.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres200.mu1.lambda1<-summary.simres(alpha,res200.mu1.lambda1,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.05

res200.mu2.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres200.mu2.lambda1<-summary.simres(alpha,res200.mu2.lambda1,sim_args) #summarizing simulations


#mu=40 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-40
sim_args$lambda<-0.08

res200.mu1.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres200.mu1.lambda2<-summary.simres(alpha,res200.mu1.lambda2,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.08

res200.mu2.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres200.mu2.lambda2<-summary.simres(alpha,res200.mu2.lambda2,sim_args) #summarizing simulations

#Sample size 500

sim_args$n<-500

#mu=40 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)

sim_args$mu<-40
sim_args$lambda<-0.05

res500.mu1.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres500.mu1.lambda1<-summary.simres(alpha,res500.mu1.lambda1,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.05

res500.mu2.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres500.mu2.lambda1<-summary.simres(alpha,res500.mu2.lambda1,sim_args) #summarizing simulations


#mu=40 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-40
sim_args$lambda<-0.08

res500.mu1.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres500.mu1.lambda2<-summary.simres(alpha,res500.mu1.lambda2,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.08

res500.mu2.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres500.mu2.lambda2<-summary.simres(alpha,res500.mu2.lambda2,sim_args) #summarizing simulations



#Sample size 1000

sim_args$n<-1000

#mu=40 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)

sim_args$mu<-40
sim_args$lambda<-0.05

res1000.mu1.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres1000.mu1.lambda1<-summary.simres(alpha,res1000.mu1.lambda1,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.05

res1000.mu2.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres1000.mu2.lambda1<-summary.simres(alpha,res1000.mu2.lambda1,sim_args) #summarizing simulations


#mu=40 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-40
sim_args$lambda<-0.08

res1000.mu1.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres1000.mu1.lambda2<-summary.simres(alpha,res1000.mu1.lambda2,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.08

res1000.mu2.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres1000.mu2.lambda2<-summary.simres(alpha,res1000.mu2.lambda2,sim_args) #summarizing simulations


#Sample size 2000

sim_args$n<-2000

#mu=40 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)

sim_args$mu<-40
sim_args$lambda<-0.05

res2000.mu1.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres2000.mu1.lambda1<-summary.simres(alpha,res2000.mu1.lambda1,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.05

res2000.mu2.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres2000.mu2.lambda1<-summary.simres(alpha,res2000.mu2.lambda1,sim_args) #summarizing simulations


#mu=40 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-40
sim_args$lambda<-0.08

res2000.mu1.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres2000.mu1.lambda2<-summary.simres(alpha,res2000.mu1.lambda2,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.08

res2000.mu2.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres2000.mu2.lambda2<-summary.simres(alpha,res2000.mu2.lambda2,sim_args) #summarizing simulations


#Sample size 3500

sim_args$n<-3500

#mu=40 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)

sim_args$mu<-40
sim_args$lambda<-0.05

res3500.mu1.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres3500.mu1.lambda1<-summary.simres(alpha,res3500.mu1.lambda1,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.05 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.05

res3500.mu2.lambda1 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres3500.mu2.lambda1<-summary.simres(alpha,res3500.mu2.lambda1,sim_args) #summarizing simulations


#mu=40 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-40
sim_args$lambda<-0.08

res3500.mu1.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres3500.mu1.lambda2<-summary.simres(alpha,res3500.mu1.lambda2,sim_args) #summarizing simulations

#mu=45 (mean of outcome among those still alive), lambda=0.08 (rate of terminal event)
sim_args$mu<-45
sim_args$lambda<-0.08

res3500.mu2.lambda2 <- lava::sim(onerun, Nsim, args = sim_args) |> as.data.frame() #running simulations
sumres3500.mu2.lambda2<-summary.simres(alpha,res3500.mu2.lambda2,sim_args) #summarizing simulations

################################################
#stitching summary results into one data.frame##
################################################
sumres200<-cbind(rbind(sumres200.mu1.lambda1,sumres200.mu2.lambda1),rbind(sumres200.mu1.lambda2,sumres200.mu2.lambda2)[,4:5])[,c(1:3,4,6,5,7)]
sumres500<-cbind(rbind(sumres500.mu1.lambda1,sumres500.mu2.lambda1),rbind(sumres500.mu1.lambda2,sumres500.mu2.lambda2)[,4:5])[,c(1:3,4,6,5,7)]
sumres1000<-cbind(rbind(sumres1000.mu1.lambda1,sumres1000.mu2.lambda1),rbind(sumres1000.mu1.lambda2,sumres1000.mu2.lambda2)[,4:5])[,c(1:3,4,6,5,7)]
sumres2000<-cbind(rbind(sumres2000.mu1.lambda1,sumres2000.mu2.lambda1),rbind(sumres2000.mu1.lambda2,sumres2000.mu2.lambda2)[,4:5])[,c(1:3,4,6,5,7)]
sumres3500<-cbind(rbind(sumres3500.mu1.lambda1,sumres3500.mu2.lambda1),rbind(sumres3500.mu1.lambda2,sumres3500.mu2.lambda2)[,4:5])[,c(1:3,4,6,5,7)]

sumres<-mets::dsort(rbind(sumres200,sumres500,sumres1000,sumres2000,sumres3500),~Hypothesis+mu+Sample_size)
rownames(sumres)<-NULL
#Creating *tex table

library(knitr)

knitr::kable(sumres, digits=4, format="latex", booktabs=TRUE)

#########################################################
#Resulting table#########################################


\begin{tabular}{lrrrrrr}
\toprule
Hypothesis & Sample\_size & mu & type1.w1 .lambda= 0.05 & type1.w1 .lambda= 0.08 & type1.w2 .lambda= 0.05 & type1.w2 .lambda= 0.08\\
\midrule
H123 & 200 & 40 & 0.0260 & 0.0256 & 0.0242 & 0.0259\\
H123 & 500 & 40 & 0.0220 & 0.0264 & 0.0239 & 0.0267\\
H123 & 1000 & 40 & 0.0238 & 0.0259 & 0.0223 & 0.0259\\
H123 & 2000 & 40 & 0.0265 & 0.0238 & 0.0257 & 0.0254\\
H123 & 3500 & 40 & 0.0266 & 0.0227 & 0.0263 & 0.0242\\
\addlinespace
H123 & 200 & 45 & 0.0255 & 0.0252 & 0.0262 & 0.0240\\
H123 & 500 & 45 & 0.0231 & 0.0246 & 0.0235 & 0.0246\\
H123 & 1000 & 45 & 0.0261 & 0.0261 & 0.0255 & 0.0277\\
H123 & 2000 & 45 & 0.0247 & 0.0240 & 0.0256 & 0.0255\\
H123 & 3500 & 45 & 0.0238 & 0.0269 & 0.0247 & 0.0272\\
\addlinespace
H12 & 200 & 40 & 0.0262 & 0.0261 & 0.0249 & 0.0268\\
H12 & 500 & 40 & 0.0218 & 0.0261 & 0.0239 & 0.0260\\
H12 & 1000 & 40 & 0.0238 & 0.0258 & 0.0218 & 0.0250\\
H12 & 2000 & 40 & 0.0260 & 0.0242 & 0.0261 & 0.0267\\
H12 & 3500 & 40 & 0.0267 & 0.0225 & 0.0265 & 0.0238\\
\addlinespace
H12 & 200 & 45 & 0.0249 & 0.0251 & 0.0274 & 0.0220\\
H12 & 500 & 45 & 0.0233 & 0.0259 & 0.0242 & 0.0264\\
H12 & 1000 & 45 & 0.0261 & 0.0262 & 0.0250 & 0.0280\\
H12 & 2000 & 45 & 0.0251 & 0.0247 & 0.0268 & 0.0264\\
H12 & 3500 & 45 & 0.0236 & 0.0276 & 0.0236 & 0.0256\\
\addlinespace
H13 & 200 & 40 & 0.0238 & 0.0253 & 0.0238 & 0.0254\\
H13 & 500 & 40 & 0.0225 & 0.0274 & 0.0241 & 0.0283\\
H13 & 1000 & 40 & 0.0249 & 0.0270 & 0.0229 & 0.0269\\
H13 & 2000 & 40 & 0.0256 & 0.0234 & 0.0246 & 0.0234\\
H13 & 3500 & 40 & 0.0260 & 0.0237 & 0.0260 & 0.0254\\
\addlinespace
H13 & 200 & 45 & 0.0238 & 0.0276 & 0.0256 & 0.0239\\
H13 & 500 & 45 & 0.0246 & 0.0254 & 0.0244 & 0.0255\\
H13 & 1000 & 45 & 0.0259 & 0.0266 & 0.0240 & 0.0286\\
H13 & 2000 & 45 & 0.0246 & 0.0256 & 0.0249 & 0.0261\\
H13 & 3500 & 45 & 0.0238 & 0.0273 & 0.0251 & 0.0269\\
\addlinespace
H23 & 200 & 40 & 0.0240 & 0.0257 & NA & NA\\
H23 & 500 & 40 & 0.0239 & 0.0271 & NA & NA\\
H23 & 1000 & 40 & 0.0219 & 0.0260 & NA & NA\\
H23 & 2000 & 40 & 0.0266 & 0.0255 & NA & NA\\
H23 & 3500 & 40 & 0.0264 & 0.0239 & NA & NA\\
\addlinespace
H23 & 200 & 45 & 0.0263 & 0.0238 & NA & NA\\
H23 & 500 & 45 & 0.0237 & 0.0254 & NA & NA\\
H23 & 1000 & 45 & 0.0254 & 0.0272 & NA & NA\\
H23 & 2000 & 45 & 0.0255 & 0.0261 & NA & NA\\
H23 & 3500 & 45 & 0.0243 & 0.0278 & NA & NA\\
\bottomrule
\end{tabular}



