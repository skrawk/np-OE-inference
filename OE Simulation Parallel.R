########################################################################
### This code executes a simulation comparing a number of methods for
### producing propensities of selection for a non-probability dataset.
### The simulation setup is similar to what was used in Chen et al (2020),
### and parts of the code is based on what was used for their simulation.
### My sincere thanks to the authors, who provided me a copy of their code.
###
### Author: Lyndon Ang
### Year: 2024
###
######################################################################

### Note: In this code, the dataset A is the non-probability sample, while
###       the dataset B represents the probability sample.

rm(list=ls())

library("xtable")
library("VIM")
library(sampling)

### Libraries to facilitate parallel computing
library(foreach)
library(doParallel)
library(doRNG)

starttime <- Sys.time()

### Set the number of core processors to use
registerDoParallel(cores=4)


#----------------------------------------------------------
#------------------self-defined functions------------------
#----------------------------------------------------------

#randomized systematic sampling

syspps=function(N_f,x){ 
  ##Population is first randomized
  U=sample(N_f,N_f)
  xx=x[U]
  cx=cumsum(xx)
  r=runif(1)
  s=numeric()
  for(i in 1:N_f){
    if(cx[i]>=r){
      s=c(s,U[i])
      r=r+1
    }
  }
  return(s)
}


#Estimate theta with Chen et al method

hat_theta_f=function(xsa,xsb,wsb){
  #xsa <- X_SA_res
  #xsb <- X_SB_res
  #wsb <- w_SB
  
  col_xsa=colSums(xsa)
  theta_0=rep(0,length(col_xsa))
  theta_1=solve(t(xsb)%*%(wsb*xsb),col_xsa-t(xsb)%*%wsb)
  while (abs(max(theta_1-theta_0))>10^(-8)){
    theta_0=theta_1
    ps=1/(1+exp(-xsb%*%theta_0))
    theta_1=theta_0+solve(t(xsb)%*%(c(wsb*ps*(1-ps))*xsb),col_xsa-t(xsb)%*%c(wsb*ps))
  }
  
  return(theta_1)
}


### Estimate theta using subset of B excluding A

hat_theta2_f <- function(xsa, xsb, wsb) {
  #xsa <- X_SA_res
  #xsb <- X_SBnA_res
  #wsb <- w_SBnA
  
  col_xsa <- colSums(xsa) ### Calculate column sums of x variables in A
  theta_0 <- rep(0.1, length(col_xsa)) ### Initial values for theta
  theta_1 <- theta_0 - 0.1
  
  while(abs(max(theta_1 - theta_0)) > 10^(-8)) {
    theta_0 <- theta_1
    ps_A <- 1/(1 + exp(-xsa %*% theta_0)) ### Calculate propensity for units in A
    ps_BnA <- 1/(1 + exp(-xsb %*% theta_0)) ### Calculate propensity for units in B not A
    
    ### Now find the next iteration of theta
    theta_1 <- theta_0 + solve(t(xsa) %*% (c(ps_A*(1 - ps_A)) * xsa) + t(xsb)%*%(c(wsb*ps_BnA*(1-ps_BnA))*xsb),
                               col_xsa - t(xsa) %*% ps_A - t(xsb)%*%c(wsb*ps_BnA))
    
  }
  
  return(theta_1)
  
}



### Estimate theta using a 'predicted' delta_i value

hat_theta_EM_f <- function(xsa, xsb, wsb, wsbna) {
  #xsa <- X_SA_res
  #xsb <- X_SBnA_res
  #wsb <- w_SBnA
  
  col_xsa <- colSums(xsa) ### Calculate column sums of x variables in A
  theta_0 <- rep(0.1, length(col_xsa)) ### Initial values for theta
  theta_1 <- theta_0 - 0.1
  
  while(abs(max(theta_1 - theta_0)) > 10^(-8)) {
    theta_0 <- theta_1
    ps_A <- 1/(1 + exp(-xsa %*% theta_0)) ### Calculate propensity for units in A
    ps_B <- 1/(1 + exp(-xsb %*% theta_0)) ### Calculate propensity for units in B
    
    ### Now find the next iteration of theta
    theta_1 <- theta_0 + solve(t(xsa) %*% (c(ps_A*(1 - ps_A)) * xsa) + t(xsb)%*%(c(wsb*ps_B*(1-ps_B))*xsb) - t(xsb)%*%(c(wsbna*ps_B*(1-ps_B))*xsb),
                               col_xsa - t(xsa) %*% ps_A - t(xsb)%*%c(wsb*ps_B) + t(xsb)%*%c(wsbna*ps_B))
    
  }
  
  return(theta_1)
  
}


### Implement EM algorithm to produce final theta values based on predicted delta_i values

EM_loop_f <- function(theta_em_init,xsa_em,xsb_em,wsb) {
  
  mod.coef.old <- theta_em_init
  mod.coef.new <- theta_em_init*2
  theta_em <- theta_em_init
  
  i <- 0
  
  ### EM iterations
  while (max(abs(mod.coef.new-mod.coef.old)) > 10^(-8)) {
    mod.coef.old <- theta_em
    i<-i+1
    
    ### E-step - estimate pi.B based on current set of theta parameters
    hat_score_SB.em <- 1/(1+exp(-xsb_em%*%theta_em))
    
    ### Derive a new "weight" value which is the design weight multiplied by hat_score_SB.em
    w_SBna_pred <- wsb * hat_score_SB.em
    
    ### M-step
    ### Calculate new hat_theta_em values using the new weights
    theta_em <- hat_theta_EM_f(xsa_em, xsb_em, wsb, w_SBna_pred)
    
    mod.coef.new <- theta_em
    
    ### To stop infinite loops
    if(i>100) break
    
  }
  
  return(theta_em)
  
}


### Estimate theta using A and B linked together (Kim and Wang method)

hat_theta_KW_f <- function(xsbpa, xsbna, wsbpa, wsbna) {
  
  #xsbpa <- X_SBpA_res
  #xsbna <- X_SBnA_res
  #wsbpa <- w_SBpA
  #wsbna <- w_SBnA
  
  col_xsbpa <- t(xsbpa) %*% wsbpa  ### Calculate weighted column sums of x variables in B
  theta_0 <- rep(0.1, length(col_xsbpa)) ### Initial values for theta
  theta_1 <- theta_0 - 0.1
  
  while(abs(max(theta_1 - theta_0)) > 10^(-8)) {
    theta_0 <- theta_1
    ps_BpA <- 1/(1 + exp(-xsbpa %*% theta_0)) ### Calculate propensity for units in B linked to A
    ps_BnA <- 1/(1 + exp(-xsbna %*% theta_0)) ### Calculate propensity for units in B not A
    
    ### Now find the next iteration of theta
    theta_1 <- theta_0 + solve(t(xsbpa) %*% (c(wsbpa*ps_BpA*(1 - ps_BpA)) * xsbpa) + t(xsbna)%*%(c(wsbna*ps_BnA*(1-ps_BnA))*xsbna),
                               col_xsbpa - t(xsbpa) %*% c(wsbpa*ps_BpA) - t(xsbna)%*%c(wsbna*ps_BnA))
    
  }
  
  return(theta_1)
  
}


### Estimate theta with Valliant and Dever (VD) method

hat_theta_VD_f <- function(xsa, xsb, wsb) {
  #xsa <- X_SA_res
  #xsb <- X_SB_res
  #wsb <- w_VD
  
  col_xsa <- colSums(xsa) ### Calculate column sums of x variables in A
  theta_0 <- rep(0.1, length(col_xsa)) ### Initial values for theta
  theta_1 <- theta_0 - 0.1
  
  while(abs(max(theta_1 - theta_0)) > 10^(-8)) {
    theta_0 <- theta_1
    ps_A <- 1/(1 + exp(-xsa %*% theta_0)) ### Calculate propensity for units in A
    ps_B <- 1/(1 + exp(-xsb %*% theta_0)) ### Calculate propensity for units in B
    
    ### Now find the next iteration of theta
    theta_1 <- theta_0 + solve(t(xsa) %*% (c(ps_A*(1 - ps_A)) * xsa) + t(xsb)%*%(c(wsb*ps_B*(1-ps_B))*xsb),
                               col_xsa - t(xsa) %*% ps_A - t(xsb)%*%c(wsb*ps_B))
    
  }
  
  return(theta_1)
  
}



### Bootstrap function for the EM algorithm approach using predicted delta_i values

boot_var_em_f=function(B_f,ysa,xsa_res,xsb_res,xsa_reg,xsb_reg,wsb,nsa,nsb){
  
  byest=numeric(0)
  
  for (bi in 1:B_f){
    ### 1 - Select the NP sample A
    
    bsa_index=sample(nsa,nsa,replace=TRUE)
    bxsa_reg=xsa_reg[bsa_index,]
    bxsa_res=xsa_res[bsa_index,]
    bysa=ysa[bsa_index]
    
    ### 2 - Select the P sample B
    bsb_index=sample(nsb,nsb,replace=TRUE)
    bxsb_reg=xsb_reg[bsb_index,]
    bxsb_res=xsb_res[bsb_index,]
    bwsb=wsb[bsb_index]
    
    ### 3a - Estimate initial theta for the bootstrap sample using CLW
    btheta_init=hat_theta_f(bxsa_res,bxsb_res,bwsb)
    
    ### 3b - Adjust initial hat_theta using EM algorithm
    btheta_em <- EM_loop_f(btheta_init,bxsa_res,bxsb_res,bwsb)
    
    ### 4 - Calculate the pi.A values for the NP sample
    bscore=1/(1+exp(-bxsa_res%*%btheta_em))
    
    ### 5 - Calculate estimate Y-hat for the bootstrap sample
    byest[bi] <- sum(bysa/bscore)/sum(1/bscore)
  }
  
  return(var(byest))
}


##################parameter settings##############



set.seed(20231)
#set.seed(2024)
B=1                            #bootstrap sample size
N=200000                             #population total
iteration=10                     #simulation iteration                        
beta=c(2,1,1,1,1)                 #True beta for "y=x beta + e"
theta=c(0.1,0.2,0.3,0.4)          #True theta for PS, no intercept
#theta=c(0.2,0.3,0.5,0.8)
#Covariate
x1=rbinom(N,1,0.5) 
x2=runif(N,min=0,max=2)+0.3*x1                                                   
x3=rexp(N,1)+0.2*(x1+x2)                                                    
x4=rchisq(N,4)+0.1*(x1+x2+x3)



#Design matrix
X=cbind(rep(1,N),x1,x2,x3,x4)
#working matrices
X_reg=cbind(rep(1,N),x1,x2,x3,x4)
X_res=cbind(rep(1,N),x1,x2,x3,x4)
e=rnorm(N,mean=0,sd=1)

rho=0.3            #cor(Y,XB)
mean_n_SA=750   #sample size of S_A 
mean_n_SB=5000  #sample size of S_B 

### Number of multiple imputes
mi.num <- 1

#find sigma
sigma=sqrt((1/rho^2-1)*var(c(X%*%beta)))

#Generate y
Y=X%*%beta+e*sigma
#Check cor(Y,XB)
cor(Y,X%*%beta)
mu=mean(Y)

### We will now put Y and X into a data frame
popn <- data.frame(Y,X)
### Add a column with a unit identifier
popn$id <- 1:nrow(popn)

###############################################

#find intercept
#mean_n_SA <- 180000
#theta=c(0.2,0.3,0.5,0.8)
etap=as.vector(X[,-1]%*%theta)
dif=1
L=-10
R=-1
while(dif>0){
  M=(L+R)/2
  ps=exp(M+etap)/(1+exp(M+etap))
  if(sum(ps)<=mean_n_SA) L=M
  if(sum(ps)>=mean_n_SA) R=M
  if(abs(sum(ps)-mean_n_SA)<=0.5) dif=0
}


#true score
score=exp(M+etap)/(1+exp(M+etap))   
#score <-0.5
(min(score))
(max(score))
(sum(score))

#set max(z1)/min(z1)=50
#find c

z=x3
c=(max(z)-10*min(z))/9
s=sum(z+c)
w=s/(z+c)/mean_n_SB
z1=(z+c)/s


### Calculate population value of b
a.pop <- t(X_res)%*%(c(score*(1-score))*X_res)
b.pop <- solve(a.pop,t(X_res)%*%c((1 - score)*(Y-mean(Y))))

### Now calculate probability of selection into S_B (probability sample)
cost.i <- 1

### Limit the spread of weights for the optimal B sample
z.opt <- (score * (sqrt(1 - score)) * abs(X_res %*% b.pop)/sqrt(cost.i))
c.opt <- (max(z.opt) - 10*min(z.opt))/9
s.opt <- sum(z.opt + c.opt)
w.opt <- s.opt/(z.opt + c.opt)/mean_n_SB
z1.opt <- (z.opt + c.opt)/s.opt

pi.B <- z1.opt * mean_n_SB
w_opt <- 1/pi.B

### Form final population level dataset including the design weights for probability sample B
popn <- cbind(popn,w,w_opt)


#----------------------------------------------------------
#------------------iteration starts here-------------------
#----------------------------------------------------------

Bias=0
MSE=0
est_var=0
c_rate=0

simresults <- foreach(g = 1:iteration, .combine=rbind) %dorng% {
  
  #draw the non-probability sample 
  sam=rbinom(N,1,score)
  n_SA=sum(sam)
  index_SA=order(-sam)[1:n_SA] 
  Y_SA=Y[index_SA]
  X_SA_reg=X_reg[index_SA,]
  X_SA_res=X_res[index_SA,]
  w_SA=w[index_SA]
  
  samp_A <- popn[index_SA,]
  
  
  ### Approach 1 - Draw probability sample with sys method (not optimal, based on x3)
  #pik <- sampling::inclusionprobabilities(z1, mean_n_SB)
  #index_SB=sampling::UPopips(pik, type="pareto")
  #s<-sampling::UPpivotal(pik)
  #index_SB <- (1:length(pik))[s==1]
  index_SB=syspps(N,z1*mean_n_SB)
  
  n_SB=mean_n_SB
  X_SB_res=X_res[index_SB,]
  X_SB_reg=X_reg[index_SB,]
  w_SB=w[index_SB]
  Y_SB <- Y[index_SB]
  
  samp_B <- popn[index_SB,]
  
  
  ### Approach 2 - Draw a probability sample with sys method split-population approach
  
  ### Find the number of observations in sample B \ A above - this will be our sample size for 
  ### the split-population B sample
  #mean_n_SC <- nrow(samp_B[!(samp_B$id %in% samp_A$id),])
  ### Or... set sample to B, so can directly compare gains
  mean_n_SC <- mean_n_SB
  
  ### Obtain the subset C = U \ A which we will base another probability sample on
  index_SC <- popn$id[!(popn$id %in% samp_A$id)]
  N_SC <- N - n_SA
  X_SC_res <- X_res[index_SC,]
  X_SC_reg <- X_reg[index_SC,]
  x3_SC <- x3[index_SC]
  
  popn.C <- popn[index_SC,]
  
  
  ### Now obtain our size measure for subpopulation C
  z.C=x3_SC
  c=(max(z.C)-10*min(z.C))/9
  s.C=sum(z.C+c)
  w.C=s.C/(z.C+c)/mean_n_SC
  z1.C=(z.C+c)/s.C
  
  pi.C <- z1.C * mean_n_SC
  w_C <- 1/pi.C
  
  ### Form final population level dataset including the design weights for probability sample C
  popn.C <- cbind(popn.C,w_C)
  
  
  ### Draw a probability sample which is based on the (N - n_A) units as the population
  
  index_SB.C=syspps(N_SC,pi.C)
  n_SB.C=mean_n_SC
  X_SB.C_res=X_SC_res[index_SB.C,]
  X_SB.C_reg=X_SC_res[index_SB.C,]
  w_SB.C=w_C[index_SB.C]
  
  samp_B.C <- popn.C[index_SB.C,]
  
  
  ### Approach 3 - Draw a probability sample with sys method - using optimal probabilities
  index_SB_opt <- syspps(N,pi.B)
  n_SB_opt <- mean_n_SB
  X_SB_opt_res <- X_res[index_SB_opt,]
  X_SB_opt_reg <- X_reg[index_SB_opt,]
  w_SB_opt <- w_opt[index_SB_opt]
  
  samp_B_opt <- popn[index_SB_opt,]
  
  ############################################################
  
  
  ### Obtain the subset B \ A which we will use for the alternative method
  
  ### Approach 1 - Assume we can directly obtain B \ A
  
  samp_B_not_A <- samp_B[!(samp_B$id %in% samp_A$id),]
  
  ### Now put together the various datasets we need for B \ A
  X_SBnA_res <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  X_SBnA_reg <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  w_SBnA <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id),c("w")])
  
  
  ### Approach 2 - Assume that we must estimate the links between A and B
  
  ### Our initial hat_theta is from CLW
  hat_theta_em_init=hat_theta_f(X_SA_res,X_SB_res,w_SB)
  
  ### Adjust initial hat_theta using EM algorithm
  hat_theta_em <- EM_loop_f(hat_theta_em_init,X_SA_res,X_SB_res,w_SB)
  
  
  ### Split population approach ###
  ### Obtain B \ A which will be used for the alternative method - for the split-population approach
  ### this will be the full sample B since there is no overlap between A and B
  
  samp_B.C_not_A <- samp_B.C[!(samp_B.C$id %in% samp_A$id),]
  X_SBCnA_res <- as.matrix(samp_B.C[!(samp_B.C$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  X_SBCnA_reg <- as.matrix(samp_B.C[!(samp_B.C$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  w_SBCnA <- as.matrix(samp_B.C[!(samp_B.C$id %in% samp_A$id),c("w_C")])
  
  ### Optimal allocation
  ### For our alternative method obtain datasets using optimally sampled B
  
  samp_B_opt_not_A <- samp_B_opt[!(samp_B_opt$id %in% samp_A$id),]
  X_SBnA_opt_res <- as.matrix(samp_B_opt[!(samp_B_opt$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  X_SBnA_opt_reg <- as.matrix(samp_B_opt[!(samp_B_opt$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  w_SBnA_opt <- as.matrix(samp_B_opt[!(samp_B_opt$id %in% samp_A$id),c("w_opt")])
  
  
  ### Also find the subset B and A which will be used for the Kim and Wang method
  samp_B_and_A <- samp_B[samp_B$id %in% samp_A$id,]
  
  ### Now put together the various datasets we need for B and A
  X_SBpA_res <- as.matrix(samp_B[samp_B$id %in% samp_A$id, c("V1","x1","x2","x3","x4")])
  X_SBpA_reg <- as.matrix(samp_B[samp_B$id %in% samp_A$id, c("V1","x1","x2","x3","x4")])
  w_SBpA <- as.matrix(samp_B[samp_B$id %in% samp_A$id,c("w")])
  
  #check sum
  N_B=sum(w_SB)
  
  ### Now calculate the weight for the VD method
  w_VD <- w_SB*(N_B - n_SA)/N_B
  
  
  #estimated theta and scores
  hat_theta=hat_theta_f(X_SA_res,X_SB_res,w_SB)
  hat_score_SA=1/(1+exp(-X_SA_res%*%hat_theta))
  hat_score_SB=1/(1+exp(-X_SB_res%*%hat_theta))
  
  ### Estimated theta for the alternative method
  hat_theta2 <- hat_theta2_f(X_SA_res, X_SBnA_res, w_SBnA)
  hat_score2_SA=1/(1+exp(-X_SA_res%*%hat_theta2))
  hat_score2_SB=1/(1+exp(-X_SB_res%*%hat_theta2))
  
  ### Estimated theta for alternative method - predicted deltas
  hat_score2_SAp=1/(1+exp(-X_SA_res%*%hat_theta_em))
  hat_score2_SBp=1/(1+exp(-X_SB_res%*%hat_theta_em))
  
  ### Estimated theta for the alternative method using split population
  hat_theta2_C <- hat_theta2_f(X_SA_res, X_SBCnA_res, w_SBCnA)
  hat_score2_SA_C <- 1/(1+exp(-X_SA_res%*%hat_theta2_C))
  hat_score2_SB_C <- 1/(1+exp(-X_SB.C_res%*%hat_theta2_C))
  
  ### Estimated theta for alternative method using optimal B sample
  hat_theta2_opt <- hat_theta2_f(X_SA_res, X_SBnA_opt_res, w_SBnA_opt)
  hat_score2_opt_SA <- 1/(1+exp(-X_SA_res%*%hat_theta2_opt))
  hat_score2_opt_SB <- 1/(1+exp(-X_SB_opt_res%*%hat_theta2_opt))
  
  ### Estimated p for WVL method - use the algorithm from our method as we have the same expression
  # Scale weights for B
  lw_SB <- w_SB/mean(w_SB)
  
  hat_theta_WVL <- hat_theta2_f(X_SA_res,X_SB_res,lw_SB)
  hat_p_SA <- 1/(1+exp(-X_SA_res%*%hat_theta_WVL))
  hat_p_SB <- 1/(1+exp(-X_SB_res%*%hat_theta_WVL))
  # Now transform to theta-hats
  hat_score_WVL_SA <- hat_p_SA/(1-hat_p_SA)
  hat_score_WVL_SB <- hat_p_SB/(1-hat_p_SB)
  
  
  ### Estimated theta for KW method
  hat_theta3 <- hat_theta_KW_f(X_SBpA_res, X_SBnA_res, w_SBpA, w_SBnA)
  hat_score3_SA <- 1/(1+exp(-X_SA_res%*%hat_theta3))
  hat_score3_SB <- 1/(1+exp(-X_SB_res%*%hat_theta3))
  
  
  ### Estimated theta for VD method
  hat_theta_VD <- hat_theta_VD_f(X_SA_res, X_SB_res, w_VD)
  hat_score_VD_SA <- 1/(1+exp(-X_SA_res%*%hat_theta_VD))
  hat_score_VD_SB <- 1/(1+exp(-X_SB_res%*%hat_theta_VD))
  
  
  ### Estimate the theta for conditional likelihood
  X_AuB <- rbind(X_SA_res, X_SBnA_res)
  X_AuB <- X_AuB[,-1] # Remove the intercept term as it will be fit in the model later
  
  delta_A <- rep(1, nrow(X_SA_res))
  delta_B <- rep(0, nrow(X_SBnA_res))
  delta_AuB <- c(delta_A, delta_B)
  os_A <- -log(1/w_SA)
  os_B <- -log(1/w_SBnA)
  os_AuB <- c(os_A, os_B)
  ### Fit a glm using X_AuB as covariates and os_AuB as offset
  hat_theta_cond <- glm(delta_AuB ~ X_AuB, offset=os_AuB, family=binomial(link="logit"))
  hat_score_cond_SA <- 1/(1+exp(-X_SA_res%*%hat_theta_cond$coefficients))
  hat_score_cond_SB <- 1/(1+exp(-X_SB_res%*%hat_theta_cond$coefficients))
  
  ### Produce estimates of population size
  
  N_A=sum(1/hat_score_SA)
  N_A2=sum(1/hat_score2_SA)
  N_A2p <- sum(1/hat_score2_SAp)
  N_A2_C <- sum(1/hat_score2_SA_C)
  N_A2_opt <- sum(1/hat_score2_opt_SA)
  N_A3 <- sum(1/hat_score3_SA)
  N_A_WVL <- sum(1/hat_score_WVL_SA)
  N_A_VD <- sum(1/hat_score_VD_SA)
  N_A_cond <- sum(1/hat_score_cond_SA)
  
  
  #-----------------------------------------------------
  #--------------------point estimates------------------
  #-----------------------------------------------------
  
  mu_naive=mean(Y_SA)
  mu_ht <- sum(Y_SB*w_SB)/N_B
  
  mu_ipw1=sum(Y_SA/hat_score_SA)/N
  mu_ipw2=sum(Y_SA/hat_score_SA)/N_A
  
  mu_new <- sum(Y_SA/hat_score2_SA)/N_A2
  mu_newp <- sum(Y_SA/hat_score2_SAp)/N_A2p
  mu_new_sp <- sum(Y_SA/hat_score2_SA_C)/N_A2_C
  mu_new_opt <- sum(Y_SA/hat_score2_opt_SA)/N_A2_opt
  
  mu_wvl <- sum(Y_SA/hat_score_WVL_SA)/N_A_WVL
  
  mu_kw <- sum(Y_SA/hat_score3_SA)/N_A3
  
  mu_vd <- sum(Y_SA/hat_score_VD_SA)/N_A_VD
  
  mu_cond <- sum(Y_SA/hat_score_cond_SA)/N_A_cond
  
  #-----------------------------------------------------
  #---------Calculate estimated optimal probs-----------
  #-----------------------------------------------------
  
  ### The estimated mean calculated above is the last term we need for
  ### the probabilities.  Now we can calculate the estimated optimal probabilities
  
  ### Use the non-probability sample to impute x onto popn
  ### We assume x2 is available on popn and A, so impute the other x variables
  ### the hot-deck method is used for imputation
  
  ### Create a copy of the popn dataset for the imputed values
  popn2 <- popn[,c("x1","x2","x3","x4","id")]
  
  ### Set x1,x2,x4 to missing if the row is not in samp_B
  popn2$x1[!(popn$id %in% samp_A$id)] <- NA
  popn2$x2[!(popn$id %in% samp_A$id)] <- NA
  popn2$x4[!(popn$id %in% samp_A$id)] <- NA
  
  ### Produce decile categories for x3
  popn2$x3.quantile <- with(popn2, factor(
    findInterval( x3, c(-Inf,
                        quantile(x3, probs=c(0.1, .2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), Inf)), 
    labels=c("Q1","Q2","Q3","Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")
  ))
  
  ### Set up an empty matrix to contain the optimal pi values
  pi.B.estopt.matrix <- matrix(NA, nrow=nrow(popn2), ncol=mi.num)
  
  ### Perform multiple hot deck imputation to estimate optimal pi.B values
  
  for (mi.loop in 1:mi.num) {
    
    ### Create a random order variable
    popn2$rand.ord <- runif(nrow(popn2))
    
    ### Now impute the variables x1, x2, x4 using the hot-deck method
    popn.imp <- VIM::hotdeck(data=popn2, 
                              ord_var="rand.ord",
                              variable=c("x1","x2","x4"),
                              domain_var="x3.quantile"
    )
    
    ### Rename x1, x2 and x4 to x1.est, x2.est, x4.est
    names(popn.imp)[names(popn.imp)=="x1"] <- "x1.est"
    names(popn.imp)[names(popn.imp)=="x2"] <- "x2.est"
    names(popn.imp)[names(popn.imp)=="x4"] <- "x4.est"
    
    ### Merge the x1, x2 and x4 variables back on to popn
    popn_sim <- merge(popn,popn.imp[,c("x1.est","x2.est","x4.est","id")],by="id")
    
    X_P_res <- as.matrix(popn_sim[,c("V1","x1.est","x2.est","x3","x4.est")])
    X_P_reg <- as.matrix(popn_sim[,c("V1","x1.est","x2.est","x3","x4.est")])
    
    ### Use the estimated theta hat_theta2 on the estimated x values in popn
    ### Assumes that we have produced hat_theta2 from a previous sample B
    hat_score2_estopt_P <- 1/(1+exp(-X_P_res%*%hat_theta2))
    
    
    ### Calculate estimated population value of b using previous sample B
    a.estpop <- t(X_SA_res)%*%(c(hat_score2_SA*(1 - hat_score2_SA))*X_SA_res)
    b.estpop <- solve(a.estpop,t(X_SA_res)%*%c((1 - hat_score2_SA)*(Y_SA - mu_new)))
    
    ### Now calculate probability of selection into S_B (probability sample)
    cost.i <- 1
    
    
    ### Limit the spread of weights for the optimal B sample
    z.estopt <- (hat_score2_estopt_P * (sqrt(1 - hat_score2_estopt_P)) * abs(X_P_res %*% b.estpop)/sqrt(cost.i))
    c.estopt <- (max(z.estopt) - 10*min(z.estopt))/9
    s.estopt <- sum(z.estopt + c.estopt)
    w.estopt <- s.estopt/(z.estopt + c.estopt)/mean_n_SB
    z1.estopt <- (z.estopt + c.estopt)/s.estopt
    
    pi.B.imp <- z1.estopt * mean_n_SB
    
    pi.B.estopt.matrix[,mi.loop] <- pi.B.imp
  }
  
  ### Calculate the mean of the estimated optimal probabilities
  pi.B.estopt <- rowMeans(pi.B.estopt.matrix)
  
  w_estopt <- 1/pi.B.estopt
  
  #-----------------------------------------------------
  #-------------------Draw new sample-------------------
  #-----------------------------------------------------
  
  ### Draw a probability sample with sys method - using estimated 
  ### popn optimal probabilities
  index_SB_estopt <- syspps(N,pi.B.estopt)
  n_SB_estopt <- mean_n_SB
  X_SB_estopt_res <- X_res[index_SB_estopt,]
  X_SB_estopt_reg <- X_reg[index_SB_estopt,]
  w_SB_estopt <- w_estopt[index_SB_estopt]
  
  popn3 <- cbind(popn,w_estopt)
  
  samp_B_estopt <- popn3[index_SB_estopt,]
  
  ### For our alternative method obtain datasets using approximately optimally sampled B
  samp_B_estopt_not_A <- samp_B_estopt[!(samp_B_estopt$id %in% samp_A$id),]
  X_SBnA_estopt_res <- as.matrix(samp_B_estopt[!(samp_B_estopt$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  X_SBnA_estopt_reg <- as.matrix(samp_B_estopt[!(samp_B_estopt$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  w_SBnA_estopt <- as.matrix(samp_B_estopt[!(samp_B_estopt$id %in% samp_A$id),c("w_estopt")])
  
  #-----------------------------------------------------
  #---Calculate estimates for approx optimal sample-----
  #-----------------------------------------------------
  
  ### Estimated theta for alternative method using approximately optimal B sample
  hat_theta2_estopt <- hat_theta2_f(X_SA_res, X_SBnA_estopt_res, w_SBnA_estopt)
  hat_score2_estopt_SA <- 1/(1+exp(-X_SA_res%*%hat_theta2_estopt))
  hat_score2_estopt_SB <- 1/(1+exp(-X_SB_estopt_res%*%hat_theta2_estopt))
  
  N_A2_estopt <- sum(1/hat_score2_estopt_SA)
  
  mu_new_estopt <- sum(Y_SA/hat_score2_estopt_SA)/N_A2_estopt
  
  
  ### Create estimates of mean
  muhat=c(mu_naive,mu_ipw1,mu_ipw2,mu_new,mu_newp,mu_new_sp,mu_new_opt,mu_new_estopt,mu_wvl,mu_kw,mu_vd,mu_cond,mu_ht)
  
  ### Calculate the overlap between the samples A and B
  overlap_AB <- sum(samp_A$id %in% samp_B$id)
  
  #-----------------------------------------------------
  #-----------------------Variance-----------------------
  #-----------------------------------------------------
  a=t(X_SB_res)%*%(c(w_SB*hat_score_SB*(1-hat_score_SB))*X_SB_res)
  anew <- t(X_SB_res)%*%(c(w_SB*hat_score2_SB*(1-hat_score2_SB))*X_SB_res)
  anew2 <- t(X_SA_res)%*%(c(1-hat_score2_SA)*X_SA_res)
  
  p_i <- hat_score_cond_SA/(hat_score_cond_SA + 1/w_SA - hat_score_cond_SA/w_SA)
  p_iA <- hat_score_cond_SB/(hat_score_cond_SB + 1/w_SB - hat_score_cond_SB/w_SB)
  a_cond <- t(X_SA_res) %*% (c(1 - p_i)*X_SA_res)
  
  hat_b1=solve(a,t(X_SA_res)%*%c((1/hat_score_SA-1)*Y_SA))
  hat_b2=solve(a,t(X_SA_res)%*%c((1/hat_score_SA-1)*(Y_SA-mu_ipw2)))

  hat_bnew <- solve(anew2,t(X_SA_res)%*%c((1/hat_score2_SA-1)*(Y_SA-mu_new)))
  
  hat_bcond <- solve(a_cond, t(X_SA_res)%*%c((1/hat_score_cond_SA-1)*(Y_SA-mu_cond)))
  
  ### Find b for the split population approach
  #Original denominator
  #a.sp <- t(X_SA_res) %*% (c(1 + hat_score2_SA_C)*c(1 - hat_score2_SA_C)*X_SA_res)
  #Revised when likelihood over U_C is used
  a.sp <- t(X_SA_res) %*% (c(hat_score2_SA_C)*c(1 - hat_score2_SA_C)*X_SA_res) + t(X_SB.C_res) %*% (c(w_SB.C*hat_score2_SB_C)*c(1 - hat_score2_SB_C)*X_SB.C_res)
  hat_bnew_sp <- solve(a.sp, t(X_SA_res)%*%c((1/hat_score2_SA_C - 1)*(Y_SA-mu_new_sp)))
  
  
  ### Estimating b for optimal probs
  anew2_opt <- t(X_SA_res)%*%(c(1-hat_score2_opt_SA)*X_SA_res)
  hat_bnew_opt <- solve(anew2_opt,t(X_SA_res)%*%c((1/hat_score2_opt_SA-1)*(Y_SA-mu_new_opt)))
  
  ### Estimating b for anticipated optimal probabilities
  anew2_estopt <- t(X_SA_res)%*%(c(1-hat_score2_estopt_SA)*X_SA_res)
  hat_bnew_estopt <- solve(anew2_estopt,t(X_SA_res)%*%c((1/hat_score2_estopt_SA-1)*(Y_SA-mu_new_estopt)))
  
  ### Variance component for IPW1 in CLW
  tt=c(X_SB_res%*%hat_b1)*c(hat_score_SB)
  st=sum(w_SB*tt)
  est_vd1=1/N^2*sum((w_SB*tt-st/n_SB)^2)
  
  ### Variance component for IPW2 in CLW
  tt=c(X_SB_res%*%hat_b2)*c(hat_score_SB)
  st=sum(w_SB*tt)
  est_vd2=1/N^2*sum((w_SB*tt-st/n_SB)^2)
  
  ### Variance component for conditional likelihood
  tt=c(X_SB_res%*%hat_bcond)*(c(p_iA)*c(1/w_SB)*(1-c(hat_score_cond_SB)))
  st=sum(w_SB*tt)
  est_vdcond=1/N^2*sum((w_SB*tt-st/n_SB)^2)
  
  ### variance component for new method
  tt=c(X_SB_res%*%hat_bnew)*(c(hat_score2_SB)*(1-c(hat_score2_SB)))
  st=sum(w_SB*tt)
  est_vdnew=1/N^2*sum((w_SB*tt-st/n_SB)^2)
  
  
  ### variance component for new method - split population
  tt <- c(X_SB.C_res %*% hat_bnew_sp) * c(hat_score2_SB_C)
  st <- sum(w_SB.C * tt)
  est_vdnew_sp <- (1/N^2) * sum((w_SB.C*tt - st/n_SB.C)^2)
  
  
  ### variance component for optimal probs
  tt <- c(X_SB_opt_res%*%hat_bnew_opt)*(c(hat_score2_opt_SB)*(1-c(hat_score2_opt_SB)))
  st <- sum(w_SB_opt*tt)
  est_vdnew_opt <- 1/N^2*sum((w_SB_opt*tt-st/n_SB_opt)^2)
  
  ### variance component for estopt probs
  tt <- c(X_SB_estopt_res%*%hat_bnew_estopt)*(c(hat_score2_estopt_SB)*(1-c(hat_score2_estopt_SB)))
  st <- sum(w_SB_estopt*tt)
  est_vdnew_estopt <- 1/N^2*sum((w_SB_estopt*tt-st/n_SB_estopt)^2)
  
  
  est_vr1=(1/N^2)*sum((1-hat_score_SA)*(Y_SA/hat_score_SA-X_SA_res%*%hat_b1)^2)
  est_vr2=(1/N^2)*sum((1-hat_score_SA)*((Y_SA-mu_ipw2)/hat_score_SA-X_SA_res%*%hat_b2)^2)

  ### estimated variance for new method
  est_vr2new <- (1/N^2)*sum((1-hat_score2_SA)*((Y_SA-mu_new)/hat_score2_SA-X_SA_res%*%hat_bnew)^2)
  
  est_vrnew <- (1/N^2) * sum((hat_score2_SA)^2*(1 - hat_score2_SA)*(w_SA - 1)*(X_SA_res%*%hat_bnew)^2)
  est_var_new <- est_vr2new + est_vdnew + est_vrnew
  
  
  ### estimated variance for conditional method
  est_vr1cond <- (1/N^2) * sum((1 - hat_score_cond_SA)*((Y_SA - mu_cond)/hat_score_cond_SA - (1 - p_i*(1 - 1/w_SA))*(X_SA_res%*%hat_bcond))^2)
  est_vr3cond <- (1/N^2) * sum( (1 - hat_score_cond_SA) * p_i * (1/w_SA)^2 * (w_SA - 1) * (X_SA_res%*%hat_bcond)^2)
  est_var_cond <- est_vr1cond + est_vdcond + est_vr3cond
  
  
  ### vr variance component for split population
  est_vr2new_sp <- (1/N^2)*sum((1-hat_score2_SA_C)*((Y_SA-mu_new_sp)/hat_score2_SA_C-X_SA_res%*%hat_bnew_sp)^2)
  #est_vrnew_sp <- (1/N^2) * sum((hat_score2_SA_C)^2*(1 - hat_score2_SA_C)*(X_SA_res%*%hat_bnew_sp)^2)
  est_vrnew_sp <- (1/N^2) * sum((hat_score2_SA_C)*(1 - hat_score2_SA_C)*(hat_score2_SA_C - 3)*(X_SA_res%*%hat_bnew_sp)^2)
  est_vynew_sp <- (1/N^2) * sum((1/hat_score2_SA_C)*2*(Y_SA - mu_new_sp) * (X_SA_res%*%hat_bnew))
  est_var_new_sp <- est_vr2new_sp + est_vdnew_sp + est_vrnew_sp + est_vynew_sp
  
  
  ### vr2 variance component for optimal probs
  est_vr2new_opt <- (1/N^2)*sum((1-hat_score2_opt_SA)*((Y_SA-mu_new_opt)/hat_score2_opt_SA-X_SA_res%*%hat_bnew_opt)^2)
  
  ### vr2 variance component for estimated optimal probs
  est_vr2new_estopt <- (1/N^2)*sum((1-hat_score2_estopt_SA)*((Y_SA-mu_new_estopt)/hat_score2_estopt_SA-X_SA_res%*%hat_bnew_estopt)^2)
  
  ### estimated variance for new method using optimal probs
  w_SA_opt <- w_opt[index_SA]
  est_vrnew_opt <- (1/N^2) * sum((hat_score2_opt_SA)^2*(1 - hat_score2_opt_SA)*(w_SA_opt - 1)*(X_SA_res%*%hat_bnew_opt)^2)
  est_var_new_opt <- est_vr2new_opt + est_vdnew_opt + est_vrnew_opt
  
  
  ###estimated variances for new method using estimated optimal probs
  w_SA_estopt <- w_estopt[index_SA]
  est_vrnew_estopt <- (1/N^2) * sum((hat_score2_estopt_SA)^2*(1 - hat_score2_estopt_SA)*(w_SA_estopt - 1)*(X_SA_res%*%hat_bnew_estopt)^2)
  est_var_new_estopt <- est_vr2new_estopt + est_vdnew_estopt + est_vrnew_estopt
  
  
  ### estimated variances for IPW1 and IPW2 in CLW
  est_var_ipw1=est_vd1+est_vr1
  est_var_ipw2=est_vd2+est_vr2
  
  
  #bootstrap variance estimator for the EM algorithm
  boot_var_em <- boot_var_em_f(B,Y_SA,X_SA_res,X_SB_res,X_SA_reg,X_SB_reg,w_SB,n_SA,n_SB)
  
  ### Estimated variance
  est_var=c(est_var_ipw1,est_var_ipw2,est_var_new,boot_var_em,est_var_new_sp,est_var_new_opt,est_var_new_estopt,est_var_cond)
  
  
  #sd
  sd_ipw1=sqrt(est_var_ipw1)
  sd_ipw2=sqrt(est_var_ipw2)
  sd_new <- sqrt(est_var_new)
  sd_newp <- sqrt(boot_var_em)
  sd_new_opt <- sqrt(est_var_new_opt)
  sd_new_estopt <- sqrt(est_var_new_estopt)
  sd_new_sp <- sqrt(est_var_new_sp)
  sd_cond <- sqrt(est_var_cond)
  
  
  #CI
  CI_ipw1=c(mu_ipw1-1.96*sd_ipw1,mu_ipw1+1.96*sd_ipw1)
  CI_ipw2=c(mu_ipw2-1.96*sd_ipw2,mu_ipw2+1.96*sd_ipw2)
  CI_new <- c(mu_new - 1.96*sd_new, mu_new + 1.96*sd_new)
  CI_newp <- c(mu_newp - 1.96*sd_newp, mu_newp + 1.96*sd_newp)
  CI_new_opt <- c(mu_new_opt - 1.96*sd_new_opt, mu_new_opt + 1.96*sd_new_opt)
  CI_new_estopt <- c(mu_new_estopt - 1.96*sd_new_estopt, mu_new_estopt + 1.96*sd_new_estopt)
  CI_new_sp <- c(mu_new_sp - 1.96*sd_new_sp, mu_new_sp + 1.96*sd_new_sp)
  CI_cond <- c(mu_cond - 1.96*sd_cond, mu_cond + 1.96*sd_cond)

  ### Average length of the confidence intervals
  al_ipw1=CI_ipw1[2]-CI_ipw1[1]
  al_ipw2=CI_ipw2[2]-CI_ipw2[1]
  al_new=CI_new[2]-CI_new[1]
  al_newp=CI_newp[2]-CI_newp[1]
  al_new_opt=CI_new_opt[2]-CI_new_opt[1]
  al_new_estopt=CI_new_estopt[2]-CI_new_estopt[1]
  al_new_sp=CI_new_sp[2]-CI_new_sp[1]
  al_cond=CI_cond[2]-CI_cond[1]
  
  avelength <- c(al_ipw1,al_ipw2,al_new,al_newp,al_new_opt,al_new_estopt,al_new_sp, al_cond)
  
  
  #cp
  
  CI_ipw1_c=(mu>=CI_ipw1[1])*(mu<=CI_ipw1[2])
  CI_ipw2_c=(mu>=CI_ipw2[1])*(mu<=CI_ipw2[2])
  CI_new_c <- (mu >= CI_new[1])*(mu <= CI_new[2])
  CI_newp_c <- (mu >= CI_newp[1])*(mu <= CI_newp[2])
  CI_new_opt_c <- (mu >= CI_new_opt[1])*(mu <= CI_new_opt[2])
  CI_new_estopt_c <- (mu >= CI_new_estopt[1])*(mu <= CI_new_estopt[2])
  CI_new_sp_c <- (mu >= CI_new_sp[1])*(mu <= CI_new_sp[2])
  CI_cond_c <- (mu >= CI_cond[1])*(mu <= CI_cond[2])
  
  
  c_rate=c(CI_ipw1_c,CI_ipw2_c,CI_new_c, CI_newp_c,CI_new_sp_c, CI_new_opt_c, CI_new_estopt_c, CI_cond_c)
  
  
  if(floor(g/10)==g/10) print(g)
  
  ### Create list of results to output
  list(ests=muhat,est_var=est_var,c_rate=c_rate,avelength=avelength,overlap=overlap_AB)
  
}



endtime <- Sys.time()
endtime-starttime


#################final result###############


sim.ests <- matrix(unlist(simresults[,1]),byrow=T,nrow=iteration)
sim.vars <- matrix(unlist(simresults[,2]), byrow=T,nrow = iteration)
sim.c <- matrix(unlist(simresults[,3]), byrow=T, nrow=iteration)
sim.al <- matrix(unlist(simresults[,4]), byrow=T, nrow=iteration)
sim.overlap <- matrix(unlist(simresults[,5]), byrow=T, nrow=iteration)

### Produce the column means of sim.ests and then bias
mean.ests <- colMeans(sim.ests)
Bias <- mean.ests - mu
RB <- 100*Bias/mu

### Produce estimated MSE using variance methods
est.var <- colMeans(sim.vars)
est.RSE <- 100*sqrt(est.var)/mu

### Produce MC MSE
MC.MSE <- apply(sim.ests,2,var) + Bias^2
MC.RRMSE <- 100*sqrt(MC.MSE)/mu
MC.variance <- MC.MSE - Bias^2
MC.RSE <- 100*sqrt(MC.variance)/mu


### Produce coverage rate
c_rate <- colMeans(sim.c)

### Average length of CIs
al <- colMeans(sim.al)

### Overlap
overlap <- mean(sim.overlap)


###################LATEX##############
library(xtable)
output1=data.frame("Estimator"=c("mu_naive","mu_ipw1","mu_ipw2","mu_new","mu_newp","mu_new_sp","mu_new_opt","mu_new_estopt","mu_wvl","mu_kw","mu_vd","mu_cond","mu_ht"),
                   "RB"=RB,"RRMSE"=MC.RRMSE)
rownames(output1) <- NULL

table_title=paste0("na=",mean_n_SA,",","nb=",mean_n_SB,", ","cor(Y,XB)=",rho*100)
print(xtable(output1, caption = table_title,digits=3),
      caption.placement = 'top',include.colnames =T, include.rownames = F)




