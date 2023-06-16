library(survival)
library(MASS)
library(designmatch)
#library(gurobi)
library(quadprog)
library(survival)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(sbw)
library(ATE)
library(ebal)
library(WeightIt)
library(CBPS)
library(glmnet)
library(cobalt)

#################################
## Beta Bissection
#################################

beta_bissection = function(U, Treatment, target, LP, lambda, eta){

  a_b = c(-50, 40)
  f_a_b = c(Inf,Inf)
  i=0
  convergence = FALSE
  max_iter = 200
  while (i<max_iter && !convergence){
    #print(a_b)
    c = (a_b[1]+a_b[2])/2
    f_a = marginal_beta(betaT=a_b[1], U=U, Treatment=Treatment,
                        LP=LP, lambda=lambda, eta=eta)
    f_b = marginal_beta(betaT=a_b[2], U=U, Treatment=Treatment,
                        LP=LP, lambda=lambda, eta=eta)
    f_c = marginal_beta(betaT=c, U=U, Treatment=Treatment,
                        LP=LP, lambda=lambda, eta=eta)
    f_a = f_a - target
    f_b = f_b - target
    f_c = f_c - target

    if (sign(f_c) == sign(f_a)){
      a_b[1] = c
    }
    if (sign(f_c) == sign(f_b)){
      a_b[2] = c
    }
    if(abs(a_b[1]-a_b[2])< 10^-6){
      convergence = TRUE
    }
    i = i+1
    if(i==max_iter-1){
      stop('Warning: beta did not converge')
    }
  }

  beta_T = c
  return (list(beta_T = beta_T, a_b = a_b))
}

#################################
## Calculates a Marginal Beta
## Given a LP and a conditional beta
#################################

marginal_beta = function(betaT, U, Treatment, LP, lambda, eta){

  n = length(U)
  T1 <- rep(1,n)
  T2 <- rep(0,n)
  Y1 <- ((- U)/(lambda*exp(betaT * T1 + LP))) ^ (1/eta)
  Y2 <- ((- U)/(lambda*exp(betaT * T2 + LP))) ^ (1/eta)
  Y1 = Y1[which(Treatment==1)]
  Y2 = Y2[which(Treatment==1)]
  T1 = T1[which(Treatment==1)]
  T2 = T2[which(Treatment==1)]
  event = rep(1,length(Y1))
  data = data.frame(Y = c(Y1,Y2), T_final = c(T1,T2), event = event)
  fmla <- paste("Surv(time = Y, event = event) ~ T_final ")
  cph <- coxph(as.formula(fmla), data = data)
  return(exp(cph$coefficients))

}

### Function generates data according to a weibull distrubution or uniform
### That can be used to compare and censor
censoring_data = function(censoring_p = 0.5, type = 'uniform',
                          betaT, Treatment, LP, eta, lambda){
  n = length(LP)
  alpha.t = eta
  theta1 <- sqrt(1/lambda)
  lambda_i <- theta1*exp(-(betaT * Treatment + LP)/alpha.t)
  max.lambda_i<-max(lambda_i)

  density.fun.lambda<-function(x){
    pred.y <- predict(y.loess, newdata=x)
    return(pred.y)
  }
  ### compute the density function ###
  h = 1.06*min(sd(lambda_i), IQR(lambda_i)/1.34)*(n^(-0.2))
  dens<-density(lambda_i,n=n,bw=h,from=0,to=max.lambda_i,na.rm=TRUE)
  x<-dens$x
  y<-dens$y
  y.loess <-loess(y~x,span=0.1, control=loess.control(surface="direct"))
  #plot(dens,lty=1,lwd=2)
  #lines(y.loess,col=2)
  #integrate(function(u){predict(y.loess, newdata=u)}, 0, max.lambda_i)

  censor.prop<-function(theta,args){
    censoring_p<-as.numeric(args[1])
    cen.P<-integrate(function(u){
      alpha.t = as.numeric(args[2])
      type = args[3]
      lambda.i<-u
      part1<-density.fun.lambda(lambda.i)
      if (type == 'weibull'){part2 = 1/(1 + ((theta/lambda.i)^alpha.t))}
      if (type == 'uniform'){
        part2 = (lambda.i/(alpha.t*theta))*pgamma((theta/lambda.i)^alpha.t,1/alpha.t)*gamma(1/alpha.t)}
      return(part1*part2)
    },0,max.lambda_i, subdivisions=1000, rel.tol = 2e-4)$value
    return(cen.P-censoring_p)
  }
  args<-c(censoring_p,alpha.t, type)
  theta<-uniroot(censor.prop,args=args,c(0.0001,10000),tol=0.00000001)$root
  if (type == 'weibull'){censor_values = rweibull(n = n, shape = eta, scale = theta)}
  if (type == 'uniform'){censor_values = runif(n = n, min = 0, max = theta)}
  return(censor_values)
}

#################################
## Data generating function
#################################

misspec_surv_data = function(n = 1000, lambda = 0.00002, eta = 2,
                             scenario = 'A',
                             m = 1, k = 1, target_beta = 0.8, logcoef = FALSE,
                             censoring_p = 0, type_censoring = 'uniform',
                             hidden = 0, b_0 = 0) {
  event = rep(1,n)
  X_1 = rbinom(n = n, size = 1, prob = 0.5)
  X_3 = rbinom(n = n, size = 1, prob = 0.5)
  X_5 = rbinom(n = n, size = 1, prob = 0.5)
  X_6 = rbinom(n = n, size = 1, prob = 0.5)
  X_8 = rbinom(n = n, size = 1, prob = 0.5)
  X_9 = rbinom(n = n, size = 1, prob = 0.5)

  X_2 = rnorm(n = n)
  X_4 = rnorm(n = n)
  X_7 = rnorm(n = n)
  X_10 = rnorm(n = n)

  X_22 = X_2^2
  X_44 = X_4^2
  X_77 = X_7^2

  X_13 = X_1*X_3
  X_24 = X_2*X_4
  X_35 = X_3*X_5
  X_46 = X_4*X_6
  X_57 = X_5*X_7
  X_16 = X_1*X_6
  X_23 = X_2*X_3
  X_34 = X_3*X_4
  X_45 = X_4*X_5
  X_56 = X_5*X_6

  b_1 = 0.8 * k
  b_2 = -0.25 * k
  b_3 = 0.6 * k
  b_4 = -0.4 * k
  b_5 = -0.8 * k
  b_6 = -0.5 * k
  b_7 = 0.7 * k

  a_1 = 0.3 * m
  a_2 = -0.36 * m
  a_3 = -0.73 * m
  a_4 = -0.2 * m
  a_5 = 0.71 * m
  a_6 = -0.19 * m
  a_7 = 0.26 * m

  a_h_1 = 0.5 * m
  a_h_2 = 0.5 * m
  b_h_1 = 0.5 * k
  b_h_2 = 0.5 * k

  X_h_1 = rnorm(n = n)
  X_h_2 = rnorm(n = n)

  X = as.matrix(cbind(X_1,X_2,X_3,X_4,X_5,X_6,X_7,X_8,X_9,X_10))

  # Scenario A
  if(scenario == 'A'){
    if (hidden == 0){
      logit_pr_T = b_0 + b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7
    }
    if (hidden == 1){
      logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
        b_h_1*X_h_1
    }
    if (hidden == 2){
      logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
        b_h_1*X_h_1 + b_h_2*X_h_2
    }
  }
  # Scenario B
  if(scenario == 'B'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_2*X_22
  }
  # Scenario C
  if(scenario == 'C'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_2*X_22 + b_4*X_44 + b_7*X_77
  }
  # Scenario D
  if(scenario == 'D'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_1*0.5*X_13 + b_2*0.7*X_24 + b_4*0.5*X_45 + b_5*0.5*X_56
  }
  # Scenario E
  if(scenario == 'E'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_2*X_22 +
      b_1*0.5*X_13 + b_2*0.7*X_24 + b_4*0.5*X_45 + b_5*0.5*X_56
  }
  # Scenario F
  if(scenario == 'F'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_1*0.5*X_13 + b_2*0.7*X_24 + b_3*0.5*X_35 + b_4*0.7*X_46 + b_5*0.5*X_57 +
      b_1*0.5*X_16 + b_2*0.7*X_23 + b_3*0.5*X_34 + b_4*0.5*X_45 + b_5*0.5*X_56
  }
  # Scenario G
  if(scenario == 'G'){
    logit_pr_T = b_1*X_1 + b_2*X_2 + b_3*X_3 + b_4*X_4 + b_5*X_5 + b_6*X_6 + b_7*X_7 +
      b_2*X_22 + b_4*X_44 + b_7*X_77 +
      b_1*0.5*X_13 + b_2*0.7*X_24 + b_3*0.5*X_35 + b_4*0.7*X_46 + b_5*0.5*X_57 +
      b_1*0.5*X_16 + b_2*0.7*X_23 + b_3*0.5*X_34 + b_4*0.5*X_45 + b_5*0.5*X_56
  }
  prob_T <- exp(logit_pr_T) / (1 + exp(logit_pr_T))
  Treatment <- rbinom(n, 1, p = prob_T)

  U <- runif(n)
  U = log(U)
  p = ncol(X)

  if(hidden == 0){
    LP = a_1*X_1 + a_2*X_2 + a_3*X_3 + a_4*X_4 + a_5*X_8 + a_6*X_9 + a_7*X_10
  }
  if(hidden == 1){
    LP = a_1*X_1 + a_2*X_2 + a_3*X_3 + a_4*X_4 + a_5*X_8 + a_6*X_9 + a_7*X_10 + a_h_1*X_h_1 + a_h_2*X_h_2
  }
  if(hidden == 2){
    LP = a_1*X_1 + a_2*X_2 + a_3*X_3 + a_4*X_4 + a_5*X_8 + a_6*X_9 + a_7*X_10 + a_h_1*X_h_1 + a_h_2*X_h_2
  }

  suppressWarnings({betaT = beta_bissection(U=U,Treatment = Treatment, target = target_beta,
                                            LP = LP, lambda = lambda, eta = eta)$beta_T})
  #latent survival time
  Y <- ((- U)/(lambda*exp((betaT * Treatment) + LP))) ^ (1/eta)
  data_out = data.frame(cbind(Y, event, Treatment, X))

  if (censoring_p>0){
    censor_values = censoring_data(censoring_p = censoring_p, type = type_censoring,
                                   betaT = betaT, Treatment = Treatment,
                                   LP = LP, eta = eta, lambda = lambda)
    Ycens <- pmin(data_out$Y, censor_values)
    event <- as.numeric(data_out$Y <= censor_values)
    data_out$event <- event
    data_out$Y <- Ycens
  }

  names(data_out) = c('Y', 'event', 'Treatment', sprintf("X%s",seq(1:p)))

  return(list(data = data_out, betaT = betaT, prob_T = prob_T, target_beta = target_beta))
}



#######################################
## Lasso estimation of propensity score
#######################################

ps_est_func <- function(data, weights = NULL, balance_tresh = 0.1,
                        type_analysis = 'mean_balance') {

  logit <- function(x){
    res <- log(x/ (1 - x))
    return(res)
  }
    col_remove <- which(colnames(data)  %in% c("stage", 'Y', 'event', 'Treatment', 'censored'))
    x_train <- model.matrix( ~. -1, data[ , -c(col_remove)]) #make into matrix
    if (!is.null(weights)) { 
      psmod <- glmnet(x = x_train, y = data$Treatment, family = "binomial",
                      standardize=FALSE, nlambda = 100, weights = weights)
    } else{ psmod <- glmnet(x = x_train, y = data$Treatment, family = "binomial",
                           standardize=FALSE, nlambda = 100)}
    lambda_seq <- psmod$lambda
    lambda_seq = c(lambda_seq,0)
    
    if (!is.null(weights)) { 
      psmod <- glmnet(x = x_train, y = data$Treatment, family = "binomial",
                      standardize=FALSE, lambda = lambda_seq, weights = weights)
    } else{ psmod <- glmnet(x = x_train, y = data$Treatment, family = "binomial",
                            standardize=FALSE, lambda = lambda_seq)}
    
    lambda_seq <- psmod$lambda
    psmod_pred <- predict(psmod, newx = x_train,type = "response")
    uni_ps <- apply(psmod_pred, 2 , function(x) length(unique(x))) #Number of unique values in a column
    if(any(uni_ps == 1)){
      psmod_pred <- psmod_pred[ ,-which(uni_ps == 1)] #Remove who had just 1 unique value
    }
    Rep <- ncol(psmod_pred)
    reslist<-vector("list", Rep)
    mean_balance_values = c()
    abs_balance = c()
    for(i in 1:Rep){
      data2 = data
      col_remove <- which(colnames(data2)  %in% c('Y', 'event', 'censored','weights'))
      data2 <- data2[ , -c(col_remove)]
      lambda_i <- lambda_seq[i]
      data2$ps = psmod_pred[,i]
      data2$weights = ifelse(data2$Treatment == 1, 1, data2$ps/(1-data2$ps))
      col_remove <- which(colnames(data2)  %in% c('ps', 'weights','Treatment'))
      bal_tab = bal.tab(data2[,-c(col_remove)], treat = data2$Treatment ,weights = data2$weights,
                        m.threshold = 0.2, v.threshold = 2, estimand = "ATT"
                        ,method = 'weighting')
      #print(balanced)
      mean_balance_values[i] = mean(abs(bal_tab$Balance$Diff.Adj))
      abs_balance[i] = sum(abs(bal_tab$Balance$Diff.Adj)>=balance_tresh)
    }
    # We choose the minimum because this is returning the
    # scaled difference between treated and non-treated.
    # So we want the case with the lowest difference
    
    if(type_analysis == 'mean_balance'){
      opt_ind = which.min(mean_balance_values)
      glmnet_ps = psmod_pred[,opt_ind]
    }
    
    if(type_analysis == 'both'){
      which_balanced = which(abs_balance == min(abs_balance))
      opt_ind = which.min(mean_balance_values[which_balanced])
      glmnet_ps = psmod_pred[,opt_ind]
    }

    #glm_res = glm(Treatment ~. -Y -event, data = data, family = 'binomial')
    #weights = ifelse(data$Treatment == 1, 1, glm_res$fitted.values/(1-glm_res$fitted.values))
    #bal_tab2 = bal.tab(data[,-c(1,2,3)], treat = data$Treatment ,weights = weights,
    #                  m.threshold = 0.2, v.threshold = 2, estimand = "ATT",
    #                  s.d.denom = 'pooled', method = 'weighting')
    #mean(abs(bal_tab2$Balance$Diff.Adj))
    #psmod <- cv.glmnet(x = x_train, y = data$Treatment, family = "binomial",
    #                   lambda = lambda_seq, standardize=FALSE)
    #best = which(psmod$lambda == psmod$lambda.min)
    #
    #model <- glmnet(x = x_train, y = data$Treatment, family = "binomial",
    #                lambda = psmod$lambda.min, standardize=FALSE)
    #lasso_pred <- predict(model, newx = x_train,type = "response")

  #return(lasso_pred)
  return(glmnet_ps)
}


#####################################################################
## Estimation of the marginal hazard ratio using the weighted sample
#####################################################################

weights_cox_reg = function(weighting_data){

  weighting_data$set <- 1:nrow(weighting_data)
  fmla <- paste("Surv(time = Y, event = event) ~ Treatment + cluster(set)")
  weights = weighting_data$weights
  weighted_coxph <- coxph(as.formula(fmla), data = weighting_data, weights = weights)
  test.ph <- cox.zph(weighted_coxph)

 
  beta <- weighted_coxph$coefficients
  exp_beta <- exp(weighted_coxph$coefficients)
  rob_se <- sqrt(weighted_coxph$var)
  names(beta) <- names(exp_beta) <- names(rob_se) <- NULL
  ciL <- exp(beta - 1.96 * rob_se)
  ciU <- exp(beta + 1.96 * rob_se)
  wald <- (beta / rob_se) ^ 2
  pval <- pchisq(wald, df = 1, lower.tail = FALSE)

  return(list(beta = beta , exp_beta = exp_beta,
              ciL = ciL, ciU = ciU, pval = pval,
              rho = test.ph$table[1], chisq = test.ph$table[2],
              p = test.ph$table[3]))
}




############################################
## Calculate weight balancing for covariates
############################################

balance = function(df, sbw = FALSE){
  col_remove <- which(colnames(df)  %in% c('Y', 'event', 'censored', 'Treatment','weights'))
  Treatment = df$Treatment
  weights = df$weights
  df <- df[ , -c(col_remove)]

  base_balance = c()
  weighted_balance = c()

  for (col in 1:ncol(df)) {
    non_w_mean1 = mean(df[,col][which(Treatment == 1)])
    non_w_mean0 = mean(df[,col][which(Treatment == 0)])
    non_w_var1 = var(df[,col][which(Treatment == 1)])
    non_w_var0 = var(df[,col][which(Treatment == 0)])
    base_balance[col] = (abs(non_w_mean1 - non_w_mean0))/sqrt((non_w_var1+non_w_var0)/2)

    weighted_mean1 = wtd.mean(df[,col][which(Treatment == 1)], weights[which(Treatment == 1)])
    weighted_mean0 = wtd.mean(df[,col][which(Treatment == 0)], weights[which(Treatment == 0)])
    weighted_var1 = wtd.var(df[,col][which(Treatment == 1)], weights[which(Treatment == 1)], normwt = sbw)
    weighted_var0 = wtd.var(df[,col][which(Treatment == 0)], weights[which(Treatment == 0)], normwt = sbw)
    weighted_balance[col] = (abs(weighted_mean1 - weighted_mean0))/sqrt((weighted_var1+weighted_var0)/2)
  }
  return(list(non_weighted = base_balance, weighted = weighted_balance))
}


##########################
## Summarize the results
##########################

results_table <- function(results){
  results_df = as.data.frame(results)
  
  ps_bias <- mean((results_df$exp_beta_PS - results_df$target))
  ps_sd <- sd(results_df$exp_beta_PS)
  ps_rmse <- sqrt((mean(results_df$exp_beta_PS - results_df$target) ^ 2) + var(results_df$exp_beta_PS))
  ps_length_ci <- mean(abs(results_df$ciU_PS- results_df$ciL_PS))
  ps_ciL <- mean(results_df$ciL_PS)
  ps_ciU <- mean(results_df$ciU_PS)
  
  glm_bias <- mean((results_df$exp_beta_GLM - results_df$target))
  glm_sd <- sd(results_df$exp_beta_GLM)
  glm_rmse <- sqrt((mean(results_df$exp_beta_GLM - results_df$target) ^ 2) + var(results_df$exp_beta_GLM))
  glm_length_ci <- mean(abs(results_df$ciU_GLM- results_df$ciL_GLM))
  glm_ciL <- mean(results_df$ciL_GLM)
  glm_ciU <- mean(results_df$ciU_GLM)
  
  raw_bias <- mean((results_df$exp_beta_raw - results_df$target))
  raw_sd <- sd(results_df$exp_beta_raw)
  raw_rmse <- sqrt((mean(results_df$exp_beta_raw - results_df$target) ^ 2) + var(results_df$exp_beta_raw))
  raw_length_ci <- mean(abs(results_df$ciU_raw- results_df$ciL_raw))
  raw_ciL <- mean(results_df$ciL_raw)
  raw_ciU <- mean(results_df$ciU_raw)
  
  ebal_bias <- mean((results_df$exp_beta_ebal - results_df$target))
  ebal_sd <- sd(results_df$exp_beta_ebal)
  ebal_rmse <- sqrt((mean(results_df$exp_beta_ebal - results_df$target) ^ 2) + var(results_df$exp_beta_ebal))
  ebal_length_ci <- mean(abs(results_df$ciU_ebal - results_df$ciL_ebal))
  ebal_ciL <- mean(results_df$ciL_ebal)
  ebal_ciU <- mean(results_df$ciU_ebal)
  
  cal_bias <- mean((results_df$exp_beta_cal - results_df$target))
  cal_sd <- sd(results_df$exp_beta_cal)
  cal_rmse <- sqrt((mean(results_df$exp_beta_cal - results_df$target) ^ 2) + var(results_df$exp_beta_cal))
  cal_length_ci <- mean(abs(results_df$ciU_cal - results_df$ciL_cal))
  cal_ciL <- mean(results_df$ciL_cal)
  cal_ciU <- mean(results_df$ciU_cal)
  
  sbw_bias <- mean((results_df$exp_beta_SBW - results_df$target))
  sbw_sd <- sd(results_df$exp_beta_SBW)
  sbw_rmse <- sqrt((mean(results_df$exp_beta_SBW - results_df$target) ^ 2) + var(results_df$exp_beta_SBW))
  sbw_length_ci <- mean(abs(results_df$ciU_SBW - results_df$ciL_SBW))
  sbw_ciL <- mean(results_df$ciL_SBW)
  sbw_ciU <- mean(results_df$ciU_SBW)
  
  cbps_bias <- mean((results_df$exp_beta_CBPS - results_df$target))
  cbps_sd <- sd(results_df$exp_beta_CBPS)
  cbps_rmse <- sqrt((mean(results_df$exp_beta_CBPS - results_df$target) ^ 2) + var(results_df$exp_beta_CBPS))
  cbps_length_ci <- mean(abs(results_df$ciU_CBPS - results_df$ciL_CBPS))
  cbps_ciL <- mean(results_df$ciL_CBPS)
  cbps_ciU <- mean(results_df$ciU_CBPS)

  

  results_table <- rbind(c(ps_bias, ps_sd, ps_rmse, ps_length_ci, ps_ciL, ps_ciU),
                         c(sbw_bias, sbw_sd, sbw_rmse, sbw_length_ci, sbw_ciL, sbw_ciU),
                         c(glm_bias, glm_sd, glm_rmse, glm_length_ci, glm_ciL, glm_ciU),
                         c(raw_bias, raw_sd, raw_rmse, raw_length_ci, raw_ciL, raw_ciU),
                         c(ebal_bias, ebal_sd, ebal_rmse, ebal_length_ci, ebal_ciL, ebal_ciU),
                         c(cal_bias, cal_sd, cal_rmse, cal_length_ci, cal_ciL, cal_ciU),
                         c(cbps_bias, cbps_sd, cbps_rmse, cbps_length_ci, cbps_ciL, cbps_ciU))
  colnames(results_table) <- c("Bias", "SD", "RMSE", "CIW", "CIL", "CIU")
  rownames(results_table) <- c("IPTW-PS", "IPTW-SBW", "IPTW-GLM", "RAW", "EBAL","CAL",'CBPS')
  return(results_table)

}


